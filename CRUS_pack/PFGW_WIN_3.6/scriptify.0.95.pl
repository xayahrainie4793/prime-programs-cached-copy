#!/usr/bin/perl -w
# 0.95 Phil Carmody, 2004/07/21
# Thinking of changing my name to P. Shermani (look it up.)
# Important contributions by Jim Fougeron.
#
# Turns surprisingly C-like syntax into PFGW scripts.

# This is the version of the back end SCRIPT you want to be compatible with
# Is that how you spell compatable?
my $scriptVersion=20031012;

# adds debugging printfs if a certain file exists
my $debug= -f 'scriptifydebugging';

# adds comments to the script file
my $comments=($scriptVersion>=20031012);


sub deEscapeString
{
    $_=$_[0];
    s/\\([\"\\\*\/])/$1/g;
    $_;
}

sub snarf
{
    open(IF, "<$_[0]");
    my @lines=(<IF>);
    close(IF);
    return @lines;
}

print "SCRIPT // autobodgified by scriptify.pl\n";
print ": pre-declare scriptify's globals\n" if($comments);
print "DIMS PCtmpString\n";


my ($isInt,$isString,$isFile,$isLabel,$isFunc)=
    ('integer', 'string', 'file', 'label','function');
my %canonical=(); # mapping of canonical names to original case

# all of these indexed by canonical names
my %vars=();

my %ints=();

my %strings=();

my %files=();
my %filestate=();

my %labels=();
my %labelsneeded=();

my %functions=();
my %functionproto=();
my %functionrets=();
my %functionsneeded=();
my $funcregexp="";
my $functionsloppyregexp='PC\w*';

my $infunction=undef;
my $localvarrewrite=undef;
my $localvardecls=undef;


my %hashhash=('integer' => \%ints,
              'string'  => \%strings,
              'file'    => \%files,
              'label'   => \%labels,
              'function'=> \%functions);

# helpers to check variable existance
sub existOrDie
{
    my ($name, $type, $message)=@_;
    my $canonical=lc($name);
    if(!defined($vars{$canonical}))
    {
        die "$.: Error: Can't $message undeclared $type $name.\n";
    }
    elsif($vars{$canonical} ne $type)
    {
        die "$.: Error: Trying to $message $vars{$canonical} $name from line ${hashhash{$vars{$canonical}}->{$canonical}}, should be type $type.";
    }
    elsif($canonical{$canonical} ne $name)
    {
        print STDERR "$.: Warning: Treating $vars{$canonical} $name as the same as $canonical{$canonical} from line ${hashhash{$vars{$canonical}}->{$canonical}}.\n";
    }
    $canonical;
}
sub existAndDie
{
    my ($name, $type)=@_;
    my $canonical=lc($name);
    if(defined($vars{$canonical}))
    {
        die "$.: Error: Can't declare $name as $type, already a $vars{$canonical} declared at ${hashhash{$vars{$canonical}}->{$canonical}}.\n";
    }
    $canonical{$canonical}=$name;
    $vars{$canonical}=$type;
    $canonical;
}


# Disable redefinition of builtin functions
for my $f('C', 'F', 'GCD', 'IF', 'L', 'len', 'LINEAR', 'LUCASU', 'LUCASV', 'P', 'PHI',
      'PRIMU', 'PRIMV', 'R', 'S', 'SM', 'SMR', 'U', 'V', 'W',
      'PRINT', 'PRP', 'WRITE', 'FACTORIZE',
      'fopen', 'fread', 'read'. 'fclose', 'puts')
{
    my $cf=existAndDie($f, $isFunc);
    $functions{$cf}='pfgw';
    $functionproto{$cf}='...';
}
for my $v('ISPRIME', 'ERRORLEVEL', 'FACTORFOUND', 'MAXF', 'MINF')
{
    my $cv=existAndDie($v, $isInt);
    $ints{$cv}='pfgw';
}
for my $v('MODF')
{
    my $cv=existAndDie($v, $isString);
    $strings{$cv}='pfgw';
}




my $counter='a';
my @context=();
my $breakage=0;
my @pending=();
my $immediatelypending=undef;


# helper function that understands integer sets
sub handleSet
{
    print STDERR "# set\n" if($debug);
    my ($r,$op,$e)=@_;
    existOrDie($r,$isInt,'set');
    # sanity-check apparant user-function calls (which are verboten)
    if($e =~ m/\b($functionsloppyregexp)\s*\(/)
    {
        print STDERR "$.: Warning: Possible use of user-function $1 in arithmetic expresssion.\n";
    }
    if(defined($op))
    {
        if($op eq '<<') { $op="*2^"; }
        elsif($op eq '>>') { $op="/2^"; }
        print "SET $r,$r$op($e)\n";
    }
    else
    {
        print "SET $r,$e\n";
    }
}



while(defined($_=($#pending>=0
          ? shift(@pending)
          : <>)))
{
    # First trim start-of-line whitespace, and EOL
    chomp;
    s/^\s+//; 
    s/^``*//; # emacs fontify gets confused on a single backtick.

    # Did we have a partial line 
    if($immediatelypending)
    { 
        $_="$immediatelypending $_";
        $immediatelypending=undef;
    }

    # Do we still have a partial line?
    if(s/\s*\\$//) 
    {
        print STDERR "# propagating ``$_''\n";
        $immediatelypending=$_;
        next;
    }

    # now trim end-of-line white-space.
    s/\s*$//;

    print STDERR "# ", join('/',@context), " : got ``$_''\n" if($debug);
# Begin pre-processing commands

    # Emit overrides _everything_, and short-circuit processing this line
    if(m/^\#\s*pragma\s+__emit\s+(.*)$/)
    {
        print STDERR "#pragma __emit\n" if($debug);
        # warning!! NO error checking.  The text is simply output to the script.
        print "$1\n";
        next;
    }
    
    # pragma warns short-circuit processing of this line
    if(m/^\#\s*pragma\s+warn\s+(.*)$/)
    {
        print STDERR "$.: \#Warning: $1\n";
        next;
    }

    # Apart from that, remove/spit out comments as the top priority
    while(m/(\/[\/\*])/)
    {
        # what's the first comment, a /* or a // ?
        if($1 eq '/*')
        {
            # block comment - is it terminated?
            if(s@\s*/\*.*?\*/\s*@ @)
            {
                s/^\s+//;
                s/\s+$//;
                print STDERR "# extracted block comment - now got ``$_''\n" if($debug);
            }
            else
            {
                # unterminated block comment - build up abbreviated version
                s@/\*.*@/* ...@;
                $immediatelypending=$_; 
                last; 
            }
        }
        else
        {
            s@\s*//\s*(.*)@@;
            print ": $1\n" if($comments);
        }
    }
    
    # If we have an unterminated block comment - get more _now_
    if($immediatelypending) { next; }

    # and short-circuit an otherwise empty line
    if(!$_) { next; }

    # Include files short-circuit processing of this line
    if(m/^\#\s*include\s*["<]([^\">]+)[>"]$/)
    {
        print STDERR "# include\n" if($debug);
        if(-r $1)
        {
            push @pending, snarf($1);
        }
        else
        {
            die "$.: Error: Can't include file $1\n";
        }
        next;
    }

    # Now process line proper.

    # rewrite local variables if we're in a function
    if($localvarrewrite)
    {
        if(s/\b($localvarrewrite)\b/PC_${infunction}_$1/gi)
        {
            print STDERR "# rewriting to $_\n" if($debug);
        }
    }

    if(!m/[{};:]$/)
    {
        print STDERR "$.: Warning: Almost certainly an error. Missing :, {, }, or ;\n";
    }

# Begin with the declarations
# 
# These typically have an 'exist and die' semantic.

    # Function definition ::= type name ( params ) {
    if(m/^(int(eger)?|string)\s+(\w+)\s*
       \(\s*
       (
        (int(eger)?|string)\s*(&|\s)\s*(\w+)
        (\s*,\s*
         (int(eger)?|string)\s*(&|\s)\s*(\w+)
         )*
        )?
       \s*\)\s*
       (\{|;)$/x)
    {
        my ($rt,$fn,$pl,$def)=($1,$3,$4||'',$14);
        print STDERR "# funcdef :$rt:$fn:$pl:$def:\n" if($debug);
        if($infunction && $def eq '{') 
        {
            die "$.: Error: Cannot nest function $fn within function $infunction.";
        }

        my $cfn=lc($fn);
        my $firsttime=1;

        # Canonicalise the function type
        $rt=~s/\bint\b/integer/g;

        $pl=~s/\bint\b/integer/g;
        $pl=~s/\s\s+/ /g;
        $pl=~s/,\s+/,/g;
        $pl=~s/\s*&\s*(\w+)\s*/ &$1/g; # canonicalise to "type &name"

        # Have we seen anything like this before?
        if($functions{$cfn})
        { 
            # Yes, so this must agree with prior prototype
            $firsttime=0;

            print STDERR "# oldproto $functionrets{$cfn} ($functionproto{$cfn})\n";
            print STDERR "# newproto $rt ($pl)\n";

            if($functionrets{$cfn} ne $rt)
            {
                die "$.: Error: Function $fn\'s return type $rt doesn't match type $functionrets{$cfn} in forward declaration on line $functions{$cfn}.";
            }
            if($functionproto{$cfn} ne $pl)
            {
                die "$.: Error: Function $fn\'s protoype ($pl) doesn't match ($functionproto{$cfn}) in forward declaration on line $functions{$cfn}.";
            }
            if($canonical{$cfn} ne $fn)
            {
                print STDERR "$.: Warning: Treating function $fn as the same as pre-declared function $canonical{$cfn} from line $functions{$cfn}}.\n";
            }
        }
        else
        {
            # This is a new function, or forward declaration
            $cfn=existAndDie($fn, $isFunc);
            $functions{$cfn}=$.;
            $functionsloppyregexp.="|$fn";
        }

        # If we've got a declaration, then we need to do variable rewrites
        if($def eq '{')
        {
            $localvarrewrite='';
            $localvardecls=$comments?": local decls for function $fn\n":'';
        }

        my $funcregexpsadd;
        if($firsttime)
        {
            $funcregexpsadd='';
            print "\n:: $rt $fn($pl) $def\n" if($comments);
            print "DIM", ($rt eq 'string'?'S':''), " PCret_$cfn\n";
        }

        my @pl=split(/\s*,\s*/,$pl);
        for my $param(@pl)
        {
            print STDERR "# - param $param\n" if($debug);
            my ($t,$r,$n)=($param =~ m/(\w+)\s*(&|\s)\s*(\w+)/);
            
            # First check for shadowing, this is not a fatal error.
            my $cn=lc($n);
            if(defined($vars{$cn}))
            {
                # die "$.: Error: Local variable $t $n must have unique name";
                print STDERR "$.: Warning: Parameter $t $n shadows global $vars{$cn} $canonical{$cn}.\n";
            }
            
            $localvarrewrite .= "|$n" if($def eq '{');

            if($firsttime)
            {
                # Now look at the _real_ local parameter, this may die.
                my $newn="PC_${fn}_$n";
                $cn=existAndDie($newn, $t) if($firsttime);
                $funcregexpsadd .= '\s*,\s*';
                my $typeletter;

                if($t eq 'string') 
                { 
                    $strings{$cn}="$.\@$fn"; 
                    $funcregexpsadd .= '(\w+|"(|.*[^\\\\](\\\\\\\\)*)")';
                    $typeletter='S';
                }
                else 
                { 
                    #!!!!!!
                    # Phil - permit parameters to be expressions.
                    $ints{$cn}="$.\@$fn"; 
                    $funcregexpsadd .= '(\w+)';
                    $typeletter='';
                }
                print "DIM$typeletter $newn\n";
            }
        }

        if($def eq '{')
        {
            $infunction=$fn;
            $localvarrewrite=substr($localvarrewrite,1) if($localvarrewrite);
            unshift @context, "PCfunc_$fn";
            print "GOTO PCfunc_${fn}_end\n";
            print ":\n" if($comments);
            print "LABEL PCfunc_$fn\n";
            print STDERR "# local var regexp = $localvarrewrite\n" if($debug);
            
            if($functionsneeded{$cfn}) { delete($functionsneeded{$cfn}); }
        }

        if($firsttime)
        {
            $functionrets{$cfn}=$rt;
            $functionproto{$cfn}=$pl;
            $funcregexp .= ($funcregexp?'|':'') 
                . $fn
                . '\s*\(\s*'
                . ($funcregexpsadd?substr($funcregexpsadd,7).'\s*':'')
                . '\)';
            print STDERR "# funcy regexp = $funcregexp\n" if($debug);
        }
        if($def eq ';' && !$firsttime)
        {
            print STDERR "$.: Warning: Forward declaration of $fn, while valid, duplicates that from line $functions{$cfn}.\n";
        }
    }
    elsif(m/^int(eger)?\s+(\w+)\s*(=\s*(.*\S))?\s*;$/)
    {
        print STDERR "#",$infunction?" local":''," decl int\n" if($debug);
        my ($n,$e)=($2,$4||'');
        if($infunction)
        {
            $ln="PC_${infunction}_$n";
            my $cn=existAndDie($ln,$isInt);
            $ints{$cn}="$.\@$infunction";
            # don't declare anything, however, initialise
            if($e) { print "SET $ln,$e\n"; }
            $localvarrewrite .= ($localvarrewrite?'|':'') . $n;
            $localvardecls .= "DIM $ln;\n";
            print STDERR "# local var regexp = $localvarrewrite\n" if($debug);
        }
        else
        {
            my $cn=existAndDie($n,$isInt);
            $ints{$cn}=$.;
            print "DIM $n",$e?",$e\n":"\n";
        }
    }
    elsif(m/^string\s+(\w+)(\s*=\s*"(|.*[^\\](\\\\)*)")?\s*;$/)
    {
        print STDERR "#",$infunction?" local":''," decl str\n" if($debug);
        my ($n,$e,$ee)=($1,$2,deEscapeString($3||''));
        if($infunction)
        {
            $ln="PC_${infunction}_$n";
            my $cn=existAndDie($ln,$isString);
            $strings{$cn}="$.\@$infunction";
            # don't declare anything, however, initialise
            if($e) { print "SETS $ln,$ee\n"; }
            $localvarrewrite .= ($localvarrewrite?'|':'') . $n;
            $localvardecls .= "DIMS $ln;\n";
            print STDERR "# local var regexp = $localvarrewrite\n" if($debug);
        }
        else
        {
            my $cn=existAndDie($n,$isString);
            $strings{$cn}=$.;
            if(defined($e)) 
            {
                if($scriptVersion>=20031012)
                {
                    print "DIMS $1,$ee\n"; 
                }
                else
                {
                    print "DIMS $1\n";
                    print "SETS $1,%s;$ee\n";
                }
            }
            else
            {
                print "DIMS $1\n"; 
            }
        }
    }
    elsif(m/^file\s+(\w+)(\s*=\s*fopen\s*\(\s*("([^ ,\"])+"|\w+)\s*,\s*"([wra])t?"\s*\))?\s*;$/)
    {
        print STDERR "# decl+fopen\n" if($debug);
	my $fv=$1;
        my $cn=existAndDie($fv,$isFile);
        if($2)
        {
	    my($fn, $m)=($3,$5);
	    if($fn=~s/"([^ ,\"]+)"/$1/) 
	    {
		# literal string
	    }
	    else
	    {
		# variable
		$fn=existOrDie($fn,$isString,'fopen');
	    }
            if($m eq 'r')
            {
                $files{$cn}=$.;
                $filestate{$cn}=0;
                print "OPENFILEIN $fv,$fn\n";
            }
            elsif($m eq 'w')
            {
                $files{$cn}=$.;
                $filestate{$cn}=1;
                print "OPENFILEOUT $fv,$fn\n";
            }
            elsif($m eq 'a')
            {
                $files{$cn}=$.;
                $filestate{$cn}=2;
                print "OPENFILEAPP $fv,$fn\n";
            }
        }
        else
        {
            $files{$cn}=$.;
            $filestate{$cn}=-1;
        }
    }
    elsif(m/^(\w+)\s*:$/)
    {
        print STDERR "# label\n" if($debug);
        my $cn=existAndDie($1,$isLabel);
        $labels{$cn}=$.;
        if(defined($labelsneeded{$cn})) { delete($labelsneeded{$cn}); }
        print "LABEL $1\n";
    }

# End of 'declarations'
  
# Start of Function calls
#
# These typically have 'exist or die' semantics, except where noted.

    elsif(m/^(\w+)\s*=\s*fopen\s*\(\s*("([^ ,\"])+"|\w+)\s*,\s*"([rwa])t?"\s*\)\s*;$/)
    {
        print STDERR "# fopen2\n" if($debug);
	my ($fv,$fn,$m)=($1,$2,$4);
        my $cn=existOrDie($fv,$isFile,'fopen');
        if($filestate{$cn}!=-1)
        {
            print STDERR "$.: Warning: possible use of already open file $fv.\n";
        }
	if($fn=~s/"([^ ,\"]+)"/$1/) 
	{
	    # literal string
	}
	else
	{
	    # variable
	    $fn=existOrDie($fn,$isString,'fopen');
	}
        if($m eq 'r')
        {
            $filestate{$cn}=0;
            print "OPENFILEIN $fv,$fn\n";
        }
        elsif($m eq 'w')
        {
            $filestate{$cn}=1;
            print "OPENFILEOUT $fv,$fn\n";
        }
        elsif($m eq 'a')
        {
            $filestate{$cn}=2;
            print "OPENFILEAPP $fv,$fn\n";
        }
    }
    elsif(m/^(\w+)\s*=\s*f?read\s*\(\s*(\w+)\s*(,\s*&\s*(\w+)\s*)?\)\s*;$/)
    {
        print STDERR "# fread\n" if($debug);
        existOrDie($1,$isInt,'read into');
        my $cfn=existOrDie($2,$isFile,'read from');
        if(defined($3))
        {
            existOrDie($4,$isString,'read into (reference to)');
        }
        # Sanity check the file's state
        if($filestate{$cfn}==-1)
        {
            print STDERR "$.: Warning: possible read of unopened file $2.\n";
        }
        elsif($filestate{$cfn}!=0)
        {
            print STDERR "$.: Warning: possible read of write/append file $2.\n";
        }
        print "GETNEXT $1,$2",$3?",$4\n":"\n";
    }
    elsif(m/^puts\((\w+)\)\s*;$/)
    {
        print STDERR "# puts\n" if($debug);
        existOrDie($1,$isString,'puts');
        print "PRINT $1\n";
    }
    elsif(m/^print\((\w+)\)\s*;$/)
    {
        print STDERR "# print\n" if($debug);
        existOrDie($1,$isInt,'print');
        print "PRINT $1\n";
    }
    elsif(m/^printf\s*\(\s*"(|.*?[^\\](\\\\)*)"\s*(,.*)?\)\s*;$/)
    {
        # !!!!!!
        # Phil - permit format string to be a string variable
        my ($fstring,$varlist)=($1,$3||'');
        print STDERR "# printf:$1:$varlist \n" if($debug);
        $varlist=join(';',map{s/"(.*)"/$1/;$_}split(/\s*,\s*/,$varlist));
        if($fstring=~s/([^\\]);/${1}£/g)
        {
            print STDERR "$.: Warning: literal ';' found in format string. Uglifying.\n";
        }
	# By default, we don't want a newline, so use ",c"
	my $modifier=',c';
	# Could replace a trailing \n with nothing, dropping the ,c 
	#if($fstring =~ s/^((|.*[^\\])(\\\\)*)\\n$/$1/)
	#{
	#    $modifier='';
	#}
        #$fstring=deEscapeString($fstring);
        print "SETS PCtmpString,$fstring$varlist\n";
        print "PRINT PCtmpString$modifier\n";
    }
    elsif(m/^([sf])printf\s*\(\s*(\w+)\s*,\s*"(|.*[^\\](\\\\)*)"\s*(,.*)?\)\s*;$/)
    {
        # !!!!!!
        # Phil - permit format string to be a string variable
        my ($dest,$varname,$fstring,$varlist)=($1,$2,$3,$5||'');
        print STDERR "# $1printf:$2:$3:$varlist\n" if($debug);
        if($dest eq 's')
        {
            existOrDie($varname, $isString, 'sprintf');
        }
        else
        {
            existOrDie($varname, $isFile, 'fprintf');
        }
        $varlist=~s/,\s*/;/g;
        if($fstring=~s/([^\\]);/${1}£/g)
        {
            print STDERR "$.: Warning: literal ';' found in format string. Uglifying.\n";
        }
        #$fstring=deEscapeString($fstring);
        if($dest eq 's')
        {
            print "SETS $varname,$fstring$varlist\n";
        }
        else
        {
            print ": synthesise fprintf\n" if($comments);
            print "SETS PCtmpString,$fstring$varlist\n";
            print "WRITE $varname,PCtmpString\n";
        }
    }
    elsif(m/^write\s*\((\w+)\s*,\s*(\w+)\s*\)\s*;$/)
    {
        print STDERR "# write\n" if($debug);
        my $cf=existOrDie($1,$isFile,'write to');
        if($filestate{$cf}==-1)
        {
            print STDERR "$.: Warning: possible write to unopened file $1.\n";
        }
        elsif($filestate{$cf}==0)
        {
            print STDERR "$.: Warning: possible write/append to file opened for reading $1.\n";
        }
        
        if(!defined($ints{lc($2)}) && !defined($strings{lc($2)}))
        {
            die "$.: Error: Can't write undeclared variable $1.";
        }
        
        print "WRITE $1,$2\n";
    }
    elsif(m/^fclose\s*\((\w+)\)\s*;$/)
    {
        print STDERR "# fclose\n" if($debug);
        my $cf=existOrDie($1,$isFile,'fclose');
        if($filestate{$cf}==-1)
        {
            print STDERR "$.: Warning: possible close of unopened file $1.\n";
        }
        $filestate{$cf}=-1;
        print "CLOSEFILE $1\n";
    }
    elsif(m/^((\w+)\s*=\s*)?system\((\w+)\)\s*;$/)
    {
        print STDERR "# system\n" if($debug);
        if($1) { existOrDie($2,$isInt,'put system return into'); }
        existOrDie($3,$isString,'system');
        print "SYSTEM $3\n";
        if($1) { print "SET $2,ERRORLEVEL\n"; }
    }
    elsif(m/^(\w+)\s*=\s*powmod\s*\(\s*([^,]+?)\s*,\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)\s*;$/)
    {
        print STDERR "# powmod\n" if($debug);
        my ($r,$b,$n,$p)=($1,$2,$3,$4);
        existOrDie($r,$isInt,'put powmod return into');
        # Very basic attempts to detect user errors. 
        # This is not an expression parser, but if it looks like
        # a variable, then we can at least check that it is one.
        if($b!~m/^\d+$/ && $b=~m/^\w+$/) 
        { existOrDie($b,$isInt, 'pass to powmod'); }
        if($n!~m/^\d+$/ && $n=~m/^\w+$/) 
        { existOrDie($n,$isInt, 'pass to powmod'); }
        if($p!~m/^\d+$/ && $p=~m/^\w+$/) 
        { existOrDie($p,$isInt, 'pass to powmod'); }

        print "POWMOD $r,$b,$n,$p\n";
    }
    elsif(m/^((\w+)\s*=\s*)?factori[zs]e\s*\(\s*([^,]+)\s*(,\s*([^,]+)\s*(,\s*([^,]+))?\s*)?\)\s*;$/)
    {
        print STDERR "# factorize\n" if($debug);
        my ($f,$n,$l,$ll,$h,$hh)=($2,$3,$4,$5,$6,$7);
        if(defined($f)) { existOrDie($f,$isInt,'put factorize return into'); }
        if($n!~m/^\d+$/ && $n=~m/^\w+$/) 
        { existOrDie($n,$isInt, 'factorize'); }
        # do we have bounds?
        if(defined($l))
        {
            if($ll!~m/^\d+$/ && $ll=~m/^\w+$/) 
            { existOrDie($ll,$isInt, 'pass to factorize'); }
            print "SET MINF,$ll\n";
            if(defined($h))
            {
                if($hh!~m/^\d+$/ && $hh=~m/^\w+$/) 
                { existOrDie($hh,$isInt, 'pass to factorize'); }
                print "SET MAXF,$hh\n";
            }
        }
        print "FACTORIZE $n\n";
        if($f) { print "SET $f,FACTORFOUND\n"; }
    }
    elsif(m/^((\w+)\s*=\s*)?PRP\s*\(\s*([^,\)]+)\s*(\s*,\s*(\w+)\s*)?\)\s*;$/)
    {
        # !!!!!!
        # Phil - Merge with following literal version
        print STDERR "# PRP\n" if($debug);
        if($1) { existOrDie($2,$isInt,'put PRP return into'); }
        if($4) { existOrDie($5,$isString, 'supply string to PRP'); }
        print "PRP $3",$4?",$5\n":"\n";
        if($1) { print "SET $2,ISPRIME\n"; }
    }
    elsif(m/^((\w+)\s*=\s*)?PRP\s*\(\s*([^,\)]+)\s*,\s*"(|.*[^\\](\\\\)*)"\s*\)\s*;$/)
    {
        print STDERR "# PRP2\n" if($debug);
        if($1) { existOrDie($2,$isInt,'put PRP return into'); }
        my $literal=deEscapeString($4);
        print "SETS PCtmpString,$literal\n";
        print "PRP $3,PCtmpString\n";
        if($1) { print "SET $2,ISPRIME\n"; }
    }
    elsif(m/^((\w+)\s*=\s*)?(atoi|eval)\s*\(\s*("([^ ,\"])+"|\w+)\s*\)\s*;$/)
    {
        print STDERR "# atoi\n" if($debug);
        if($1) { existOrDie($2,$isInt,'atoi into'); }
        if(substr($4,0,1) eq '"')
        {
            my $literal=deEscapeString(substr($4,1,-1));
            print "SETS PCtmpString,$literal\n";
            print "STRTOINT $2,PCtmpString\n";
        }
        else
        {
            existOrDie($4,$isString,'atoi from');
            print "STRTOINT $2,$4\n";
        }
    }
    elsif($funcregexp && 
          m/^((\w+)\s*=\s*)?($funcregexp)\s*;$/) # reevaluate it each time
    {
        print STDERR "# function call\n" if($debug);
        my ($rv,$call)=($1?$2:undef,$3);
        my ($fn,$params)=($call=~m/^(\w+)\s*\(\s*(.*\S|)\s*\)$/);
        my $cfn=existOrDie($fn,$isFunc,'call user-function');
        if($rv) { existOrDie($rv,$functionrets{$cfn},"put $fn return into"); }

        my $refsout='';
        if(defined($params))
        {
            my @params=split(/\s*,\s*/,$params);
            my @args=split(/\s*,\s*/,$functionproto{$cfn});
            for my $i (0..$#params)
            {
                my $pn=$params[$i];
                my ($at,$ar,$an)=($args[$i] =~ m/(\w+)\s*(&|\s)\s*(\w+)/);
                # print STDERR "# - is $pn of type $at ($ar) for $an\n" if($debug);
                if($at eq 'integer')
                {
                    # literal numeric tokens permitted
                    if($pn=~m/^\d+$/)
                    {
                        if($ar eq '&')
                        {
                            die "$. Error: $fn parameter $an is a reference, can't be passed literal $pn.\n";
                        }
                    }
                    else
                    {
                        existOrDie($pn,$isInt,'pass as integer parameter');
                        if($ar eq '&')
                        {
                            $refsout .= "SET $pn,PC_${fn}_$an\n";
                        }
                    }
                    print "SET PC_${fn}_$an,$pn\n";
                }
                else
                {
                    #print STDERR "# - $an is string\n";
                    if($pn=~m/^"(|.*[^\\])"$/)
                    {
                        if($ar eq '&')
                        {
                            die "$. Error: $fn parameter $an is a reference, can't be passed literal $pn.\n";
                        }
                        print "SETS PC_${fn}_$an,$1\n";
                    }   
                    else
                    {
                        existOrDie($pn,$isString,'pass as string parameter');
                        print "SETS PC_${fn}_$an,%s;$pn\n";
                        if($ar eq '&')
                        {
                            $refsout .= "SETS $pn,%s;PC_${fn}_$an\n";
                        }
                    }
                }
            }
        }
        print "GOSUB PCfunc_$fn\n";
        if(defined($rv))
        {
            print "SET", ($functionrets{$cfn} eq 'string'?'S':''), " $rv,PCret_$fn\n";
        }
        if($refsout) { print $refsout; }
    }
    elsif(m/^(\w+)\s*(\+|-|<<|>>|\*|\/|\%)?=\s*([^\"]+)\s*;$/)
    {
        handleSet($1,$2,$3);
    }
    elsif(m/^(\w+)\s*=\s*"(|.*[^\\](\\\\)*)"\s*$/)
    {
        print STDERR "# set str\n" if($debug);
        existOrDie($1,$isString,'set from literal string');
        my $literal=deEscapeString($2);
        print "SETS $1,$literal\n";
    }

# End of function calls

# Begin flow control
#
# These typically have 'exists or die' semantics

    elsif(m/^goto\s*(\w+)\s*;$/)
    {
        print STDERR "# goto $1\n" if($debug);
        # Labels are permitted to not pre-exist
        my $cn=lc($1);
        # if we don't know the label, remember to check it's eventually declared
        if(defined($labels{$cn}) || defined($labelsneeded{$cn}))
        {
            if($1 ne $canonical{$cn})
            {
                print STDERR "$.: Warning: Treating label $1 as the same as label $canonical{$cn} from line $labels{$cn}.\n";
            }
        }
        else
        {
            print STDERR "# - will later need label $1\n" if($debug);
            $labelsneeded{$cn}=$.;
            $canonical{$cn}=$1;
        }
        print "GOTO $1\n";
    }
    elsif(m/^return\s*\((.*)\)\s*;$/ ||
          m/return\s+(.*)\s*;$/)
    {
        if(!$infunction)
        {
            print STDERR "# global return\n" if($debug);
            print "END\n";
        }
        else
        {
            print STDERR "# function return\n" if($debug);
            if(defined($1))
            {
                print 'SET', ($functionrets{lc($infunction)} eq 'string' ? 'S':''), " PCret_$infunction,$1\n";
            }
            print "RETURN\n";
        }
    }
    elsif(m/^exit\s*\(.*\)\s*;$/)
    {
        print STDERR "# exit\n" if($debug);
        print "END\n";
    }
    elsif(m/^\}$/)
    {
        print STDERR "# close brace\n" if($debug);
        # must end an if/else or a while or a do.
        if($#context<0)
        {
            die "$.: Error: \} with no context.";
        }
        my $location=shift @context;
        if($location =~ m/^PCfunc_/)
        {
            print "RETURN\n";
            print "LABEL ${location}_end\n";
            if($localvardecls) { print $localvardecls; }
            print "::\n\n" if($comments);
            my $pl=$functionproto{lc($infunction)};
            # !!!!!! 
            # Phil - why are we removing the local variables?
            if(0)
            {
                my @pl=split(/\s*,\s*/, $pl);
                for$param(@pl)
                {
                    my ($t,$r,$n)=($param =~ m/(\w+)\s*(&|\s)\s*(\w+)/);
                    $newn="PC_${infunction}_$n";
                    print STDERR "# forgetting locals $newn = $t$r$n\n" if($debug);
                    delete($vars{$n});
                    delete($canonical{$n});
                    if($t eq 'string') { delete($strings{$n}); }
                    else { delete($ints{$n}); }
                }
            }
            $infunction=undef;
            $localvarrewrite=undef;
            $localvardecls=undef;
        }
        elsif($location =~ m/^PCwhile_/)
        {
            # While:
            # always jump to the top condition.
            # always need an end label, for when that condition fails
            print "GOTO $location\n";
            print "LABEL ${location}_end\n";
            $breakage>>=1;
        }
        elsif($location =~ m/^PCfor_/)
        {
            # for:
            # perform the post op (if there is one)
            # always jump to the top condition.
            # always need an end label, for when that condition fails
            unshift @pending, "goto $location;\n", "${location}_end:\n";
            my $forop = shift @fortailcmd;
            if($forop) { unshift @pending, "$forop;\n"; }
            $breakage>>=1;
        }
        else
        {
            # This is the end of an if or an else.
            # they both always needs an end label, 
            # However, they're already name munged
            print "LABEL $location\n";
        }
    }
    elsif(m/^if\s*(\(\s*[^\)]+\))\s*\{$/)
    {
        print STDERR "# if\n" if($debug);
        unshift @context, "PCnotif_$counter";
        print "IF !$1 THEN GOTO PCnotif_$counter\n";
        ++$counter;
    }
    elsif(m/^}\s*else\s*{$/)
    {
        print STDERR "# else\n" if($debug);
        if($#context<0)
        {
            die "$.: Error: \} else \{ with no context.";
        }
        my $location=shift @context;
        unshift @context, "PCendelse_$counter\n";
        print "GOTO PCendelse_$counter\n";
        ++$counter;
        print "LABEL $location\n";
    }
    elsif(m/^do\s*\{$/)
    {
        print STDERR "# do\n" if($debug);
        unshift @context, "PCdo_$counter";
        print "LABEL PCdo_$counter\n";
        ++$counter;
        $breakage<<=1;
    }
    elsif(m/^\}\s*while\s*(\(\s*[^\)]+\))\s*;$/)
    {
        print STDERR "# do's while\n" if($debug);
        if($#context<0)
        {
            die "$.: Error: \} while with no context.";
        }
        my $location=shift @context;
        print "IF $1 THEN GOTO $location\n";
        # As do's are bottom-conditioned, they only need a label
        # if there's been someone trying to break out of the loop.
        if($breakage&1) 
        {
            # print STDERR "$.: Info: Must break $location\n";
            print "LABEL ${location}_end\n";
        }
        $breakage>>=1;
    }
    elsif(m/^for\s*\(\s*([^\;]*?)\s*;\s*([^\;]*?)\s*;\s*([^\)]*)\s*\)\s*\{$/)
    {
        print STDERR "# for ($1;$2;$3) {\n" if ($debug);
        print ": for ($1;$2;$3) {\n" if($comments);
        unshift @context, "PCfor_$counter";
        # Create the label so that we can use the engine to call a goto.
        my $cl=existAndDie("PCfor_$counter",$isLabel,"set up for loop");
        $labels{$cl}=$.;

        unshift @fortailcmd, $3;
        my $my2 = $2;

        if($1) # is first part non-empty?
        {
            if($1=~m/^(\w+)\s*(\+|-|<<|>>|\*|\/|\%)?=\s*([^\"]+)\s*/)
            {
                handleSet($1,$2,$3);
            }
            else
            {
                die "$.:  Invalid first param in for() loop\n";
            }
        }
        else
        {
            print STDERR "# blank first param in for() loop\n" if($debug);
        }

        print "LABEL PCfor_$counter\n";
        if (!$my2) # is 2nd part empty
        {
           print ": empty compare section of the for loop\n" if($comments);
        }
        else
        {
           print "IF !($my2) THEN GOTO PCfor_${counter}_end\n";
        }
        ++$counter;
        $breakage<<=1;
    }
    elsif(m/^while\s*(\(\s*[^\)]+\))\s*\{$/)
    {
        print STDERR "# while\n" if($debug);
        unshift @context, "PCwhile_$counter";
        print "LABEL PCwhile_$counter\n";
        print "IF !$1 THEN GOTO PCwhile_${counter}_end\n";
        ++$counter;
        $breakage<<=1;
    }
    elsif(m/^break(\s*\((\d+)\))?\s*;$/)
    {
        print STDERR "# break\n" if($debug);
        my $find=$2?$3:1;
        $breakage|=1<<($find-1);
        my $loop=0;
        my $location=undef;

        while(defined($context[$loop]))
        {
            if($context[$loop]=~m/func/) { last; }
            if($context[$loop]=~m/do|while/ && !--$find)
            {
                $location=$context[$loop];
                last;
            }
            ++$loop;
        }
        if($location)
        {
            print "GOTO ${location}_end\n";
        }
        else
        {
            die "$.: Error: Can't break requested number of loops.\n";
        }
    }
# End of flow control

    elsif($_) # were we trying to parse anything at all?
    {
        die "$.: Error: Can't parse line.";
    }
}

foreach(keys(%labelsneeded))
{
    # can only die once.
    die "$labelsneeded{$_}: Error: goto target ``$_'' not found.";
}
foreach(keys(%functionsneeded))
{
    # can only die once.
    die "$functionsneeded{$_}: Error: function ``$_'' not found.";
}
if($#context>=0)
{
    die 'EOF: Error: Context=('.join('|',@context).').';
}
if($breakage)
{
    die "EOF: Error: Breakage=$breakage.";
}

# finish the script with an explicit end.
print "END\n";


# 0.3
# Added \" handling to string initialise and assign
# Added \" handling to sprintf
# Increased variable/label name sanity-checking
# Various other minor bugfixes

# 0.4
# Added functions
# Increases label/variable name sanity checking
# Added unresolved goto checking
# use DIMS s,val syntax
# stuck debugging prints in
# Make tokens as well as vars legal as parameters.

# 0.5
# Added reference parameters.
# Added local scoping of parameters
# Added atoi -> STRTOINT

# 0.6
# Factored out exist(And|Or)Die code for variable checking
# Leave ';' at end of each line, add to each regexp that needs it
# Jim's __emit concept added
# Case insensitivity of all variables/labels added
# pragma warn added
# Silly bugs - "_end" labels.
# powmod and factorize added

# 0.7
# Local variables for functions using unique renaming.
# Jim added for(;;) with almost full C semantics!
# I shrank it by:
# - factoring out the a=b handling code into a helper function
# - using @pending to inject goto-creation code at the }
# Multi-line comments added, using $immediatelypending
# Silly bugs - SETS syntax wrong, $label{}, \<\< not needed in regexp

# 0.8
# Forward declarations
# Improved literal string handling to accept \\ and \" escapes
# Silly bugs - \s* added to a few regexps. printf handles literal string params.

# 0.9
# Silly bugs - comparing strings against '', as !$s is true when $s='0'.
