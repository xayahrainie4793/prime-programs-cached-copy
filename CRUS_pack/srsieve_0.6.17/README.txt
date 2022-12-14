The latest source and some binaries can be found at:

  http://www.geocities.com/g_w_reynolds/srsieve/

For discussion check the Sierpinski/Riesel Base 5 forum at

  http://www.mersenneforum.org/

To report bugs, email me (Geoff): g_w_reynolds at yahoo.co.nz.


SRSIEVE
=======

Srsieve was written to speed up sieving for the Sierpinski/Riesel Base 5
projects which look for primes of the form k*5^n+1 and k*5^b-1 for certain
even values of k. However it uses a fairly general algorithm and should work
for any integer sequence in n of the form k*b^n+c, subject to the following
limits:

  1 <=  k  < 2^64, k relatively prime to c.
  2 <=  b  < 2^32, b relatively prime to c.
  0 <=  n  < 2^32-1
  1 <= |c| < 2^63, c relatively prime to b and k.

Note that for sequences of the form k*2^n+1 or k*2^n-1 there are other
programs (NewPGen, Proth-sieve, jjsieve that I know of) which may be faster
than srsieve.

An example to get started:

  $ srsieve --newpgen --nmin 5000 --nmax 10000 --pmax 400000 "24*7^n+1"

This will sieve the single sequence 24*7^n+1 with n in the range 5,000 -
10,000 for all factors up to 400,000 and write the remaining terms in newpgen
format to the file t16_b7_k24.npg. It is equivalent (and should produce an
identical file) to the following NewPGen invocation:

  $ newpgen -wp=t16_k7_b24.npg -t=16 -base=7 -k=24 -nmin=5000 -nmax=10000 \
    -osp=400000


Short forms of the command line switches and integers in exponent notation
are accepted, and any number of sequences (of the same base) can be given on
the command line. The following will sieve the current smallest candidate
sequence from each of the base 5 projects over the same range as above,
writing the remaining terms in the srsieve format file srsieve.out

  $ srsieve -n 5e3 -N 10e3 -P 4e6 "1396*5^n-1" "5114*5^n+1"


To resume sieving where the previous job ended, for primes up to 20,000,000:

  $ srsieve --pmax 20e6 srsieve.out


If you want to sieve many more sequences you can write them all in a text
file one per line and invoke:

  $ srsieve --nmax 1e6 myfile.txt


If --nmin or --pmin are not specified they default to 0. If --pmax is not
specified, sieving will continue until p=pmin+4e12 or interrupted by ctrl-c
or shutdown.

The information srsieve generates can be output in two main forms: a sieve
file (or files) listing the remaining unfactored terms; or a factors file
listing each factor found. A factors file plus checkpoint file plus starting
sieve file (if any) together contain all the information in a final sieve
file. Thus if a factors file is generated, it is not necessary to save a
final sieve file (which could be quite large), and a simple checkpoint file
containing the current sieve prime is sufficient.

The --checkpoint switch will cause a checkpoint file to be written
periodically, and if the --factors switch is also used then a final sieve
file will not be written unless explicitly requested. This sounds a bit
complicated, so here are a few examples, assume sieve.in is a previously
saved sieve file:

  $ srsieve sieve.in

The above will resume sieving from sieve.in, write periodic snapshots of the
sieve and a final sieve (when stopped) to srsieve.out.

  $ srsieve --newpgen sieve.in

The above will resume sieving from sieve.in, write periodic snapshots and a
final sieve to (possibly more than one) NewPGen format file t*_b*_k*.npg.

  $ srsieve --checkpoint --factors sieve.in

The above will resume sieving sieve.in at the prime in checkpoint.txt, if
that file exists, or from sieve.in if not. Periodic checkpoint.txt files
will be written, and each found factor will be appended to srfactors.txt,
but no sieve file will be written. To resume sieving if this job is stopped,
simply repeat the above command. Optionally, the factors in srfactors.txt
may be removed from sieve.in before resuming using the srfile utility below.

  $ srsieve --checkpoint --factors --pfgw sieve.in

The above will act the same as the previous example except that a final
sieve file will also be written in a format suitable for input to PFGW.


[The following option should not normally be needed, as of version 0.50
srsieve takes into account the quadratic character of sequences to avoid
unnecessarily checking primes that cannot be factors].

If you know that the factors p for all terms of all sequences being sieved
satisfy the congruence p = a,b,c (mod X) then you can avoid checking for
other factors with the command line option `--mod=X,a,b,c'. For example the
following will sieve the sequence 5*36^n-1, which has factors that all end
in 1 or 9:

  $ srsieve --mod=10,1,9 --nmax=10e3 "5*36^n-1"


To get a full list of currently supported options, type:

  $ srsieve --help


Beware that recent operating systems might reduce the CPU speed when only
low priority programs are running. If this is not what you want then use the
-Z switch to run sr5sieve at normal priority, or select 'performance mode'
in the operating system or BIOS configuration.



SRFILE
======

The srfile utility can be used to convert sieve files to different formats.
The following will read sieves from myfile1 and myfile2, which can be in any
supported format, and write them to a file sr_5.pfgw suitable for input to
PFGW, i.e. sorted in order of increasing n:

  $ srfile --pfgw myfile1 myfile2


To distribute sieving over a number of different machines, run srsieve with
the --factors switch and with different --pmin and --pmax parameters on each
machine. Then take the sieve file output by one machine and remove all the
factors found by the other machines (in srfactors.txt) with this invocation:

  $ srfile --known-factors srfactors.txt srsieve.out


To remove all terms of the sequence 1396*5^n-1 from the sieve file sieve.out
(used for example when a prime is found by prp testing):

  $ srfile --delete "1396*5^n-1" srsieve.out


To create a Prime95 worktodo.ini file for P-1 factoring all terms in
srsieve.out, with parameters chosen on the assumption that 1 PRP test will
be saved if a factor is found:

  $ srfile --pfactor 1 srsieve.out > worktodo.ini


To get a full list of currently supported options, type:

  $ srfile --help
