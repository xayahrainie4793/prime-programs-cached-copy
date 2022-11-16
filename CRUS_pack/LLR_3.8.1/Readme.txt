
		Welcome to LLR program Version 3.8.1 !

N.B. : Three bugs from version 3.8.0 are now fixed in this version :

- The PRP tests of k*b^n+c numbers with abs(c) != 1 got false negative results
due to a misusage of gwnum's "gwsetaddin" function, they are now correct!

- One "pathological" k value caused a false negative result in some very rare
cases, when using gwnum V25.13 ; this bug has disappeared with V25.14

- Round off error recovery has been improved, and infinite error loops are
now avoided (always, I hope...) by restarting with next larger FFT length.

- This version is identical to the development one at the date of July 1 2010
It will not be further changed without changing the version number.


1) Main features :

  LLR is a primality proving program for numbers of the form N = k*b^n +/- 1,
  (k < b^n), or numbers which can be rewritten in this form, like 
  Gaussian-Mersenne norms.
  It can also do strong and efficient PRP tests on more general k*b^n+c forms,
  on Wagstaff numbers (2^p+1)/3, repunits (10^p-1)/9 and generalized repunits
  (b^p-1)/(b-1), b!=2.

  This version uses the most recent release (25.14) of George Woltman's Gwnum
  library, to do fast multiplications and squarings of large integers modulo N.

  Since version 25.11, gwnum library is no longer restricted to base two for
  efficient computing modulo k*b^n+c numbers (but SSE2 is required if b != 2), 
  and LLR greatly takes advantage of this improvement.

  LLR can run on all machines where gwnum code can run, so, on all Intel x86
  compatible machines.

2) User interfaces of LLR :

--> The command mode applications can run on Windows, Linux, FreeBSD and
    Mac OS X platforms. All have exactly the same directions for use.

--> The "llr.exe" GUI application is only available on Windows.

--> "cllr.exe" is the Windows console application.

--> "llr" is the dynamically linked application on Linux, BSD and Mac OS X.

--> "sllr" is the same one, but statically linked.

3) How to use LLR :

  - Windows GUI application (llr.exe) :

    You may simply double-click on the icone of the application, and then 
    use the scrolling menues. Then, the options, and results and .ini 
    file names have the default values.

    Note : This text can be displayed in the GUI by selecting "Help Contents" !

    In order to get better flexibility, you rather should launch LLR from a
    command prompt window :
    
    >llr [-a<nnnn>] [-w<directory>] [-o<keyword>=<value>]...

    - With -a, you choose another set of .ini and results file.
    - With -w, you choose another working directory.
    - with -o, you set one or more user options(10 max.) in the .ini file.
    - and then, use the GUI...

  - Command line application (cllr.exe, llr or sllr)

    - First, you may type : llr -h to get the online help information...

    - To use the program in BATCH MODE, type (for example, with cllr) :
    
    >cllr [-d] [-a<nnnn>] [-w<dir.>] [-o<keyw.>=<val.>]... <input file name>

    (If you ommit -d, the program will work silently!)

    - To test a SINGLE k*b^n+c number, type :

    >cllr [-d] [-a<nnnn>] [-w<dir.>] [-o<keyw.>=<val.>]... -q"expression"

    - To use the program in INTERACTIVE MODE, type :

    >cllr -m [-a<nnnn>] [-w<dir.>] [-o<keyw.>=<val.>]...

    then, you get the main menu, and continue according to your choice(s).

4) Input, output and intermediate data :

    LLR can take its INPUT DATA from Newpgen output files, and also from
    some particular ABC format files. The file name is the user's choice.

    Except when testing a single number in command mode, using -q"k*b^n+c",
    the successful results (prime or PRP) are always registered in an
    OUTPUT FILE ; again, the file name is the user's choice, but the
    format of this file is the same as the input file's one.
    So, it can be used as input for another program.

    In all cases, the results are all registered in a RESULT FILE,
    (not to be confused with the output file). Its default name is 
    "lresults.txt". When using the command flag -a<nnnn>, its name is 
    "lresu<nnnn>.txt"

    Default, working and user's OPTIONS are registered in a .INI FILE
    Its name is "llr.ini" (default), or "llr<nnnn>.ini" if using -a<nnnn>.
    These data are needed when stopping and resuming a job.

    To avoid the need to restart a test from the beginning, after a
    crash or an user's stop, intermediate data are registered at 
    regular intervals, or when stopping, in a temporary file,
    which is automatically named by the program, and removed when
    the test is completed. Its name is a letter followed by 7 digits
    computed after the candidate value currently tested.

5) Main user options (not set by default) :

    Verbose=1 : Get more details in the results file (default : 1 line/result).
    StopOnSuccess=1 : Stop the job when a prime or PRP is found.
    BeepOnSuccess=1 : Make noise at a positive result,
      if both Stop and Beep are set, make noise until stopped by the user!
    StopOnPrimedK=<number> : after <number> sucesses with this k value,
      skip further pairs having the same k value (usually, <number> = 1).
    StopOnPrimedN=<number> : Same thing, involving the value of n.
    StopOnPrimedB=<number> : Same thing, involving the base value.
    Verify=1  : Suppress prefactoring or previous PRP test.
    NoPrefactoring=1 : Suppress prefactoring (Gaussian Mersenne or Wagstaff).
    ErrorCheck=1 : Check errors on each iteration (it's time consuming!).
    Testdiff=1 : Check sum inputs != sum outputs (only for real FFT's, c<0).
    FacTo=<bits> : Used to launch a prefactoring only job (Wagstaff or
    Gaussian-Mersenne norms candidates only).

6) Options used to change a default value :

    OutputIterations=<number> : Nb. of iters between outputs (def. 10,000).
    DiskWriteTime=<number> : Time elapsed between disk savings (def. 30mn.).
    FBase=<number> : The base for the Fermat PRP test (default is 3).
    MaxRestarts=<number> : Max. restarts of an N+1 or N-1 test (default 10).
    - There are several other values you have almost no reason to change...

7) More special options :

    ForcePRP=1 : Do only a PRP test, even if a deterministic one is possible.
    LucasPRPtest=1 : For a PRP test, use only the Lucas+Frobenius algorithm.
    FermatPRPtest=1 : For a PRP test, use only the Fermat SPRP algorithm.
    (The default is Fermat SPRP, followed by Frobenius on positive results)
    TestGM=1 : Register the Gaussian-Mersenne norm if it is prime (default).
    TestGQ=1 : Register the associated (2^p+-2^((p+1)/2)+1)/5 if it is PRP.
    VrbaReixTest=1 : Test a Wagstaff number using the Vrba-Reix algorithm.
    (the default is to do a strong Fermat PRP test)
    DualTest=1 : Test again a Wagstaff PRP with the alternate algorithm.

8) Input file formats :

    - In previous versions, the format of input data was known by LLR
    according to the very first (header) line. In this new version, you
    may now have MULTIPLE DATA FORMATS in the same input file, because
    data format descriptors can be inserted anywhere in the file. In
    return of that, if an invalid descriptor is found, input lines are
    flushed until finding the next valid one...And, indeed, the first
    input line must be a valid descriptor!

    - The NEWPGEN DESCRIPTOR has five fields separated by colons :

      <sieved to>:<letter code>:<chain length>:<base>:<mask>

      for example 1:P:1:2:1 for a Proth test, 1:M:1:2:2 for a Riesel one.
    - The second and last field describe the expression to be tested.
      (yes, it is redundant, the mask overrides the letter and should
      be preferred)
    - <base> is the value of b in k*b^n+c
    - <chain length> should be 1, excepted for Cunnigham chains.
    - <sieved to> integer is ignored by LLR (only copied in output file).
    - All NewPgen file formats are accepted, except the Primorial ones.
    - For more details, consult in-line help of NeWPgen or Appendix below.

    - Moreover, LLR accepts these ABC FORMAT DESCRIPTORS :

      - Two numbers per data line formats :

	- Fixed k and c : ABC%d*$a^$b+%d or ABC%d*$a^$b-%d
	- Fixed b and c : ABC$a*%d^$b+%d or ABC$a*%d^$b-%d
	- Fixed n and c : ABC$a*$b^%d+%d or ABC$a*$b^%d-%d

      - Three numbers per data line formats :

	- Fixed k : ABC%d*$a^$b$c
	- Fixed b : ABC$a*%d^$b$c
	- Fixed n : ABC$a*$b^%d$c
	- Fixed c : ABC$a*$b^$c+%d or ABC$a*$b^$c-%d

      - General k*b^n+c format (four numbers per data line) :

	- ABC$a*$b^$c$d

      - Some special ABC formats :

	- ABC4^$a+1 : Gaussian-Mersenne norm candidates
	- ABC(2^$a+1)/3 : Wagstaff PRP candidates
	- ABC(10^$a-1)/9 : Repunits PRP candidates
	- ABC($a^$b-1)/($a-1) : Generalized Repunits PRP candidates
	- ABC$a*$b^$a$c : (Generalized) Cullen/Woodall candidates
	- ABC(2^$a$b)^2-2 : near-square (ex-Carol/Kynea) candidates

9) Basic Algorithms :

    - The base two numbers (with k<2^n) are the fastest to test :
     k*2^n+1 numbers are tested using the Proth algorithm.
     k*2^n-1 numbers are tested using the Lucas-Lehmer-Riesel algorithm.

    - Non base two numbers (with k<b^n) :
     k*b^n+1 numbers are tested using the N-1 Pocklington algorithm.
     k*b^n-1 numbers are tested using the N+1 Morrison algorithm.

    - K*b^n+c numbers with |c| <> 1 or k > b^n can only be PRP tested.

10) Special algorithms :

    - GAUSSIAN-MERSENNE NORMS are tested for primality by using the
    factorization  : 4^p+1 = (2^p+2^((p+1)/2)+1)(2^p-2^((p+1)/2)+1)
    if p is prime, one factor may be prime, and then, is the norm of
    a prime complex Gauss integer of the form (1+or-i)^p-1. 5 divides
    always the second factor, but the quotient by 5 may be also prime.
    The algorithm is fast, because squarings are done modulo 4^p+1.
    Moreover, the primality test of the GM factor, and the PRP test of
    the quotient GQ are done in the same loop.

    - The test of WAGSTAFF NUMBERS W = (2^p+1)/3 can be done with two
    algorithms : a strong Fermat PRP, and/or the Vrba-Reix algorithm.
    Both are fast, because the squarings are done modulo 2^p+1.
    The Vrba-Reix algorithm, proposed by Tony Reix and Anton Vrba, is
    similar to Lucas-Lehmer for Mersenne or Fermat numbers,
    but with 3/2 modulo W as a seed.
    Both tests are known as a necessary condition for primality.
    For now, the Vrba-Reix test is only a PRP test. It has not yet
    been proved to be a primality test for Wagstaff numbers.

11) Error checking and recovery :

    - Error checking is done on the first and last 50 iterations, before 
    writing an intermediate file (either user-requested stop or a
    30 minute interval expired), and every 128th iteration.
    - After an excessive (> 0.40) and reproducible round off error,
    the iteration is redone using a slower and more reliable method.
    - If this error was not reproducible, or if the iteration fails again,
    the test is restarted from the beginning, using the next larger
    FFT length...

Appendix :

    - Complements about NewPGen descriptors :

    - The letter is a one character code as follows :

      P : k.b^n+1 (Plus)
      M : k.b^n-1 (Minus)
      T: k.b^n+-1 (Twin)
      S: k.b^n-1; k.b^(n+1)-1 (Sophie Germain (CC 1st kind len 2))
      C: k.b^n+1; k.b^(n+1)+1 (CC 2nd kind len 2)
      B: k.b^n+-1; k.b^(n+1)+-1 (BiTwin)
      J: k.b^n+-1; k.b^(n+1)-1 (Twin/SG)
      K: k.b^n+-1; k.b^(n+1)+1 (Twin/CC)
      Y : k.b^n+1 + others (Lucky Plus)
      Z : k.b^n-1 + others (Lucky Minus)
      1: CC 1st kind chain
      2: CC 2nd kind chain
      3: BiTwin chain

    - NEWPGEN output files use the mask as defined below :

      0x01 :   k.b^n+1
      0x02 :   k.b^n-1
      0x04 :   k.b^(n+1)+1
      0x08 :   k.b^(n+1)-1
      0x10 :   k.b^(n+2)+1
      0x20 :   k.b^(n+2)-1
      0x40 :   PRIMORIAL - can't handle this
      0x80 :   k.b^n+5
      0x200 :  2^n+2k-1
      0x400 :  MODE_NOTGENERALISED, so :
      0x404 :  2k.b^n+1
      0x408 :  2k.b^n-1
      0x410 :  4k.b^n+1
      0x420 :  4k.b^n-1
      0x800 :  k.b^n+7
      0x1000 : 2k.b^n+3
      0x8000 : MODE_DUAL, so :
      0x8001 : b^n+k
      0x8002 : b^n-k

