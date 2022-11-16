
		Welcome to LLR program Version 4.0.3 !

0) What is new in this version :
  I added two new ABC formats, principally to help PRP searchers.
  - k*b^n+c format with k, b, c fixed, for example : ABC 22*17^n+13
  - (k*b^n+c)/d format with k, b, c, d fixed, for example : ABC (1*16^n+619)/5
  - In Proth or LLR tests, even values of k yield a false result...
  These bugs are now fixed.
  The Affinity managing was not really implemented on LLR...
  This issue is now fixed on Linux and WIN32 versions.
  The option -oAffinity=2 allows the progam to run on logical core 2.
  You may choose a list of cores by setting -oAffinity="2,3,5" for example.
  In previous Version 4.0.0, one call to free() function was missing in
  Gerbicz error checking code ; this caused an important memory leak...
  This issue is now fixed here!
  No much new feature, but some improvements related to reliability and speed.
  - I implemented the patch sent to me Serge Batalov, and the new ABC Format :
  ABC DivPhi($a*$b^$c+1) may be used to do these tests.
  - By default, all tests on base two numbers use Gerbicz error checking.
  This is the case for PRP Fermat and SPRP tests as in Prime95 or Mprime,
  but also for the deterministic prime tests of Proth numbers.
  LLR tests on Riesel numbers are only done after a positive Fermat PRP result.
  Also, if b==2, k==+1 and abs(c)==1, a random shift on the PRP base is done.
  It is especially interesting for the prime test of Gaussian Mersenne norms.
  - The hwloc library is now linked with the program and may be used.
  - -oCpuSupportsAVX512F=0 option is supported, as required by some users.
  - I compiled GNU gmp6.1.0 on Windows 32 bits VC6.0, so, this library is now
  linked on this Windows platform, which allows to do APRCL tests without the
  need of an external binary... 
  

1) Main features :

  LLR is a primality proving program for numbers of the form N = k*b^n +/- 1,
  (k < b^n), or numbers which can be rewritten in this form, like 
  Gaussian-Mersenne norms or b^n-b^m +/- 1 with n>m (new feature).
  The identity Phi(3,-X) = X^2-X+1 is now used with X=b^n to search for
  Generalized Unique Primes.

  It can also do strong and efficient PRP tests on more general k*b^n+c forms,
  on Wagstaff numbers (2^p+1)/3, repunits (10^p-1)/9 and generalized repunits
  (b^p-1)/(b-1), b!=2.

  In this new version, the base b can now be an arbitrary large integer 
  (however restricted to 64 bits for fixed b formats). This is especially
  useful while testing Generalized Fermat numbers.

  This new version can be built on 64bit platforms, but then, the prefactoring
  code is no more available. This affects only the Gaussian-Mersenne norm and
  Wagstaff tests, for which the prefactoring must be done using a 32bit program.

  This version uses the last release version (30.6) of George Woltman's Gwnum
  library, to do fast multiplications and squarings of large integers modulo N.
  A multithreading AVX bug due to misuse of the POSTFFT feature by LLR is
  fixed (28/01/18).
  A bug that affected tests using multithreading and FMA3 has been fixed.
  A bug has been fixed in 28.6 version, and a related new issue in 28.7 : 
  They affected only tests on a CPU which supports AVX code, and, indeed, 
  if this feature was activated. This bug existed in gwnum versions 27.1
  and higher, and so, in LLR from version 3.8.9. The AVX bug seems to be really
  fixed in this version;
  A bug that affected tests using AVX512 code on Windows64 machines is
  fixed.

  The main advantage of this gwnum version is better performances on 64bit
  machines ; also, several internal bugs in V26 have been corrected.

  Since version 25.11, gwnum library is no longer restricted to base two for
  efficient computing modulo k*b^n+c numbers (but SSE2 is required if b != 2), 
  and LLR greatly takes advantage of this improvement.

  LLR can run on all machines where gwnum code can run, so, on all Intel x86
  compatible machines.

  The MULTITHREADING can now be used by setting -t<number> or 
  -oThreadsPerTest=<number> in the command line. Many thanks to
  Serge Batalov who showed me how simple it was to implement this feature!

  11/02/16 : Three new features has been added in the small primes ranges
  (< 1000 digits), all of them using the GNU-MP library and the APR-CL codes :

  - The search for Wieferich primes in the range and with the base stated
  by the user.
  - The test of Wieferich prime candidates with the base stated by user.
  - The APR-CL primality test of general digit strings.
  Please, see below for the ABC formats required for these features.
  
  - The option -oNoSaveFile=1 has been added at the request of an user.
  Indeed, if it is set, any test is restarted at beginning if stopped...

  - 05/05/15 : The Roundoff error recovery code has been completely rewritten.
  The logic is now very similar to Prime95 / Mprime's one, to avoid intensive
  usage of slow careful iterations. Continuing the test using the next FFT
  length is forced if the error is not reproducible, if it occured in a careful
  iteration, or if too much errors were encountered (MAX_ERROR_COUNT=5).

  - 11/02/16 : To avoid endless retries, AbortOnRoundoff option is forced when
  testing with the next FFT length has been done more than MaxFFTinc times
  (5 by default). Then, the default is continuing the test with the next term
  in the input file. But the user can override this behavior by setting
  -oStopOnAbort=1, and then, manually deciding to continue.

  - To improve reliability, error checking may now be forced, if the program
  is working near the current FFT limit. This feature may be adjusted by using
  the option -oPercentFFTLimit=dd.d, the default value being 0.5 ; note that
  setting echk to one is generally not much time consuming : typically 5% more.
  Also, this feature may be wiped out by setting PercentFFTLimit to 0.0!

  - For those wo do not like to force error checking, I implemented a new
  option : -oNextFFTifNearLimit=1 (default is zero). If activated, and if the
  default FFT length at setting is too near the limit, then, FFT_Increment is
  incremented by one, a message is displayed and the test is immediatly 
  restarted ; indeed, in this case, echk can no more be forced...

  - In all versions of LLR, a simple trial division test was done for candidates
  not larger than 32 bits ; now, an APR-CL test as been added as a new feature
  in all primality or PRP tests, for all candidates not larger than 100 decimal
  digits. On Windows platforms, this code is implemented in a binary application
  called as a child process.

  - The new option -oBPSW=1 allows to replace the standard Lucas PRP test by
  the Baillie-PSW one (see http://www.trnicely.net/misc/bpsw.html for details)
  the Frobenius test then follows as usual.

  - The two (mutually exclusive) options -oRising_ns=1 and -oRising_ks=1 have
  been added in this version. They allow to process an input file
  while it is sieved by another process of the same machine (LLR closes the file
  during the test of the current candidate, to allow its updating by the sieve).

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
    use the scrolling menues. Then, the options, results and .ini 
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

    - To test a SINGLE k*b^n+c or b^n-b^m+c number, type :

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

    ThreadsPerTest=<number> or -t<number> : (default is 1).
    OutputIterations=<number> : Nb. of iters between outputs (def. 10,000).
    DiskWriteTime=<number> : Time elapsed between disk savings (def. 30mn.).
    FBase=<number> : The base for the Fermat PRP test (default is 3).
    PBase=<number> : The starting P value for a Lucas test (default is 3).
    MaxRestarts=<number> : Max. restarts of an N+1 or N-1 test (default 10).
    MaxN=<number> : Stop the batch when this exponent value is reached...
    NoLresultFile=1 : Suppress the recording of the result file...
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

	- ABC$a^$b+1 : Generalized Fermat candidates
	- ABC4^$a+1 : Gaussian-Mersenne norm candidates
	- ABC$a^$b-$a^$c-1 : a^b-a^c-1 candidates
	- ABC$a^$b-$a^$c+1 : a^b-a^c+1 candidates
	- ABC(2^$a+1)/3 : Wagstaff PRP candidates
	- ABC(10^$a-1)/9 : Repunits PRP candidates
	- ABC($a^$b-1)/($a-1) : Generalized Repunits PRP candidates
	- ABC$a*$b^$a$c : (Generalized) Cullen/Woodall candidates
	- ABC(2^$a$b)^2-2 : near-square (ex-Carol/Kynea) candidates
	- ABC$a$b$c : Used to launch a Wieferich prime search, 
	  the range being $b to $b and the base $c (new feature!)
	- ABC$a$b : Used to test a Wieferich prime candidate $a, base $b
	- ABC$a : General APR-CL primality test of number $a

9) Basic Algorithms :

    - The base two numbers (with k<2^n) are the fastest to test :
     k*2^n+1 numbers are tested using the Proth algorithm.
     k*2^n-1 numbers are tested using the Lucas-Lehmer-Riesel algorithm.

    - Non base two numbers (with k<b^n) :
     k*b^n+1 numbers are tested using the N-1 Pocklington algorithm.
     k*b^n-1 numbers are tested using the N+1 Morrison algorithm.

    - K*b^n+c numbers with |c| <> 1 or k > b^n can only be PRP tested.
     If the number is found PRP, the % of factorization is then shown,
     but note that it is relevant only if c == +1 or -1...

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
    the test is restarted from the last save file, using the next larger
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

