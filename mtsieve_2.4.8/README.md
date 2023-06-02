The Multi-threaded Sieve Framework

The Multi-threaded Sieve Framework (called mtsieve) is a framework that I have been using for many of the sieving programs that I have written over the years.  Although those programs reference the framework, it was never promoted for other people to take advantage of.  The main idea of this framework is to help optimize the search for factors of a collection of numbers that follow a well-defined pattern.  Some examples are:

factorials - where each successive term can only be calculated by the previous term
number in the form k*b^n+c for a fixed k, b, and c and a range of n - a Discrete Log can be used to find an n such that k*b^n+1 is divisible by a given p
numbers in the form k*b^n+c for a fixed b, n, and c and a range of k - a simple function can compute k such that k*b^n+1 is divisible by a given p

This framework can also be used to do other functions with primes, such as determining if each prime is a Wieferich Prime.

You can d/l mtsieve sources and all sieving programs using that framework from here.
Open and Portable
I cannot tell you how many times I wanted to search for primes on my Mac only to find that before I could do that I would have to sieve on a Windows OS.  Two of the most notorious examples are newpgen and fermfact.  Although the newpgen source is available, that source is a complete mess.  I understand that it was great during its heyday on 32-bit CPUs, but there was nothing in it to take advantage of 64-bit CPUs.  On the other hand fermfact source is not available despite my requests to the original author.  This also has 32-bit limitations and the command line arguments are confusing at best.  And even though I am the author of MultiSieve, it was 32-bit and Windows GUI only.  gfndsieve built upon the mtsieve framework replaces fermfact.  Other sieves built upon the framework only cover some of the functionality of newpgen and MultiSieve.  If I get requests, I will provide work on extending mtsieve to replace other sieves from those programs.

Since mtsieve is also 64-bit the memory limitations of the old 32-bit programs are gone.

With mtsieve, we now have a framework that can built and run on Windows, Linux, and Mac OS X.  For Windows I build with msys2so that I can use many of the standard libraries that many "old time programmers" are familiar with rather than doing something the "Visual Studio" way.  For Mac and Linux, the software should build out of the box with the given makefile.

What about ARM?
Most of the programs in the framework can be built on ARM.  Exceptions are afsieve, pixsieve, and xyyxsieve which rely on x86 ASM.  This might change in the future, but nobody is asking, therefore it is not high priority.

What about the GPU?
Many of the programs in the framework support GPUs.  Since I use the OpenCL SDK, they can run on both AMD and NVIDIA GPUs.

What about Apple Metal?
I am working on it, but it is on the back burner right now.  I'm close, but Apple's framework for using Metal can be challenging.
What about srsieve/srfile/sr1sieve/sr2sieve?
These are all integrated into srsieve2/srsieve2cl.  srsieve2 is faster than srsieve.  sr1sieve is faster than srsieve2.  sr2sieve is faster than srsieve2, when using Legendre tables.  srsieve2cl is faster than all of them.  In other words, if you have a GPU, use srsieve2cl.

Why isn't srsieve2 always faster than sr1sieve/sr2sieve?
This is due to hand-tuned x86 ASM code in sr1sieve/sr2sieve.  Since I care about portability and because I use srsieve2cl exclusively I did not port those routines.  I also support srsieve, sr1sieve, and sr2sieve and will continue to support them.

Why isn't the mtsieve program faster than another program that does the same thing?
The likely have hand-tuned x86 ASM or are using other algorithms to get that speed bump.  In some cases I do not have access to their source code.  In some cases the code is only intelligible to the original developer.  In some cases it is far too specialized.  I do what I can to keep my code fairly clean so that people with average math and coding skills can understand it.

What about sieving for <xxx>?
It depends upon what <xxx> is.  If it is something I wrote, I'll get around to it.  If you are a developer, take a look at the bottom of this page.  If you are still stuck I'll help you use the framework to create your application.

Can I contribute?
Absolutely.  Others have contributed performance improvements either to the framework or to individual sieves.  I am easy to get hold of via e-mail via rogue _at__ wi.rr.com.

Can you add features to <xxx>?
Maybe.  It depends upon how well they fit into that program.  If it requires a new algorithm in the worker class, then probably not.

Can you write a program to sieve for me?
Maybe.  It would have to be of interest to me.  I am more likely to help you optimize the worker than do the rest.  I strongly encourage you to learn C++.  You will thank me later.

What about 32-bit?
All of the current programs using the framework create 64-bit builds.  In theory most will compile as 32-bit applications, but I haven't tried.  Nobody has asked me for a 32-bit build in years, so I have no intention of supporting one.
What the Framework Provides
The main goal of the framework is to abstract Windows, Linux, and OS X functionality from the developer using the framework.  For example, if you are a Windows software developer, you don't need to know how to create a mutex on Linux.  Likewise if you are a Mac software developer, you don't need to know how Windows threading works.  As a developer you only need to know the interface (the .h file), and not the execution (the .cpp file).  The sources are there so that you can see how it is done on other platforms, but you won't need to understand the details, you just have to call the correct functions.

The framework also tries to abstract the GPU so you don't need to know the gory details of OpenCL or Metal.

The Source
Classes in core
The classes in this directory are the ones shared by all of the program using the framework.  These classes include:

main.cpp - this class contains the main() function and the handler for signals
Clock.cpp - this class contains static functions for getting the current clock and processor times, in microseconds
Parser.cpp - this class contains functions for parsing command line options.  Many command line options support scientific notation, i.e. 1e3 --> 1000, and 1e6 --> 1000000
App.cpp - this is an abstract class that implements the core features shared by all applications built using mtsieve.  There is only one instance of this class at any time when your program is running.  This thread will sleep whenever there are no WorkerThreads available.FactorApp.cpp - this class is used for capturing runtime information so that factor rates and be shown while programs are running
FactorApp.cpp - this is an abstract class that implements features shared by applications that use mtsieve for factoring.  It extends App.cpp.
AlgebraicFactorApp.cpp - this is an abstract class that implements features shared by applications that use mtsieve for factoring but which also have algebraic factorizations.  It extends FactorApp.cpp.
Worker.cpp - this is an abstract class that implements the core functions for the process of each list of primes.  There will be one instance of this class for each worker thread when your program is running.
HashTable.cpp - this class provides a hash table, which is used programs that need one
SharedMemoryItem.cpp - this class provides a mutex for variables that can be read and modified by multiple threads.
GpuDevice.cpp - this abstract class manages GPU devices.
GpuKernel.cpp - this abstract class manages kernels that run in the GPU.
MpArith.h and MpArithVec.h - these contain the inline mulmod/powmod logic used by most of the sieves.
Please refer to App.h and Worker.h to see which method are abstract.  These are the methods that must be coded in any concrete class that extends App or Worker.

Classes in gpu_opencl
The classes in this directory are the ones shared by all of the program using the framework that rely on OpenCL kernels.  These classes include:

OpenCLDevice.cpp - this class contains functions for identifying the available GPU devices
OpenCLKernel.cpp - this class contains functions for each kernel to be run in the GPU
OpenCLKernelArgument.cpp - this class contains functions to define arguments for each Kernel as well as allocating GPU memory for those arguments
OpenCLErrorChecker.cpp - this class contains static functions to check the status of each call to an OpenCL function
Classes in sieve
This is the sieveing source from primesieve, a library whose sole purpose is to generate a list of primes as quickly as possible.  Per the primsieve license, no changes should be made to this code.

Runtime Options
mtsieve has a few command line options that are available to all programs using the framework.  These are:

-h - to print help for the command line options.  For GPU enabled sievers this will list the available platforms and devices.
-p - the minimum prime returned by the sieve
-P - the maximum prime returned by the sieve.  Most sieves have a limit of 2^62.  Some have a limit of 2^52.  The limit is shown at runtime.
-w - to specify how many primes each thread should work on at a time before asking for more primes.  The default is 1e6, but each program in the framework can change the default.  The application will adjust this value to create chunks of work that take between 1 and 5 seconds per chunk.
-W - to specify the number of CPU worker threads.  The default is 1.  It cannot be set to 0 even if using GPU workers.
-A - to apply factors or to reformat a candidate file without sieving.
Programs with GPU support have some additional options available to them:

-d - device to use on the platform
-D - GPU platform to use
-g - the number of blocks of primes per worker thread.  This number of primes per block dependent upon the GPU and the kernel.  The number of primes per worker will be output when the program starts.
-G - the number of GPU worker threads.  The default is 0.
Programs extending FactorApp or AlgebraicFactorApp (instead of App) have some additional options available to them:

-i - the name of an input file with input terms to factor.  When resuming sieving, the candidate list will be built from this file.  It will override any other options used when starting a new sieve.
-I - the name of an input file with factors to apply to terms. 
-o - the name of an output file to write terms to.  When the program is terminated, any numbers that have not been factored are written to this file.
-O - the name of an input file with any new factors that have been found.  This file can be used as input to pfgw to verify factors found by individual applications.
Each program using this framework, whether it extends App or FactorApp, can add its own command line options, but it cannot replace any that are built into the framework.

Software Using the mtsieve Framework
Here is a list of programs I've written the use this framework.  All are bundled with mtsieve and are available from the download link.

afsieve/afsievecl
This program searches for factors of Alternating Factorials.  This program supports these additional parameters:

-n - the minimum n
-N - the maximum n
-S - the amount n is iterated by per GPU kernel execution
-M - the maximum number of factors allowed per GPU kernel execution
cksieve/cksievecl
This program searches for factors of Carol / Kynea numbers. These numbers a form of Near Square numbers with the form (b^n-1)^2-2 and (b^n+1)^2-2.  This program supports these additional parameters:

-b - the base to search
-n - the minimum n
-N - the maximum n
-M - the maximum number of factors allowed per GPU kernel execution
dmdsieve
This program searches for factors of number of the form 2*k*(2^p-1)+1.  Numbers that are not removed from the sieve are potential divisors of Double Mersenne numbers.  This program supports these additional parameters:

-k - the minimum k to search
-K - the maximum k to search
-n - the n to search
-x - test remaining terms for Double Mersenne divisiblity
-X - when using -x, the number of k to sieve at a time
fbncsieve
This program searches for factors of numbers in the form k*b^n+1 and k*b^n-1.  It is specifically designed to replace similar functionality in newpgen.  This program supports these additional parameters:
-k - the minimum k to search
-K - the maximum k to search
-s - the sequence to find factors of.  The sequence must be of the form k*b^n+c where b, n and c take decimal values.
-f - the format of the output file (A = ABC, D = ABCD, N = NEWPGEN)
-r - remove k where k % base = 0
fkbnsieve
This program searches for factors of the form k*b^n+c for fixed k, b, and n and variable c.   This program supports these additional parameters:
-c - the minimum c to search
-C - the maximum c to search
-s - the sequence to find factors of.  The sequence must be of the form k*b^n+c where k, b and n take decimal values.
gcwsieve/gcwsievecl
This program searches for factors of Cullen and Woodall numbers.  This program supports these additional parameters:

-b - the base to search
-n - the minimum n to search
-N - the maximum n to search
-a - use AVX routines (only on x86)
-s - sign to sieve for (+ = Cullen, - = Woodall, b = both)
-f - the format of the output file (A = ABC, L = LLR)
-S - the number of steps iterated per GPU kernel execution
-M - the maximum number of factors allowed per GPU kernel execution
gfndsieve/gfndsievecl
This program searches for factors of numbers in the form k*2^n+1 for a range of k and a range of k.  It is specifically designed to replace fermfact.  The output from this sieve should be used with pfgw and the -gxo switch to find GFN divisors. This program supports these additional parameters:

-k - the minimum k to search
-K - the maximum k to search
-n - the minimum n to search
-N - the maximum n to search
-T - the number of n per output file
k1b2sieve
This program searches for factors of numbers of the form b^n+c for fixed b and variable c.  This program supports these additional parameters:

-c - the minimum c to search
-C - the maximum c to search
-n - the minimum n to search
-N - the maximum n to search
kbbsieve
This program searches for factors of numbers of the form k*b^b+1 or k*b^b-1 for fixed k and variable b.  This program supports these additional parameters:

-k - the k value to search
-b - the minimum b to search
-B - the maximum b to search
mfsieve/mfsievecl
This program searches for factors of MultiFactorials.  This program supports these additional parameters:

-n - the minimum n
-N - the maximum n
-m - multifactoral (ie x!m where m = 1 -> x!, m = 2 -> x!!, m = 3 -> x!!!, etc.)
-S - the number of steps iterated per GPU kernel execution
-M - the maximum number of factors allowed per GPU kernel execution
pixsieve/pixsievecl
This program searches for factors of numbers that are a substring of a long decimal string where each successive term adds on decimal digit to the end of the previous decimal term.  This program supports these additional parameters:

-l - the minimum length to search
-L - the maximum length to search
-s - the file which contains a decimal representation of a number (eg pi, e, the Champernowne constant, etc.)
-S - the starting point of the substring
-N - the number of steps iterated per GPU kernel execution
-M - the maximum number of factors allowed per GPU kernel execution
psieve/psievecl
This program searches for factors of primorials.  This program supports these additional parameters:

-n - the minimum n to search
-N - the maximum n to search
-S - the number of steps iterated per GPU kernel execution
-M - the maximum number of factors allowed per GPU kernel execution
sgsieve
A program to sieve for Sophie-Germain primes of the form k*b^n-1 with variable k, fixed b, and fixed n.

-k - the minimum k to search -K - the maximum k to search
-b - the b to search
-n - the n to search
-g - multiply the second term by b instead of 2 for the generalized form
-f - the format of the output file (D = ABCD, N = NEWPGEN)
smsieve/smsievecl
A program to sieve for Smarandache primes.

-n - the minimum n to search
-N - the maximum n to search
-S - the number of steps iterated per GPU kernel execution
-M - the maximum number of factors allowed per GPU kernel execution
srsieve2/srsieve2cl
A program to sieve Sierpinski/Riesel sequences of the form k*b^n+1 or k*b^n-1 for fixed b, variable n, and multiple k.  If you need to use -K > 1, then consider adjusting -b as it can reduce the number of kernels you need.  -U, -V, and -X and advanced options to play around with if you are trying to maximize the sieving rate.

-n - the minimum n to search
-N - the maximum n to search
-s - a sequence or a file of sequences
-f - the format of the output file (A = ABC, D = ABCD, B = BOINC, P = ABC with number_primes)
-l - bytes to use for the Legendre tables, only supported if abs(c) = 1 for all sequences
-L - input/output directory to hold the Legendre tables.  No files are created if -L is not specified.
-M - the maximum number of factors per 1e6 terms per GPU kernel execution
-K - the number of kernels when splitting large numbers of sequences for the GPU
-C - the number of chunks of primes per GPU worker
-R - remove specified sequence
-b = used when calculating the number of baby steps and giant steps.  As be increases so does the number of baby steps.
-U - multiplied by 2 to compute BASE_MULTIPLE
-V - multiplied by BASE_MULTIPLE to compute POWER_RESIDUE_LCM
-X - mulitplied by POWER_RESIDUE_LCM to compute LIMIT_BASE
twinsieve
This program searches for factors of twin numbers of the form k*b^n+1 and k*b^n-1.  This program supports these additional parameters:

-k - the minimum k to search
-K - the maximum k to search
-b - the base to search
-n - the n to search
-f - the format of the output file (A = ABC, D = ABCD, N = NEWPGEN)
-r - remove k where k % base = 0
-s  to sieve +1 and -1 sides independently
xyyxsieve/xyyxsievecl
This program searches for factors of x^y+y^x and x^y-y^x numbers.  This program supports these additional parameters:

-x - the minimum x
-X - the maximum x
-y - the minimum y
-Y - the maximum y
-D - disable AVX routines
-s - the sign (+, - or b)
-S - the number of steps iterated per GPU kernel execution
-M - the maximum number of factors allowed per GPU kernel execution

How To Write Your Own Sieve Using the Framework
The best way to create your own sieve is to start with an existing sieve that is close to what you need.  This will save you a lot of time.

1)  Copy the folder for a sieve similar to what you want.
2)  Rename the folder to something meaningful (no spaces).
3)  Inside that folder, rename the .h and .cpp files, but keep the App.h, App.cpp, Worker.h, Worker.cpp portions of the name.
4)  Using NotePad++ open those files.
5)  Use the "find and replace all open files" to rename the class names to the same as the file names.
6)  Update the makefile and create a new object list for your sieve.
7)  Update the makefile and add an entry for your program near the end.  This will tell make what objects are needed for the executable.
For example, let's say that you want a new sieve similar to cksieve.  Let's call it mcsieve, short for "my custom sieve".
1)  Copy the carol_kynea folder and rename as my_custom
2)  Rename CarolKyneaApp.h to MyCustomApp.h.
3)  Rename CarolKyneaApp.cpp to MyCustomApp.cpp.
4)  Rename CarolKyneaWorker.h to MyCustomWorker.h.
5)  Rename CarolKyneaWorker.cpp to MyCustomWorker.cpp.
6)  Edit those four files in NotePad++.
7)  Use "find and replace in all files" to  change "CarolKynea" to "MyCustom" and save.
8)  Add MC_OBJS to makefile setting to "my_custom/MyCustomApp.o my_custom/MyCustomWorker.o"
9)  Copy the entry for "cksieve" and name as "mcsieve".
10)  For mcsieve, change CK_OBJS to MC_OBJS
11)  Type "make mcsieve" from the command line and it should build without errors.

Now for the fun part, writing your custom code.
MyApp constructor
Call SetBanner() to set the banner, which is printed when the program runs.  Call SetLogFileName() to set the name of the log file used by this program.  Set the variables specific to the program that can.
MyApp::Help()
Call FactorApp::ParentHelp() first, then print help for each of the program specific command line options.
MyApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
Call FactorApp::ParentAddCommandLineOptions first, then append shortOpts as necessary. For each short option, add a colon after the option to signify that the option takes a parameter.  Before returning call AppendLongOpt to add long options that will be supported by the command line.
MyApp::ParseOption(int opt, char *arg, const char *source)
Call FactorApp::ParentParseOption() first and if it returns P_UNSUPPORTED, then parse the argument for the commnad line option.  For numeric options, use Parser::Parse passing the argument along with the min and min values for the option and the variable to hold the value that is parsed. Return the status for parsing the option.
MyApp::ValidateOptions(void)
Use this method to apply additional validations to the input numbers, such as verifying that one input is less than another input.  It can also set or adjust some command line options if they were not specified on the command line.  Call FatalError() if the value for a parameter is invalid.
This method is responsible for creating the bit map representing candidates used throughout the application.
Call FactorApp::ParentValidateOptions() before returning.
MyApp::CreateWorker(uint32_t id, bool gpuWorker)
This method will create an instance of the correct worker class and return it to the caller.  The gpuWorker flag indicates if this method should create a GPU worker instead of a CPU worker.  GPU worker can only be true if ib_SupportsGPU was set to true in the constructor.
MyApp::ProcessInputTermsFile(bool haveBitMap) 
This method is called twice.  The first time it is called to determine the min/max of the candidates in the source file.  The second time it will apply the candidates from the source file to the bit map created by ValidateOptions().
MyApp::WriteOutputTermsFile(uint64_t largestPrime)
This method will create a file for the program that does the PRP testing.  This method must lock ip_FactorAppLock while accessing the list of candidates then release ip_FactorAppLock before returning.  The largestPrime is typically written to the first line of the output file and will be used as the starting prime if sieving is resumed from the file.
MyApp::ApplyFactor(const char *term)
This method is called at start up to apply factors from an input file.  The input term represents a candidate number, i.e. the actual string representing the number tested by the PRP program.  Parse this string to find the candidate then remove that candidate from the from the bitmap.  Decrement il_TermCount for each candidate removed.  All factors should be verified before updating the bitmap.
MyApp::ReportFactor
Each sieving application will need a version of this method.  In other words each application will call this method with a different list of parameters.  This method will take the inputs and turn of the bit in the bitmap.  For each found factor, it must decrement il_TermCount and increment il_FactorCount.  This method must lock ip_FactorAppLock while accessing the list of candidates then release ip_FactorAppLock before returning.
This method should return a boolean indicating if the term was removed from the bitmap.  All factors should be verified before updating the bitmap.
MyWorker constructor
Allocate memory or resourced needed by this worker.  This can also be used to set instance variables from MyApp, such as the min/max for the candidates being tested.
MyWorker::TestPrimeChunk(uint64_t &largestPrimeTested, uint64_t &primesTested)
This is where the heavy lifting goes.  The code in this method will iterate thru the it_Primes vector to find factors for the candidates being sieved.  It must set largestPrimeTested and primesTested before returning to the caller.
MyWorker::CleanUp(void)
This method will free memory and other resources allocated or created by the constructor.
