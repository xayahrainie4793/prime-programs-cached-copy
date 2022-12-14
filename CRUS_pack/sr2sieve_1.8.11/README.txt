Latest binaries can be found at:

  http://www.geocities.com/g_w_reynolds/sr2sieve/

Source is at:

  http://www.geocities.com/g_w_reynolds/sr5sieve/

To report bugs, email me (Geoff): g_w_reynolds at yahoo.co.nz.


SR2SIEVE
========

Sr2sieve is Srsieve specialised for sequences k*b^n+/-1 and b^n+/-k.

Command line options:

 -s --sierpinski       Sieve sequences k*2^n+1 in dat format file 'SoB.dat'
 -r --riesel           Sieve sequences k*2^n-1 in dat format file 'riesel.dat'
 -i --input FILE       Read sieve from FILE instead of `sr2data.txt'.
                        Instead of `SoB.dat',`riesel.dat' when used with -s,-r.
 -p --pmin P0          Sieve for factors p in the range P0 <= p <= P1 instead
 -P --pmax P1           of reading ranges from work file `sr2work.txt'.
 -f --factors FILE     Append found factors to FILE instead of factors.txt.
 -u --uid STRING       Append -STRING to base of per-process file names.
 -c --cache            Save Legendre symbol tables to `sr2cache.bin'.
 -C --cache-file FILE  Load (or save) Legendre symbol tables from (or to) FILE.
 -l --L1-cache SIZE    Assume L1 data cache is SIZE Kb.
 -L --L2-cache SIZE    Assume L2 cache is SIZE Kb.
 -B --baby METHOD      Use METHOD for baby step mulmods.
 -G --giant METHOD     Use METHOD for giant step mulmods.
    --ladder METHOD    Use METHOD for ladder mulmods.
 -H --hashtable SIZE   Force use of a SIZE Kb hashtable.
 -Q --subseq Q         Force sieving k*b^n+c as subsequences (k*b^d)*(b^Q)^m+c.
    --scale-giant X    Scale the number of giant steps by X (default 1.0).
    --min-giant NUM    Always perform at least NUM giant steps (default 1).
 -D --duplicates FILE  Append duplicate factors to FILE.
 -S --save TIME        Write checkpoint every TIME seconds. (default 300).
 -z --lower-priority   Run at low priority. (-zz lower).
 -Z --raise-priority   Run at high priority. (-ZZ higher).
 -A --affinity N       Set affinity to CPU N.
 -j --sobistrator      Run in Sobistrator compatibility mode.
 -X --skip-cubic       Skip cubic and higher power residue tests.
 -x --no-lookup        Don't pre-compute Legendre symbol lookup tables.
 -d --dual             Sieve b^n+/-k instead of k*b^n+/-1.
 -q --quiet            Don't print found factors.
 -v --verbose          Print some extra messages. -vv prints more.
 -h --help             Print this help.

Some of the following additional options may also be available:
    --amd              Use CMOV optimisations.
    --intel            Don't use CMOV optimisations.
    --sse2             Use SSE2 vector optimisations.
    --no-sse2          Don't use SSE2 vector optimisations.
 -t --threads NUM      Start NUM child threads. (Default 0).

If -p and -P are not provided, then sr2sieve will read the sieve range from
a file called sr2work.txt in the current directory. Each line of this file
should contain a pair of integers A,B seperated by a comma, space, or dash.
A and B are in billions, the actual sieving range is A*10^9 to B*10^9.
E.g. invoking `sr2sieve' with the line `1000,1001' in sr2work.txt is
equivalent to invoking `sr2sieve -p 1000e9 -P 1001e9'.

For sequences k*b^n+/-1 and prime factor p, the following limits apply:

  1 < k < 2^32.
  1 < b < 2^32.
  core(k)*core(b) < 2^31. (core(k) is the squarefree part of k).
  0 < n < 2^32.
  k,b < p < 2^52. (or 2^62 for x86/x86-64 builds, or 2^63 for ppc64).

For sequences b^n+/-k, the following additional limit applies:

  1 < k < 2^31.

Only factors larger than k and b will be found, so the sieve must be started
and sieved until p >= max(k,b) with another program (NewPGen, Srsieve).

If core(k)*core(b) is large, then a lot of initialization time and memory
may be required.  To avoid this the -x switch can be used, at some cost to
performance. Alternatively the -c or -C switches can be used to write
initialization data to file, which will be used to speed up initialization
the next time the sieve is started.

If neither -s nor -r switches are used then the input file must be in ABCD
format. The srfile utility can be used to convert files from NewPGen format
to ABCD format and vice versa, and can also be used to remove found factors
from the ABCD file (sr2sieve doesn't update the input file itself). Srfile
is included in srsieve at http://www.geocities.com/g_w_reynolds/srsieve/

If sr2sieve is interrupted (e.g. pressing ctrl-c) or terminated nicely
(e.g. by running: kill `pidof sr2sieve`) then any output files will be
updated before it finishes. Beware that this doesn't happen in Windows when
closing the sr2sieve window by clicking the close button, so press ctrl-c.

When comparing the output with Proth sieve, be aware that sr2sieve only
reports the first factor it finds for a number and all subsequent factors
for the same number are considered duplicates. Proth sieve however doesn't
count these as duplicates, and so fact.txt may contain some extra factors.
Use 'sort -k3,3 fact.txt | uniq -f2 -d' to find them.

If no command line arguments are given but `sr2sieve-command-line.txt'
exists in the current directory, then the command line will be as if the
first line of this file had been used to invoke sr2sieve. This may be useful
on some GUI machines where the command shell and batch files have been
disabled for security reasons.

If the `-j --sobistrator' switch is given then sr2sieve will behave in a
similar way to JJsieve or proth_sieve, for compatibility with Sobistrator.
Instead of reading ranges from `sr2work.txt' and writing checkpoints to
`checkpoint.txt', checkpoints will be written to `SoBStatus.dat' (or
`RieselStatus.dat' if the -r switch is used) and subsequent ranges read from
`nextrange.txt'. In addition factors will be written to `fact.txt' and
duplicates to `factexcl.txt' (these file names can be overridden with the
-f and -D switches.)

To compile sr2sieve from the sr5sieve source, just change the definition of
BASE in sr5sieve.h from 5 to 0, run make as usual (see INSTALL), and rename
the executable from sr5sieve to sr2sieve.


ABCD file
=========

While sr2sieve can read SoB.dat and riesel.dat files when sieving sequences
of the form k*2^n+/-1, for more general use it is necessary to use ABCD
format.

The ABCD format expected by sr2sieve (which is only a subset of the possible
ABCD format files) has this form:

  ABCD <K>*<B>^$a<+/-C> [<N0>] // Sieved to <P0>
  <D1>
  <D2>
  <D3>
  ...


<N0> is the lowest N in the sequence K*B^N+C, <D1> is the offset to the
next N, etc.  <P0> is optional.


Example: An ABCD file for the terms 254*5^76-1, 254*5^162-1, 254*5^186-1:

  ABCD 254*5^$a-1 [76] // Sieved to 6703249519 with srsieve
  86
  24


Example: An ABCD file for the terms 2^1044-31859, 2^1092-31859, 2^1140-31859
and 2^1075+10223, 2^1159+10223, 2^1531+10223:

  ABCD 1*2^$a-31859 [1044] // Sieved to 1000000 with srsieve
  48
  48
  ABCD 1*2^$a+10223 [1075]
  84
  372


When using srsieve to start the sieving, use the -a switch to save the sieve
in ABCD format for later use with sr2sieve. E.g.

  srsieve -a -n 1e6 -N 2e6 -P 1e9 "2^n-31859" "2^n+10223"


To convert a bunch of NewPGen sieve files into ABCD format, use the srfile
utility included with srsieve. E.g.

  srfile -a *.npg


To convert an ABCD file called myfile.abcd back into a bunch of NewPGen
files for testing with PRP or LLR:

  srfile -g myfile.abcd
