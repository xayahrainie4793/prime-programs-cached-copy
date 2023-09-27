#include <stdint.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#include <math.h>
#include <vector>

#ifdef WIN32
#include <windows.h>
#else
#include <inttypes.h>
#endif

#include "primegen.h"

#define NO_PHASE  999
#define MIN_N_FOR_RESIDUES 25000

typedef struct
{
   int32_t minN;
   int32_t maxN;
} phase_t;

int32_t     gi_base = 0;
int32_t     gi_c = 0;
int32_t     gi_recoveryPhase;
int32_t     gi_recoverySubPhase;
int32_t     gi_recoveryStep = 0;
uint32_t    gt_start_time = 0;
uint64_t    gl_minK = 0;
uint64_t    gl_maxK = 0;
uint64_t    gl_remainingK = 0;
uint32_t    gi_maxfbncprimes = 0;
uint32_t    gi_maxNfbncsieve = 0;
uint32_t    gi_maxKsrsieve2 = 0;
phase_t     gs_phases[105];
int32_t     gi_phaseCount = 0;
int32_t     gi_currentPhase = 0;
int32_t     gi_currentSubPhase = 0;
bool        gb_recovery = false;
bool        gb_quitting = false;

std::vector<bool> gb_KTerms;

#define BIT(k)       ((k) - gl_minK)

#ifdef WIN32
#define atoll  _atoi64

#define PRIu64 "I64u"
#define PRId64 "I64d"
#define PRIx64 "I64x"
#define SCNu64 "I64u"
#define SCNd64 "I64d"
#define SCNx64 "I64x"
#endif

#define STEP_SIEVING                      21
#define STEP_PRP_TESTING                  31
#define STEP_PRIMALITY_TESTING            41
#define STEP_MERGING_RESULTS              91
#define STEP_PHASE_DONE                   99

void setQuitting(int sig);
void getInputs(void);
void checkRecovery(void);
bool buildRemainingBitMap(void);
void buildNewBitMap(void);
void removeTrivial(void);
void removeRAlgeb(void);
void removeRAlgebMods(uint32_t &removed);
void removeSAlgeb(void);
void removeGFN(void);
void removeMOB(void);
void runfbncsieve(uint32_t n);
void processFbncAbcdFile(uint32_t n, char *fileName);
void doPhase(int phase);
void do_sieving(uint64_t max_k);
void do_prp_testing(void);
void do_primality_testing(void);
void process_results(const char *prpFileName, const char *primeFileName);
void merge_results(void);
void delete_temp_files(void);
void output_remain(const char *message);
void prepare_recovery(uint32_t step);
void delete_file(const char *fileName);
void checkForProgram(const char *programName);
double getAverageTestTime(uint64_t k, uint32_t n);
void verifyFbncsieveRanToCompletion(uint32_t n);
void verifySrsieve2RanToCompletion(void);
void verifyPfgwRanToCompletion(bool prpTest);

#ifdef _WIN32
#include <windows.h>

// Some old MinGW/CYGWIN distributions don't define this:
#ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
#define ENABLE_VIRTUAL_TERMINAL_PROCESSING  0x0004
#endif

void setQuitting(int sig)
{
   printf("CTRL-C accepted\n");
   gb_quitting = true;
}

static HANDLE stdoutHandle;
static DWORD outModeInit;

void setupConsole(void) {
	DWORD outMode = 0;
	stdoutHandle = GetStdHandle(STD_OUTPUT_HANDLE);

	if(stdoutHandle == INVALID_HANDLE_VALUE) {
		return;
	}
	
	if(!GetConsoleMode(stdoutHandle, &outMode)) {
		return;
	}

	outModeInit = outMode;
	
    // Enable ANSI escape codes
	outMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;

	if(!SetConsoleMode(stdoutHandle, outMode)) {
		return;
	}	
}

void restoreConsole(void) {
    // Reset colors
    printf("\x1b[0m");	
	
    // Reset console mode
	if(!SetConsoleMode(stdoutHandle, outModeInit)) {
		return;
	}
}
#else
void setupConsole(void) {}

void restoreConsole(void) {
    // Reset colors
    printf("\x1b[0m");
}
#endif

void error(const char *fmt, ...)
{
   va_list args;
   char    buffer1[500], buffer2[500];

   setupConsole();

   sprintf(buffer1, "\nError: ");

   va_start(args, fmt);
   vsprintf(buffer2, fmt, args);
   va_end(args);

   printf("\x1b[31m%s %s\n", buffer1, buffer2);

   restoreConsole();

   exit(0);
}

void reportWithDateTime(const char *fmt, ...)
{
   FILE *fPtr = fopen("srbsieve.log", "a+");

   time_t report_time;
   time(&report_time);
   struct tm *stm = localtime(&report_time);
   char    buffer1[500], buffer2[500];

   va_list args;

   va_start(args, fmt);
   vsprintf(buffer2, fmt, args);
   va_end(args);

   sprintf(buffer1, "Status (%04d-%02d-%02d %02d:%02d:%02d): ", stm->tm_year+1900, stm->tm_mon+1, stm->tm_mday, stm->tm_hour, stm->tm_min, stm->tm_sec);

   setupConsole();

   // Output the text in green
   printf("\x1b[32m%s %s\n", buffer1, buffer2);

   restoreConsole();

   fprintf(fPtr, "%s %s\n", buffer1, buffer2);

   fclose(fPtr);
}

int main(int argc, char **argv)
{
   int32_t  phase;

   gi_recoveryPhase = NO_PHASE;
   gi_recoverySubPhase = 0;
   
   gt_start_time = (uint32_t) time(NULL);

   getInputs();

   if (argc > 1)
   {
      if (atol(argv[1]) < 1 || atol(argv[1]) > gi_phaseCount)
         error("Invalid phase specified");
      
      phase = gi_recoveryPhase = atol(argv[1]);
      FILE    *ckpt = fopen("srbsieve.ckpt", "r");
      if (ckpt)
         error("srbsieve.ckpt was found.  Either delete that file or do not specify the starting phase on the command line");
   }
   else
      checkRecovery();

   if (gi_recoveryPhase != NO_PHASE)
   {
      if (!buildRemainingBitMap()) 
         error("pl_remain.txt not found.  File must exist when specifying starting phase on the command line");
      else
      {
         if (gi_recoverySubPhase == 0)
         {
            if (gl_remainingK == 0)
               error("%c%u phase %u k = %u to %u has 0 k's remaining.  mink or maxk may be invalid", (gi_c == 1 ? 'S' : 'R'), gi_base, gi_recoveryPhase, gl_minK, gl_maxK);
            else
               reportWithDateTime("Started %c%u phase %u.  %" PRIu64" %s remaining", (gi_c == 1 ? 'S' : 'R'), gi_base, gi_recoveryPhase, gl_remainingK, (gl_remainingK == 1 ? "k" : "k's"));
         }
         else
            reportWithDateTime("Continuing %c%u phase %u subphase %u.  %" PRIu64" %s remaining at start of phase", (gi_c == 1 ? 'S' : 'R'), gi_base, gi_recoveryPhase, gi_recoverySubPhase, gl_remainingK, (gl_remainingK == 1 ? "k" : "k's"));
      }
      
      phase = gi_recoveryPhase;
      gb_recovery = true;
   }

#ifdef SIGHUP
   signal(SIGHUP, SIG_IGN);
   signal(SIGQUIT, setQuitting);
#endif

   signal(SIGINT, setQuitting);
   signal(SIGTERM, setQuitting);

   if (!gb_recovery)
   {         
      if (buildRemainingBitMap()) 
         error("pl_remain.txt was found.  Either delete or specify the phase to start with");

      delete_file("pl_prime.txt");
      delete_file("pl_unknown.txt");
      delete_file("pl_MOB.txt");
      delete_file("pl_algeb.txt");
      delete_file("pl_GFN.txt");
      delete_file("pl_remain.txt");
      delete_file("pl_trivial.txt");

      buildNewBitMap();

      removeTrivial();

      if (gi_c == -1)
         removeRAlgeb();
      else
         removeSAlgeb();

      removeGFN();
      removeMOB();

      if (gl_remainingK == 0)
      {
         reportWithDateTime("Finished preliminary steps. 0 k's remaining");
         exit(0);
      }
      else
         reportWithDateTime("Finished preliminary steps. n = 1 is next: %d %s remaining", gl_remainingK, (gl_remainingK == 1 ? "k" : "k's"));

      for (uint32_t n=1; n<=gi_maxNfbncsieve; n++)
      {
         runfbncsieve(n);

         if (gl_remainingK == 0)
            exit(0);
      }

      if (gi_phaseCount == 0)
      {
         output_remain("Finished all initial steps");
         exit(0);
      }

      output_remain("Finished all initial steps.  Phase 1 is next");
      phase = 1;
   }

#ifdef WIN32
   checkForProgram("srsieve2.exe");
   checkForProgram("pfgw64.exe");
#else
   checkForProgram("srsieve2");
   checkForProgram("pfgw64");
#endif

   for ( ; phase<=gi_phaseCount; phase++)
   {
      doPhase(phase);

      if (gl_remainingK == 0)
         break;
   }

   delete_file("srbsieve.ckpt");
   exit(0);
}

void getInputs(void)
{
   char  buffer[100], *pos;
   FILE *ini = fopen("srbsieve.ini", "r");
   int   index;

   if (!ini)
      error("srbsieve.ini not found");

   gi_base = 0;
   gi_c = 0;
   gl_minK = 0;
   gl_maxK = 0;
   gi_phaseCount = 0;
   gi_maxNfbncsieve = 0;

   while (fgets(buffer, 100, ini) != NULL)
   {
      if (!memcmp(buffer, "base=", 5))
         gi_base = atoi(buffer+5);
      if (!memcmp(buffer, "c=", 2))
         gi_c = atoi(buffer+2);
      if (!memcmp(buffer, "mink=", 5))
         gl_minK = atoll(buffer+5);
      if (!memcmp(buffer, "maxk=", 5))
         gl_maxK = atoll(buffer+5);
      if (!memcmp(buffer, "maxNfbncsieve=", 14))
         gi_maxNfbncsieve = atoll(buffer+14);
      if (!memcmp(buffer, "maxKsrsieve2=", 13))
         gi_maxKsrsieve2 = atoll(buffer+13);
      if (!memcmp(buffer, "maxfbncprimes=", 14))
         gi_maxfbncprimes = atoll(buffer+14);

      if (!memcmp(buffer, "phase=", 6))
      {
         gi_phaseCount++;
         index = gi_phaseCount;

         if (gi_phaseCount > 100)
            error("Reached limit of 100 phase entries");

         pos = buffer+6;

         // gs_phases[0] will be empty

         if (sscanf(pos, "%d", &gs_phases[index].maxN) != 1)
            error("Could not process line %s", buffer);

         if (gs_phases[index].maxN <= gs_phases[index-1].maxN)
            error("Phases need to have increasing n");

         if (gi_phaseCount == 1)
         {
            if (gs_phases[1].maxN < gi_maxNfbncsieve + 10)
               error("Phase 1 must be at least 10 greater than maxNfbncsieve");

            gs_phases[1].minN = gi_maxNfbncsieve + 1;
         }
         else
            if (gs_phases[index].maxN < gs_phases[index-1].maxN + 10)
               error("Phases need to increment by at least 10");

         // Set minN for next phase
         gs_phases[index+1].minN = gs_phases[index].maxN + 1;
      }
   }

   fclose(ini);

   if (gi_base == 0)
      error("base must be specified in srbsieve.ini");
   if (gi_c == 0)
      error("c must be specified in srbsieve.ini");
   if (gl_minK == 0)
      error("mink must be specified in srbsieve.ini");
   if (gl_maxK == 0)
      error("maxk must be specified in srbsieve.ini");

   if (gi_base < 2)
      error("base cannot be less than 2");

   if (gi_base > 1000000000)
      error("base cannot be greater than 1e9");

   if (gi_c != 1 && gi_c != -1)
      error("c must be -1 or 1");

   if (gl_minK >= gl_maxK)
      error("mink must be less than maxk");

   if (gl_maxK > 1000000000000000000)
      error("maxk cannot be greater than 1e18");

   if (gi_maxNfbncsieve < 2)
   {
      if (sqrt(gl_minK) * gi_base < 1000000000)
         error("maxNfbncsieve must be at least 2 in srbsieve.ini");
      else
         if (gi_maxNfbncsieve < 1)
         {
            if (sqrt(gl_minK) * gi_base < 100000000000)
               error("maxNfbncsieve must be at least 1 in srbsieve.ini");
         }
   }

   if (gi_maxKsrsieve2 < 1000)
      error("maxKsrsieve2 must be at least 1000 in srbsieve.ini");
}

void checkRecovery(void)
{
   char     buffer[1000];
   FILE    *ckpt = fopen("srbsieve.ckpt", "r");
   int32_t  minN;
   int32_t  maxN;
   int32_t  index;

   if (!ckpt)
      return;

   while (fgets(buffer, 100, ckpt) != NULL)
   {
      if (!memcmp(buffer, "phaseInProgress=", 16))
         gi_recoveryPhase = atoi(buffer+16);
      if (!memcmp(buffer, "subphaseInProgress=", 19)) 
         gi_recoverySubPhase = atoi(buffer+19);
      if (!memcmp(buffer, "phase=", 6))
         if (sscanf(buffer, "phase=%d,%d", &minN, &maxN) != 2)
            error("Could not process line %s", buffer);
      if (!memcmp(buffer, "currentStep=", 12))
         gi_recoveryStep = atoi(buffer+12);
   }

   fclose(ckpt);

   if (gi_recoveryPhase < 1 || gi_recoveryPhase > gi_phaseCount)
      error("Invalid phase for recovery");
   if (gi_recoverySubPhase < 1)
      error("Invalid subphase for recovery");
   if (gi_recoveryStep < 9 || gi_recoveryStep > 99)
      error("Invalid recovery step");
   if (minN != gs_phases[gi_recoveryPhase].minN || maxN != gs_phases[gi_recoveryPhase].maxN)
      error("checkpoint phase does not match input phase");
}

void buildNewBitMap(void)
{
   uint64_t k = gl_minK;
   uint32_t removed = 0;

   gl_remainingK = gl_maxK - gl_minK + 1;

   reportWithDateTime("Started %c%u with %" PRIu64" k's", (gi_c == 1 ? 'S' : 'R'), gi_base, gl_remainingK);

   gb_KTerms.resize(gl_remainingK);
   std::fill(gb_KTerms.begin(), gb_KTerms.end(), true);

   if (gi_base & 1)
   {
      // For odd bases remove odd k so start with odd k
      if (!(k & 1))
         k++;

      for ( ; k<=gl_maxK; k+=2)
      {
         gb_KTerms[BIT(k)] = false;
         gl_remainingK--;
         removed++;
      }
      reportWithDateTime("Removed %u %s due to odd k on odd base", removed, (removed == 1 ? "k" : "k's"));
   }
}

bool buildRemainingBitMap(void)
{
   FILE *remain = fopen("pl_remain.txt", "r");
   char     buffer[1000];
   uint64_t k;

   if (!remain)
      return false;

   gl_remainingK = gl_maxK - gl_minK + 1;

   gb_KTerms.resize(gl_remainingK);
   std::fill(gb_KTerms.begin(), gb_KTerms.end(), false);

   gl_remainingK = 0;

   while (fgets(buffer, sizeof(buffer), remain) != NULL)
   {
      if (buffer[0] == '\n')
         continue;

      char *pos = strchr(buffer, '*');
      if (pos != NULL)
         *pos = 0;

      k = atoll(buffer);

      if (k < gl_minK || k > gl_maxK)
         continue;
      
      gb_KTerms[BIT(k)] = true;
      gl_remainingK++;
   }

   fclose(remain);
      
   return true;
}

void removeTrivial(void)
{
   int32_t  left = gi_base - 1;
   int32_t  p;
   int32_t  factors[10];
   int32_t  factorCount = 0;
   uint64_t k;
   uint32_t removed = 0;
   primegen pg;
   FILE    *trivials = NULL;

   // Find all prime divisors of the base - 1
   primegen_init(&pg);
   do
   {
      p = primegen_next(&pg);

      if (left % p == 0)
      {
         factors[factorCount++] = p;
         while (left % p == 0 && left > 0)
            left /= p;
      }
   } while (left > p);

   // If base is odd, skip factor 2
   int f = (gi_base & 1 ? 1 : 0);

   for ( ; f<factorCount; f++)
   {
      // Find the smallest k >= minK where (k + c) is divisible by factors[f]
      k = gl_minK + factors[f] - 1 - ((gl_minK + gi_c + factors[f] - 1) % factors[f]);

      for ( ; k <= gl_maxK; k += factors[f])
      {
         if (gb_KTerms[BIT(k)])
         {
            gb_KTerms[BIT(k)] = false;
            gl_remainingK--;
            removed++;
            
            if (!trivials)
               trivials = fopen("pl_trivial.txt", "w");
            
            fprintf(trivials, "%" PRId64"\n", k);
         }
      }
   }

   if (trivials)
      fclose(trivials);

   if (removed > 0)
      reportWithDateTime("Removed %u %s due to trivial factorization", removed, (removed == 1 ? "k" : "k's"));
}

void removeRAlgeb(void)
{
   FILE    *algeb = NULL;
   uint32_t pwr[10] = {0};
   int32_t  pwr_cnt = 0;
   uint32_t removed = 0;
   uint32_t work;
   uint64_t k;
   uint64_t k_pwr;

   // check if base is a prime power
   // set max power to floor (log(base) / log(2))
   uint32_t max_pwr = round( (double)log(gi_base) / (double)log(2) );
   if (pow(2, max_pwr) > gi_base)
      max_pwr--;

   uint32_t prime_tbl[10] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
   for (int i=0; prime_tbl[i]<=max_pwr; i++)
   {
      work = round( (double)pow(gi_base, 1.0/prime_tbl[i]) );
      work = pow(work, prime_tbl[i]);

      if (gi_base == work)
         pwr[pwr_cnt++] = prime_tbl[i];
   }

   // Find all k's that have a same power as the base.
   // They will be removed as they cannot have primes.
   // First set k to the smallest pwr[i] power >= minK.
   for (int i=0; i<pwr_cnt; i++)
   {
      k_pwr = round( (double)pow(gl_minK, 1.0/pwr[i]) );
      if (pow(k_pwr, pwr[i]) < gl_minK)
         k_pwr++;

      k = pow(k_pwr, pwr[i]);
      while (k <= gl_maxK)
      {
         if (gb_KTerms[BIT(k)])
         {
            gb_KTerms[BIT(k)] = false;
            gl_remainingK--;
            removed++;

            if (!algeb)
               algeb = fopen("pl_algeb.txt", "w");

            fprintf(algeb, "%" PRId64"\n", k);
         }
         k_pwr++;
         k = pow(k_pwr, pwr[i]);
      }
   }

   if (algeb)
      fclose(algeb);

   if (pwr[0] != 2)
      removeRAlgebMods(removed);

   if (removed > 0)
      reportWithDateTime("Removed %u %s due to algebraic factorization", removed, (removed == 1 ? "k" : "k's"));
}

void removeRAlgebMods(uint32_t &removed)
{
   FILE    *algeb = NULL;
   int32_t  left = gi_base + 1;
   int32_t  factors[30];
   int32_t  factorCount = 0;
   int32_t  p;
   primegen pg;
   uint32_t work;
   uint64_t k;
   uint64_t k_root;
   uint64_t k_cmpr;
   uint64_t start_root;
   uint32_t alg_fact[10];
   uint32_t alg_mod[10][2];
   int32_t  alg_cnt = 0;
   int32_t  k_cnt = 0;
   int32_t  size = round( sqrt(gl_maxK - gl_minK) * 1.5 ) + 10;
   uint64_t *k_tbl = new uint64_t[size];

   // For unsquared bases, find k's that have partial algebraic factors to
   // make a full covering set. Where k = m^2, if base + 1 has prime factor
   // f % 4 = 1 and sqrt(f-1 mod f) = m, remove the k.
   // Example base 337: 338 has prime factors f = 2 & 13. 13 % 4 = 1
   // so sqrt(12 mod 13) = both 5 and 8 since 5^2 and 8^2 % 13 = 12.
   // Hence k = 5^2, 8^2, 18^2, 21^2, 31^2, 34^2, etc. as well as
   // k = base * 5^2, 8^2, 18^2, 21^2, 31^2, 34^2, etc. will be dropped.
   // There can be multiple prime factors f % 4 = 1 that must be checked.

   // find all prime divisors of base + 1
   primegen_init(&pg);
   do
   {
      p = primegen_next(&pg);

      if (left % p == 0)
      {
         factors[factorCount++] = p;
         while (left % p == 0)
            left /= p;
      }
   } while (left > p);

   // find algebraic modulos if factor % 4 = 1
   for (int i=0; i<factorCount; i++)
   {
      if (factors[i] % 4 == 1)
      {
         for (uint64_t j=1; j<factors[i]; j++)
         {
            if (j*j % factors[i] == factors[i] - 1)
            {
               alg_fact[alg_cnt] = factors[i];
               alg_mod[alg_cnt][0] = j;
               alg_mod[alg_cnt][1] = factors[i] - j;
               alg_cnt++;
               break;
            }
         }
      }
   }

   // remove k's that are algebraic modulos squared
   // start at first squared k >= minK
   for (int i=0; i<alg_cnt; i++)
   {
      k_root = round( (double)sqrt(gl_minK) );
      if (k_root * k_root < gl_minK)
         k_root++;

      work = k_root % alg_fact[i];
      start_root = k_root;
      for (int j=0; j<2; j++)
      {
         k_root = start_root - work + alg_mod[i][j];
         if (work > alg_mod[i][j])
            k_root += alg_fact[i];

         k = k_root * k_root;
         while (k <= gl_maxK)
         {
            k_tbl[k_cnt] = k;
            k_cnt++;
            k_root += alg_fact[i];
            k = k_root * k_root;
         }
      }
   }

   // If a base has k's that contain algebraic factors per the above condition,
   // it may have additional k's that contain algebraic factors for k's
   // that are divisible by the base or if the base contains squares.
   if (alg_cnt > 0)
   {
      // Find all prime divisors of the base.  Store all
      // occurrences of factors to later find squares.
      left = gi_base;
      factorCount = 0;
      primegen_init(&pg);
      do
      {
         p = primegen_next(&pg);

         while (left % p == 0)
         {
            factors[factorCount++] = p;
            left /= p;
         }
      } while (left >= p);

      // See if base contains a square by looking
      // for factors that occur more than once.
      uint32_t sqr_fact = 1;
      for (int i=1; i<factorCount; i++)
      {
         if (factors[i] == factors[i-1])
         {
            sqr_fact *= factors[i] * factors[i];
            i++;
         }
      }

      if (sqr_fact == 1)
      {
         // Base does not contain a square.  Remove k's that are base * algebraic
         // modulos squared.  Start at first base * k squared >= minK.
         for (int i=0; i<alg_cnt; i++)
         {
            k_root = round( (double)sqrt(gl_minK / gi_base) );
            k_cmpr = gi_base * k_root * k_root;
            if (k_cmpr < gl_minK)
               k_root++;

            work = k_root % alg_fact[i];
            start_root = k_root;
            for (int j=0; j<2; j++)
            {
               k_root = start_root - work + alg_mod[i][j];
               if (work > alg_mod[i][j])
                  k_root += alg_fact[i];				

               k = gi_base * k_root * k_root;
               while (k <= gl_maxK)
               {
                  k_tbl[k_cnt] = k;
                  k_cnt++;
                  k_root += alg_fact[i];
                  k = gi_base * k_root * k_root;
               }
            }
         }
      }
      else
      {
         // For bases that contain squares, find k's that have partial algebraic
         // factors to make a full covering set. Where s = square-free part of
         // base and k = s * m^2, if base + 1 has prime factor f % 4 = 1 and
         // sqrt(1 mod f) = m, remove the k. Example for base 414:
         // 415 has prime factors f = 5 & 83. 5 % 4 = 1 so sqrt(1 mod 5) =
         // both 1 and 4 since 1^2 and 4^2 % 5 = 1.  Here s = 414 / 9 = 46.
         // Hence k = 46*1^2, 46*4^2, 46*6^2, 46*9^2, etc. will be dropped.
         uint32_t base_nosqr = gi_base / sqr_fact;
         int32 alg_cnt_save = alg_cnt;
         for (int i=0; i<alg_cnt_save; i++)
         {
            for (uint64_t j=1; j<alg_fact[i]; j++)
            {
               if (base_nosqr*j*j % alg_fact[i] == 1)
               {
                  alg_fact[alg_cnt] = alg_fact[i];
                  alg_mod[alg_cnt][0] = j;
                  alg_mod[alg_cnt][1] = alg_fact[i] - j;
                  alg_cnt++;
                  break;
               }
            }
         }
         // remove k's that are a square-free part of
         // the base * algebraic modulos squared
         // start at first square-free part * squared k >= minK
         for (int i=alg_cnt_save; i<alg_cnt; i++)
         {
            k_root = round( (double)sqrt(1.0 * gl_minK / base_nosqr) );
            k_cmpr = base_nosqr * k_root * k_root;
            if (k_cmpr < gl_minK)
               k_root++;

            work = k_root % alg_fact[i];
            start_root = k_root;
            for (int j=0; j<2; j++)
            {
               k_root = start_root - work + alg_mod[i][j];
               if (work > alg_mod[i][j])
                  k_root += alg_fact[i];

               k = base_nosqr * k_root * k_root;
               while (k <= gl_maxK)
               {
                  k_tbl[k_cnt] = k;
                  k_cnt++;
                  k_root += alg_fact[i];
                  k = base_nosqr * k_root * k_root;
               }
            }
         }
      }
   }

   for (int i=0; i<k_cnt; i++)
   {
      if (gb_KTerms[BIT(k_tbl[i])])
      {
         gb_KTerms[BIT(k_tbl[i])] = false;
         gl_remainingK--;
         removed++;

         if (!algeb)
            algeb = fopen("pl_algeb.txt", "a");

         fprintf(algeb, "%" PRId64"\n", k_tbl[i]);
      }
   }

   delete[] k_tbl;

   if (algeb)
      fclose(algeb);
}

void removeSAlgeb(void)
{
   FILE *algeb = NULL;
   uint32_t pwr[10] = {0};
   int32_t  pwr_cnt = 0;
   uint32_t removed = 0;
   uint32_t work;
   uint64_t k;
   uint64_t k_pwr;
   int32_t  k_cnt = 0;
   int32_t  size = round( pow((gl_maxK - gl_minK), 1.0/3) * 1.5 ) + 10;
   uint64_t *k_tbl = new uint64_t[size];

   // check if base is a prime power > 2
   // set max power to floor (log(base) / log(2))
   uint32_t max_pwr = round( (double)log(gi_base) / (double)log(2) );
   if (pow(2, max_pwr) > gi_base)
      max_pwr--;

   uint32_t prime_tbl[9] = {3, 5, 7, 11, 13, 17, 19, 23, 29};
   for (int i=0; prime_tbl[i]<=max_pwr; i++)
   {
      work = round( (double)pow(gi_base, 1.0/prime_tbl[i]) );
      work = pow(work, prime_tbl[i]);

      if (gi_base == work)
         pwr[pwr_cnt++] = prime_tbl[i];
   }

   // Find all k's that have a same power as the base.
   // They will be removed as they cannot have primes.
   // First set k to the smallest pwr[i] power >= minK.
   for (int i=0; i<pwr_cnt; i++)
   {
      k_pwr = round( (double)pow(gl_minK, 1.0/pwr[i]) );
      if (pow(k_pwr, pwr[i]) < gl_minK)
         k_pwr++;

      k = pow(k_pwr, pwr[i]);

      while (k <= gl_maxK)
      {
         k_tbl[k_cnt] = k;
         k_pwr++;
         k_cnt++;
         k = pow(k_pwr, pwr[i]);
      }
   }

   // Look for bases that are a 4th power and the k is 4 times a 4th power.
   // The k's will be removed as they cannot have primes.
   work = round( (double)pow(gi_base, 0.25) );
   work = pow(work, 4);

   if (gi_base == work)
   {
      // set k to the smallest 4 * 4th power >= minK
      k_pwr = round( (double)pow(gl_minK / 4.0, 0.25) );

      if (4 * pow(k_pwr, 4) < gl_minK)
         k_pwr++;

      k = 4 * pow(k_pwr, 4);

      while (k <= gl_maxK)
      {
         k_tbl[k_cnt] = k;
         k_pwr++;
         k_cnt++;
         k = 4 * pow(k_pwr, 4);
      }
   }

   for (int i=0; i<k_cnt; i++)
   {
      if (gb_KTerms[BIT(k_tbl[i])])
      {
         gb_KTerms[BIT(k_tbl[i])] = false;
         gl_remainingK--;
         removed++;

         if (!algeb)
            algeb = fopen("pl_algeb.txt", "w");

         fprintf(algeb, "%" PRId64"\n", k_tbl[i]);
      }
   }

   delete[] k_tbl;

   if (algeb)
      fclose(algeb);

   if (removed > 0)
      reportWithDateTime("Removed %u %s due to algebraic factorization", removed, (removed == 1 ? "k" : "k's"));
}

void removeGFN(void)
{
   // Only b^n+1 can be GFNs and b must be even, else b^n+1 is always even.
   if ((gi_base & 1) || gi_c == -1)
      return;

   FILE *gfns = NULL;
   uint32_t max_pwr = round( (double)log(gi_base) / (double)log(2) );
   uint32_t root_pwr;
   int32_t  root, work = 0;
   uint32_t removed = 0;
   uint64_t k;

   // Find the lowest root for the base, e.g. if base = 216, then root = 6 (216 = 6^3).
   // Work from largest to smallest power.
   for (root_pwr = max_pwr; work != gi_base; root_pwr--)
   {
      root = round( (double)pow(gi_base, 1.0/root_pwr) );
      work = pow(root, root_pwr);
   }

   // set k to the smallest root ^ root_pwr >= minK
   root_pwr = round( (double)log(gl_minK) / (double)log(root) );
   if (pow(root, root_pwr) < gl_minK)
      root_pwr++;

   k = pow(root, root_pwr);

   // Find all k where k*b^n+1 = b^m+1 for some integer m.  They will be
   // removed since they are GFNs.  GFNs can only be prime when m=2^q
   // so are not considered in the conjectures.
   for ( ; k <= gl_maxK; k *= root)
   {
      if (gb_KTerms[BIT(k)])
      {
         gb_KTerms[BIT(k)] = false;
         gl_remainingK--;
         removed++;

         if (!gfns)
            gfns = fopen("pl_GFN.txt", "w");

         fprintf(gfns, "%" PRId64"\n", k);
      }
   }

   if (gfns)
      fclose(gfns);

   if (removed > 0)
      reportWithDateTime("Removed %u %s due to GFN", removed, (removed == 1 ? "k" : "k's"));
}

void removeMOB(void)
{
   uint32_t maxFact = round(sqrt(gl_maxK + gi_c)) + 1;
   uint32_t estPrimes = maxFact / log(maxFact) * 1.15 + 10;
   primegen pg;
   uint64_t *factors;
   uint64_t k, k_init;
   uint32_t removed = 0, adder;
   uint32_t bit, byte;
   uint8_t *sieveBitMap;
   FILE    *mobs = NULL;

   factors = (uint64_t *) malloc(estPrimes * sizeof(uint64_t));

   primegen_init(&pg);
   factors[0] = primegen_next(&pg);
   for (int pi=1; factors[pi-1]<=maxFact; pi++)
      factors[pi] = primegen_next(&pg);

   sieveBitMap = (uint8_t *) malloc(1 + (size_t) 1 + (gl_maxK - gl_minK + 1 / 8));
   memset(sieveBitMap, 0xff, 1 + (size_t) 1 + (gl_maxK - gl_minK + 1 / 8));

   // Initialize k to the smallest k >= minK that is divisible by the base
   k_init = gl_minK + gi_base - 1 - ((gl_minK - 1) % gi_base);

   // Sieve through all k between minK and maxK that are divisible
   // by the base to find those where (k + c) is composite.
   for (int i=0; factors[i] < maxFact; i++)
   {
      if (gi_base % factors[i] == 0)
         continue;

      // Use the largest of factors[i] and (base + c) for loop incrementing.
      if (factors[i] >= gi_base + gi_c)
      {
         // Find the smallest k >= k_init where (k + c) is divisible by factors[i].
         k = k_init + factors[i] - 1 - ((k_init + gi_c - 1) % factors[i]);
         
         if (k + gi_c == factors[i])
            k += factors[i];

         // (k + c) is now > factors[i] and divisible by factors[i].
         // Increment by factors[i] until k is divisible by the base.
         while (k % gi_base != 0 && k <= gl_maxK)
            k += factors[i];
      }
      else
      {
         // (k + c) is now > factors[i] and k is divisible by the base.
         // Increment by the base until (k + c) is divisible by factors[i].
         k = k_init;
         while ((k + gi_c) % factors[i] != 0 && k <= gl_maxK)
            k += gi_base;
      }

      for ( ; k <= gl_maxK; k += (factors[i] * gi_base))
      {
         byte = (uint32_t) ((k - gl_minK) >> 3);
         bit = (uint8_t) (1 << ((k - gl_minK) & 7));

         // Mark as composite
         sieveBitMap[byte] &= ~bit;
      }
   }

   // Reset k to the smallest k >= minK that is divisible by the base
   k = k_init;

   // If base and k are odd, make k even
   if ((gi_base & 1) && (k & 1))
      k += gi_base;

   // If the base is odd, then add base*2 to get the next k
   adder = (gi_base & 1 ? (gi_base*2) : gi_base);

   // For each k between minK and maxK that is divisible by the base, if (k + c) is
   // composite then k will have the same prime as k / base so it can be removed.
   for ( ; k <= gl_maxK; k += adder)
   {
      byte = (uint32_t) ((k - gl_minK) >> 3);
      bit = (uint8_t) (1 << ((k - gl_minK) & 7));

      if (!(sieveBitMap[byte] & bit))
      {
         if (gb_KTerms[BIT(k)])
         {
            gb_KTerms[BIT(k)] = false;
            gl_remainingK--;
            removed++;
            
            if (!mobs)
               mobs = fopen("pl_MOB.txt", "w");
            
            fprintf(mobs, "%" PRId64"\n", k);
         }
      }
   }

   free(factors);
   free(sieveBitMap);

   if (mobs)
      fclose(mobs);

   if (removed > 0)
      reportWithDateTime("Removed %u %s due to MOB", removed, (removed == 1 ? "k" : "k's"));
}

void runfbncsieve(uint32_t n)
{
   char     program[50];
   char     sequence[50];
   char     command[200];
   char     fileName[100];
   FILE    *fPtr;

   sprintf(fileName, "n%u.abcd", n);
   
#ifdef WIN32
   fPtr = fopen("fbncsieve.exe", "r");
   
   sprintf(program, "fbncsieve");
   sprintf(sequence, "k*%u^^%u%+d", gi_base, n, gi_c);
#else
   fPtr = fopen("fbncsieve", "r");

   sprintf(program, "./fbncsieve");
   sprintf(sequence, "k*%u^%u%+d", gi_base, n, gi_c);
#endif

   if (!fPtr)
      error("fbncsieve not found");
   
   fclose(fPtr);

   delete_file("fbncsieve.log");
   
   if (gi_maxfbncprimes > 0)
      sprintf(command, "%s -fD -k%llu -K%llu -s%s -w%uf -o %s", program, gl_minK, gl_maxK, sequence, gi_maxfbncprimes, fileName);
   else
      sprintf(command, "%s -fD -k%llu -K%llu -s%s -o %s", program, gl_minK, gl_maxK, sequence, fileName);

   printf("command: %s\n", command);
   system(command);
   
   // This will stop the program if fbncsieve was stopped before completing
   verifyFbncsieveRanToCompletion(n);
   
   processFbncAbcdFile(n, fileName);
   
   delete_file(fileName);
}

void processFbncAbcdFile(uint32_t n, char *fileName)
{
   char     buffer[100], message[200];
   char    *pos;
   uint64_t k, diff, pmax;
   uint32_t removed = 0;
   uint32_t base;
   int32_t  c;
   
   if (n > gs_phases[0].minN)
      gs_phases[0].minN = n + 1;
   
   FILE *fptr = fopen(fileName, "r");
   FILE *primes = fopen("pl_prime.txt", "a");

   if (!fptr)
   {
      reportWithDateTime("fbncsieve did not generate an abcd file.  Assuming all terms had factors and none are prime.  Continuing", fileName);
      return;
   }

   // Read the first line
   if (fgets(buffer, 100, fptr) == NULL)
   {
      fclose(fptr);
      fclose(primes);

      reportWithDateTime("Removed 0 k's from file for n = %u: %" PRIu64" %s remaining", n, gl_remainingK, (gl_remainingK == 1 ? "k" : "k's"));
      return;      
   }

   if (sscanf(buffer, "ABCD $a*%u^%d%d  [%" SCNu64"] // Sieved to %" SCNu64"", &base, &n, &c, &k, &pmax) != 5)
      error("Invalid first line of file %s", fileName);

   if (base != gi_base)
      error("Read base %d from %s", base, fileName);

   if (c != +1 && c != -1)
      error("Read c %+d from %s", c, fileName);

   // Note that fbncsieve will only output even k if the base is odd.
   // Fortunately sbrsieve will already have set odd k to false in the vector.
   if (gb_KTerms[BIT(k)])
   {
      gb_KTerms[BIT(k)] = false;
      gl_remainingK--;
      removed++;
      
      fprintf(primes, "%" PRId64"*%d^%d%+d\n", k, gi_base, n, gi_c);
   }

   while (fgets(buffer, sizeof(buffer), fptr) != NULL)
   {
      if (sscanf(buffer, "%" SCNu64 , &diff) != 1)
         error("Line %s is malformed", buffer);
      
      k += diff;

      if (gb_KTerms[BIT(k)])
      {
         gb_KTerms[BIT(k)] = false;
         gl_remainingK--;
         removed++;
         
         fprintf(primes, "%" PRId64"*%d^%d%+d\n", k, gi_base, n, gi_c);
      }
   }
   
   fclose(fptr);
   fclose(primes);

   sprintf(message, "Removed %u %s from file for n = %u", removed, (removed == 1 ? "k" : "k's"), n);
   reportWithDateTime("%s: %u %s remaining", message, gl_remainingK, (gl_remainingK == 1 ? "k" : "k's"));
}

void doPhase(int32_t phase)
{
   FILE     *sieve_in = NULL;
   uint32_t  count = 0;
   uint64_t  startingK = 0, completedK = 0;
   int32_t   subPhase = 0;
   uint64_t  k, max_k, min_k;
   int64_t   index;
   uint32_t  xOfY, chunks;

   gi_currentPhase = phase;
   gi_currentSubPhase = 0;

   gt_start_time = (uint32_t) time(NULL);

   k = gl_minK;

   // Skip over all k tested so far
   if (gb_recovery && gi_recoveryPhase > 0)
   {
      if (gi_recoverySubPhase == 0)
      {
         printf("Recovering phase %d with %" PRIu64" %s\n", gi_recoveryPhase, gl_remainingK, (gl_remainingK == 1 ? "k" : "k's"));
         gi_recoveryPhase = NO_PHASE;
      }
      else 
      {
         printf("Recovering from phase %d, subphase %d with %" PRIu64" %s\n", gi_recoveryPhase, gi_recoverySubPhase, gl_remainingK, (gl_remainingK == 1 ? "k" : "k's"));

         index = (gi_recoverySubPhase - 1) * gi_maxKsrsieve2 + 1;
         max_k = 0;

         for (; ; k++)
         {         
            if (gb_KTerms[BIT(k)])
            {
               max_k++;
               
               if (index == max_k)
                  break;
            }
         }

         gi_currentSubPhase = gi_recoverySubPhase - 1;
         gi_recoverySubPhase = 0;
         completedK = (index - 1);
      }
   }

   startingK = gl_remainingK;
   xOfY = 0;
   chunks = 1 + (gl_remainingK / gi_maxKsrsieve2);

   for ( ; k<=gl_maxK; k++)
   {
      if (gb_KTerms[BIT(k)])
      {
         max_k = k;

         if (count == 0)
         {
            gi_currentSubPhase++;
            sieve_in = fopen("sieve.in", "w");
            min_k = k;
         }

         fprintf(sieve_in, "%" PRId64"*%d^n%+d\n", k, gi_base, gi_c);
         count++;

         if (count == gi_maxKsrsieve2)
         {
            fclose(sieve_in);
            sieve_in = NULL;

            xOfY++;

            printf("Phase %d:  Processing k from %" PRIu64" to %" PRIu64". Chunk %u of %u.\n", gi_currentPhase, min_k, max_k, xOfY, chunks);

            do_sieving(max_k);
            do_prp_testing();
            do_primality_testing();
            delete_temp_files();

            completedK += count;
            count = 0;

            reportWithDateTime("Completed %.02f pct of phase %d.", 100.0 * ((double) completedK) / ((double) startingK), phase);
            printf("\n\n");
         }
      }
   }

   if (sieve_in)
      fclose(sieve_in);

   if (count > 0)
   {
      printf("Phase %d:  Processing k from %" PRIu64" to %" PRIu64"\n", gi_currentPhase, min_k, max_k);

      do_sieving(max_k);
      do_prp_testing();
      do_primality_testing();
      delete_temp_files();

      completedK += count;
   }

   merge_results();

   if (gl_remainingK == 0)
      reportWithDateTime("Completed phase %u.  0 k's remaining", phase);
   else
      reportWithDateTime("Completed phase %u to n = %u.  %" PRIu64" %s remaining", phase, gs_phases[gi_currentPhase].maxN, gl_remainingK, (gl_remainingK == 1 ? "k" : "k's"));

   printf("\n\n");

   prepare_recovery(STEP_PHASE_DONE);

   gb_recovery = false;
}

void do_sieving(uint64_t max_k)
{
   char      command[500];
   char      program[100];
   char      fileName[100];
   char      rate[20];
   char      worksize[10];
   FILE     *fPtr;
   uint32_t  n_min = gs_phases[gi_currentPhase].minN;
   uint32_t  n_max = gs_phases[gi_currentPhase].maxN;

   if (gi_recoveryPhase != NO_PHASE)
      if (gi_recoveryStep != STEP_SIEVING)
         return;

   double averageTestTime = getAverageTestTime(max_k, n_min + 2*((n_max - n_min) / 3));

   prepare_recovery(STEP_SIEVING);

   delete_file("srsieve2.log");

   sprintf(fileName, "sr_%c%u.pfgw", (gi_c == 1 ? 'S' : 'R'), gi_base);

#ifdef WIN32
   strcpy(program, "srsieve2");
   checkForProgram("srsieve2.exe");
#else
   strcpy(program, "./srsieve2");
   checkForProgram("srsieve2");
#endif

   // If average test time >= 1 second, use average test time / 2 to determine the last number of minutes that srsieve2
   // uses for the average removal rate.  This will allow for >= ~30 factors to be found to determine the removal rate.
   // Example: 30 seconds / factor would cause srsieve2 to use the last 30 / 2 = 15 minutes for the removal rate.

   if (averageTestTime < 1.0)
      sprintf(rate, "-4 %f", 1.0 / averageTestTime);
   else
   {
      uint32_t sieve_mins = (uint32_t) round(averageTestTime / 2.0);
      sprintf(rate, "-5 %f -6 %u", averageTestTime, sieve_mins);
   }

   // If many k are remaining for a small test depth, use smaller work size
   // to allow srsieve2 to stop sooner when removal rate is reached.
   if (n_max <= 200 && gl_remainingK >= 50000)
      sprintf(worksize, "1e3");
   else
      sprintf(worksize, "1e4");

   fPtr = fopen(fileName, "r");
   if (fPtr == NULL)
      sprintf(command, "%s -fP -l0 %s -o%s -w%s -n%d -N%d -s sieve.in ", program, rate, fileName, worksize, n_min, n_max);
   else
   {
      sprintf(command, "%s -fP -l0 %s -o%s -w%s -i%s ", program, rate, fileName, worksize, fileName);
      fclose(fPtr);
   }

   printf("command: %s\n", command);
   system(command);

   verifySrsieve2RanToCompletion();
}

double getAverageTestTime(uint64_t k, uint32_t n)
{
   char command[500];
   char buffer[1000];

   delete_file("pfgw.ini");
   delete_file("rate.out");
   delete_file("rate.in");

   if (k < 10)
      k = 10;

   FILE *rate_in = fopen("rate.in", "w");

   // Take the average of 10 tests with varying k because each k could use
   // a different FFT size.
   fprintf(rate_in, "%" PRIu64"*%u^%u%+d\n", k-9, gi_base, n, gi_c);
   fprintf(rate_in, "%" PRIu64"*%u^%u%+d\n", k-8, gi_base, n, gi_c);
   fprintf(rate_in, "%" PRIu64"*%u^%u%+d\n", k-7, gi_base, n, gi_c);
   fprintf(rate_in, "%" PRIu64"*%u^%u%+d\n", k-6, gi_base, n, gi_c);
   fprintf(rate_in, "%" PRIu64"*%u^%u%+d\n", k-5, gi_base, n, gi_c);
   fprintf(rate_in, "%" PRIu64"*%u^%u%+d\n", k-4, gi_base, n, gi_c);
   fprintf(rate_in, "%" PRIu64"*%u^%u%+d\n", k-3, gi_base, n, gi_c);
   fprintf(rate_in, "%" PRIu64"*%u^%u%+d\n", k-2, gi_base, n, gi_c);
   fprintf(rate_in, "%" PRIu64"*%u^%u%+d\n", k-1, gi_base, n, gi_c);
   fprintf(rate_in, "%" PRIu64"*%u^%u%+d\n", k-0, gi_base, n, gi_c);

   fclose(rate_in);

#ifdef WIN32
   sprintf(command, "pfgw64 -Cquiet -lrate.out -f0 rate.in > nul");
   checkForProgram("pfgw64.exe");
#else
   sprintf(command, ".\pfgw64 -Cquiet -lrate.out -f0 rate.in > nul");
   checkForProgram("pfgw64");
#endif

   system(command);

   FILE *rate_out = fopen("rate.out", "r");
   double testTime, totalTime = 0.0;
   double count = 0.0;
   while (fgets(buffer, 1000, rate_out) != NULL)
   {
      if (strstr(buffer, "inf") != NULL)
         continue;

      char *pos1 = strchr(buffer, '(');
      if (pos1 == NULL)
         error("could not get rate from rate.out (cannot find start of time)");

      char *pos2 = strchr(pos1, 's');
      if (pos2 == NULL)
         error("could not get rate from rate.out (cannot find end of time)");

      *pos1 = *pos2 = 0;
      if (sscanf(pos1+1, "%lf", &testTime) != 1)
         error("could not get rate from rate.out (failed scan)");

      totalTime += testTime;
      count += 1.0;
   }

   fclose(rate_out);

   delete_file("pfgw.ini");
   delete_file("rate.out");
   delete_file("rate.in");
   delete_file("pfgw.log");

   if (count == 0.0 || totalTime == 0.0)
      return 0.00001;

   return totalTime / count;
}

void do_prp_testing(void)
{
   FILE     *pfgw_in;
   FILE     *pfgw_out;
   char      fileName[300];
   char      buffer[1000];
   char      command[500];
   char     *pos;
   uint32_t  n_max = gs_phases[gi_currentPhase].maxN;

   if (gi_recoveryPhase != NO_PHASE)
   {
      if (gi_recoveryStep != STEP_PRP_TESTING)
         return;
   }

   sprintf(fileName, "sr_%c%u.pfgw", (gi_c == 1 ? 'S' : 'R'), gi_base);

   gi_recoveryPhase = NO_PHASE;   
   prepare_recovery(STEP_PRP_TESTING);

#ifdef WIN32
   sprintf(command, "pfgw64 -Cquiet -f0%s %s > nul", (n_max > MIN_N_FOR_RESIDUES ? " -l" : ""), fileName);
   checkForProgram("pfgw64.exe");
#else
   sprintf(command, "./pfgw64 -Cquiet -f0%s %s > /dev/null", (n_max > MIN_N_FOR_RESIDUES ? " -l" : ""), fileName);
   checkForProgram("pfgw64");
#endif

   printf("command: %s\n", command);
   system(command);

   verifyPfgwRanToCompletion(true);

   if (n_max > MIN_N_FOR_RESIDUES)
   {
      sprintf(fileName, "results-%c%u-%u.txt", (gi_c == 1 ? 'S' : 'R'), gi_base, n_max);
      delete_file(fileName);
      rename("pfgw.out", fileName);
   }
}

void do_primality_testing(void)
{
   char      command[500];
   char      fileName[300];
   FILE     *prpFile;
   uint32_t  n_max = gs_phases[gi_currentPhase].maxN;

   if (gi_recoveryPhase != NO_PHASE)
   {
      if (gi_recoveryStep != STEP_PRIMALITY_TESTING)
         return;
   }
   else
   {
      prpFile = fopen("sr_b.prp", "r");

      if (!prpFile)
      {
         prpFile = fopen("pfgw.log", "r");

         // This could happen if no PRPs were found
         if (!prpFile) {
            printf("Neither sr_b.prp nor pfgw.log were found.  Skipping primality testing step.\n");
            return;
         }

         fclose(prpFile);

         rename("pfgw.log", "sr_b.prp");
      }
      else
         fclose(prpFile);
      
      // Yes, this will cause pfgw to start from the beginning of the file
      // but it shouldn't cost too much time to retest whatever is in the file.
      delete_file("pfgw.ini");
      delete_file("pfgw-prime.log");
   }

   prepare_recovery(STEP_PRIMALITY_TESTING);

#ifdef WIN32
   sprintf(command, "pfgw64 -Cquiet -f0%s -t%c sr_b.prp > nul", (n_max > MIN_N_FOR_RESIDUES ? " -l" : ""), (gi_c == 1 ? 'm' : 'p'));
#else
   system("./pfgw64 -Cquiet -f0%s -t%c sr_b.prp > /dev/null", (n_max > MIN_N_FOR_RESIDUES ? " -l" : ""), (gi_c == 1 ? 'm' : 'p'));
#endif

   printf("command: %s\n", command);
   system(command);

   verifyPfgwRanToCompletion(false);

   if (n_max > MIN_N_FOR_RESIDUES)
   {
      sprintf(fileName, "results-prime-%c%u-%u.txt", (gi_c == 1 ? 'S' : 'R'), gi_base, n_max);
      delete_file(fileName);
      rename("pfgw.out", fileName);
   }

   process_results("sr_b.prp", "pfgw-prime.log");

   delete_file("sr_b.prp");
}

void process_results(const char *prpFileName, const char *primeFileName)
{
   char buffer[500];
   uint64_t k;
   FILE *fPtr;
   FILE *out;
   
   out = fopen("results.tmp", "a");
   
   std::vector<bool> ksWithVerifiedPrime; 
   ksWithVerifiedPrime.resize(gl_maxK - gl_minK + 1);
   
   std::fill(ksWithVerifiedPrime.begin(), ksWithVerifiedPrime.end(), false);
   
   fPtr = fopen(primeFileName, "r");
   
   if (fPtr)
   {
      while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
      {
         if (strlen(buffer) < 3)
            continue;
         
         if (sscanf(buffer, "%" SCNu64"*", &k) != 1)
            error("missing k from %s line %s\n", primeFileName, buffer);
         
         ksWithVerifiedPrime[k - gl_minK] = true;
      }
      
      fclose(fPtr);
   }

   fPtr = fopen(prpFileName, "r");
   
   if (fPtr)
   {
      while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
      {
         if (strlen(buffer) < 3)
            continue;
         
         if (sscanf(buffer, "%" SCNu64"*", &k) != 1)
            error("missing k from %s line %s\n", prpFileName, buffer);
         
         if (ksWithVerifiedPrime[k - gl_minK])
            fprintf(out, "1,%s", buffer);
         else
            fprintf(out, "0,%s", buffer);
      }
      
      fclose(fPtr);
   }

   fclose(out);
}

void merge_results(void)
{
   FILE    *in;
   FILE    *out1, *out2, *out3 = NULL;
   char     buffer[1000], *pos;
   uint64_t k;

   in = fopen("results.tmp", "r");
   if (!in) return;

   out2 = fopen("pl_prime.txt", "a");
   
   while (fgets(buffer, 1000, in) != NULL)
   {
      if (buffer[0] == '1')
         fprintf(out2, "%s", buffer+2);
      else
      {
         if (out3 == NULL)
            out3 = fopen("pl_unknown.txt", "a");
         
         fprintf(out3, "%s", buffer+2);
      }

      // Remove k from bitmap
      pos = strchr(buffer, '*');
      if (!pos)
         error("Missing multiplier");

      *pos = 0;
      sscanf(buffer+2, "%" SCNu64"", &k);

      gb_KTerms[BIT(k)] = false;
   }
   
   fclose(in);
   fclose(out2);
   
   if (out3 != NULL)
      fclose(out3);

   gl_remainingK = 0;
   out1 = fopen("pl_remain.tmp", "w");

   for (uint64_t k=gl_minK; k<=gl_maxK ; k++)
   {
      if (gb_KTerms[BIT(k)])
      {
         fprintf(out1, "%" PRIu64"*%u^n%+d\n", k, gi_base, gi_c);
         gl_remainingK++;
      }
   }

   fclose(out1);

   delete_file("results.tmp");
   delete_file("pl_remain.txt");
   
   if (gl_remainingK > 0)
      rename("pl_remain.tmp", "pl_remain.txt");
   else
      delete_file("pl_remain.tmp");
}

void output_remain(const char *message)
{
   FILE *remain = fopen("pl_remain.txt", "w");

   for (uint64_t k=gl_minK; k<=gl_maxK; k++)
   {
      if (gb_KTerms[BIT(k)])
         fprintf(remain, "%" PRIu64"*%u^n%+d\n", k, gi_base, gi_c);
   }

   fclose(remain);

   reportWithDateTime("%s: %u %s remaining", message, gl_remainingK, (gl_remainingK == 1 ? "k" : "k's"));
}

void delete_temp_files()
{
   char  command[200];
   char fileName[100];

   sprintf(fileName, "sr_%c%u.pfgw", (gi_c == 1 ? 'S' : 'R'), gi_base);

   delete_file("sieve.in");
   delete_file(fileName);
   delete_file("pfgw.ini");
   delete_file("pfgw-prime.log");
   delete_file("pfgw.log");

#ifdef WIN32
   sprintf(command, "if exist alg_* del alg_*");
#else
   sprintf(command, "rm alg_* 2>/dev/null");
#endif

   system(command);
}

void prepare_recovery(uint32_t currentStep)
{
   FILE    *ckpt = fopen("srbsieve.ckpt_tmp", "w");

   fprintf(ckpt, "phaseInProgress=%d\n", gi_currentPhase);
   fprintf(ckpt, "subphaseInProgress=%d\n", gi_currentSubPhase);
   fprintf(ckpt, "currentStep=%u\n", currentStep);
   fprintf(ckpt, "phase=%d,%d\n", gs_phases[gi_currentPhase].minN, gs_phases[gi_currentPhase].maxN);

   fclose(ckpt);

   delete_file("srbsieve.ckpt");
   rename("srbsieve.ckpt_tmp", "srbsieve.ckpt");

   gi_recoveryPhase = NO_PHASE;
}

void delete_file(const char *fileName)
{   
   FILE *fPtr = fopen(fileName, "r");
   if (!fPtr) return;
   fclose(fPtr);

   remove(fileName);
}

void checkForProgram(const char *programName)
{
   FILE *fPtr = fopen(programName, "r");
   
   if (fPtr)
   {
      fclose(fPtr);
      return;
   }
   
   error("Program %s was not found", programName);
}

void verifyFbncsieveRanToCompletion(uint32_t n)
{
   char buffer[1000];
   bool runCompleted = false;
   
   if (gb_quitting)
      error("Ending fbncphase with n=%u due to ^C", n);
   
   FILE *fPtr = fopen("fbncsieve.log", "r");
   
   if (!fPtr)
   {
      error("Ending fbncphase with n=%u due to missing fbncsieve.log", n);
      return;
   }

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (strlen(buffer) < 5)
         continue;
      
      if (strstr(buffer, "Sieve completed") != NULL)
      {
         runCompleted = true;
         break;
      }
   }

   fclose(fPtr);
   
   delete_file("fbncsieve.log");
   
   if (!runCompleted)
      error("Ending fbncphase with n=%u as fbncsieve did not complete normally", n);
}

void verifySrsieve2RanToCompletion()
{
   char buffer[1000];
   bool runCompleted = false;
   
   if (gb_quitting)
      error("Ending sieving in phase %u due to ^C", gi_currentPhase);
   
   FILE *fPtr = fopen("srsieve2.log", "r");
   
   if (!fPtr)
   {
      error("Ending sieving in phase %u due to missing srsieve2.log", gi_currentPhase);
      return;
   }

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (strlen(buffer) < 5)
         continue;
      
      if (strstr(buffer, "Sieve interrupted") != NULL)
      {
         runCompleted = true;
         break;
      }
   }

   fclose(fPtr);
   
   delete_file("fbncsieve.log");
   
   if (!runCompleted)
      error("Ending sieving in phase %u as srsieve2 did not complete normally", gi_currentPhase);
}

void verifyPfgwRanToCompletion(bool prpTest)
{
   char buffer[1000];
   bool runCompleted = false;
   
   if (gb_quitting)
      error("Ending %s testing in phase %u due to ^C", (prpTest ? "PRP" : "primality"), gi_currentPhase);
   
   FILE *fPtr = fopen("pfgw.ini", "r");
   
   if (!fPtr)
   {
      error("Ending %s testing in phase %u due to missing pfgw.ini", (prpTest ? "PRP" : "primality"), gi_currentPhase);
      return;
   }

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (strlen(buffer) < 5)
         continue;
      
      if (strstr(buffer, "CurFileProcessing=false") != NULL)
      {
         runCompleted = true;
         break;
      }
   }

   fclose(fPtr);
   
   delete_file("pfgw.ini");
   
   if (!runCompleted)
      error("Ending %s testing in phase %u as pfgw did not complete normally", (prpTest ? "PRP" : "primality"), gi_currentPhase);
}
