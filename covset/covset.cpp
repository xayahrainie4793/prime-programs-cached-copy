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

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "primegen.h"

int32_t     gi_base = -999;
int64_t     gi_minK = -999;
int64_t     gi_maxK = -999;
int32_t     gi_c = -999;
int32_t     gi_period = -999;
int32_t     gi_maxFact = -999;
int32_t     gi_facDisp = -999;
int64_t     gi_kDisp = -999;
bool        gb_quitting = false;

#ifdef WIN32
#define atoll  _atoi64
#define PRId64 "I64d"
#endif

void setQuitting(int sig);
void getInputs(void);
void findCovset(void);
void delete_file(const char *fileName);

#ifdef _WIN32
#include <windows.h>

// Some old MinGW/CYGWIN distributions don't define this:
#ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
#define ENABLE_VIRTUAL_TERMINAL_PROCESSING  0x0004
#endif

void setQuitting(int sig)
{
   printf("\nCtrl-C accepted. ");
   gb_quitting = true;
}

static HANDLE stdoutHandle;
static DWORD outModeInit;

void setupConsole(void) {
	DWORD outMode = 0;
	stdoutHandle = GetStdHandle(STD_OUTPUT_HANDLE);

	if(stdoutHandle == INVALID_HANDLE_VALUE) {
		exit(GetLastError());
	}
	
	if(!GetConsoleMode(stdoutHandle, &outMode)) {
		exit(GetLastError());
	}

	outModeInit = outMode;
	
    // Enable ANSI escape codes
	outMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;

	if(!SetConsoleMode(stdoutHandle, outMode)) {
		exit(GetLastError());
	}	
}

void restoreConsole(void) {
    // Reset colors
    printf("\x1b[0m");	
	
    // Reset console mode
	if(!SetConsoleMode(stdoutHandle, outModeInit)) {
		exit(GetLastError());
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
   FILE *fPtr = fopen("covset.log", "a+");

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
   fprintf(fPtr, "\n");

   fclose(fPtr);
}

int gcd(int num1, int num2)
{
    if (num2 == 0)
       return num1;
    else
       return gcd(num2, num1 % num2);
}
 
int lcm_of_array(std::vector<int> arr)
{
    int lcm = arr[0];
    for (int i = 1; i < 2; i++)
    {
        int num1 = lcm;
        int num2 = arr[i];
        int gcd_val = gcd(num1, num2);
        lcm = (lcm * arr[i]) / gcd_val;
    }
    return lcm;
}

int main()
{
   getInputs();
   
#ifdef SIGHUP
   signal(SIGHUP, SIG_IGN);
   signal(SIGQUIT, setQuitting);
#endif

   signal(SIGINT, setQuitting);
   signal(SIGTERM, setQuitting);
   delete_file("covering.txt");

   printf("\n");
   reportWithDateTime("Covering sets for %c%d k=%" PRId64" to %" PRId64" with factor<=%d and period=%d",
   (gi_c == 1 ? 'S' : 'R'), gi_base, gi_minK, gi_maxK, gi_maxFact, gi_period);
   printf("\n");

   findCovset();
   exit(0);
}

void  getInputs(void)
{
   char  buffer[100], *pos;
   FILE *ini = fopen("covset.ini", "r");

   if (!ini)
      error("covset.ini not found");

   while (fgets(buffer, 100, ini) != NULL)
   {
      if (!memcmp(buffer, "base=", 5))
         gi_base = atoi(buffer+5);
      if (!memcmp(buffer, "mink=", 5))
         gi_minK = atoll(buffer+5);
      if (!memcmp(buffer, "maxk=", 5))
         gi_maxK = atoll(buffer+5);
      if (!memcmp(buffer, "c=", 2))
         gi_c = atoi(buffer+2);
      if (!memcmp(buffer, "period=", 7))
         gi_period = atoll(buffer+7);
      if (!memcmp(buffer, "maxfact=", 8))
         gi_maxFact = atoll(buffer+8);
      if (!memcmp(buffer, "facdisp=", 8))
         gi_facDisp = atoll(buffer+8);
      if (!memcmp(buffer, "kdisp=", 6))
         gi_kDisp = atoi(buffer+6);
   }

   fclose(ini);

   if (gi_base == -999)
      error("base must be specified in covset.ini");
   if (gi_minK == -999)
      error("mink must be specified in covset.ini");
   if (gi_maxK == -999)
      error("maxk must be specified in covset.ini");
   if (gi_c == -999)
      error("c must be specified in covset.ini");

   if (gi_period == -999 || gi_period == 0)
      gi_period = 144;

   if (gi_maxFact == -999 || gi_maxFact == 0)
      gi_maxFact = 1000000;

   if (gi_facDisp == -999)
      gi_facDisp = 0;

   if (gi_kDisp == -999)
      gi_kDisp = 0;

   if (gi_base < 2)
      error("base cannot be less than 2");

   if (gi_base > 1000000000)
      error("base cannot be greater than 1e9");

   if (gi_minK > gi_maxK)
      error("mink cannot be greater than maxk");

   if (gi_minK < 1)
      error("mink cannot be less than 1");

   if (gi_maxK > 1000000000000000000)
      error("maxk cannot be greater than 1e18");

   if (gi_c != 1 && gi_c != -1)
      error("c must be -1 or 1");

   if (gi_period < 2)
      error("period cannot be less than 2");

   if (gi_period > 5040)
      error("period cannot be greater than 5040");

   if (gi_maxFact < 5)
      error("maxfact cannot be less than 5");

   if (gi_maxFact > 1000000000)
      error("maxfact cannot be greater than 1e9");

   if (gi_facDisp < 0)
      error("facdisp cannot be negative");

   if (gi_facDisp > 1000000000)
      error("facdisp cannot be greater than 1e9");

   if (gi_kDisp < 0)
      error("kdisp cannot be negative");

   if (gi_kDisp > 1000000000)
      error("kdisp cannot be greater than 1e9");
}

void  findCovset(void)
{
   uint32_t estPrimes = gi_maxFact / log(gi_maxFact) * 1.15 + 10;
   primegen pg;
   uint32_t *factors;
   int32_t  left = gi_base - 1;
   int32_t  p;
   int32_t  base_facs[10];
   int32_t  base_num_facs = 0;
   uint32_t covsets = 0;
   uint32_t n, n_cntr, fac_cntr;
   uint32_t fac_tbl_cntr = 0;
   uint32_t period_chk = gi_period * 2;
   uint32_t fac_tbl[200][2];
   uint32_t fac_per_tbl[200][2];
   uint32_t cov_set_fac[30][2];
   int64_t  kbnc, base_mod;
   int64_t  k, k_mod;

   FILE    *covering = NULL;
   FILE    *fPtr = fopen("covset.log", "a");

   std::vector<std::vector<uint32_t> > n_tbl_facs;
   n_tbl_facs.resize(period_chk + 1);
   for (int i = 0; i <= period_chk; ++i)
      n_tbl_facs[i].resize(10);

   factors = (uint32_t *) malloc(estPrimes * sizeof(uint32_t));

   primegen_init(&pg);
   factors[0] = primegen_next(&pg);
   for (int pi=1; factors[pi-1] <= gi_maxFact; pi++)
      factors[pi] = primegen_next(&pg);

   // find all prime divisors of the base - 1
   primegen_init(&pg);
   do
   {
      p = primegen_next(&pg);

      if (left % p == 0)
      {
         base_facs[base_num_facs++] = p;
         while (left % p == 0 && left > 0)
            left /= p;
      }
   } while (left > p);

   for (fac_cntr=0; factors[fac_cntr] < gi_maxFact; fac_cntr++)
   {
      if (gi_facDisp > 0)
      {
         // display status of factor load
         if (factors[fac_cntr + 1] / gi_facDisp > factors[fac_cntr] / gi_facDisp &&
         gi_maxFact / gi_facDisp > factors[fac_cntr] / gi_facDisp)
         {
            printf("Loading at factor = %u.", factors[fac_cntr + 1] / gi_facDisp * gi_facDisp);

            for (int i=1; i<=35; i++)
               printf("\b");
         }
      }

      if (gb_quitting)
      {
         printf("Ending factor load at factor=%u.\n", factors[fac_cntr]);
         return;
      }

      // skip factors where the base or (base - 1) is divisible by the factor
      if (gi_base % factors[fac_cntr] == 0 || (gi_base - 1) % factors[fac_cntr] == 0)
         continue;

      // loop thru stored prime factors to find first n that
      // contains each factor; n is the period of the factor
      base_mod = 1;
      for (n=1; n <= gi_period; n++)
      {
         base_mod = (base_mod * gi_base) % factors[fac_cntr];
         kbnc = base_mod - 1 + factors[fac_cntr];

         if (kbnc % factors[fac_cntr] == 0)
         {
            // place all factors in a table where the specified
            // period is divisible by the factor period
            if (gi_period % n == 0)
            {
               fac_tbl[fac_tbl_cntr][0] = factors[fac_cntr];
               fac_tbl[fac_tbl_cntr][1] = n;
               fac_tbl_cntr++;
            }
            break;
         }
      }
   }

   free(factors);

   if (fac_tbl_cntr > 0)
   {
      fprintf(fPtr, "factor and period list:\n");

      for (fac_cntr=0; fac_cntr < fac_tbl_cntr; fac_cntr++)
         fprintf(fPtr, "%u  %u\n", fac_tbl[fac_cntr][0], fac_tbl[fac_cntr][1]);

      fprintf(fPtr, "\n");
   }

   fclose(fPtr);

   if (fac_tbl_cntr < 2)
   {
      reportWithDateTime("0 covering sets found for %c%d", (gi_c == 1 ? 'S' : 'R'), gi_base);
      fprintf(fPtr, "\n");
      return;
   }

   uint32_t fac_per_cntr = 0;

   // build factor period table
   for (n_cntr=2; n_cntr <= gi_period; n_cntr++)
   {
      for (fac_cntr=0; fac_cntr < fac_tbl_cntr; fac_cntr++)
      {
         if (fac_tbl[fac_cntr][1] == n_cntr)
         {
            fac_per_tbl[fac_per_cntr][0] = fac_tbl[fac_cntr][0];
            fac_per_tbl[fac_per_cntr][1] = fac_tbl[fac_cntr][1];
            fac_per_cntr++;
         }
      }
   }

   for (k=gi_minK; k <= gi_maxK; k++)
   {
      if (gb_quitting)
      {
         printf("Ending k search at k=%" PRId64".\n", k);
         return;
      }

      if (gi_kDisp > 0)
      {
         // display status of k search
         if (k % gi_kDisp == 0)
         {
            printf("Searching at k = %" PRId64". %u covering set%s found. ", k, covsets, (covsets == 1 ? "" : "s"));

            for (int i=1; i<=65; i++)
               printf("\b");
         }
      }

      // skip k's that are a multiple of the base
      if (k % gi_base == 0)
         continue;

      // skip k's where (k + c) contains a factor of (base - 1)
      for (fac_cntr=0; fac_cntr < base_num_facs; fac_cntr++)
      {
         if ((k + gi_c) % base_facs[fac_cntr] == 0)
            break;
      }
      if (fac_cntr < base_num_facs)
         continue;

      for (n_cntr=1; n_cntr <= period_chk; n_cntr++)
         for (int i=0; i < 10; i++)
            n_tbl_facs[n_cntr][i] = 0;

      // loop thru factor period table to mark all n for the k that contain a factor
      for (fac_cntr=0; fac_cntr < fac_tbl_cntr; fac_cntr++)
      {
         k_mod = k % fac_per_tbl[fac_cntr][0];
         base_mod = 1;
         for (n=1; n <= period_chk; n++)
         {
            base_mod = (base_mod * gi_base) % fac_per_tbl[fac_cntr][0];
            kbnc = k_mod * base_mod + gi_c + fac_per_tbl[fac_cntr][0];

            if (kbnc % fac_per_tbl[fac_cntr][0] == 0)
            {
               int i;
               for (n_cntr=n; n_cntr <= period_chk; n_cntr += fac_per_tbl[fac_cntr][1])
               {
                  // find first available place in n's for factor
                  // and mark applicable n's as covered by the factor
                  for (i=0; n_tbl_facs[n_cntr][i] != 0; i++)
                     ;
                  n_tbl_facs[n_cntr][i] = fac_per_tbl[fac_cntr][0];
               }
               break;
            }
         }

         if (fac_cntr < 2)
            n_cntr = 1;
         else
         {
            // check if all n's are covered
            for (n_cntr=1; n_cntr <= period_chk; n_cntr++)
            {
               if (n_tbl_facs[n_cntr][0] == 0)
                  break;
            }
            if (n_cntr > period_chk)
               break;
         }
      }

      if (n_cntr > period_chk)
      {
         // all n's were covered
         int cov_set_cntr = 0;
         int fac_tbl_per = 0;
         int cov_set_per = 1;
         int cv_fac;

         // build covering set table
         for (fac_cntr=0; fac_cntr < fac_tbl_cntr; fac_cntr++)
         {
            for (n_cntr=1; n_cntr <= period_chk; n_cntr++)
            {
               if (fac_tbl[fac_cntr][0] == n_tbl_facs[n_cntr][0])
               {
                  cov_set_fac[cov_set_cntr][0] = fac_tbl[fac_cntr][0];
                  cov_set_fac[cov_set_cntr][1] = fac_tbl[fac_cntr][1];
                  cov_set_cntr++;

                  if (fac_tbl[fac_cntr][1] > fac_tbl_per)
                     fac_tbl_per = fac_tbl[fac_cntr][1];

                  // determine period by finding the least common multiple of all factor periods
                  cv_fac = fac_tbl[fac_cntr][1];
                  std::vector<int> arr = { cov_set_per, cv_fac };
                  cov_set_per = lcm_of_array(arr);
                  break;
               }
            }
         }

         // check if all factors are necessary by dropping factors with periods
         // that do not divide into the maximum period of all factors and then
         // check if all n's are still covered
         if (cov_set_per > fac_tbl_per)
         {
            for (fac_cntr=0; fac_cntr < cov_set_cntr; fac_cntr++)
            {
               if (fac_tbl_per % cov_set_fac[fac_cntr][1] != 0)
               {
                  // zero out covering set period for later use if factor is dropped
                  cov_set_fac[fac_cntr][1] = 0;
                  for (n_cntr=1; n_cntr <= period_chk; n_cntr++)
                  {
                     for (int i=0; i < 10; i++)
                     {
                        // zero out factor with disallowed period in n table
                        if (n_tbl_facs[n_cntr][i] == cov_set_fac[fac_cntr][0])
                           n_tbl_facs[n_cntr][i] = 0;

                        if (n_tbl_facs[n_cntr][i] == 0)
                           break;
                     }
                  }
               }
            }

            // check if all n's are still covered by other factors
            for (n_cntr=1; n_cntr <= period_chk; n_cntr++)
            {
               int i;
               for (i=0; i < 10; i++)
               {
                  if (n_tbl_facs[n_cntr][i] != 0)
                     break;
               }
               if (i >= 10)
                  // n is not covered
                  break;
            }

            if (n_cntr > period_chk)
            {
               // all n's are all still covered so update the period
               // and zero out unneeded factors in covering set table
               cov_set_per = fac_tbl_per;

               for (fac_cntr=0; fac_cntr < cov_set_cntr; fac_cntr++)
               {
                  if (cov_set_fac[fac_cntr][1] == 0)
                     cov_set_fac[fac_cntr][0] = 0;
               }
            }
         }

         // write k, period, and covering set to the covering file
         covering = fopen("covering.txt", "a");

         fprintf(covering, "%" PRId64"  %u  [", k, cov_set_per);

         for (fac_cntr=0; cov_set_fac[fac_cntr][0] == 0; fac_cntr++)
            ;

         fprintf(covering, "%u", cov_set_fac[fac_cntr][0]);

         for (fac_cntr=fac_cntr+1; fac_cntr < cov_set_cntr; fac_cntr++)
            if (cov_set_fac[fac_cntr][0] > 0)
               fprintf(covering, ", %u", cov_set_fac[fac_cntr][0]);

         fprintf(covering, "]\n");

         fclose(covering);
         covsets++;
      }
   }

   reportWithDateTime("%u covering set%s found for %c%d", covsets, (covsets == 1 ? "" : "s"), (gi_c == 1 ? 'S' : 'R'), gi_base);
   fprintf(fPtr, "\n");
}

void delete_file(const char *fileName)
{   
   FILE *fPtr = fopen(fileName, "r");

   if (!fPtr)
      return;

   fclose(fPtr);
   remove(fileName);
}
