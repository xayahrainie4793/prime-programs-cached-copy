#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <vector>

#ifdef WIN32
#include <windows.h>
#else
#include <inttypes.h>
#endif

uint64_t    gl_minK = 0;
uint64_t    gl_maxK = 0;
int32_t     gi_base = 0;
int32_t     gi_c = 0;

std::vector<bool> gb_k;
std::vector<bool> gb_file_k;

#define BIT(k)    ((k) - gl_minK)

#ifdef WIN32
#define atoll  _atoi64

#define PRIu64 "I64u"
#define PRId64 "I64d"
#define PRIx64 "I64x"
#define SCNu64 "I64u"
#define SCNd64 "I64d"
#define SCNx64 "I64x"
#endif

void getInputs(void);
void apply(const char *fileName, bool kOnly, bool havePrime);

uint32_t sortFile(const char *fileName);
uint32_t sortPrimes(void);
uint32_t sortRemain(void);

void error(const char *fmt, ...)
{
   va_list args;

   fprintf(stderr, "Error: ");
   va_start(args, fmt);
   vfprintf(stderr, fmt, args);
   va_end(args);
   fprintf(stderr, "\n");
   exit(0);
}

int main(int argc, char **argv)
{
   getInputs();

   gb_k.resize(gl_maxK - gl_minK + 1);
   gb_file_k.resize(gl_maxK - gl_minK + 1);
   std::fill(gb_k.begin(), gb_k.end(), false);

   uint32_t oddK = 0;
   uint64_t k = gl_minK;

   if (gi_base & 1)
   {
      // For odd bases odd k were already excluded by srbsieve so they
      // won't appear in any of the files.
      
      // Start with odd k
      if (!(k & 1))
         k++;

      for ( ; k<=gl_maxK; k+=2)
      {
         gb_k[BIT(k)] = true;
         oddK++;
      }
   }
   
   apply("pl_prime.txt", false, true);
   apply("pl_remain.txt", false, false);
   apply("pl_trivial.txt", true, false);
   apply("pl_algeb.txt", true, false);
   apply("pl_GFN.txt", true, false);
   apply("pl_MOB.txt", true, false);
   
   FILE *mPtr = NULL;
   uint32_t missingCount = 0;
   for (k=gl_minK; k<=gl_maxK; k++)
   {
      if (!gb_k[BIT(k)])
      {
         if (!mPtr)
            mPtr = fopen("pl_missing.txt", "w");
         
         fprintf(mPtr, "%" PRIu64"*%u^n%+d\n", k, gi_base, gi_c);
         missingCount++;
      }
   }
   
   uint32_t trivials = sortFile("pl_trivial.txt");
   uint32_t algeb = sortFile("pl_algeb.txt");
   uint32_t gfns = sortFile("pl_GFN.txt");
   uint32_t mobs = sortFile("pl_MOB.txt");
   uint32_t primes = sortPrimes();
   uint32_t remain = sortRemain();

   if (mPtr)
   {
      fclose(mPtr);
      error("%u %s missing.  %s in pl_missing.txt", missingCount, (missingCount == 1 ? "k is" : "k's are"), (missingCount == 1 ? "It is" : "They are"));
   }

   printf("Woohoo!  All %" PRIu64" k's are accounted for.\n", gl_maxK - gl_minK + 1);
   
   printf("Final tally for %c%u with k range from %" PRIu64" to %" PRIu64":\n", (gi_c == 1 ? 'S' : 'R'), gi_base, gl_minK, gl_maxK);
   
   if (gi_base & 1)
      printf("  Odd k and base:  %10u\n", oddK);
   
   printf("  Trivials:        %10u\n", trivials);
   printf("  Algebraic:       %10u\n", algeb);
   printf("  GFNs:            %10u\n", gfns);
   printf("  MOBs:            %10u\n", mobs);
   printf("  Primes:          %10u\n", primes);
   printf("  Remaining:       %10u\n", remain);
   
   printf("All results have been sorted by ascending k in their respective files\n");
   
   exit(0);
}

void  getInputs(void)
{
   char  buffer[100], *pos;
   FILE *ini = fopen("srbsieve.ini", "r");
   int   index;

   if (!ini)
      error("srbsieve.ini not found");

   gl_minK = 0;
   gl_maxK = 0;

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
   }

   fclose(ini);

   if (gl_minK >= gl_maxK)
      error("mink must be less than maxk");
}

void  apply(const char *fileName, bool kOnly, bool havePrime)
{
   char buffer[50];
   uint64_t k;
   int32_t  c;
   uint32_t b, n;
   uint32_t primeCount = 0;

   FILE *fPtr = fopen(fileName, "r");

   if (!fPtr)
      return;

   std::fill(gb_file_k.begin(), gb_file_k.end(), false);

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (buffer[0] == '\n')
         continue;

      if (kOnly)
      {
         if (sscanf(buffer, "%" SCNu64"\n", &k) != 1)
            error("Could not parse k from line %s in %s", buffer, fileName);
      }
      else
      {
         if (havePrime)
         {
            if (sscanf(buffer, "%" SCNu64"*%u^%u%d\n", &k, &b, &n, &c) != 4)
               error("Could not parse k/b/n/c from line %s in %s", buffer, fileName);

            if (n == 0)
               error("Prime n=0 is invalid from line %s in %s", buffer, fileName);
         }
         else
         {
            if (sscanf(buffer, "%" SCNu64"*%u^n%d\n", &k, &b, &c) != 3)
               error("Could not parse k/b/c from line %s in %s", buffer, fileName);

         }

         if (b != gi_base)
            error("Found b=%u, but expected b=%u from line %s in %s", b, gi_base, buffer, fileName);

         if (c != gi_c)
            error("Found c=%+d, but expected c=%+d from line %s in %s", c, gi_c, buffer, fileName);
      }

      if (k < gl_minK)
         error("k = %" PRIu64" is less than mink.  Look in file %s", k, fileName);

      if (k > gl_maxK)
         error("k = %" PRIu64" is greater than maxk.  Look in file %s", k, fileName);

      if ((gi_base & 1) && (k & 1))
         error("k = %" PRIu64" is an odd k on an odd base.  Look in file %s", k, fileName);

      // if k has previously been found and was not in the current file, then it is a duplicate.
      if ( gb_k[BIT(k)] && !gb_file_k[BIT(k)] )
         error("k = %" PRIu64" found twice.  Look in file %s", k, fileName);

      gb_k[BIT(k)] = true;
      gb_file_k[BIT(k)] = true;
   }

   fclose(fPtr);
}

uint32_t  sortFile(const char *fileName)
{
   char buffer[50];
   uint64_t k;
   int32_t  c;
   uint32_t b, n;
   uint32_t kCount = 0;
   uint32_t file_kcnt = 0;
   
   FILE *fPtr = fopen(fileName, "r");

   if (!fPtr)
      return 0;

   std::fill(gb_k.begin(), gb_k.end(), false);

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (buffer[0] == '\n')
         continue;
      
      if (sscanf(buffer, "%" SCNu64"\n", &k) != 1)
         error("Could not parse k from line %s in %s", buffer, fileName);

      gb_k[BIT(k)] = true;
      file_kcnt++;
   }

   fclose(fPtr);

   fPtr = fopen(fileName, "w");
   
   for (k=gl_minK; k<=gl_maxK; k++)
   {      
      if (gb_k[BIT(k)])
      {
         kCount++;
         fprintf(fPtr, "%" PRIu64"\n", k);
      }
   }
   
   fclose(fPtr);

   if (file_kcnt > kCount)
      printf("%u duplicate %s removed from %s.\n", file_kcnt - kCount, (file_kcnt - kCount == 1 ? "k" : "k's"), fileName);

   return kCount;
}

uint32_t  sortPrimes(void)
{
   char buffer[50];
   uint64_t k;
   int32_t  c;
   uint32_t b, n;
   uint32_t primeCount = 0;
   uint32_t file_kcnt = 0;

   FILE *fPtr = fopen("pl_prime.txt", "r");

   if (!fPtr)
      return 0;

   uint16_t *ns = (uint16_t *) malloc((gl_maxK - gl_minK + 1) * sizeof(uint16_t));

   memset(ns, 0x00, (gl_maxK - gl_minK + 1) * sizeof(uint16_t));

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (buffer[0] == '\n')
         continue;

      if (sscanf(buffer, "%" SCNu64"*%u^%u%d\n", &k, &b, &n, &c) != 4)
         error("Could not parse k/b/n/c from line %s in pl_prime.txt", buffer);

      if (b != gi_base)
         error("Found b=%u, but expected b=%u from line %s in pl_prime.txt", b, gi_base);

      if (c != gi_c)
         error("Found c=%+d, but expected c=%+d from line %s in pl_prime.txt", c, gi_c);

      // Store n for prime.  If prime previously exists for k, keep only the smallest n.
      if (ns[k - gl_minK] == 0 || ns[k - gl_minK] > n)
         ns[k - gl_minK] = n;

      file_kcnt++;
   }

   fclose(fPtr);

   fPtr = fopen("pl_prime.txt", "w");

   for (k=gl_minK; k<=gl_maxK; k++)
   {
      n = ns[k - gl_minK];
      
      if (n == 0)
         continue;

      primeCount++;
      fprintf(fPtr, "%" PRIu64"*%u^%u%+d\n", k, gi_base, n, gi_c);
   }

   fclose(fPtr);

   free(ns);

   if (file_kcnt > primeCount)
      printf("%u duplicate %s removed from pl_prime.txt.\n", file_kcnt - primeCount, (file_kcnt - primeCount == 1 ? "k" : "k's"));

   return primeCount;
}

uint32_t  sortRemain(void)
{
   const char *fileName = "pl_remain.txt";
   char buffer[50];
   uint64_t k;
   int32_t  c;
   uint32_t b, n;
   uint32_t kCount = 0;
   uint32_t file_kcnt = 0;

   FILE *fPtr = fopen(fileName, "r");

   if (!fPtr)
      return 0;

   std::fill(gb_k.begin(), gb_k.end(), false);

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (buffer[0] == '\n')
         continue;

      if (sscanf(buffer, "%" SCNu64"*%u^n%d\n", &k, &b, &c) != 3)
         error("Could not parse k/b/c from line %s in %s", buffer, fileName);

      if (b != gi_base)
         error("Found b=%u, but expected b=%u from line %s in %s", b, gi_base, fileName);

      if (c != gi_c)
         error("Found c=%+d, but expected c=%+d from line %s in %s", c, gi_c, fileName);

      gb_k[BIT(k)] = true;
      file_kcnt++;
   }

   fclose(fPtr);

   fPtr = fopen(fileName, "w");

   for (k=gl_minK; k<=gl_maxK; k++)
   {      
      if (gb_k[BIT(k)])
      {
         kCount++;
         fprintf(fPtr, "%" PRIu64"*%u^n%+d\n", k, gi_base, gi_c);
      }
   }

   fclose(fPtr);

   if (file_kcnt > kCount)
      printf("%u duplicate %s removed from %s.\n", file_kcnt - kCount, (file_kcnt - kCount == 1 ? "k" : "k's"), fileName);

   return kCount;
}
