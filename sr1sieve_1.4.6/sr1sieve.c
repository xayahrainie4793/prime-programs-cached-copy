/* sr1sieve.c -- (C) Geoffrey Reynolds, June 2006.

   A version of srsieve specialised for a single sequence k*b^n+/-1.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <errno.h>
#include "sr1sieve.h"
#include "version.h"
#include "arithmetic.h"

#if HAVE_MALLOPT
#include <malloc.h>
#endif

#define NAME "sr1sieve"
#define DESC "A sieve for one sequence k*b^n+/-1"

const char *output_file_name = NULL;
const char *factors_file_name = NULL;
uint32_t save_period = DEFAULT_SAVE_PERIOD;
int verbose_opt = 0;
int quiet_opt = 0;
uint32_t L1_cache_size = 0;
uint32_t L2_cache_size = 0;
const char *baby_method_opt = NULL;
const char *giant_method_opt = NULL;
#if HASHTABLE_SIZE_OPT
uint32_t hashtable_size = 0;
#endif

static void banner(void)
{
  printf(NAME " " XSTR(MAJOR_VER) "." XSTR(MINOR_VER) "." XSTR(PATCH_VER)
         " -- " DESC ".\n");
#ifdef __GNUC__
  if (verbose_opt)
    printf("Compiled on " __DATE__ " with GCC " __VERSION__ ".\n");
#endif
}

static void help(void)
{
  printf("Usage: %s -P P1 -i FILE <-o FILE | -f FILE> [OPTION ...]\n", NAME);
  printf(" -p --pmin P0\n");
  printf(" -P --pmax P1          "
         "Sieve for factors p in the range P0 <= p <= P1\n");
  printf(" -i --input FILE       "
         "Read sieve from input file FILE (can be ABC, ABCD, or NEWPGEN format).\n");
  printf(" -o --output FILE      "
         "Write sieve to NewPGen format file FILE.\n");
  printf(" -f --factors FILE     "
         "Append new factors to file FILE.\n");
#if LEGENDRE_CACHE
  printf(" -C --cache-file FILE  "
         "Load (or save) Legendre symbol tables from (or to) FILE.\n");
#endif
  printf(" -s --save TIME        "
         "Update output file every TIME (default %u) minutes.\n",
         (unsigned int)DEFAULT_SAVE_PERIOD/60);
#if MULTI_PATH
#ifdef __i386__
  printf("    --amd              "
         "Use CMOV optimisations. (Default for AMDs if supported).\n");
  printf("    --intel            "
         "Don't use CMOV optimisations. (Default for Intels).\n");
#endif
  printf("    --sse2             "
         "Use SSE2 vector optimisations. (Default if supported).\n");
  printf("    --no-sse2          "
         "Don't use SSE2 vector optimisations.\n");
#endif
#if HAVE_FORK
  printf(" -t --threads NUM      "
         "Start NUM child threads. (Default 0).\n");
#endif
  printf(" -l --L1-cache SIZE    "
         "Assume L1 data cache is SIZE Kb.\n");
  printf(" -L --L2-cache SIZE    "
         "Assume L2 cache is SIZE Kb.\n");
  printf(" -B --baby METHOD      "
         "Use METHOD for baby step mulmods.\n");
  printf(" -G --giant METHOD     "
         "Use METHOD for giant step mulmods.\n");
#if HASHTABLE_SIZE_OPT
  printf(" -H --hashtable SIZE   "
         "Force use of a SIZE Kb hashtable.\n");
#endif
#if SUBSEQ_Q_OPT
  printf(" -Q --subseq Q         "
         "Force sieving k*b^n+c as subsequences (k*b^d)*(b^Q)^m+c.\n");
#endif
  printf(" -F --format=f         "
         "Format of output file (A=ABC, D=ABCD (default), N=NEWPGEN)");
  printf(" -x --no-lookup        "
         "Don't pre-compute Legendre symbol lookup tables.\n");
  printf(" -z --lower-priority   "
         "Run at low priority. (-zz lower).\n");
  printf(" -Z --raise-priority   "
         "Run at high priority. (-ZZ higher).\n");
  printf(" -A --affinity N       "
         "Set affinity to CPU number N.\n");
  printf(" -d --duplicates       "
         "Report factors that don't eliminate any composite.\n");
  printf(" -q --quiet            "
         "Don't print found factors.\n");
  printf(" -v --verbose          "
         "Print some extra messages. -vv prints more.\n");
  printf(" -h --help             "
         "Print this help.\n\n");

  exit(EXIT_SUCCESS);
}

static const char *short_opts =
#if LEGENDRE_CACHE
  "C:"
#endif
#if HAVE_FORK
  "t:"
#endif
#if HASHTABLE_SIZE_OPT
  "H:"
#endif
#if SUBSEQ_Q_OPT
  "Q:"
#endif
  "p:P:i:o:f:s:l:L:B:G:F:xzZA:dqvh";

static const struct option long_opts[] =
  {
    {"pmin",       required_argument, 0, 'p'},
    {"pmax",       required_argument, 0, 'P'},
    {"input",      required_argument, 0, 'i'},
    {"output",     required_argument, 0, 'o'},
    {"factors",    required_argument, 0, 'f'},
#if LEGENDRE_CACHE
    {"cache-file", required_argument, 0, 'C'},
#endif
    {"save",       required_argument, 0, 's'},
#if MULTI_PATH
    {"no-sse2",    no_argument,       0, '1'},
    {"sse2",       no_argument,       0, '2'},
#ifdef __i386__
    {"amd",        no_argument,       0, '3'},
    {"intel",      no_argument,       0, '4'},
#endif
#endif /* MULTI_PATH */
#if HAVE_FORK
    {"threads",    required_argument, 0, 't'},
#endif
    {"L1-cache",   required_argument, 0, 'l'},
    {"L2-cache",   required_argument, 0, 'L'},
    {"baby",       required_argument, 0, 'B'},
    {"giant",      required_argument, 0, 'G'},
#if HASHTABLE_SIZE_OPT
    {"hashtable",  required_argument, 0, 'H'},
#endif
#if SUBSEQ_Q_OPT
    {"subseq",     required_argument, 0, 'Q'},
#endif
    {"format",     required_argument, 0, 'F'},
    {"no-lookup",  no_argument,       0, 'x'},
    {"lower-priority", no_argument,   0, 'z'},
    {"raise-priority", no_argument,   0, 'Z'},
    {"affinity",   required_argument, 0, 'A'},
    {"duplicates", no_argument,       0, 'd'},
    {"quiet",      no_argument,       0, 'q'},
    {"verbose",    no_argument,       0, 'v'},
    {"help",       no_argument,       0, 'h'},
    {0, 0, 0, 0}
  };


static time_t start_date;
static int opt_ind, opt_c;

static void attribute ((noreturn)) argerror(const char *str)
{
  if (long_opts[opt_ind].val == opt_c)
    error("--%s %s: argument %s.", long_opts[opt_ind].name, optarg, str);
  else
    error("-%c %s: argument %s.", opt_c, optarg, str);
}

static uint64_t parse_uint(uint64_t limit)
{
  uint64_t num;
  uint32_t expt;
  char *tail;

  errno = 0;
  num = strtoull(optarg,&tail,0);

  if (errno == 0 && num <= limit)
  {
    switch (*tail)
    {
      case 'e':
      case 'E':
        expt = strtoul(tail+1,&tail,0);
        if (errno != 0)
          goto range_error;
        if (*tail != '\0')
          break;
        while (expt-- > 0)
          if (num > limit/10)
            goto range_error;
          else
            num *= 10;
      case '\0':
        return num;
    }

    argerror("is malformed");
  }

 range_error:
  argerror("is out of range");
}

int main(int argc, char **argv)
{
  const char *input_file_name = NULL;
  int priority_opt =  -2;
#if MULTI_PATH
  int code_path = 0;
  int sse2_opt = 0;
#ifdef __i386__
  int cmov_opt = 0;
#endif
#endif /* MULTI_PATH */
  uint32_t i, L1_opt = 0, L2_opt = 0;
  int want_help = 0;

  set_accumulated_time(0.0);

#if HAVE_MALLOPT
  /* All memory allocation takes place during initialization, so it doesn't
     have to be fast. Reducing this threshold allows more memory to be
     returned to the system when free() is called.
  */
  mallopt(M_MMAP_THRESHOLD,16000);
#endif

#if USE_COMMAND_LINE_FILE
  /* If no comand line arguments are given, attempt to read them from file.
   */
  if (argc == 1)
    read_argc_argv(&argc,&argv,NAME COMMAND_LINE_FILE_NAME_SUFFIX);
#endif

  while ((opt_c = getopt_long(argc,argv,short_opts,long_opts,&opt_ind)) != -1)
    switch (opt_c)
    {
      case 'p':
        p_min = parse_uint(UINT64_MAX);
        break;
      case 'P':
        p_max = parse_uint(UINT64_MAX);
        break;
      case 'i':
        input_file_name = optarg;
        break;
      case 'o':
        output_file_name = optarg;
        break;
      case 'f':
        factors_file_name = optarg;
        break;
#if LEGENDRE_CACHE
      case 'C':
        cache_file_name = optarg;
        break;
#endif
      case 's':
        save_period = strtoul(optarg,NULL,0) * 60;
        break;
#if MULTI_PATH
      case '1':
        sse2_opt = -1; /* --no-sse2 */
        break;
      case '2':
        sse2_opt = 1; /* --sse2 */
        break;
#ifdef __i386__
      case '3':
        cmov_opt = 1; /* --amd */
        break;
      case '4':
        cmov_opt = -1; /* --intel */
        break;
#endif
#endif /* MULTI_PATH */
#if HAVE_FORK
      case 't':
        num_children = strtol(optarg,NULL,0);
        if (num_children < 0 || num_children > MAX_CHILDREN)
          argerror("out of range");
        break;
#endif
      case 'l':
        L1_opt = strtoul(optarg,NULL,0);
        break;
      case 'L':
        L2_opt = strtoul(optarg,NULL,0);
        break;
      case 'B':
        baby_method_opt = optarg;
        break;
      case 'G':
        giant_method_opt = optarg;
        break;
#if HASHTABLE_SIZE_OPT
      case 'H':
        i = strtoul(optarg,NULL,0) * 1024;
        for (hashtable_size = 1024; hashtable_size < i; )
          hashtable_size <<= 1;
        break;
#endif
#if SUBSEQ_Q_OPT
      case 'Q':
        subseq_Q = strtoul(optarg,NULL,0);
        break;
#endif
      case 'x':
        no_lookup_opt = 1;
        break;
      case 'z':
        priority_opt--;
        break;
      case 'Z':
        priority_opt++;
        break;
      case 'A':
#if HAVE_FORK
        if (affinity_opt < MAX_CHILDREN)
          child_affinity[affinity_opt++] = strtol(optarg,NULL,0);
#else
        set_cpu_affinity(strtol(optarg,NULL,0));
#endif
        break;
      case 'd':
        duplicates_opt = 1;
        break;
      case 'q':
        quiet_opt++;
        break;
      case 'v':
        verbose_opt++;
        break;
      case 'h':
        want_help = 1;
        break;
      case 'F':
         if (optarg[0] == 'A')
            outputFileFormat = FF_ABC;
         if (optarg[0] == 'D')
            outputFileFormat = FF_ABCD;
         if (optarg[0] == 'N')
            outputFileFormat = FF_NEWPGEN;
         break;
      default:
        return 1;
    }

  banner();

  if (want_help || argc < 2)
    help();


  /* We install these signal handlers early because the default handler
     terminates the program. SIGTERM and SIGINT handlers are not installed
     until sieving begins, but termination is the right default for them.
  */
#ifdef SIGUSR1
  signal(SIGUSR1,handle_signal);
#endif
#ifdef SIGUSR2
  signal(SIGUSR2,handle_signal);
#endif

#if SUBSEQ_Q_OPT
  if (subseq_Q % BASE_MULTIPLE != 0)
    error("Subsequence base exponent Q=%"PRIu32" is not a multiple of %u.",
          subseq_Q, (unsigned int)BASE_MULTIPLE);
  if (subseq_Q != 0 && LIMIT_BASE % subseq_Q != 0)
    error("Subsequence base exponent Q=%"PRIu32" is not a divisor of %u.",
          subseq_Q, (unsigned int)LIMIT_BASE);
#endif


#if MULTI_PATH
  /* Choose code path based on switches --amd, --intel, --sse2, --no-sse2.
   */
#if __i386__
  if (sse2_opt < 0)
  {
    if (cmov_opt > 0)
      code_path = 3;
    else if (cmov_opt < 0)
      code_path = 1;
  }
  else if (sse2_opt > 0)
    code_path = 2;
#elif __x86_64__
  if (sse2_opt < 0)
    code_path = 1;  /* Force use of x86 FPU */
#endif
  /* If code path was not fully determined above then test hardware features.
   */
  if (code_path == 0)
  {
#if __x86_64__
    /* x86-64 always has SSE2. Choose default based on p_max instead. */
    code_path = (p_max < (UINT64_C(1) << 51)) ? 2 : 1;
#elif __i386__
    if (sse2_opt >= 0 && have_sse2()) /* Prefer SSE2 if available. */
      code_path = 2;
    else if (have_cmov() && is_amd()) /* Prefer CMOV on AMD CPUs */
      code_path = 3;
    else
      code_path = 1;
#endif
  }
  /* Report code path selection.
   */
  if (verbose_opt)
  {
    if(code_path == 2)
      printf("SSE2 code path, ");
#if __x86_64__
    else
      printf("x87 code path, ");
#elif __i386__
    else if (code_path == 3)
      printf("AMD code path, ");
    else
      printf("Intel code path, ");
#endif
  }
#endif /* MULTI_PATH */

  set_cache_sizes(L1_opt,L2_opt);

#if LEGENDRE_CACHE
  if (no_lookup_opt != 0 && cache_file_name != NULL)
    error("--no-lookup and --cache switches are mutually exclusive.");
#endif

  if (priority_opt)
    set_process_priority(priority_opt);

  if (optind < argc)
    error("Unhandled argument %s.", argv[optind]);
  if (p_max == 0)
    error("--pmax is a required argument.");
  if (input_file_name == NULL)
    error("--input is a required argument.");
  if (output_file_name == NULL && factors_file_name == NULL)
    error("At least one of --output or --factors is a required argument.");

  read_input_file(input_file_name);

  if (p_min <= b_term || p_min>=p_max || p_max >= UINT64_C(1)<<MOD64_MAX_BITS)
    error("Sieve range P0 <= p <= P1 must be in %"PRIu32" < P0 < P1 < 2^%u.\n",
          b_term, (unsigned int)MOD64_MAX_BITS);

#if HAVE_FORK && MOD64_MAX_BITS > 62
  /* When multithreading with fork we need to encode p and parity (2 bits)
     into a 64-bit message. */
  if (num_children > 0 && p_max > UINT64_C(1) << 62)
    error("--pmax is limited to 2^62 when used with --threads.");
#endif

  finish_candidate_seqs();

#if MULTI_PATH
  switch (code_path)
  {
    default:
    case 1:
      sieve();
      break;
    case 2:
      sieve_sse2();
      break;
#ifdef __i386__
    case 3:
      sieve_cmov();
      break;
#endif
  }
#else
  sieve();
#endif

  return EXIT_SUCCESS;
}

/* Return the fraction of the current range that has been done.
 */
double frac_done(uint64_t p)
{
  return (double)(p-p_min)/(p_max-p_min);
}

#define STRFTIME_FORMAT "ETA %d %b %H:%M"
void print_status(uint64_t p, uint32_t p_per_sec, uint32_t secs_per_factor)
{
  static int toggle = 0;
  char buf1[32], buf2[32];

  toggle = !toggle; /* toggle between reporting ETA and factor rate. */

  //if (toggle && factor_count && p_per_sec)
  {
    // uint64_t p_per_factor = (p-work_pmin)/factor_count;
    snprintf(buf1,31,"%d sec/factor", secs_per_factor); //p_per_factor/p_per_sec);
  }
  //else
  {
    double done = (double)(p-p_min)/(p_max-p_min);
    buf2[0] = '\0';
    if (done > 0.0) /* Avoid division by zero */
    {
      time_t finish_date = start_date+(time(NULL)-start_date)/done;
      struct tm *finish_tm = localtime(&finish_date);
      if (!finish_tm || !strftime(buf2,sizeof(buf2),STRFTIME_FORMAT,finish_tm))
        buf2[0] = '\0';
    }
  }
  buf2[31] = '\0';

  report(0,"p=%"PRIu64", %"PRIu32" p/sec, %"PRIu32" factor%s, %.1f%% done, %s, %s",
         p,p_per_sec,factor_count,plural(factor_count),100.0*frac_done(p),buf1,buf2);
}

static double expected_factors(uint32_t n, uint64_t p0, uint64_t p1)
{
  /* TODO: Use a more accurate formula. This one is only reasonable when
     p0/p1 is close to 1.
   */

  return n*(1-log(p0)/log(p1));
}

void start_srsieve(void)
{
  if (verbose_opt)
    report(1,"Expecting to find factors for about %.2f terms.",
           expected_factors(ncount,p_min,p_max));

  logger(1,"%s started: %"PRIu32" <= n <= %"PRIu32", %"PRIu64" <= p <= %"PRIu64,
         NAME " " XSTR(MAJOR_VER) "." XSTR(MINOR_VER) "." XSTR(PATCH_VER),
         n_min, n_max, p_min, p_max);

  start_date = time(NULL);
}

void finish_srsieve(const char *reason, uint64_t p)
{
  logger(1,"%s stopped: at p=%"PRIu64" because %s.",
         NAME " " XSTR(MAJOR_VER) "." XSTR(MINOR_VER) "." XSTR(PATCH_VER),
         p, reason);
  write_output_file(1,p,output_file_name);

  logger(verbose_opt,"Found factors for %"PRIu32" term%s in %.3f sec."
         " (expected about %.2f)",factor_count,plural(factor_count),
         get_accumulated_time(),expected_factors(ncount,p_min,p));
}
