/* sr1sieve.h -- (C) Geoffrey Reynolds, April 2006.

   Srsieve specialised for the Riesel Prime Search project.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _SR1SIEVE_H
#define _SR1SIEVE_H

#include <limits.h>
#include <stdint.h>
#include "config.h"

#define XSTR(ARG) STR(ARG)
#define STR(ARG) #ARG

/* Also check config.h for compiler settings.
 */

/* Set DEFAULT_L1_CACHE_SIZE,DEFAULT_L2_CACHE_SIZE to the L1,L2 data cache
   size, in Kb. (If not set in the Makefile).
*/
#ifndef DEFAULT_L1_CACHE_SIZE
#define DEFAULT_L1_CACHE_SIZE 16
#endif
#ifndef DEFAULT_L2_CACHE_SIZE
#define DEFAULT_L2_CACHE_SIZE 256
#endif

/* Set BASE_MULTIPLE to the smallest exponent Q for which sieving in
   subsequence base b^Q will be considered. Must be a multiple of 2.
*/
#define BASE_MULTIPLE 30

/* Allow sieving in base b^Q for Q chosen from the divisors of LIMIT_BASE.
   Must be a multiple of BASE_MULTIPLE.
*/
#define LIMIT_BASE 720

/* For a prime p that satisfies p=1 (mod r), an "r-th power residue test"
   checks whether a subsequence of k*b^n+c can possibly contain any terms of
   the form x^r (mod p). If there are none then that subsequence can be
   omitted from the BSGS step.

   To conduct r-th power residue tests for each r in a set R of prime
   powers, set POWER_RESIDUE_LCM to lcm(R), and set POWER_RESIDUE_DIVISORS
   to the number of divisors of POWER_RESIDUE_LCM. POWER_RESIDUE_LCM must be
   a multiple of BASE_MULTIPLE, a divisor of LIMIT_BASE, and must be less
   than 2^15. E.g.

   R={2,3,4,5}: POWER_RESIDUE_LCM=60, POWER_RESIDUE_DIVISORS=12
   R={2,3,4,5,8}: POWER_RESIDUE_LCM=120, POWER_RESIDUE_DIVISORS=16
   R={2,3,4,5,8,9}: POWER_RESIDUE_LCM=360, POWER_RESIDUE_DIVISORS=24
   R={2,3,4,5,8,9,16}: POWER_RESIDUE_LCM=720, POWER_RESIDUE_DIVISORS=30
   R={2,3,4,5,7,8,9,16}: POWER_RESIDUE_LCM=5040, POWER_RESIDUE_DIVISORS=60

   Memory use is proportional to POWER_RESIDUE_LCM*POWER_RESIDUE_DIVISORS
*/
#define POWER_RESIDUE_LCM 720
#define POWER_RESIDUE_DIVISORS 30

/* For a hash table expected to hold M elements, use a main table of at
   least M/HASH_MAX_DENSITY and at most M/HASH_MIN_DENSITY elements. The
   size will be increased further to minimise density within this range,
   if L1_cache_size is high enough.
*/
#define HASH_MAX_DENSITY 0.60
#define HASH_MIN_DENSITY 0.10

/* Files to write log reports and report factors to.
 */
#define LOG_FILE_NAME "sr1sieve.log"

/* Date format used for log file entries (see date --help)
 */
#define LOG_STRFTIME_FORMAT "%c "

/* Report period in seconds.
 */
#define REPORT_PERIOD 60

/* Save period in seconds.
 */
#define DEFAULT_SAVE_PERIOD 3600

/* Set CHECK_FACTORS=1 to double check found factors.
 */
#define CHECK_FACTORS 1

/* Set HASHTABLE_SIZE_OPT=1 to provide the -H command-line switch, which
   forces use of a specific size hashtable.
*/
#define HASHTABLE_SIZE_OPT 1

/* Set SUBSEQ_Q_OPT=1 to provide the -Q command-line switch, which
   forces use of a specific subsequence base B^Q.
*/
#define SUBSEQ_Q_OPT 1

/* Set LEGENDRE_CACHE=1 to read Legendre lookup tables from a cache file if
   one is found, and to enable the --cache command line switch.
*/
#define LEGENDRE_CACHE 1

/* Set CHECK_FOR_GFN=1 to check whether a sequence consists of Generalised
   Fermat numbers when generating the Legendre symbol tables. A better
   table can be used for these sequences.
*/
#define CHECK_FOR_GFN 1

/* When USE_COMMAND_LINE_FILE=1 and no command line arguments are given,
   read the command line from a file called NAME-command-line.txt, where
   NAME is the name under which the program was invoked. This can be useful
   where the use of the shell and batch files has been disabled.
*/
#define USE_COMMAND_LINE_FILE 1
#define COMMAND_LINE_FILE_NAME_SUFFIX "-command-line.txt"


/* Nothing below here should normally need adjustment. */

#ifdef __GNUC__
/* macros that evaluate their arguments only once. From the GCC manual.
 */
#define MIN(a,b) \
   ({ typeof (a) _a = (a); \
      typeof (b) _b = (b); \
      _a < _b ? _a : _b; })
#define MAX(a,b) \
   ({ typeof (a) _a = (a); \
      typeof (b) _b = (b); \
      _a > _b ? _a : _b; })
#define ABS(a) \
   ({ typeof (a) _a = (a); \
      _a < 0 ? -_a : _a; })
#else
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define ABS(a) ((a) < 0 ? -(a) : (a))
#endif

#if (MULTI_PATH && CODE_PATH==2)
# define CODE_PATH_NAME(name) name##_sse2
#elif (MULTI_PATH && CODE_PATH==3)
# define CODE_PATH_NAME(name) name##_cmov
#else
# define CODE_PATH_NAME(name) name
#endif


/* bsgs.c */

void sieve(void);

#if MULTI_PATH
void sieve_sse2(void);
void sieve_cmov(void);
#endif


/* choose.c */

uint32_t find_best_Q(uint32_t *subseqs);


/* cpu.c */

void set_cache_sizes(uint32_t L1_opt, uint32_t L2_opt);
#if defined(__i386__) || defined(__x86_64__)
int have_sse2(void);
int have_cmov(void);
int is_amd(void);
#endif

/* clock.c */

double get_accumulated_cpu(void);
void set_accumulated_cpu(double seconds);
uint32_t millisec_elapsed_time(void);
double get_accumulated_time(void);
void set_accumulated_time(double seconds);
uint64_t timestamp(void);


/* fork.c */
#if HAVE_FORK
#define MAX_CHILDREN 8
extern int num_children;
extern int child_num;
extern int affinity_opt;
extern int child_affinity[MAX_CHILDREN];

uint64_t get_lowtide(void);
void child_eliminate_term(uint32_t subseq, uint32_t m, uint64_t p);
void parent_thread(uint64_t p, int parity);
void init_threads(uint64_t pmin, void (*fun)(uint64_t,int));
int fini_threads(uint64_t pmax);
#endif


/* events.c */

typedef enum
{
  initialise_events,
  received_sigterm,
  received_sigint,
  received_sighup,
#if HAVE_FORK
  received_sigpipe,
  received_sigchld,
#endif
  report_due,
  save_due,

  factor_found,

  /* add more events here */

  last_event
} event_t;

static inline void check_events(uint64_t current_prime)
{
  extern volatile int event_happened;
  extern void process_events(uint64_t);

  if (event_happened)
  {
#if HAVE_FORK
    /* When multithreading, the current prime is just the one most recently
       dispatched to the children threads, there is no guarantee that it has
       actually been tested yet. get_lowtide() is a lower bound on all
       untested primes.
    */
    if (num_children > 0)
      current_prime = get_lowtide();
#endif
    process_events(current_prime);
  }
}
void handle_signal(int signum);
void notify_event(event_t event);
void check_progress(void);
void notify_factor(void);

/* factors.c */

void save_factor(uint32_t n, uint64_t p);
#if CHECK_FACTORS
int is_factor(uint32_t n, uint64_t p);
#endif


/* files.c */

typedef enum { FF_UNKNOWN = 1, FF_ABCD, FF_ABC, FF_NEWPGEN } format_t;

extern format_t outputFileFormat;

#ifdef EOF
FILE *xfopen(const char *fn, const char *mode, void (*fun)(const char *,...));
void xfclose(FILE *f, const char *fn);
#endif
void read_input_file(const char *file_name);
void write_output_file(int scroll, uint64_t p, const char *file_name);
#if USE_COMMAND_LINE_FILE
void read_argc_argv(int *argc, char ***argv, const char *file_name);
#endif


/* legendre.c */

extern uint32_t seq_mod;
extern const uint_fast32_t *seq_map[2];
extern int seq_parity;
extern int no_lookup_opt;
extern int64_t kc_core;
#if LEGENDRE_CACHE
extern const char *cache_file_name;
#endif
#if CHECK_FOR_GFN
extern int seq_gfn;
#endif

void generate_legendre_lookup_table(void);


/* primes.c */

void init_prime_sieve(uint64_t pmax);
void prime_sieve(uint64_t low_prime, uint64_t high_prime,
                 void (*bsgs_fun)(uint64_t,int));
void fini_prime_sieve(void);

#if CHECK_FOR_GFN
void prime_sieve_gfn(uint64_t low_prime, uint64_t high_prime,
                 void (*bsgs_fun)(uint64_t,int), int y);
#endif


/* priority.c */

void set_process_priority(int level);
void set_cpu_affinity(int cpu_number);


/* sequences.c */

extern uint32_t b_term;
extern int32_t c_term;
extern uint64_t k_term;

extern uint32_t *N;        /* list of remaining n */
extern uint32_t ncount;  /* number of remaining n */

extern uint64_t p_min;
extern uint64_t p_max;
extern uint32_t n_min;
extern uint32_t n_max;

extern int16_t div_shift[POWER_RESIDUE_LCM/2];
extern uint8_t divisor_index[POWER_RESIDUE_LCM+1];

const char *kbnc_str(uint32_t n);
const char *kbc_str(void);
void add_seq_n(uint32_t n);
void finish_candidate_seqs(void);


/* sr1sieve.c */

extern const char *factors_file_name;
extern const char *output_file_name;
extern uint32_t save_period;
extern int verbose_opt;
extern int quiet_opt;
extern uint32_t L1_cache_size;
extern uint32_t L2_cache_size;
extern const char *baby_method_opt;
extern const char *giant_method_opt;
#if HASHTABLE_SIZE_OPT
extern uint32_t hashtable_size;
#endif

void print_status(uint64_t p, uint32_t p_per_sec, uint32_t secs_per_factor);
void start_srsieve(void);
void finish_srsieve(const char *reason, uint64_t p);


/* subseq.c */

typedef struct
{
  uint32_t d;           /* Each term k*b^n+c satisfies n=d (mod Q) */
  uint_fast32_t *M;     /* Bitmap of terms, bit i corresponds to n=iQ+d */
#ifndef NDEBUG
  uint32_t mcount;      /* Number of remaining terms (popcount of M) */
#endif
#if CHECK_FOR_GFN
  uint32_t a, b;
#endif
} subseq_t;

extern subseq_t *SUBSEQ;
extern uint32_t subseq_count;
extern uint32_t subseq_Q;
extern uint32_t factor_count;
extern int duplicates_opt;
extern int benchmarking;

#define SUBSEQ_MAX UINT16_MAX
extern const uint16_t **sc_lists[3][POWER_RESIDUE_DIVISORS];
extern const uint16_t **sc_ladders[3][POWER_RESIDUE_DIVISORS];

void make_subseqs(void);
void make_subseq_congruence_tables(void);
void eliminate_term(uint32_t subseq, uint32_t m, uint64_t p);
uint32_t get_first_term(void);
uint32_t for_each_term(void (*fun)(uint32_t,void *), void *arg);


/* util.c */

void error(const char *fmt, ...) attribute ((noreturn,format(printf,1,2)));
void warning(const char *fmt, ...) attribute ((format(printf,1,2)));
void logger(int print, const char *fmt, ...) attribute ((format(printf,2,3)));
void *xmalloc(uint32_t sz) attribute ((malloc));
void *xrealloc(void *d, uint32_t sz);
void *xmemalign(uint32_t align, uint32_t size) attribute ((malloc));
void xfreealign(void *mem);
void report(int scroll, const char *fmt, ...) attribute ((format(printf,2,3)));
uint32_t gcd32(uint32_t a, uint32_t b) attribute ((const));
static inline const char *plural(int n)
{
  return (n == 1) ? "" : "s";
}

#endif /* _SR1SIEVE_H */
