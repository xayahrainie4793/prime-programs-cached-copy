/* bsgs.c -- (C) Geoffrey Reynolds, April 2006.

   Implementation of a baby step giant step algorithm for finding all
   n in the range nmin <= n <= nmax satisfying b^n=d_i (mod p) where b
   and each d_i are relatively prime to p.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sr1sieve.h"
#include "arithmetic.h"
#include "bitmap.h"
#include "hashtable.h"

/*
  Giant step baby step algorithm, from Wikipaedea:

  input: A cyclic group G of order n, having a generator a and an element b.

  Output: A value x satisfying a^x = b (mod n).

  1. m <-- ceiling(sqrt(n))
  2. For all j where 0 <= j < m:
     1. Compute a^j mod n and store the pair (j, a^j) in a table.
  3. Compute a^(-m).
  4. g <-- b.
  5. For i = 0 to (m - 1):
     1. Check if g is the second component (a^j) of any pair in the table.
     2. If so, return im + j.
     3. If not, g <-- ga^(-m) mod n.
*/


static uint32_t m;     /* Number of baby steps */
static uint32_t M;     /* Number of giant steps */
static uint32_t sieve_low;

typedef struct {uint32_t m; uint32_t M;} steps_t;
static steps_t *steps;

/* Use functions defined in hash-i386.S or hash-x86-64.S.
 */
#if SHORT_HASHTABLE && !CONST_EMPTY_SLOT && USE_ASM && __i386__
uint32_t build_hashtable_i386(uint32_t m) attribute ((regparm(1)));
#define build_hashtable _build_hashtable_i386
#elif SHORT_HASHTABLE && !CONST_EMPTY_SLOT && USE_ASM && __x86_64__
uint32_t build_hashtable_x86_64(uint32_t m);
#define build_hashtable build_hashtable_x86_64
#else
static uint32_t build_hashtable(uint32_t m)
{
  uint64_t bj0;
  uint32_t j;

  bj0 = BJ64[0];
  BJ64[m] = bj0;
  clear_hashtable(hsize);
  j = 0;
  do { post_insert64(j); }
  while (BJ64[++j] != bj0);

  return (j < m) ? j : 0;
}
#endif

static int baby_vec_len;
static void
#if USE_ASM && defined(__i386__) && defined(__GNUC__)
__attribute__((regparm(3)))
#endif
     (*baby_vec_mulmod64)(const uint64_t *,uint64_t *,int);

static uint32_t baby_steps(uint64_t inv_b, uint64_t p)
{
  uint64_t b;

  /* b <- inv_b^Q (mod p) */
  b = powmod64(inv_b,subseq_Q,p);

  BJ64[0] = powmod64(b,sieve_low,p);
  BJ64[1] = mulmod64(BJ64[0],b,p);
  b = sqrmod64(b,p);
  if (baby_vec_len > 2)
  {
    vec_mulmod64_initb(b);
#ifdef HAVE_vector_sse2
    vec2_mulmod64_sse2(BJ64,BJ64+2,baby_vec_len-2);
#else
    vec2_mulmod64(BJ64,BJ64+2,baby_vec_len-2);
#endif
    vec_mulmod64_finib();
    b = powmod64(b,baby_vec_len/2,p);
  }

  /* Don't bother to check whether m > baby_vec_len. This could result in a
     total overrun of up to 2*baby_vec_len-1.
   */
  vec_mulmod64_initb(b);
  baby_vec_mulmod64(BJ64, BJ64+baby_vec_len, m-baby_vec_len);
  vec_mulmod64_finib();

  return build_hashtable(m);
}

static uint64_t bQ;           /* b^Q (mod p) */
static const uint16_t *C16;   /* List of candidate subsequences for BSGS */
static uint64_t *D64;         /* D64[j] holds -c/(k*b^(im+d)) (mod p) for
                                 subsequence C16[j] at giant step i. */

/* Use functions defined in hash-i386.S or hash-x86-64.S.
 */
#if SHORT_HASHTABLE && !CONST_EMPTY_SLOT && USE_ASM && __i386__
uint32_t search_hashtable_i386(const uint64_t *D64, uint32_t cc)
     attribute ((pure));
#define search_hashtable search_hashtable_i386
#elif SHORT_HASHTABLE && !CONST_EMPTY_SLOT && USE_ASM && __x86_64__
uint32_t search_hashtable_x86_64(const uint64_t *D64, uint32_t cc)
     attribute ((pure));
#define search_hashtable search_hashtable_x86_64
#else
static uint32_t search_hashtable(const uint64_t *D64, uint32_t cc)
{
  uint32_t k;

  for (k = 0; k < cc; k++)
    if (lookup64(D64[k]) != HASH_NOT_FOUND)
      break;

  return k;
}
#endif

#if SHORT_HASHTABLE && !CONST_EMPTY_SLOT && USE_ASM && __i386__ && 0
void giant_i386(const uint64_t *D64, uint32_t cc, uint32_t M,
                void (*fun)(uint32_t,uint32_t,uint32_t,uint64_t),
                uint64_t b, uint64_t p);
#define giant_steps giant_i386
#define HAVE_new_giant_steps
#elif SHORT_HASHTABLE && !CONST_EMPTY_SLOT && USE_ASM && __x86_64__
#if USE_FPU_MULMOD
void giant4_x87_64(const uint64_t *D64, uint32_t cc, uint32_t M,
                   void (*fun)(uint32_t,uint32_t,uint32_t,uint64_t),
                   uint64_t b, uint64_t p);
#define giant_steps giant4_x87_64
#else
void giant4_x86_64(const uint64_t *D64, uint32_t cc, uint32_t M,
                    void (*fun)(uint32_t,uint32_t,uint32_t,uint64_t),
                    uint64_t b, uint64_t p);
#define giant_steps giant4_x86_64
#endif
#define HAVE_new_giant_steps
#endif

static void eliminate_1(uint32_t i, uint32_t j, uint32_t k, uint64_t p)
{
  eliminate_term(C16[k],sieve_low+i*m+j,p);
}

static void eliminate_hashtable_terms(uint32_t i, uint32_t cc, uint64_t p)
{
  uint32_t j, k;

  for (k = 0; k < cc; k++)
    if ((j = lookup64(D64[k])) != HASH_NOT_FOUND)
      eliminate_1(i,j,k,p);
}

#ifndef HAVE_new_giant_steps
static void
#if USE_ASM && defined(__i386__) && defined(__GNUC__)
__attribute__((regparm(3)))
#endif
     (*giant_vec_mulmod64)(const uint64_t *,uint64_t *,int);

static void giant_steps(uint64_t p, uint32_t cc)
{
  uint64_t b;
  uint32_t i;

  assert (M > 1);
  assert (cc > 0);

  b = powmod64(bQ,m,p); /* b <- b^Qm (mod p) */
  vec_mulmod64_initb(b);
  i = 1;
  do
  {
    giant_vec_mulmod64(D64,D64,cc);
    if (search_hashtable(D64,cc) < cc)
      eliminate_hashtable_terms(i,cc,p);
  } while (++i < M);
  vec_mulmod64_finib();
}
#endif


/* Assign BJ64[i] = b^i (mod p) for each i in the ladder.
   Return b^Q (mod p).
 */
static uint64_t climb_ladder(const uint16_t *ladder, uint64_t b, uint64_t p)
{
  uint32_t i;

  BJ64[0] = 1;
  BJ64[1] = b;
  BJ64[2] = sqrmod64(b,p);
  for (i = 2; *ladder > 0; i += *ladder++)
    BJ64[i+*ladder] = mulmod64(BJ64[i],BJ64[*ladder],p);

  return BJ64[subseq_Q];
}

/* This function builds the list C16[] of subsequences (k*b^d)*(b^Q)^m+c for
   which p may be a factor (-ckb^d is a quadratic/cubic/quartic/quintic
   residue with respect to p) and initialises the table D64[] with the
   values -c/(k*b^d) (mod p). As a side effect, bQ is set to the value b^Q
   (mod p) for use later in bsgs64(), and m,M are set to the best choice of
   baby,giant steps. Returns the number of subsequences listed in C16[].
*/
static
uint32_t attribute ((noinline)) setup64(uint64_t inv_b, uint64_t p, int parity)
{
  uint64_t neg_ck, p_s, X[POWER_RESIDUE_LCM+1];
  const uint16_t *list;
  uint32_t r, h, j;
  int16_t s;

  /* neg_ck <-- -k/c (mod p) == -ck (mod p) */
  neg_ck = k_term;
  if (p <= neg_ck)
  {
    neg_ck %= p;
    if (neg_ck == 0)
      return 0;
  }
  if (c_term > 0)
    neg_ck = p - neg_ck;

  s = div_shift[(p/2)%(POWER_RESIDUE_LCM/2)];
  if (s == 0)
  {
    /* p = 1 (mod 2) is all we know, check for quadratic residues only.
     */

    /* Precompute b^d (mod p) for 0 <= d <= Q, as necessary. */
    bQ = climb_ladder(sc_ladders[parity+1][0][0],b_term,p);

    /* For each subsequence (k*b^d)*(b^Q)^(n/Q)+c, compute
       -c/(k*b^d) (mod p) and add the subsequence to the bsgs list.
    */
    PRE2_MULMOD64_INIT(neg_ck);
    C16 = sc_lists[parity+1][0][0];
    for (j = 0, list = C16; (h = *list++) < SUBSEQ_MAX; j++)
      D64[j] = PRE2_MULMOD64(BJ64[SUBSEQ[h].d],neg_ck,p); /* -c/(k*b^d) */
    PRE2_MULMOD64_FINI();
    m = steps[j].m;
    M = steps[j].M;
    return j;
  }
  else if (s > 0)
  {
    /* p = 1 (mod s), where s is not a power of 2. Check for r-th power
       residues for each prime power divisor r of s.
    */
    p_s = p/s;
  }
  else /* s < 0 */
  {
    /* p = 1 (mod 2^s), where s > 1. Check for r-th power residues for each
       divisor r of s. We handle this case seperately to avoid computing p/s
       using plain division.
    */
    p_s = p >> (-s);
#ifndef NDEBUG
    s = 1 << (-s);
#endif
  }

  /* For 0 <= r < s, X[r] <- 1/(b^r)^((p-1)/s) */
  X[0] = 1;
  X[1] = powmod64(inv_b,p_s,p);
  PRE2_MULMOD64_INIT(X[1]);
  for (r = 1; X[r] != 1; r++)
    X[r+1] = PRE2_MULMOD64(X[r],X[1],p);
  PRE2_MULMOD64_FINI();
  assert(s%r == 0);
  /* 1/(b^r)^((p-1)/s)=1 (mod p) therefore (1/(b^r)^((p-1)/s))^y=1 (mod p)
     for 0 <= y < s/r. (Could we do more with this?)
  */

  /* X[r] <- (-ck)^((p-1)/s) */
  X[r] = powmod64(neg_ck,p_s,p);

  /* Find h such that X[h]=X[r], i.e. (-ckb^h)^((p-1)/r)=1 (mod p),
     or h=r if not found. */
  for (h = 0; X[h] != X[r]; h++)
    ;

  if (h < r && (C16 = sc_lists[parity+1][divisor_index[r]][h]) != NULL)
  {
    /* -ckb^n is an r-power residue for at least one term k*b^n+c of
       this sequence.
    */

    /* Precompute b^d (mod p) for 0 <= d <= Q, as necessary. */
    bQ = climb_ladder(sc_ladders[parity+1][divisor_index[r]][h],b_term,p);
    PRE2_MULMOD64_INIT(neg_ck);
    for (j = 0, list = C16; (h = *list++) < SUBSEQ_MAX; j++)
    {
      /* -ckb^d is an r-th power residue for at least one term
         (k*b^d)*(b^Q)^(n/Q)+c of this subsequence.
      */
      D64[j] = PRE2_MULMOD64(BJ64[SUBSEQ[h].d],neg_ck,p); /* -c/(k*b^d) */
    }
    PRE2_MULMOD64_FINI();
    m = steps[j].m;
    M = steps[j].M;
    return j;
  }

  return 0;
}

static void bsgs64(uint64_t p, int parity)
{
  uint32_t i, j, k, cc;
  uint64_t inv_b;

  /* inv_b <-- 1/base (mod p) */
  inv_b = (b_term == 2) ? (p+1)/2 : invmod32_64(b_term,p);

  vec_mulmod64_initp(p);

  if ((cc = setup64(inv_b,p,parity)) > 0)
  {
    /* Baby steps. */
    if ((i = baby_steps(inv_b,p)) > 0) /* Unlikely */
    {
      /* i is the order of b (mod p). This is all the information we need to
         determine every solution for this p, so no giant steps are needed.
      */
      for (k = 0; k < cc; k++)
        for (j = lookup64(D64[k]); j < m*M; j += i)
          eliminate_term(C16[k],sieve_low+j,p);
    }
    else
    {
      /* First giant step. */
      if (search_hashtable(D64,cc) < cc)
        eliminate_hashtable_terms(0,cc,p);

      /* Remaining giant steps. */
      if (M > 1)
#ifdef HAVE_new_giant_steps
        giant_steps(D64,cc,M,eliminate_1,powmod64(bQ,m,p),p);
#else
        giant_steps(p,cc);
#endif
    }
  }

  vec_mulmod64_finip();
}

static void choose_steps(uint32_t r, uint32_t s)
{
  uint32_t m, M;

  /*
    r = range of n, s = number of subsequences.
    In the worst case we will do do one table insertion and one mulmod
    for m baby steps, then s table lookups and s mulmods for M giant
    steps. The average case depends on how many solutions are found
    and how early in the loop they are found, which I don't know how
    to analyse. However for the worst case we just want to minimise
    m + s*M subject to m*M >= r, which is when m = sqrt(s*r).
  */

  M = MAX(1,rint(sqrt((double)r/s)));
  m = MIN(r,ceil((double)r/M));

  if (m > HASH_MAX_ELTS)
  {
    /* There are three ways to handle this undesirable case:
       1. Allow m > HASH_MAX_ELTS (undersize hash table).
       2. Restrict m <= HASH_MAX_ELTS (suboptimal baby/giant-step ratio).
       3. Recompile with SHORT_HASHTABLE=0.
    */
    M = ceil((double)r/HASH_MAX_ELTS);
    m = ceil((double)r/M);
  }

  assert(m <= HASH_MAX_ELTS);

  steps[s].m = m;
  steps[s].M = M;
}

/* Increase baby steps as far as the next multiple of vec_len if doing so
   will reduce the number of giant steps. This minimises the extra work due
   to doing mulmods in blocks of vec_len.
*/
static void trim_steps(uint32_t r, int vec_len, uint32_t max_baby_steps)
{
  uint32_t i, m, M;

  for (i = 1; i <= subseq_count; i++)
  {
    m = steps[i].m;
    M = steps[i].M;
    if (m % vec_len) /* m is not a multiple of vec_len */
    {
      m += (vec_len - m%vec_len); /* next multiple of vec_len. */
      while (M > 1 && m*(M-1) >= r) /* minimise M */
        M--;
      while ((m-1)*M >= r) /* minimise m */
        m--;
      if (m <= max_baby_steps)
      {
        steps[i].m = m;
        steps[i].M = M;
      }
    }
  }
}


/* Choose which mulmod scheme to use for bsgs.
 */
typedef struct {
#if USE_ASM && defined(__i386__) && defined(__GNUC__)
  void __attribute__((regparm(3))) (*fun)(const uint64_t *,uint64_t *,int);
#else
  void (*fun)(const uint64_t *,uint64_t *,int);
#endif
  const char *desc;
  int vec_len;
} vec_fun_t;

static const vec_fun_t vec_funs[] = {
#ifdef HAVE_vector_sse2
  { vec2_mulmod64_sse2, "sse2/2", 2 },
  { vec4_mulmod64_sse2, "sse2/4", 4 },
  { vec8_mulmod64_sse2, "sse2/8", 8 },
  { vec16_mulmod64_sse2, "sse2/16", 16 },
#endif
  { vec2_mulmod64, "gen/2", 2 },
  { vec4_mulmod64, "gen/4", 4 },
#if USE_ASM && defined(__x86_64__) && defined(__GNUC__)
  { vec6_mulmod64, "gen/6", 6 },
#endif
  { vec8_mulmod64, "gen/8", 8 }
};

#define NUM_TEST_PRIMES 10
static const uint64_t test_primes[NUM_TEST_PRIMES] = {
  UINT64_C(2250000000000023),
  UINT64_C(2250000000000043),
  UINT64_C(2250000000000059),
  UINT64_C(2250000000000061),
  UINT64_C(2250000000000079),
  UINT64_C(2250000000000089),
  UINT64_C(2250000000000113),
  UINT64_C(2250000000000163),
  UINT64_C(2250000000000191),
  UINT64_C(2250000000000209) };

#define NUM_VEC_FUN sizeof(vec_funs)/sizeof(vec_fun_t)
static int choose_baby_method(int print)
{
  uint64_t t0, t1, best[NUM_VEC_FUN];
  uint32_t cc;
  int i, j;

  /* How many subsequences pass the power residue tests? Probably 1/4 or
     less pass, but the speed of the bsgs routine is more important when a
     larger number pass since then BSGS accounts for a larger proportion of
     the total work done. Benchmark assuming 1/2 have passed.
  */
  cc = (subseq_count+1)/2;
  m = steps[cc].m;
  for (i = 0; i < NUM_VEC_FUN; i++)
  {
    best[i] = UINT64_MAX;
    baby_vec_len = vec_funs[i].vec_len;
    baby_vec_mulmod64 = vec_funs[i].fun;
    for (j = 0; j < NUM_TEST_PRIMES; j++)
    {
      uint64_t inv_b, p;

      p = test_primes[j];
      inv_b = (b_term == 2) ? (p+1)/2 : invmod32_64(b_term,p);

      vec_mulmod64_initp(p);
      t0 = timestamp();
      baby_steps(inv_b,p);
      t1 = timestamp();
      best[i] = MIN(best[i],t1-t0);
      vec_mulmod64_finip();
    }
  }

  if (print >= 2)
    for (i = 0; i < NUM_VEC_FUN; i++)
      report(1,"Best time for baby step method %s: %"PRIu64".",
             vec_funs[i].desc,best[i]);

#if USE_ASM && defined(__x86_64__) && defined(__GNUC__)
  /* Prevent gen/6 method being selected unless running on AMD. */
  if (!is_amd())
    for (i = 0; i < NUM_VEC_FUN; i++)
      if (vec_funs[i].vec_len == 6)
        best[i] = UINT64_MAX;
#endif

  for (i = 1, j = 0; i < NUM_VEC_FUN; i++)
    if (best[i] < best[j])
      j = i;

  return j;
}

#ifndef HAVE_new_giant_steps
static int choose_giant_method(int print)
{
  uint64_t t0, t1, best[NUM_VEC_FUN];
  uint32_t cc;
  int i, j;

  cc = (subseq_count+1)/2;
  m = steps[cc].m;
  M = MAX(2,steps[cc].M);

  mod64_init(test_primes[0]);
  D64[0] = 23;
  for (i = 1; i < cc; i++)
    D64[i] = mulmod64(D64[i],D64[0],test_primes[0]);
  mod64_fini();

  for (i = 0; i < NUM_VEC_FUN; i++)
  {
    best[i] = UINT64_MAX;
    giant_vec_mulmod64 = vec_funs[i].fun;
    for (j = 0; j < NUM_TEST_PRIMES; j++)
    {
      uint64_t inv_b, p;

      p = test_primes[j];
      inv_b = (b_term == 2) ? (p+1)/2 : invmod32_64(b_term,p);

      vec_mulmod64_initp(p);
      baby_steps(inv_b,p);
      t0 = timestamp();
      giant_steps(p,cc);
      t1 = timestamp();
      best[i] = MIN(best[i],t1-t0);
      vec_mulmod64_finip();
    }
  }

  if (print >= 2)
    for (i = 0; i < NUM_VEC_FUN; i++)
      report(1,"Best time for giant step method %s: %"PRIu64".",
             vec_funs[i].desc,best[i]);

#if USE_ASM && defined(__x86_64__) && defined(__GNUC__)
  /* Prevent gen/6 method being selected unless running on AMD. */
  if (!is_amd())
    for (i = 0; i < NUM_VEC_FUN; i++)
      if (vec_funs[i].vec_len == 6)
        best[i] = UINT64_MAX;
#endif

  for (i = 1, j = 0; i < NUM_VEC_FUN; i++)
    if (best[i] < best[j])
      j = i;

  return j;
}
#endif

/* Set baby_steps, giant_steps, climb_ladder functions.
   Return baby step vector size.
*/
static void choose_bsgs_methods(int print)
{
  benchmarking = 1; /* Don't check or report any factors. */

  int baby_method;
  if (baby_method_opt == NULL)
    baby_method = choose_baby_method(print);
  else
  {
    for (baby_method = 0; baby_method < NUM_VEC_FUN; baby_method++)
      if (strcmp(baby_method_opt,vec_funs[baby_method].desc)==0)
        break;
    if (baby_method >= NUM_VEC_FUN)
      error("Unknown baby step method `%s'.",baby_method_opt);
  }
  baby_vec_len = vec_funs[baby_method].vec_len;
  baby_vec_mulmod64 = vec_funs[baby_method].fun;

#ifndef HAVE_new_giant_steps
  int giant_method;
  if (giant_method_opt == NULL)
    giant_method = choose_giant_method(print);
  else
  {
    for (giant_method = 0; giant_method < NUM_VEC_FUN; giant_method++)
      if (strcmp(giant_method_opt,vec_funs[giant_method].desc)==0)
        break;
    if (giant_method >= NUM_VEC_FUN)
      error("Unknown giant step method `%s'.",giant_method_opt);
  }
  giant_vec_mulmod64 = vec_funs[giant_method].fun;
#endif

#ifndef HAVE_new_giant_steps
  if (print)
    report(1,"Baby step method %s, giant step method %s.",
           vec_funs[baby_method].desc,vec_funs[giant_method].desc);
#else
  if (print)
    report(1,"Baby step method %s, giant step method new/4.",
           vec_funs[baby_method].desc);
#endif

  benchmarking = 0;
}

static void init_sieve(void)
{
  uint32_t r, s, sieve_high, max_baby_steps;

  sieve_low = n_min / subseq_Q;
  sieve_high = n_max / subseq_Q;

  r = sieve_high - sieve_low + 1;
  steps = xmalloc((subseq_count+1)*sizeof(steps_t));

  for (s = 1, max_baby_steps = 0; s <= subseq_count; s++)
  {
    choose_steps(r,s);
    assert(sieve_high < sieve_low+steps[s].m*steps[s].M);
    max_baby_steps = MAX(max_baby_steps,steps[s].m);
  }

  if (verbose_opt >= 2)
   report(1,"BSGS range: %"PRIu32"*%"PRIu32" - %"PRIu32"*%"PRIu32".",
          steps[1].m,steps[1].M,steps[subseq_count].m,steps[subseq_count].M);

  /* Allow room for vector operations to overrun the end of the array.
   */
  D64 = xmemalign(64,(subseq_count+7)*sizeof(uint64_t));
  memset(D64,0,(subseq_count+7)*sizeof(uint64_t));
  init_hashtable(max_baby_steps);

  choose_bsgs_methods(0); /* Dummy run to prep cache. */
  choose_bsgs_methods(verbose_opt);
  trim_steps(r,baby_vec_len,max_baby_steps);
}

static void fini_sieve(void)
{
  fini_hashtable();
  xfreealign(D64);
  free(steps);
}

void CODE_PATH_NAME(sieve)(void)
{
  void (*fun)(uint64_t,int);

#if HAVE_FORK
  if (num_children > 0)
    fun = parent_thread;
  else
#endif
    fun = bsgs64;

  init_sieve();
#if HAVE_FORK
  /* Start children before initializing the Sieve of Eratosthenes, to
     minimise forked image size. */
  if (num_children > 0)
    init_threads(p_min,bsgs64);
#endif
  init_prime_sieve(p_max);
  start_srsieve();

#if CHECK_FOR_GFN
  if (seq_gfn)
    prime_sieve_gfn(p_min,p_max,fun,seq_gfn);
  else
    prime_sieve(p_min,p_max,fun);
#else
  prime_sieve(p_min,p_max,fun);
#endif

#if HAVE_FORK
  /* Check that all children terminated normally before assuming that
     the range is complete. */
  if (num_children > 0 && fini_threads(p_max))
    finish_srsieve("range is incomplete",get_lowtide());
  else
#endif
    finish_srsieve("range is complete",p_max);

  fini_prime_sieve();
  fini_sieve();
}
