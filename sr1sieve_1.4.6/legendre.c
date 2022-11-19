/* legendre.c -- (C) Geoffrey Reynolds, September 2006.

   Generate lookup tables to replace the Legendre(-ckb^n,p) function.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include "sr1sieve.h"
#include "bitmap.h"

uint32_t seq_mod;
const uint_fast32_t *seq_map[2];
int seq_parity;         /* odd = -1, even = 1, both = 0 */
int no_lookup_opt = 0;  /* Set to disable lookup table generation. */
int64_t kc_core;        /* -c*core(k) for k*b^n+c */

#if LEGENDRE_CACHE
const char *cache_file_name = NULL;
#if HAVE_MMAP
#include <sys/mman.h>
#endif
#endif

#if CHECK_FOR_GFN
int seq_gfn = 0;
#endif

/*
  Return the value of the Legendre symbol (a/p) if gcd(a,p)==1, 0 otherwise.
 */
static int jacobi32(int32_t a, uint32_t p)
{
  uint32_t x, y, t;
  int sign;

  if (a < 0)
    x = -a, sign = (p % 4 == 1) ? 1 : -1;
  else
    x = a, sign = 1;

  for (y = p; x > 0; x %= y)
  {
    for ( ; x % 2 == 0; x /= 2)
      if (y % 8 == 3 || y % 8 == 5)
        sign = -sign;

    t = x, x = y, y = t;

    if (x % 4 == 3 && y % 4 == 3)
      sign = -sign;
  }

  return (y == 1) ? sign : 0;
}

/* Return the least factor d of n such that n/d is a square.
 */
static uint64_t core64(uint64_t n)
{
  /* TODO: Improve the trial division algorithm (use sieve). */

  uint64_t c = 1;
  uint_fast32_t d, q, r;

  assert(n > 0);

  while (n % 2 == 0)
  {
    n /= 2;
    if (n % 2 != 0)
      c *= 2;
    else
      n /= 2;
  }
  while (n % 3 == 0)
  {
    n /= 3;
    if (n % 3 != 0)
      c *= 3;
    else
      n /= 3;
  }
  r = sqrtl(n);
  if ((uint64_t)r*r == n)
    return c;
  for (q = 5, d = 4; q <= r; d = 6-d, q += d)
    if (n % q == 0)
    {
      do
      {
        n /= q;
        if (n % q != 0)
          c *= q;
        else
          n /= q;
      }
      while (n % q == 0);
      r = sqrtl(n);
      if ((uint64_t)r*r == n)
        return c;
    }

  return c * n;
}

#if CHECK_FOR_GFN
/* Return the greatest value y < 6 such that x = A^(2^y).
 */
static uint32_t squares(uint64_t x)
{
  uint64_t r;
  uint32_t y;

  for (y = 0; y < 5; y++)
  {
    r = sqrtl(x);
    if (r*r != x)
      break;
    x = r;
  }

  return y;
}

/* Return the greatest value y < 6 such that x = M*(2^y).
 */
static uint32_t twos(uint32_t x)
{
  uint32_t y;

  for (y = 0; y < 5; y++)
  {
    if (x % 2)
      break;
    x /= 2;
  }

  return y;
}

/* Find the greatest value y < 6 such that every term is of the form A^(2^y)+1.
*/
static uint32_t gen_fermat_y(void)
{
  uint32_t i, y;

  /* k*b^n+1 must satisfy: */

  /* 1.  n = M*2^y.  */
  for (i = 0, y = 5; i < subseq_count; i++)
    y = MIN(y,MIN(twos(SUBSEQ[i].a),twos(SUBSEQ[i].b)));
  y += squares(b_term);

  /* 2.  k = A^(2^y) */
  y = MIN(y,squares(k_term));

  return y;
}
#endif /* CHECK_FOR_GFN */

/* Sets bit i of map iff r is a quadratic residue with respect to all primes
   p = 2*i+1 (mod 2*m).
*/
#define PROGRESS_STEP 1000000
static void
build_legendre_table(int32_t r, uint32_t m, int parity, uint_fast32_t *map)
{
  uint32_t i, j;

  for (i = 0; i < m; )
  {
    /* Always indicate progress when building a large table, even if
       verbose output is not requested, to avoid giving the impression
       that the program is hung.
     */
    if (verbose_opt || m > PROGRESS_STEP)
      report(0,"Building %s Legendre symbol table: %.1f%%",
             (parity == -1) ? "odd" : "even", 100.0*i/m);

    for (j = MIN(i+PROGRESS_STEP,m); i < j; i++)
      if (jacobi32(r,2*i+1) == 1)
        set_bit(map,i);
  }
}

#if LEGENDRE_CACHE
#include <stdio.h>
#include <stdlib.h>

static FILE *cache_file;

#define CACHE_VERSION 2

#define THIS_BYTE_ORDER  0x44332211
#define OTHER_BYTE_ORDER1 0x11223344
#define OTHER_BYTE_ORDER2 0x33441122
#define OTHER_BYTE_ORDER3 0x22114433

static void xread(void *data, size_t size, size_t count)
{
  if (fread(data,size,count,cache_file) < count)
    error("Failed read from cache file `%s'.",cache_file_name);
}

static uint32_t xread32(void)
{
  uint32_t tmp;

  xread(&tmp,sizeof(tmp),1);

  return tmp;
}

static void attribute((noreturn)) cache_file_error(const char *msg)
{
  fclose(cache_file);
  error("Cache file `%s': %s.", cache_file_name, msg);
}

static int read_legendre_cache(int32_t r, uint32_t m, int parity,
                               uint_fast32_t **map0, uint_fast32_t **map1)
{
  uint_fast32_t *map, check, check0, check1;
  uint32_t i, size0, size1;

  assert(cache_file_name != NULL);

  if ((cache_file = fopen(cache_file_name,"rb")) == NULL)
    return 0;

  switch(xread32())
  {
    case THIS_BYTE_ORDER:
      break;
    case OTHER_BYTE_ORDER1:
    case OTHER_BYTE_ORDER2:
    case OTHER_BYTE_ORDER3:
      cache_file_error("Wrong byte order");
    default:
      cache_file_error("Unrecognised file type");
  }

  if (xread32() != CACHE_VERSION)
    cache_file_error("Wrong version");
  if (xread32() != sizeof(uint_fast32_t))
    cache_file_error("Wrong register width");
  if (xread32() != r || xread32() != m)
    cache_file_error("Wrong sequence");
  if (xread32() != parity)
    cache_file_error("Wrong parity");

  report(0,"Loading Legendre symbol tables from version %u cache file `%s'...",
         (unsigned int)CACHE_VERSION, cache_file_name);

  size0 = xread32();
  size1 = xread32();
  xread(&check0,sizeof(check0),1);
  xread(&check1,sizeof(check1),1);

#if HAVE_MMAP
  map = mmap(NULL,8*sizeof(uint32_t)+(size0+size1+2)*sizeof(uint_fast32_t),
             PROT_READ,MAP_PRIVATE,fileno(cache_file),0);
  if (map != (void *)-1)
  {
    if (size0 > 0)
    {
      uint32_t off = 8*sizeof(uint32_t)/sizeof(uint_fast32_t)+2;
      for (i = 0, check = 0; i < size0; i++)
        check += map[i+off];
      if (check != check0)
        cache_file_error("Bad checksum");
      if (map0 == NULL)
        cache_file_error("Bad parity");
      *map0 = map+off;
    }
    else
      cache_file_error("Bad parity");

    if (size1 > 0)
    {
      uint32_t off = 8*sizeof(uint32_t)/sizeof(uint_fast32_t)+2+size0;
      for (i = 0, check = 0; i < size1; i++)
        check += map[i+off];
      if (check != check1)
        cache_file_error("Bad checksum");
      if (map1 == NULL)
        cache_file_error("Bad parity");
      *map1 = map+off;
    }
    else if (parity == 0)
      cache_file_error("Bad parity");
  }
  else
  {
    warning("Failed to mmap() cache file `%s', will use malloc() instead.",
            cache_file_name);
#endif
    if (size0 > 0)
    {
      map = xmalloc(size0*sizeof(uint_fast32_t));
      xread(map,sizeof(uint_fast32_t),size0);
      for (i = 0, check = 0; i < size0; i++)
        check += map[i];
      if (check != check0)
        cache_file_error("Bad checksum");
      if (map0 == NULL)
        cache_file_error("Bad parity");
      *map0 = map;
    }
    else
      cache_file_error("Bad parity");

    if (size1 > 0)
    {
      map = xmalloc(size1*sizeof(uint_fast32_t));
      xread(map,sizeof(uint_fast32_t),size1);
      for (i = 0, check = 0; i < size1; i++)
        check += map[i];
      if (check != check1)
        cache_file_error("Bad checksum");
      if (map1 == NULL)
        cache_file_error("Bad parity");
      *map1 = map;
    }
    else if (parity == 0)
      cache_file_error("Bad parity");
#if HAVE_MMAP
  }
#endif

  fclose(cache_file);

  if (verbose_opt)
    report(1,"Loaded Legendre symbol tables from version %u cache file `%s'.",
           (unsigned int)CACHE_VERSION, cache_file_name);

  return 1;
}

static void xwrite(const void *data, size_t size, size_t count)
{
  if (fwrite(data,size,count,cache_file) < count)
    error("Failed write to cache file `%s'.",cache_file_name);
}

static void xwrite32(uint32_t data)
{
  xwrite(&data,sizeof(data),1);
}

static void write_legendre_cache(int32_t r, uint32_t m, int parity,
                        const uint_fast32_t *map0, const uint_fast32_t *map1)
{
  uint_fast32_t check;
  uint32_t i, size0, size1;

  assert(cache_file_name != NULL);

  if ((cache_file = fopen(cache_file_name,"wb")) == NULL)
    error("Cannot create cache file `%s'.",cache_file_name);

  report(0,"Saving Legendre symbol tables to version %u cache file `%s'...",
         (unsigned int)CACHE_VERSION, cache_file_name);

  xwrite32(THIS_BYTE_ORDER);
  xwrite32(CACHE_VERSION);
  xwrite32(sizeof(uint_fast32_t));
  xwrite32(r);
  xwrite32(m);
  xwrite32(parity);
  size0 = (map0 == NULL) ? 0 : bitmap_size(m);
  size1 = (map1 == NULL) ? 0 : bitmap_size(m);
  xwrite32(size0);
  xwrite32(size1);
  /* 16-aligned */
  for (i = 0, check = 0; i < size0; i++)
    check += map0[i];
  xwrite(&check,sizeof(check),1);
  for (i = 0, check = 0; i < size1; i++)
    check += map1[i];
  xwrite(&check,sizeof(check),1);
  if (size0 > 0)
    xwrite(map0,sizeof(uint_fast32_t),size0);
  if (size1 > 0)
    xwrite(map1,sizeof(uint_fast32_t),size1);
  fclose(cache_file);

  if (verbose_opt)
    report(1,"Saved Legendre symbol tables to version %u cache file `%s'.",
           (unsigned int)CACHE_VERSION, cache_file_name);
}
#endif

/* Set kc_core = -c*core(k). Then, unless no_lookup_opt is set:
 
   For sequences with single parity terms (parity=+/-1):
   Set seq_mod and seq_map[0] for sequence k*b^n+c so that bit (p/2)%mod
   of seq_map[0] is set if and only if
   (-ck/p)=1 for sequences with all n even,
   (-bck/p)=1 for sequences with all n odd.

   For sequences with mixed parity terms (parity=0):
   Set seq_mod, seq_map[0] and seq_map[1] for sequence k*b^n+c so that
   bit (p/2)%mod of seq_map[0] is set if and only if (-ck/p)=1,
   bit (p/2)%mod of seq_map[1] is set if and only if (-bck/p)=1.

   In the worst case each table for k*b^n+c could be 4*b*k bits long.

   TODO: Consider the case where gcd(k,b) > 1.
*/
void generate_legendre_lookup_table(void)
{
  uint64_t tmp;
  uint_fast32_t *map0, *map1;
  int32_t r;
  uint32_t i, m, b;

  for (i = 0, r = SUBSEQ[i].d % 2; i < subseq_count; i++)
    if (SUBSEQ[i].d % 2 != r)
      break;
  seq_parity = (i < subseq_count) ? 0 : r ? -1 : 1;

  tmp = core64(k_term);
  if (tmp > INT64_MAX)
    error("Square-free part %"PRIu64" of k in %s is too large.",tmp,kbc_str());

  kc_core = -c_term*tmp;
  seq_mod = 0;

  if (no_lookup_opt || tmp > UINT32_MAX)
    return;

  m = tmp;
#if CHECK_FOR_GFN
  if (c_term == 1 && m == 1 && (i = gen_fermat_y()) > 0)
  {
    /* Only need to check primes p = 1 (mod 2^(i+1)) */

    seq_mod = 1 << i;
    map0 = make_bitmap(seq_mod,"Legendre symbol table");
    set_bit(map0,0);
    if (verbose_opt)
      report(1,"Recognised Generalised Fermat sequence A^%"PRIu32"+1",seq_mod);
    seq_map[0] = map0;
    seq_gfn = i+1;
    return;
  }
#endif

  b = core64(b_term);
  /* We may need the signed product c*b*m, and the unsigned product m*b*2.
  */
  if (m >= INT32_MAX/b)
    return;
  r = -c_term*m;

  switch (seq_parity)
  {
    default:
    case -1: /* odd n, test for (-bck/p)==1 */
      m *= b;
      r *= b;
      /* Fall through */

    case 1: /* even n, test for (-ck/p)==1 */
      if ((r < 0 && (-r) % 4 != 3) || (r > 0 && r % 4 != 1))
        m *= 2;
#if LEGENDRE_CACHE
      if (cache_file_name != NULL)
      {
        if (read_legendre_cache(r,m,seq_parity,&map0,NULL) == 0)
          if ((map0 = make_bitmap(m,NULL)) != NULL)
          {
            build_legendre_table(r,m,seq_parity,map0);
            write_legendre_cache(r,m,seq_parity,map0,NULL);
#if HAVE_MMAP
            free(map0);
            if (read_legendre_cache(r,m,seq_parity,&map0,NULL) == 0)
              error("Failed to mmap new cache file `%s'.",cache_file_name);
#endif
          }
        seq_map[0] = map0;
      }
      else
#endif
        if ((map0 = make_bitmap(m,NULL)) != NULL)
        {
          build_legendre_table(r,m,seq_parity,map0);
          seq_map[0] = map0;
        }
      i = m;
      break;

    case 0: /* odd and even n, test for (-ck/p)==1 and (-bck/p)==1 */
      m = m*2*b;
#if LEGENDRE_CACHE
      if (cache_file_name != NULL)
      {
        if (read_legendre_cache(r,m,0,&map0,&map1) == 0)
          if ((map0 = make_bitmap(m,NULL)) != NULL)
          {
            if ((map1 = make_bitmap(m,NULL)) != NULL)
            {
              build_legendre_table(r,m,1,map0);
              build_legendre_table(r*b,m,-1,map1);
              write_legendre_cache(r,m,0,map0,map1);
#if HAVE_MMAP
              free(map0);
              free(map1);
              if (read_legendre_cache(r,m,0,&map0,&map1) == 0)
                error("Failed to mmap new cache file `%s'.",cache_file_name);
#endif
            }
            else
            {
              free(map0);
              map0 = NULL;
            }
          }
        seq_map[0] = map0;
        seq_map[1] = map1;
      }
      else
#endif
        if ((map0 = make_bitmap(m,NULL)) != NULL)
        {
          if ((map1 = make_bitmap(m,NULL)) != NULL)
          {
            build_legendre_table(r,m,1,map0);
            build_legendre_table(r*b,m,-1,map1);
            seq_map[0] = map0;
            seq_map[1] = map1;
          }
          else
          {
            free(map0);
            map0 = NULL;
          }
        }
      i = 2*m;
      break;
  }

  if (map0 == NULL)
  {
    if (verbose_opt)
      report(1,"Not enough memory (%"PRIu32" Kb) for Legendre symbol tables.",
             i/8/1024);
    return;
  }

  if (verbose_opt)
    report(1,"Using %"PRIu32" Kb for Legendre symbol tables.",
           i/8/1024);

  seq_mod = m;
}
