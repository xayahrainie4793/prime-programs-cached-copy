/* subseq.c -- (C) Geoffrey Reynolds, June 2006.

   Base b sequences in n=qm+r where n0 < n < n1 are represented as a number
   of base b^Q subsequences in m where n0/Q < m < n1/Q.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>
#include "sr1sieve.h"
#include "bitmap.h"

subseq_t *SUBSEQ = NULL;
uint32_t subseq_count = 0;
uint32_t subseq_Q = 0;
uint32_t factor_count;        /* candidate n eliminated so far */
int duplicates_opt;

const uint16_t **sc_lists[3][POWER_RESIDUE_DIVISORS];
const uint16_t **sc_ladders[3][POWER_RESIDUE_DIVISORS];

static uint32_t m_low;
static uint32_t m_high;

static uint_fast32_t subseq_test_m(uint32_t subseq, uint32_t m)
{
  assert(subseq < subseq_count);

  if (m < m_low || m > m_high)
    return 0;

  return test_bit(SUBSEQ[subseq].M, m-m_low);
}

static uint_fast32_t subseq_clear_m(uint32_t subseq, uint32_t m)
{
  assert(subseq < subseq_count);

  if (m < m_low || m > m_high)
    return 0;

  return clear_bit(SUBSEQ[subseq].M, m-m_low);
}

/* Return 1 iff subsequence h of k*b^n+c has any terms with n%a==b.
 */
static int congruent_terms(uint32_t h, uint32_t a, uint32_t b)
{
  uint_fast32_t i, j;
  uint32_t g;

  assert(b < a);

  g = gcd32(a,subseq_Q);
  if (b % g == SUBSEQ[h].d % g)
  {
    j = m_high-m_low;
    for (i = first_bit(SUBSEQ[h].M); i <= j; i = next_bit(SUBSEQ[h].M,i+1))
      if (((i+m_low)*subseq_Q+SUBSEQ[h].d) % a == b)
        return 1;
  }

  return 0;
}

static uint16_t *make_ladder(const uint16_t *sc_list, int len)
{
  uint32_t i, j, k, a;
  uint16_t *ladder;
  uint8_t subseq_d[LIMIT_BASE+1];

  for (i = 0; i < subseq_Q; i++)
    subseq_d[i] = 0;
  subseq_d[subseq_Q] = 1;
  for (i = 0, a = 1; i < len; i++)
    subseq_d[SUBSEQ[sc_list[i]].d] = 1, a++;

  for (i = 0; i < 3; i++)
  {
    if (subseq_d[i] == 1)
      a--;
    subseq_d[i] = 2;
  }

  while (a > 0)
  {
    for (i = 3, j = 2; i <= subseq_Q; i++)
    {
      if (subseq_d[i] == 2)
        j = i;
      else if (subseq_d[i] == 1)
        break;
    }
    assert(i <= subseq_Q);

    if (subseq_d[i-j] == 2)
      subseq_d[i] = 2, a--; /* We can use an existing rung */
    else
    {
      k = MIN(i-j,(i+1)/2); /* Need to create a new rung */
      assert(subseq_d[k]==0);
      subseq_d[k] = 1;
      a++;
      for (k++; k <= j; k++) /* Need to re-check rungs above the new one */
        if (subseq_d[k] == 2)
          subseq_d[k] = 1, a++;
    }
  }

  for (i = 3, a = 2; i <= subseq_Q; i++)
    if (subseq_d[i] == 2)
      a++;

  ladder = xmalloc((a-1)*sizeof(uint16_t));
  for (i = 3, j = 2, k = 0; i <= subseq_Q; i++)
    if (subseq_d[i] == 2)
    {
      assert(subseq_d[i-j]==2);
      ladder[k] = i-j;
      j = i;
      k++;
    }
  assert(k+2 == a);
  ladder[k] = 0;

  return ladder;
}

/* Build tables sc_lists[x][i][j] of pointers to lists of parity x
   subsequences whose terms k*b^m+c satisfy m = j (mod r) for i =
   divisor_index[r].
*/
void make_subseq_congruence_tables(void)
{
  uint32_t h, i, j, k, r, len[3];
  uint16_t subseq_list[3][POWER_RESIDUE_LCM];

  for (j = 1, r = 0; j <= POWER_RESIDUE_LCM; j++)
    if (POWER_RESIDUE_LCM % j == 0)
    {
      for (i = 0; i < 3; i++)
        if (seq_parity+1 == i || seq_parity == 0)
        {
          sc_lists[i][r] = xmalloc(j*sizeof(uint16_t *));
          sc_ladders[i][r] = xmalloc(j*sizeof(uint16_t *));
        }
      for (k = 0; k < j; k++)
      {
        for (len[0] = len[1] = len[2] = h = 0; h < subseq_count; h++)
          if (congruent_terms(h,j,k))
          {
            if (seq_parity != 1)
              if (SUBSEQ[h].d%2 == 1)
                subseq_list[0][len[0]++] = h; /* odd */
            if (seq_parity == 0)
              subseq_list[1][len[1]++] = h; /* odd and even */
            if (seq_parity != -1)
              if (SUBSEQ[h].d%2 == 0)
                subseq_list[2][len[2]++] = h; /* even */
          }
        for (i = 0; i < 3; i++)
        {
          if (len[i] == 0)
          {
            if (seq_parity+1 == i || seq_parity == 0)
            {
              sc_lists[i][r][k] = NULL;
              sc_ladders[i][r][k] = NULL;
            }
          }
          else
          {
            uint16_t *tmp = xmalloc((len[i]+1)*sizeof(uint16_t));
            for (h = 0; h < len[i]; h++)
              tmp[h] = subseq_list[i][h];
            tmp[h] = SUBSEQ_MAX;
            sc_lists[i][r][k] = tmp;
            sc_ladders[i][r][k] = make_ladder(subseq_list[i],len[i]);
          }
        }
      }
      r++;
    }

  assert(r == POWER_RESIDUE_DIVISORS);
}

void make_subseqs(void)
{
  uint32_t nmin[LIMIT_BASE], subseq[LIMIT_BASE];
#if CHECK_FOR_GFN
  uint32_t g[LIMIT_BASE];
#endif
  uint32_t j, r, n, s, s_count;
  subseq_t *tmp;

#ifndef NDEBUG
  uint32_t t_count = 0;
#endif

  subseq_Q = find_best_Q(&s_count);
  tmp = xmalloc(s_count*sizeof(subseq_t));
  m_low = n_min/subseq_Q;
  m_high = n_max/subseq_Q;
  s = 0;

  for (j = 0; j < subseq_Q; j++)
  {
#if CHECK_FOR_GFN
    g[j] = 0;
#endif
    nmin[j] = UINT32_MAX;
  }

  for (j = 0; j < ncount; j++)
  {
    n = N[j];
    r = n % subseq_Q;
    nmin[r] = MIN(nmin[r],n);
#if CHECK_FOR_GFN
    g[r] = gcd32(g[r],n-nmin[r]);
#endif
  }

  for (r = 0; r < subseq_Q; r++)
    if (nmin[r] != UINT32_MAX)
    {
      assert (s < s_count);
      tmp[s].d = r;
      tmp[s].M = make_bitmap(m_high-m_low+1,"subsequence bitmap");
#ifndef NDEBUG
      tmp[s].mcount = 0;
#endif
#if CHECK_FOR_GFN
      tmp[s].a = MAX(1,g[r]);
      tmp[s].b = nmin[r]%tmp[s].a;
#endif
      subseq[r] = s++;
    }

  for (j = 0; j < ncount; j++)
  {
    n = N[j];
    r = n % subseq_Q;
    set_bit(tmp[subseq[r]].M, n/subseq_Q-m_low);
#ifndef NDEBUG
    tmp[subseq[r]].mcount++;
#endif
  }

  free(N);
  assert (s == s_count);

  if (verbose_opt)
    report(1,"Split 1 base %"PRIu32" sequence into %"PRIu32" base %"PRIu32"^%"
           PRIu32" subsequence%s.",b_term,s_count,b_term,subseq_Q,plural(s_count));

  SUBSEQ = tmp;
  subseq_count = s_count;

#ifndef NDEBUG
  if (verbose_opt)
    report(1,"Used %"PRIu32" Kb for subsequence bitmaps.",
           subseq_count*(m_high-m_low)/8/1024);

  for (s = 0, t_count = 0; s < subseq_count; s++)
    t_count += SUBSEQ[s].mcount;
  assert(t_count == ncount);
#endif
}

/* Return 1 iff p == k*b^n+c.
 */
static int is_equal(uint32_t n, uint64_t p)
{
  p -= c_term; /* Assumes |c| < p and p-c < 2^64. */

  if (p % k_term)
    return 0;
  p /= k_term;

  for ( ; p % b_term == 0; p /= b_term)
    n--;

  return (n == 0 && p == 1);
}

int benchmarking = 0;

/* eliminate a single term n=Qm+d for this subsequence.
*/
void eliminate_term(uint32_t subseq, uint32_t m, uint64_t p)
{
  uint32_t n;

  assert(subseq < subseq_count);

  if (!benchmarking)
  {
#if HAVE_FORK
    if (child_num >= 0)
    {
      child_eliminate_term(subseq,m,p);
      return;
    }
#endif

    n = m*subseq_Q+SUBSEQ[subseq].d;

#if CHECK_FACTORS
    if (!is_factor(n,p))
      error("%"PRIu64" DOES NOT DIVIDE %s.",p,kbnc_str(n));
#endif

    if (subseq_test_m(subseq,m))
    {
      if (n <= 64 && is_equal(n,p))
      {
        /* p is an improper factor, p == k*b^n+c */
        logger(1,"%s is prime.",kbnc_str(n));
      }
      else
      {
        /* p is a proper new factor. */
        factor_count++;
        notify_factor();
        save_factor(n,p);
        if (!quiet_opt)
        {
          report(1,"%"PRIu64" | %s",p,kbnc_str(n));
          notify_event(factor_found);
        }
        subseq_clear_m(subseq,m);
      }
    }
    else if (duplicates_opt)
    {
      /* p is a duplicate factor. */
      if (!quiet_opt)
      {
        report(1,"%"PRIu64" | %s (duplicate)",p,kbnc_str(n));
        notify_event(factor_found);
      }
    }
  }
}

uint32_t get_first_term()
{
  uint32_t i, j;

  for (i = 0; i <= m_high-m_low; i++)
    for (j = 0; j < subseq_count; j++)
      if (test_bit(SUBSEQ[j].M,i))
        return ((i+m_low)*subseq_Q+SUBSEQ[j].d);

  return 0;
}

uint32_t for_each_term(void (*fun)(uint32_t,void *), void *arg)
{
  uint32_t i, j, count;

  for (count = 0, i = 0; i <= m_high-m_low; i++)
    for (j = 0; j < subseq_count; j++)
      if (test_bit(SUBSEQ[j].M,i))
        fun((i+m_low)*subseq_Q+SUBSEQ[j].d,arg), count++;

  return count;
}
