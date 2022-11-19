/* choose.c -- (C) Geoffrey Reynolds, August 2006.

   Choose the best Q for sieving in subsequence base b^Q.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <assert.h>
#ifndef NDEBUG
#include <inttypes.h>
#endif
#include <math.h>
#include <stdlib.h>
#include "sr1sieve.h"

#define BABY_WORK    1.0    /* 1 mulmod, 1 insert */
#define GIANT_WORK   1.0    /* 1 mulmod, 1 lookup */
#define EXP_WORK     0.3    /* 1 mulmod at most */
#define SUBSEQ_WORK  1.4    /* 1 mulmod, 1 lookup (giant step 0), plus some
                               unknown overhead (array references etc?). */

/* Set the number of baby/giant steps used in BSGS on s subsequences.
 */
static void choose_steps(uint32_t *baby,uint32_t *giant,uint32_t Q,uint32_t s)
{
  uint32_t r = n_max/Q-n_min/Q+1;

  *giant = MAX(1,rint(sqrt((double)r/s)));
  *baby = MIN(r,ceil((double)r/(*giant)));
}

/* Return an estimate of the work needed to do one BSGS iteration on s
   subsequences in base b^Q.
*/
static uint32_t estimate_work(uint32_t Q, uint32_t s)
{
  uint32_t baby, giant;

  choose_steps(&baby,&giant,Q,s);

  return baby*BABY_WORK + s*(giant-1)*GIANT_WORK + Q*EXP_WORK + s*SUBSEQ_WORK;
}

/* Try to estimate the work needed to sieve s subsequences in base b^Q. This
   is a bit rough.
*/
static uint32_t rate_Q(uint32_t Q, uint32_t s)
{
  uint32_t work, i, W[POWER_RESIDUE_LCM+1];

  assert (Q % 2 == 0);
  assert (Q % BASE_MULTIPLE == 0);
  assert (LIMIT_BASE % Q == 0);

  if (s >= SUBSEQ_MAX)
    return UINT32_MAX;

#if SUBSEQ_Q_OPT
  if (subseq_Q == Q)
    return 1;
#endif

  for (i = 2, work = 0; i <= POWER_RESIDUE_LCM; i += 2)
  {
    if (POWER_RESIDUE_LCM % i == 0)
      W[i] = estimate_work(Q,(s+i-1)/i);

    if (gcd32(i+1,POWER_RESIDUE_LCM) == 1)
      work += W[gcd32(i,POWER_RESIDUE_LCM)];
  }

#ifndef NDEBUG
  report(1,"Q=%"PRIu32": work=%"PRIu32".",Q,work);
#endif

  return work;
}

static uint32_t count_residue_classes(uint32_t d, uint32_t Q, const uint8_t *R)
{
  uint32_t i, count;
  uint8_t R0[LIMIT_BASE];

  assert (Q % d == 0);

  for (i = 0; i < d; i++)
    R0[i] = 0;

  for (i = 0; i < Q; i++)
    if (R[i])
      R0[i%d] = 1;

  for (i = 0, count = 0; i < d; i++)
    if (R0[i])
      count++;

  return count;
}

#define NDIVISORS (LIMIT_BASE/BASE_MULTIPLE)
uint32_t find_best_Q(uint32_t *subseqs)
{
  uint32_t i, j;
  uint32_t S[NDIVISORS], W[NDIVISORS];
  uint8_t R[LIMIT_BASE];

  for (j = 0; j < LIMIT_BASE; j++)
    R[j] = 0;
  for (j = 0; j < ncount; j++)
    R[N[j]%LIMIT_BASE] = 1;
  for (i = 0, j = 0; j < NDIVISORS; j++)
    if (NDIVISORS % (j+1) == 0)
    {
      S[j] = count_residue_classes((j+1)*BASE_MULTIPLE,LIMIT_BASE,R);
      W[j] = rate_Q((j+1)*BASE_MULTIPLE,S[j]);
      if (W[j] < W[i])
        i = j;
    }

  *subseqs = S[i];
  return (i+1)*BASE_MULTIPLE;
}
