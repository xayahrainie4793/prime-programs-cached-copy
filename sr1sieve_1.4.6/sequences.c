/* sequences.c -- (C) Geoffrey Reynolds, April 2006.

   Routines for creating and manipulating candidates sequences.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include "sr1sieve.h"

uint64_t p_min;
uint64_t p_max;
uint32_t *N;        /* list of remaining n */
uint32_t ncount;  /* number of remaining n */
uint32_t b_term;
int32_t c_term;
uint64_t k_term;
uint32_t n_min;
uint32_t n_max;

int16_t div_shift[POWER_RESIDUE_LCM/2];
uint8_t divisor_index[POWER_RESIDUE_LCM+1];

static uint32_t nsize;   /* N has room for nsize n */

static char seq_buf[48];
const char *kbnc_str(uint32_t n)
{
  sprintf(seq_buf,"%"PRIu64"*%"PRIu32"^%"PRIu32"%+"PRId32,
          k_term,b_term,n,c_term);

  return seq_buf;
}

const char *kbc_str(void)
{
  sprintf(seq_buf,"%"PRIu64"*%"PRIu32"^n%+"PRId32,k_term,b_term,c_term);

  return seq_buf;
}

#define SEQ_N_GROW_SIZE 4096
void add_seq_n(uint32_t n)
{
  uint32_t count;

  count = ncount;
  if (count == nsize)
  {
    nsize += SEQ_N_GROW_SIZE;
    N = xrealloc(N,nsize*sizeof(uint32_t));
  }

  N[count] = n;
  ncount++;
}

void finish_candidate_seqs(void)
{
  if (ncount < 1)
    error("Empty sequence.");

  n_min = N[0];
  n_max = N[ncount-1];

  make_subseqs();
  generate_legendre_lookup_table();
  make_subseq_congruence_tables();

  /* Precompute the div_shift table. Powers 2^x are represented as -x. */ 
  {
    uint32_t i, r, divide, shift;

    for (i = 0; i < POWER_RESIDUE_LCM/2; i++)
    {
      r = gcd32(2*i,POWER_RESIDUE_LCM);
      for (shift = 0, divide = r; divide % 2 == 0; shift++)
        divide /= 2;
      if (divide > 1)
        div_shift[i] = r;
      else if (shift > 1)
        div_shift[i] = -shift;
      else
        div_shift[i] = 0;
    }

    /* Build table divisor_index[r], for divisors r of POWER_RESIDUE_LCM.
     */
    for (i = 1, r = 0; i <= POWER_RESIDUE_LCM; i++)
      if (POWER_RESIDUE_LCM % i == 0)
        divisor_index[i] = r++;

    assert(r == POWER_RESIDUE_DIVISORS);
  }
}
