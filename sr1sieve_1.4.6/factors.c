/* factors.c -- (C) Geoffrey Reynolds, May 2006.

   Factors file routines and misc factoring related functions.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include "sr1sieve.h"
#if CHECK_FACTORS
#include "arithmetic.h"
#endif

void save_factor(uint32_t n, uint64_t p)
{
  FILE *factors_file;

  if (factors_file_name != NULL)
  {
    factors_file = xfopen(factors_file_name,"a",error);

    if (fprintf(factors_file,"%"PRIu64" | %s\n",p,kbnc_str(n)) < 0)
      error("Could not write to factors file `%s'.",factors_file_name);

    xfclose(factors_file,factors_file_name);
  }
}

#if CHECK_FACTORS
int is_factor(uint32_t n, uint64_t p)
{
  uint64_t res;

  assert(p > 0);

#if (USE_ASM && __GNUC__ && (__i386__ || (__x86_64__ && USE_FPU_MULMOD)))
  /* We don't know whether b/p or 1.0/p is on the FPU stack, or even whether
     the FPU is being used at all, so save the old mode and initialise.
  */
  uint16_t rnd_save = mod64_rnd;
  mod64_init(p);
#elif HAVE_FORK
  /* Parent thread also reports factors found by children, in which case the
     mod64_init() will not have been called.
   */
# if (USE_ASM && __GNUC__ && (__i386__ || __x86_64__))
  uint16_t rnd_save = mod64_rnd;
# endif
  if (num_children > 0)
    mod64_init(p);
#endif

  res = mulmod64(powmod64(b_term,n,p),k_term%p,p);

#if (USE_ASM && __GNUC__ && (__i386__ || (__x86_64__ && USE_FPU_MULMOD)))
  mod64_fini();
  mod64_rnd = rnd_save;
#elif HAVE_FORK
# if (USE_ASM && __GNUC__ && (__i386__ || __x86_64__))
  mod64_rnd = rnd_save;
# endif
  if (num_children > 0)
    mod64_fini();
#endif

  return (c_term == -1) ? (res == 1) : (res == p-1);
}
#endif
