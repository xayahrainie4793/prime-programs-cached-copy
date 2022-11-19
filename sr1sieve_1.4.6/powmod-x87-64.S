/* powmod-x87-64.S -- (C) Geoffrey Reynolds, January 2008.

   uint64_t powmod64_x87_64(uint64_t b, uint64_t n, uint64_t p);

     Returns b^n (mod p), where 0 <= a < p < 2^62.

     Assumes FPU is set to double extended precision and round to zero.
     Assumes %st(0) contains 1.0/p computed with above settings.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/


#include "config.h"

/* Set USE_CMOV=1 to use a conditional move instead of a branch, even when
   the branch is predictable. Slower on Intel, but maybe faster on AMD.
*/
#ifndef USE_CMOV
# define USE_CMOV 0
#endif

#ifdef _WIN64
#define ARG1 %rcx
#define ARG2 %rdx
#define ARG3 %r8
#define REG4 %r9
#define REG4l %r9d
#define REG5 %r10
#define REG6 %r11
#define TMP 16(%rsp) /* Shadow space */
#else
#define ARG1 %rdi
#define ARG2 %rsi
#define ARG3 %rdx
#define REG4 %rcx
#define REG4l %ecx
#define REG5 %r8
#define REG6 %r9
#define TMP -32(%rsp) /* Red zone */
#endif


  .text
  .globl _powmod64_x87_64
  .globl powmod64_x87_64

  .p2align 4,,15

_powmod64_x87_64:
powmod64_x87_64:
  push  %rbp
  lea TMP, %rbp
  mov $1, %eax
  mov $1, REG4l

  .p2align 4,,15

mulsqr_loop:
  /* %rax = REG4 = a, ARG1 = b, ARG2 = n, ARG3 = p, %st(0) = 1.0/p */

  mov ARG1, (%rbp)
  mov %rax, 8(%rbp)
  fildll  (%rbp)
  fildll  8(%rbp)
  fmul  %st(2), %st(0)
  fmul  %st(1), %st(0)
  fistpll 16(%rbp)
  fmul  %st(0), %st(0)
  fmul  %st(1), %st(0)
  fistpll 24(%rbp)
  imul  ARG1, %rax
  imul  ARG1, ARG1
  mov 16(%rbp), REG5
  mov 24(%rbp), REG6
  imul  ARG3, REG5
  imul  ARG3, REG6
  sub REG5, %rax
  sub REG6, ARG1
  mov %rax, REG5
  mov ARG1, REG6

correct_mul:
  sub     ARG3, REG5  /* CF=0 predicted */
#if USE_CMOV
  cmovnc  REG5, %rax
#else
  jc      correct_sqr
  mov     REG5, %rax
#endif

correct_sqr:
  sub     ARG3, REG6  /* CF=0 predicted */
#if USE_CMOV
  cmovnc  REG6, ARG1
#else
  jc  shift
  mov     REG6, ARG1
#endif

shift:
  /* REG4 = a, %rax = a*b, ARG1 = b^2, ARG2 = n */

  shr ARG2    /* CF unpredictable */
  cmovnc  REG4, %rax  /* Discard multiply if CF=0 */
  mov %rax, REG4
  jnz mulsqr_loop

  pop %rbp
  ret
