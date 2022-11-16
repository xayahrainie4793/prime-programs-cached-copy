#include <assert.h>
#include <inttypes.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "srsieve.h"
#include "bitmap.h"
#include "arithmetic.h"

// These variables should never be accessed globally.
uint32_t *smallPrimes, smallPrimesCount = 0;
FILE     *algebraicFactorFile = 0;

void     checkForSpecialForm(uint64_t k_seq, int32_t c_seq);
int32_t  removeAlgebraicSimpleTerm(uint_fast32_t *B, uint64_t k_seq, int32_t c_seq);
int32_t  removeAlgebraicSimpleRoot(uint_fast32_t *B, uint64_t k_seq, int32_t c_seq);
int32_t  removeAlgebraicComplexRoot(uint_fast32_t *B, uint64_t k_seq, int32_t c_seq);
int32_t  checkBase2(uint_fast32_t *B, uint64_t k_seq, int32_t c_seq);
int32_t  checkPower4(uint_fast32_t *B, uint64_t k_seq, int32_t c_seq);
int32_t  checkAndLogFactor(uint_fast32_t *B, uint64_t k_seq, uint32_t n, int32_t c_seq, const char *fmt, ...);
void     getRoot(uint64_t number, uint32_t *root, uint32_t *power);
uint32_t getFactorList(uint64_t the_number, uint32_t *factor_list, uint32_t *power_list);

// This function, called from subseq.c, will find algebraic factors for k*b^n+/-1
//
// Note that if k_seq and base are both odd, then this function is never called
uint32_t find_algebraic(uint_fast32_t *B, uint64_t k_seq, int32_t c_seq)
{
   uint32_t  startingCount, removedCount = 0;
   
   startingCount = 1 + n_max - n_min;

   checkForSpecialForm(k_seq, c_seq);

   removedCount = removeAlgebraicSimpleTerm(B, k_seq, c_seq);
   
   if (startingCount - removedCount > 0)
      removedCount += removeAlgebraicSimpleRoot(B, k_seq, c_seq);
   
   if (startingCount - removedCount > 0)
      removedCount += removeAlgebraicComplexRoot(B, k_seq, c_seq);
         
   if (startingCount - removedCount > 0)
      removedCount += checkPower4(B, k_seq, c_seq);

   if (startingCount - removedCount > 0)
      removedCount += checkBase2(B, k_seq, c_seq);
   
   return removedCount;
}

// If c = -1 and k=2^f and b=2^g for any f and g, then this is a Mersenne number.
// If c = 1 and k=x^f and b=x^g then we have a Generalized Fermat Number.
// 
// Note that this will only give a warning, but will not remove terms since it
// is not looking for algebraic factors.
void    checkForSpecialForm(uint64_t k_seq, int32_t c_seq)
{
   uint32_t  broot, bpower;
   uint32_t  kroot, kpower;
   
   // Needs to be k*b^n-1 or k*b^n+1
   if (c_seq != -1 && c_seq != 1)
      return;
   
   // Now we have base = broot^bpower
   getRoot(base, &broot, &bpower);

   // If c = -1, then b must be 2^g for some g
   if (c_seq == -1 && broot != 2)
      return;

   if (k_seq > 1)
   {
      // Now we have base = broot^bpower
      getRoot(k_seq, &kroot, &kpower);
      
      // If c = -1, then k must be 2^f for some f
      if (c_seq == -1 && kroot != 2)
         return;
      
      // If c = +1, then k and b must have the same root
      if (c_seq == +1 && kroot != broot)
         return;
   }
   else
      kpower = 0;

   if (c_seq == -1)
      warning("(sf) Sequence %"PRIu64"*%u^n-1 is the form of a Mersenne number", k_seq, base);
   else
   {
      if (kpower == 0)
      {
         if (bpower == 1)
            warning("(sf) Sequence %"PRIu64"*%u^n+1 as it is a GFN", k_seq, base);
         else
            warning("(sf) Sequence %"PRIu64"*%u^n+1 as it is a GFN --> %u^(%u*n)+1", k_seq, base, broot, bpower);
      }
      else
      {
         if (bpower == 1)
            warning("(sf) Sequence %"PRIu64"*%u^n+1 as it is a GFN --> %u^(n+%u)+1", k_seq, base, broot, kpower);
         else
            warning("(sf) Sequence %"PRIu64"*%u^n+1 as it is a GFN --> %u^(%u*n+%u)+1", k_seq, base, broot, bpower, kpower);
      }
   }
}

// If k=x^f and b=x^g then we have a special form
//
// If c=-1, then x-1 will factor all terms
// If c=+1, then x+1 will factor terms where (f+n*g) is odd
int32_t  removeAlgebraicSimpleTerm(uint_fast32_t *B, uint64_t k_seq, int32_t c_seq)
{
   uint32_t  removedCount = 0;
   uint32_t  n, broot, bpower;
   uint32_t  kroot, kpower = 0;
   
   // Now we have base = broot^bpower
   getRoot(base, &broot, &bpower);
   
   // For broot = 2 and c_seq = -1, we have 1 as a divisor, ignore it.
   if (broot + c_seq == 1)
      return 0;

   if (k_seq == 1)
   {
      kpower = 0;
      
      report("(st) For sequence %"PRIu64"*%u^n%+d -> b = %u^%u", k_seq, base, c_seq,
               broot, bpower);
      }
      else 
      {
      // Now we have base = kroot^kpower
      getRoot(k_seq, &kroot, &kpower);

      if (kroot != broot)
         return 0;
      
      report("(st) For sequence %"PRIu64"*%u^n%+d -> k = %u^%u and b = %u^%u", k_seq, base, c_seq,
               kroot, kpower, broot, bpower);
   }   

   for (n=n_min; n<=n_max; n++)
   {
      // For c = +1, the simple factor will only divide the term if kpower + n*bpower is odd
      if (c_seq == +1)
         if ((kpower + n*bpower) % 2 == 0)
            continue;
      
      removedCount += checkAndLogFactor(B, k_seq, n, c_seq, "%d", broot + c_seq);
   }
   
   report("(st) Sequence %"PRIu64"*%u^n%+d has %d terms removed has they have the factor %d", 
         k_seq, base, c_seq, removedCount, broot + c_seq);

   return removedCount;
}

// If k=x^f and b=y^g then look for algebraic factors.
// Note that all n will be covered by these factors.
//
// If c=-1, z divides f, z divides n*g, any n*g, then x^(f/z)*y^((n*g)/z)-1 will be a factor.
// If c=+1, z divides f, z divides n*g, odd n*g, then x^(f/z)*y^((n*g)/z)+1 will be a factor.
int32_t  removeAlgebraicSimpleRoot(uint_fast32_t *B, uint64_t k_seq, int32_t c_seq)
{
   uint32_t  removedCount = 0, loopRemovedCount = 0;
   uint32_t  n, idx;
   uint32_t  broot, bpower;
   uint32_t  kroot, kpower;
   
   // If k = 1, then removeAlgebraicSimpleTerm() will have handled it
   if (k_seq == 1)
      return 0;

      // Now we have base = broot^bpower
   getRoot(base, &broot, &bpower);

   // Now we have k = kroot^kpower
   getRoot(k_seq, &kroot, &kpower);

   // If kroot == broot, then removeAlgebraicSimpleTerm() will have handled it
   if (kroot == broot)
      return 0;

   // k_seq must be x^f for some x and f must be greater than 1
   if (kpower == 1)
      return 0;
   
   for (idx=2; idx<=kpower; idx++)
   {
      // Given k=kroot^kpower, find all idx where kpower%idx = 0.
      // If kpower == 6, then we look for algebraic factors with the
      // forms kroot^(6/2), kroot^(6/2), and kroot^(6/6).
      if (kpower % idx != 0)
         continue;
      
      if (c_seq == +1)
      {
         // x^1 is not a divisor of x^n+1, so skip this idx
         if (idx == kpower)
            continue;
         
         // x^y for even y is not a divisor of x^(y*n)+1, so skip this idx.
         if ((kpower/idx)%2 ==0)
            continue;
      }
      
      if (bpower == 1)
         report("(sr) For sequence %"PRIu64"*%u^n%+d -> (%u^%u)*%u^n%+d", k_seq, base, c_seq,
                  kroot, kpower, base, c_seq);
      else
         report("(sr) For sequence %"PRIu64"*%u^n%+d -> (%u^%u)*%u^(%u*n)%+d", k_seq, base, c_seq,
                  kroot, kpower, broot, bpower, c_seq);
               
      loopRemovedCount = 0;

      for (n=n_min; n<=n_max; n++)
      {
         if ((n*bpower) % idx != 0)
            continue;
         
            // For c = +1, x^y does not divide x^(y*n)+1 when y*n is even
         if (c_seq == +1)
            if ((n*bpower) % 2 == 0)
               continue;
         
         if (idx == kpower)
            loopRemovedCount += checkAndLogFactor(B, k_seq, n, c_seq, "(%u*%u^%u%+d)", kroot, broot, (bpower*n)/idx, c_seq);
         else
            loopRemovedCount += checkAndLogFactor(B, k_seq, n, c_seq, "((%u^%u)*%u^%u%+d)", kroot, kpower/idx, broot, (bpower*n)/idx, c_seq);
      }

      if (idx == kpower)
      {
         if (bpower == 1)
            report("(sr) Sequence %"PRIu64"*%u^n%+d has %d terms removed due to algebraic factors of the form %u*%u^(n/%u)%+d", 
                     k_seq, base, c_seq, loopRemovedCount, kroot, base, idx, c_seq);
         else
            report("(sr) Sequence %"PRIu64"*%u^n%+d has %d terms removed due to algebraic factors of the form %u*%u^((%u*n)/%u)%+d", 
                     k_seq, base, c_seq, loopRemovedCount, kroot, broot, bpower, idx, c_seq);

      }
      else
      {
         if (bpower == 1)
            report("(sr) Sequence %"PRIu64"*%u^n%+d has %d terms removed due to algebraic factors of the form (%u^%u)*%u^(n/%u)%+d", 
                     k_seq, base, c_seq, loopRemovedCount, kroot, kpower/idx, base, idx, c_seq);
         else
            report("(sr) Sequence %"PRIu64"*%u^n%+d has %d terms removed due to algebraic factors of the form (%u^%u)*%u^((%u*n)/%u)%+d", 
                     k_seq, base, c_seq, loopRemovedCount, kroot, kpower/idx, broot, bpower, idx, c_seq);
      }

      removedCount += loopRemovedCount;
   }

   return removedCount;
}

// If k_seq*b^ii = root^m for root > 1 and m > 1, then we can remove algebraic factors
// because we can find a simple root for some of the terms.
//
//   1)  base^n = b^(bpow*n)
//
//   2)  k_seq -> k*b^y
// 
//   3)  k_seq*base^n -> (k*b^y)*(b^(bpow*n))
//                    -> k*b^(bpow*n+y)
//
//   4)  For each ii where k*b^ii = root^m for root > 1 and m > 1:
//                    -> (k*b^ii)*b^(bpow*n+y-ii)
//                    -> root^m*b^(bpow*n+y-ii)
//
//   5)  Find prime factors of m -> [f1,f2,...,fn]
//
//   6)  For each prime factor f of m, if (bpow*n+y-ii)%f = 0, then:
//          if c = +1 and f is odd: root^f*b^(bpow*n+y-ii)%f+1 has a factor of root*b^(bpow*n+y-ii)/f+1
//          if c = -1: root^f*b^(bpow*n+y-ii)%f-1 has a factor of root*b^(bpow*n+y-ii)/f-1
int32_t  removeAlgebraicComplexRoot(uint_fast32_t *B, uint64_t k_seq, int32_t c_seq)
{
   uint32_t  removedCount = 0, loopRemovedCount;
   uint32_t  ii, m, n, z, z1, z2, idx;
   uint32_t  kf_count, bf_count, kbf_count;
   uint32_t  broot, bpower, kroot, kpower;
   uint32_t  bf_factor[50], bf_power[50];
   uint32_t  kf_factor[50], kf_power[50];
   uint32_t  kbf_factor[50], kbf_power[50];
   char      root[50], part[50], addParenthesis;

   // We want c_seq = +1 or -1
   if (c_seq != 1 && c_seq != -1)
      return 0;

   // base and k_seq must share a divisor
   if (gcd64(base, k_seq) == 1)
      return 0;
   
      // Now we have base = broot^bpower
   getRoot(base, &broot, &bpower);

      // Now we have k = kroot^kpower
   getRoot(k_seq, &kroot, &kpower);

   // If k is r^x and b is r^y for some r, x, and y, then we have simple roots that
   // are handled elsewhere.
   if (broot == kroot)
      return 0;
   
   // Step 1:  Factorize b
   bf_count = getFactorList(base, bf_factor, bf_power);
   
   // Step 2:  Factorize k
   kf_count = getFactorList(k_seq, kf_factor, kf_power);
      
      for (z1=0; z1<kf_count; z1++)
      {
         kbf_factor[z1] = kf_factor[z1];
         kbf_power[z1] = kf_power[z1];
      }
      
      kbf_count = kf_count;
      
   for (ii=1; ii<10; ii++)
   {
      // We could multiply out k_seq*base^ii and factor that value, but
      // k_seq*base^ii will likely exceed what can be stored in a 64-bit integer.
      // We'll just mulitply k*b^(ii-1) by b
      for (z2=0; z2<bf_count; z2++)
      {
         idx = 99999;

         for (z1=0; z1<kbf_count; z1++)
            if (kbf_factor[z1] == bf_factor[z2])
            {
               kbf_power[z1] += bf_power[z2];
               idx = z1;
            }
         
         if (idx == 99999)
         {
            kbf_factor[kbf_count] = bf_factor[z2];
            kbf_power[kbf_count] = bf_power[z2];
            kbf_count++;            
         }
      }

      // We now have k_seq*base^ii = f1^p1 * f2^p2 ... fn^pn
      // Step 4:  Determine m such that m = gcd(p1, p2, ..., pn)
      m = kbf_power[0];

      for (idx=1; idx<kbf_count; idx++)
         m = gcd32(m, kbf_power[idx]);
      
      // If m then k_seq*base^ii is not a perfect power
      if (m == 1)
         continue;

      for (idx=2; idx<m; idx++)
      {         
         // Given k*b^ii=root^m, find all idx where m%idx = 0.
         // If m == 6, then we look for algebraic factors with the
         // forms root^2 and root^3.
         if (m % idx != 0)
            continue;
         
         z = m / idx;
                  
      root[0] = 0;

         for (z1=0; z1<kbf_count; z1++)
      {
            if (kbf_power[z1] == 0)
            continue;
         
            if (z1 == 0)
         {
               if (kbf_power[z1] == z)
                  sprintf(root, "%u", kbf_factor[z1]);
            else
            {
               addParenthesis = 1;
                  sprintf(root, "%u^%u", kbf_factor[z1], kbf_power[z1] / z);
            }
         }
         else
         {
            addParenthesis = 1;
               sprintf(part, root);
            
               if (kbf_power[z1] == z)
                  sprintf(root, "%s*%u", part, kbf_factor[z1]);
            else
                  sprintf(root, "%s*%u^%u", part, kbf_factor[z1], kbf_power[z1] / z);
         }
      }
      
         if (addParenthesis)
      {
            if (ii == 1)
               report("(cr) For sequence %"PRIu64"*%u^n%+d -> %"PRIu64"*%u = (%s)^%u", k_seq, base, c_seq,
                        k_seq, base, root, z);
            else
               report("(cr) For sequence %"PRIu64"*%u^n%+d -> %"PRIu64"*%u^%u = (%s)^%u", k_seq, base, c_seq,
                        k_seq, base, ii, root, z);
         }
         else
         {
            if (ii == 1)
               report("(cr) For sequence %"PRIu64"*%u^n%+d -> %"PRIu64"*%u = %s^%u", k_seq, base, c_seq,
                        k_seq, base, root, z);
            else
               report("(cr) For sequence %"PRIu64"*%u^n%+d -> %"PRIu64"*%u^%u = %s^%u", k_seq, base, c_seq,
                        k_seq, base, ii, root, z);
         }

         // Don't allow n to go negative
         if (n_min < ii)
            n = n_min;
         else
            n = n_min - ii;

         // n-ii must be greater than 0
         if (n < ii+1) n = ii + 1;
         
         // To avoid 2^1-1 divisors
         if (!strcmp(root, "2") && c_seq == -1 && n-ii==1) n++;
         
         // For c = +1, if z is odd, then double it to make it even as we can only
         // provide factors for odd n.
         if (c_seq == +1 && z%2 == 1)
            z *= 2;

         // Find smallest n >= n_min where (n-ii)%z = 0
         while ((n - ii) % z != 0) n++;
         
         // For c = +1, n must be odd in order for root^(n/f)+1 to be a factor
         if (c_seq == +1 && (n - ii)%2 == 0)
         {
            // If z and n are even, then n%z will never be 0 so we can skip this z.
            if (z%2 == 0)
               continue;
            
            n += z;
         }
                     
         loopRemovedCount = 0;

         for (; n<=n_max; n+=z)
            if (addParenthesis)
               loopRemovedCount += checkAndLogFactor(B, k_seq, n, c_seq, "((%s)*%u^%u%+d)", root, base, (n-ii)/z, c_seq);
            else
               loopRemovedCount += checkAndLogFactor(B, k_seq, n, c_seq, "(%s*%u^%u%+d)", root, base, (n-ii)/z, c_seq);
         
         if (addParenthesis)
            report("(cr) Sequence %"PRIu64"*%u^n%+d has %d terms removed due to algebraic factors of the form (%s)*%u^((n-%d)/%d)%+d", 
                     k_seq, base, c_seq, loopRemovedCount, root, base, ii, z, c_seq);
         else
            report("(cr) Sequence %"PRIu64"*%u^n%+d has %d terms removed due to algebraic factors of the form %s*%u^((n-%d)/%d)%+d", 
                     k_seq, base, c_seq, loopRemovedCount, root, base, ii, z, c_seq);
         
         removedCount += loopRemovedCount;
      }     
   }

   return removedCount;
}

int32_t  checkBase2(uint_fast32_t *B, uint64_t k_seq, int32_t c_seq)
{
   uint64_t  k;
   uint32_t  removedCount = 0;
   uint32_t  m, n, root, y, idx;
   uint32_t  kf_count;
   uint32_t  b, bpow;
   uint32_t  kf_factor[50], kf_power[50];

   // We want c = 1
   if (c_seq != 1)
      return 0;
   
   bpow = 0;
   b = base;
   while (b > 1)
   {
      // if b is odd and greater than 1, then we can't use it
      if (b & 1)
         return 0;
      
      b >>= 1;
      bpow++;
   }

   y = 0;
   k = k_seq;
   while (k % 2 == 0)
   {
      k >>= 1;
      y++;
   }

   // Now that the base is removed, refactorize k
   kf_count = getFactorList(k, kf_factor, kf_power);
   
   // We now have k = f1^p1 * f2^p2 ... fn^pn
   // Now determine m such that m = gcd(p1, p2, ..., pn)
   m = kf_power[0];
      
   for (idx=1; idx<kf_count; idx++)
      m = gcd32(m, kf_power[idx]);

   // We want k where k=root^(4*m)
   if (m % 4 != 0)
      return 0;

   root = 1;
   for (idx=0; idx<kf_count; idx++)
      root *= (uint32_t) pow((double) kf_factor[idx], (double) (kf_power[idx] / 4));

   n = n_min;
   
   // Exit the function right away if n*bpow+y is not divisible
   // by 4 and the gcd of bpow and 4 is not 1 because we'll get
   // an infinite loop in the while statement right after this.
   if ((((n*bpow+y) % 4) != 2) && gcd32(bpow, 4) > 1)
      return 0;
         
   while ((n*bpow+y)% 4 != 2) n++;
   while (n > n_min && n > 4) n -= 4;
   while (n < n_min) n += 4;

   for ( ; n<=n_max; n+=4)
      removedCount += checkAndLogFactor(B, k_seq, n, c_seq, "(%u^2*2^%u-%u*2^%u+1)",
            root, (n*bpow+y)/2, root, (n*bpow+y+2)/4);

   if (removedCount > 0)
   {
      report("Removed %d algebraic factors for %"PRIu64"*%u^n%+d of the form (%u^2)*2^(n/2)-%u*2^((n+2)/4))+1 when n%%4=2",
              removedCount, k_seq, base, c_seq, root, root);
   }

   return removedCount;
}

int32_t  checkPower4(uint_fast32_t *B, uint64_t k_seq, int32_t c_seq)
{
   uint32_t  removedCount = 0;
   uint32_t  n, ninc, broot, bpower, kroot, kpower;
   uint32_t  b1, bexp;

   // We want c = +1
   if (c_seq != 1)
      return 0;

   // We want k to be divisible by 4
   if (k_seq % 4 != 0)
      return 0;
   
   // Now we have base = broot^bpower
   getRoot(base, &broot, &bpower);

   if (bpower % 4 == 0)
      ninc = 1;
   else if (bpower % 4 == 2)
   {
      ninc = 2;
      
      if (bpower > 2)
      {
         // If base = x^(y*2) then compute broot as x^y
         // In other words set broot = sqrt(base)
         bexp = bpower / 2;
         
         b1 = 1;
         while (bexp > 0)
         {
            b1 *= broot;
            bexp--;
         }
      
         bpower = 2;
         broot = b1;
      }
   }
   else
      ninc = 4;
   
   if (k_seq == 4)
      kroot = 1;
   else
   {
      // Now we have k_seq = 4*kroot^kpower
      getRoot(k_seq/4, &kroot, &kpower);
      
      // We want kpower to be a multiple of 4
      if (kpower % 4 != 0)
         return 0;
   }

   n = n_min;
   while (n % ninc > 0) n++;
   
   for (; n<=n_max; n+=ninc)
   {
      if (kroot == 1)
      {
         if (ninc == 1)
            removedCount += checkAndLogFactor(B, k_seq, n, c_seq, "(2*(%u^%u)^%u+2*(%u^%u)^%u+1)",
                  broot, bpower/2, n, broot, bpower/4, n);
         else if (ninc == 2)
            removedCount += checkAndLogFactor(B, k_seq, n, c_seq, "(2*%u^%u+2*%u^%u+1)",
                  broot, n, broot, n/2);
         else
            removedCount += checkAndLogFactor(B, k_seq, n, c_seq, "(2*%u^%u+2*%u^%u+1)",
                  base, n/2, broot, n/4);
      }
      else
      {
         if (ninc == 1)
            removedCount += checkAndLogFactor(B, k_seq, n, c_seq, "(2*%u^%u*(%u^%u)^%u+2*%u^%u*(%u^%u)^%u+1)",
                  kroot, kpower/2, broot, bpower/2, n, kroot, kpower/4, broot, bpower/4, n);
         else if (ninc == 2)
            removedCount += checkAndLogFactor(B, k_seq, n, c_seq, "(2*%u^%u*%u^%u+2*%u^%u*%u^%u+1)",
                  kroot, kpower/2, broot, n, kroot, kpower/4, broot, n/2);
         else
            removedCount += checkAndLogFactor(B, k_seq, n, c_seq, "(2*%u^%u*%u^%u+2*%u^%u*%u^%u+1)",
                  kroot, kpower/2, base, n/2, kroot, kpower/4, base, n/4);
      }
   }

   if (removedCount)
   {
      if (kroot == 1)
      {
         if (ninc == 1)
            report("Removed %d algebraic factors for %"PRIu64"*%u^n%+d of the form (2*(%u^%u)^n+2*(%u^%u)^n+1)",
                  removedCount, k_seq, base, c_seq, broot, bpower/2, broot, bpower/4);
         else if (ninc == 2)
            report("Removed %d algebraic factors for %"PRIu64"*%u^n%+d of the form (2*(%u^%u)^n+2*(%u^%u)^n+1)",
                  removedCount, k_seq, base, n, c_seq, broot, broot, n/2);
         else
            report("Removed %d algebraic factors for %"PRIu64"*%u^n%+d of the form (2*%u^(n/2)+2*%u^(n/4)+1)",
                  removedCount, k_seq, base, c_seq, broot, broot);
      }
      else
      {
         if (ninc == 1)
            report("Removed %d algebraic factors for %"PRIu64"*%u^n%+d of the form (2*%u^%u*(%u^%u)^n+2*%u^%u*(%u^%u)^n+1)",
                  removedCount, k_seq, base, c_seq, kroot, kpower/2, broot, bpower/2, kroot, kpower/4, broot, bpower/4);
         else if (ninc == 2)
            report("Removed %d algebraic factors for %"PRIu64"*%u^n%+d of the form (2*%u^%u*(%u^%u)^n+2*%u^%u*(%u^%u)^n+1)",
                  removedCount, k_seq, base, c_seq, kroot, kpower/2, broot, n, kroot, kpower/4, broot, n/2);
         else
            report("Removed %d algebraic factors for %"PRIu64"*%u^n%+d of the form (2*%u^%u*%u^(n/2)+2*%u^%u*%u^(n/4)+1)",
                  removedCount, k_seq, base, c_seq, kroot, kpower/2, broot, kroot, kpower/4, broot);
      }
   }

   return removedCount;
}

// Build list of primes below 100000
void get_small_primes(void)
{
  uint32_t minp, p;

  // There are less than 100000 primes below 1000000
  smallPrimes = xmalloc(100000 * sizeof(uint32_t));

  smallPrimes[0] = 2;
  smallPrimes[1] = 3;
  smallPrimesCount = 2;
  for (p = 5; p < 1000000; p += 2)
    for (minp = 0; minp <= smallPrimesCount; minp++)
    {
      if (smallPrimes[minp] * smallPrimes[minp] > p)
      {
        smallPrimes[smallPrimesCount] = p;
        smallPrimesCount++;
        break;
      }
      if (p % smallPrimes[minp] == 0)
        break;
    }
}

void free_small_primes(void)
{
   free(smallPrimes);
}

int32_t  checkAndLogFactor(uint_fast32_t *B, uint64_t k_seq, uint32_t n, int32_t c_seq, const char *fmt, ...)
{   
   va_list args;
   char  fName[50];

   if (!test_bit(B, n-n_min))
      return 0;
   
   clear_bit(B, n-n_min);
   
   sprintf(fName, "alg_%"PRIu64"_%u_%+d.log", k_seq, base, c_seq);

   if (!algebraicFactorFile)
      algebraicFactorFile = fopen(fName, "a");
   
   va_start(args,fmt);
   vfprintf(algebraicFactorFile, fmt, args);
   va_end(args);
   
   fprintf(algebraicFactorFile, " | %"PRIu64"*%u^%d%+d\n", k_seq, base, n, c_seq);
   
   return 1;
}

// Find root and power such that root^power = number
void  getRoot(uint64_t number, uint32_t *root, uint32_t *power)
{
   uint32_t  idx, r, rpow, rf_count;
   uint32_t  rf_factor[50], rf_power[50];

   rf_count = getFactorList(number, rf_factor, rf_power);
   
   rpow = rf_power[0];
   for (idx=1; idx<rf_count; idx++)
      rpow = gcd32(rpow, rf_power[idx]);
      
   r = 1;
   for (idx=0; idx<rf_count; idx++)
      r *= (uint32_t) pow((double) rf_factor[idx], (double) (rf_power[idx] / rpow));
   
   *root = r;
   *power = rpow;
}

uint32_t getFactorList(uint64_t theNumber, uint32_t *factorList, uint32_t *powerList)
{
   uint32_t  distinctFactors = 0;
   uint32_t  ii, power;

   // Get the prime factorization of the input number.  Note that since seq_primes is
   // limited to 1000000, that very large numbers might not get fully factored.
   for (ii=0; ii<smallPrimesCount; ii++)
   {
      power = 0;
      while (theNumber % smallPrimes[ii] == 0)
      {
         theNumber /= smallPrimes[ii];
         power++;
      }

      if (power > 0)
      {
         factorList[distinctFactors] = smallPrimes[ii];
         powerList[distinctFactors] = power;
         distinctFactors++;
         
         if (theNumber < smallPrimes[ii] * smallPrimes[ii])
            break;
      }
   }

   // If then number wasn't fully factored, that's fine.
   if (theNumber > 1)
   {
      factorList[distinctFactors] = theNumber;
      powerList[distinctFactors] = 1;
      distinctFactors++;
   }

   return distinctFactors;
}
