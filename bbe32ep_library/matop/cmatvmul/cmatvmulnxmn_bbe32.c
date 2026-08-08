/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs ('Cadence    */
/* Libraries') are the copyrighted works of Cadence Design Systems Inc.     */
/* Cadence IP is licensed for use with Cadence processor cores only and     */
/* must not be used for any other processors and platforms. Your use of the */
/* Cadence Libraries is subject to the terms of the license agreement you   */
/* have entered into with Cadence Design Systems, or a sublicense granted   */
/* to you by a direct Cadence licensee.                                     */
/* ------------------------------------------------------------------------ */
/*  IntegrIT, Ltd.   www.integrIT.com, info@integrIT.com                    */
/*                                                                          */
/* NatureDSP_Baseband Library                                               */
/*                                                                          */
/* This library contains copyrighted materials, trade secrets and other     */
/* proprietary information of IntegrIT, Ltd. This software is licensed for  */
/* use with Cadence processor cores only and must not be used for any other */
/* processors and platforms. The license to use these sources was given to  */
/* Cadence, Inc. under Terms and Condition of a Software License Agreement  */
/* between Cadence, Inc. and IntegrIT, Ltd.                                 */
/* ------------------------------------------------------------------------ */
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
  NatureDSP_Baseband library. Matrix Operations
    Complex Matrix-Matrix/Matrix-Vector Multiply
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#include "cmatvmulnxmn_common.h"
#include <string.h>
/*-------------------------------------------------------------------------
Complex Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of complex matrices or vectors. Both the block order and 
streaming order are allowed for input/output matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand complex matrices
y[L*Sy]     Sequence of right-hand complex matrices
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices 
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of complex result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/
#if !(HAVE_MULPC && 1)
DISCARD_FUN(void, cmatvmulnxmn,( void * pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q ))
size_t cmatvmulnxmn_getScratchSize(int N, int M, int L) { (void)N; (void)M; (void)L; return 0; }
#else

/* Return allocated space per one matrix (for block format), in 16-bit words. */
static int getSpace(int S, int isCplx)
{
  int m,maxM;
  /* compute multiple of next degree of 2 (max multiple is 16 for real, 
   8 for complex, 8 for real+dbl, 4 for complex+dbl ) */
  m = 30 - XT_NSA( S );
  maxM = 4 - ( isCplx != 0 );
  if (m>maxM) m=maxM;
  /* round up to the next multiple of 16 or smaller degree of 2 */
  S=(((S-1)>>m)+1)<<m;
  return S;
}


/* Block Order, MxN*Nx1->Mx1, Sx=MxN, Sy=N, Sz=M
   Restrictions:
     L must be even
     N,M must be multiples of 4
*/
void cmatvmulnxmn ( void * pScr,
                    complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int N, int M, int L, int Q )
{
  int MN = M*N;;
  NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);

  NASSERT(!(L & 1) && !(M & 3) && !(N & 3));
  if (L <= 0 || M <= 0) return;
  if (N <= 0)
  {
    memset(z, 0, 2 * M*L*sizeof(int16_t));
    return;
  }
  if (M==4&&N==4)
  {
    cmatvmul4x4n(z, x, y, L, Q);
  }
  else 
  {
    if (M>=8)
    {
      if (N>=8)
      {
        cmatvmulnxmn_Mgte8_Ngte8_L8x(pScr, z, x, y, N, M, (L&~7), Q);
        cmatvmulnxmn_Mgte8_Ngte8_L2x(pScr, z + (L&~7) * getSpace(M, 1),
                                      x + (L&~7) *  MN, 
                                      y + (L&~7) * getSpace(N, 1),
                                      N, M, (L & 7), Q);
      }
      else
      {
        cmatvmulnxmn_Mgte8_Nlt8_L8x(pScr, z, x, y, N, M, (L&~7), Q);
        cmatvmulnxmn_Mgte8_Nlt8_L2x(pScr, z + (L&~7) * getSpace(M, 1),
                                     x + (L&~7) * MN,
                                     y + (L&~7) * getSpace(N, 1),
                                     N, M, (L & 7), Q);
      }
    }
    else
    {
      cmatvmulnxmn_Mlt8_Ngte8_L8x(pScr, z, x, y, N, M, (L&~7), Q);
      cmatvmulnxmn_Mlt8_Ngte8_L2x(pScr, z + (L&~7) * getSpace(M, 1),
                                   x + (L&~7) * MN,
                                   y + (L&~7) * getSpace(N, 1),
                                   N, M, (L & 7), Q);
    }
  }
} /* cmatvmulnxmn() */

#define MAX(x,y) ((x)>(y)?(x):(y))
/* Return the scratch area size, in bytes. */
size_t cmatvmulnxmn_getScratchSize ( int N, int M, int L )
{
  size_t  sz = 0, scr_size;
  if (M <= 0 || N <= 0) return 0;
  if (M == 4 && N == 4) return 0;
  (void)L;
  if (M >= 8)
  {
    if (N >= 8)
    {
      scr_size = cmatvmulnxmn_Mgte8_Ngte8_L8x_getScratchSize(N, M); sz = MAX(scr_size, sz);
      scr_size = cmatvmulnxmn_Mgte8_Ngte8_L2x_getScratchSize(N, M); sz = MAX(scr_size, sz);
    }
    else
    {
      scr_size = cmatvmulnxmn_Mgte8_Nlt8_L8x_getScratchSize(N, M); sz = MAX(scr_size, sz);
      scr_size = cmatvmulnxmn_Mgte8_Nlt8_L2x_getScratchSize(N, M); sz = MAX(scr_size, sz);
    }
  }
  else
  {
    scr_size = cmatvmulnxmn_Mlt8_Ngte8_L8x_getScratchSize(N, M); sz = MAX(scr_size, sz);
    scr_size = cmatvmulnxmn_Mlt8_Ngte8_L2x_getScratchSize(N, M); sz = MAX(scr_size, sz);
  }
  return sz;
} /* cmatvmulnxmn_getScratchSize() */
#endif
