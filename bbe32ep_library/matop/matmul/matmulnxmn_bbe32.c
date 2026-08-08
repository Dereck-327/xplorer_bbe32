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
    Real Matrix-Matrix/Matrix-Vector Multiply
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"

/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#include "matmulnxmn_common.h"
#include <string.h>

/*-------------------------------------------------------------------------
Real Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of real matrices or vectors. Both the block order and streaming 
order are allowed for input/output matrix sequences.

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
x[L*Sx]     Sequence of left-hand input matrices
y[L*Sy]     Sequence of right-hand input matrices
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

#define MAX(x,y) ((x)>(y)?(x):(y))
/* get allocated space per one matrix (real) */
static int getSpaceR(int S)
{
  int m;
  /* compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl ) */
  m=30-XT_NSA(S);
  m=XT_MIN(m,LOG2_BBE_SIMD_WIDTH);
  /* round up to the  next multiple of 32 or lesser degree of 2 */
  S=(((S-1)>>m)+1)<<m;
  return S;
} /* getSpaceR() */

/* Block Order, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:  
     N,M must be multiples of 4
*/
void matmulnxmn ( void * pScr,
                  int16_t * restrict z,
            const int16_t * restrict x,
            const int16_t * restrict y,
            int N, int M, int L, int Q )
{
  int Sx, Sz;

  /* check restrictions */
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(M % 4 == 0);
  NASSERT(N % 4 == 0);

  if (L <= 0 || M <= 0) return;
  if (N <= 0)
  {
    memset(z, 0, M*M*L*sizeof(int16_t));
    return;
  }
  Sx = getSpaceR(N*M);
  Sz = getSpaceR(M*M);

  matmulnxmn_L16(pScr, z, x, y, N, M, L&~15, Q);
  x += Sx*(L&~15);
  y += Sx*(L&~15);
  z += Sz*(L&~15);
  L &= 15;
  if (L)
  {
    matmulnxmn_L4(pScr, z, x, y, N, M, L & 12, Q);
    x += Sx*(L & 12);
    y += Sx*(L & 12);
    z += Sz*(L & 12);
    L &= 3;
    if (L)
    {
      matmulnxmn_tail_L4(pScr, z, x, y, N, M, L, Q);
    }
  }
} /* matmulnxmn() */



/* Return the scratch area size, in bytes. */
size_t matmulnxmn_getScratchSize ( int N, int M, int L )
{
  size_t tail_scratch, L4_scratch, L16_scratch, sz = 0;
  if (M <= 0 || N <= 0) return 0;
  (void)L;
  if (L >= 16) { L16_scratch = matmulnxmn_L16_getScratchSize(N, M); sz = MAX(sz, L16_scratch); }
  L &= 15;
  if (L >= 4) { L4_scratch = matmulnxmn_L4_getScratchSize(N, M); sz = MAX(sz, L4_scratch); }
  L &= 3;
  if (L>0) { tail_scratch = matmulnxmn_tail_L4_getScratchSize(N, M); sz = MAX(sz, tail_scratch); }
  return sz;

} /* matmulnxmn_getScratchSize() */
