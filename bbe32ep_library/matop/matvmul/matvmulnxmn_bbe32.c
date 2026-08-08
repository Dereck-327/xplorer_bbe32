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
#include "matvmulnxmn_common.h"
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
#if !(HAVE_MULPC && HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, matvmulnxmn,( void *             pScr,
                  int16_t * restrict z,
            const int16_t * restrict x,
            const int16_t * restrict y,
            int N, int M, int L, int Q ))
size_t matvmulnxmn_getScratchSize ( int N, int M, int L ) {(void)N; (void)M; (void)L; return 0;};
#else
#define MAX(x,y) ((x)>(y)?(x):(y))


/* Block Order, MxN*Nx1->Mx1, Sx=MxN, Sy=(N>4)?((N+4)&~4):N, Sz=(M>4)?((M+4)&~4):M
   Restrictions:
     L must be a multiple of 4
     N, M must be multiples of 4
*/
void matvmulnxmn ( void * pScr,
                   int16_t * restrict z,
             const int16_t * restrict x,
             const int16_t * restrict y,
             int N, int M, int L, int Q )
{
  NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);

  NASSERT(!(L & 1) && !(M & 3) && !(N & 3));
  if (L <= 0 || M <= 0) return;
  if (N <= 0)
  {
    memset(z, 0, M*L*sizeof(int16_t));
    return;
  }
  if (M == 4)
  {
    if (N == 4)
    {
      matvmul4x4n(z, x, y, L, Q);
    }
    else if (N == 8)
    {
      matvmulnxmn_M4_N8_L4x(pScr, z, x, y, N, M, L, Q);
    }
    else
    {
      matvmulnxmn_M4_Ngte12_L4x(pScr, z, x, y, N, M, L, Q);
    }
  }
  else if (M == 8)
  {
    if (N == 4)
    {
      matvmulnxmn_M8_N4_L4x(pScr, z, x, y, N, M, L, Q);
    }
    else if (N == 8)
    {
      matvmul8x8n(z, x, y, L, Q);
    }
    else
    {
      matvmulnxmn_M8_Ngte12_L4x(pScr, z, x, y, N, M, L, Q);
    }
  }
  else
  {
    if (N == 4)
    {
      matvmulnxmn_Mgte12_N4_L4x(pScr, z, x, y, N, M, L, Q);
    }
    else if (N == 8)
    {
      matvmulnxmn_Mgte12_N8_L4x(pScr, z, x, y, N, M, L, Q);
    }
    else
    {
      int _N, _M, _L;

      _N = ((N + BBE_SIMD_WIDTH - 1) & ~(BBE_SIMD_WIDTH - 1));
      _M = ((M + BBE_SIMD_WIDTH - 1) & ~(BBE_SIMD_WIDTH - 1));

      if ((_L = (L&~15)))
      {
        matvmulnxmn_Mgte12_Ngte12_L16x(pScr,
          z, x, y,
          N, M, _L, Q);
      }

      if ((_L = (L & 15)))
      {
        matvmulnxmn_Mgte12_Ngte12_L4x(pScr,
          z + (L&~15)*_M,
          x + (L&~15)*M*N,
          y + (L&~15)*_N,
          N, M, _L, Q);
      }
    }
  }
} /* matvmulnxmn() */

/* Return the scratch area size, in bytes. */
size_t matvmulnxmn_getScratchSize ( int N, int M, int L )
{
  size_t M8_Ngte12_scratch, M4_Ngte12_scratch, Mgte12_N4_scratch, Mgte12_N8_scratch, sz = 0;
  size_t Mgte12_Ngte12_L16_scratch, Mgte12_Ngte12_L4_scratch;
  if (M <= 0 || N <= 0) return 0;
  if ((M == 4 && N == 4) || (M == 8 && N == 8)) return 0;
  if (M == 4 && N == 8) return 0;
  if (M == 8 && N == 4) return 0;
  (void)L;
  if (M == 4)
  {
    M4_Ngte12_scratch = matvmulnxmn_M4_Ngte12_getScratchSize(N, M); sz = MAX(sz, M4_Ngte12_scratch);
  }
  else if (M == 8)
  {
    M8_Ngte12_scratch = matvmulnxmn_M8_Ngte12_getScratchSize(N, M); sz = MAX(sz, M8_Ngte12_scratch);
  }
  else
  {
    if (N == 4)
    {
      Mgte12_N4_scratch = matvmulnxmn_Mgte12_N4_getScratchSize(N, M); sz = MAX(sz, Mgte12_N4_scratch);
    }
    else if (N == 8)
    {
      Mgte12_N8_scratch = matvmulnxmn_Mgte12_N8_getScratchSize(N, M); sz = MAX(sz, Mgte12_N8_scratch);
    }
    else
    {
      int _N, _M, _L;

      _N = ((N + BBE_SIMD_WIDTH - 1) & ~(BBE_SIMD_WIDTH - 1));
      _M = ((M + BBE_SIMD_WIDTH - 1) & ~(BBE_SIMD_WIDTH - 1));
      (void)_N,(void)_M;
      if ((_L = (L&~15)))
      {
        Mgte12_Ngte12_L16_scratch = matvmulnxmn_Mgte12_Ngte12_L16_getScratchSize(N, M); sz = MAX(sz, Mgte12_Ngte12_L16_scratch);
      }
      if ((_L = (L & 15)))
      {
        Mgte12_Ngte12_L4_scratch = matvmulnxmn_Mgte12_Ngte12_L4_getScratchSize(N, M); sz = MAX(sz, Mgte12_Ngte12_L4_scratch);
      }
    }
  }
  return sz;
} /* matvmulnxmn_getScratchSize() */
#endif
