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
    Real Matrix by Complex Matrix/Vector Multiply
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#include "rcmatmulnxmn_common.h"
#include <string.h>

#define MAX(x,y) ((x)>(y)?(x):(y))
/*-------------------------------------------------------------------------
Real Matrix by Complex Matrix/Vector Multiply 

Description: These functions perform pairwise multiplication of left-hand
real matrices by right-hand complex matrices or vectors. Both the block order
and streaming order are allowed for input/output matrix sequences.

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand real matrices
y[L*Sy]     Sequence of right-hand complex matrices or vectors
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
#if !(HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, rcmatmulnxmn,   (void* pScr,
                  complex_fract16 * restrict z, 
            const int16_t * restrict x, 
            const complex_fract16 * restrict y, 
            int N, int M, int L, int Q))
size_t rcmatmulnxmn_getScratchSize(int N, int M, int L) { (void)N; (void)M;  (void)L; return 0; }
#else


/* Block Order, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
      N,M must be multiples of 4
*/
void rcmatmulnxmn ( void * pScr,
                    complex_fract16 * restrict z, 
              const int16_t * restrict x, 
              const complex_fract16 * restrict y, 
              int N, int M, int L, int Q )
{
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(N % 4 == 0);
  NASSERT(M % 4 == 0);
  if (L <= 0 || M <= 0) return;
  if (N <= 0)
  {
    memset(z, 0, 2 * M*M*L*sizeof(int16_t));
    return;
  }
  if (N == 4) { rcmatmulnxmn_N4(pScr, z, x, y, M, L, Q); return; }
  if (M == 4) { rcmatmulnxmn_M4(pScr, z, x, y, N, L, Q); return; }
  rcmatmulnxmn_gen(pScr, z, x, y, N, M, L, Q);
} /* rcmatmulnxmn() */

/* Return the scratch area size, in bytes. */
size_t rcmatmulnxmn_getScratchSize ( int N, int M, int L )
{
  size_t gen_scratch, M4_scratch, N4_scratch, sz = 0;
  if (M <= 0 || N <= 0) return 0;
  (void)L;
  if (N == 4) { N4_scratch = rcmatmulnxmn_N4_getScratchSize(M); sz = MAX(sz, N4_scratch); }
  else
  {
    if (M == 4) { M4_scratch = rcmatmulnxmn_M4_getScratchSize(N); sz = MAX(sz, M4_scratch); }
    else
    {
      gen_scratch = rcmatmulnxmn_gen_getScratchSize(N, M); sz = MAX(sz, gen_scratch);
    }
  }
  return sz;
} /* rcmatmulnxmn_getScratchSize() */
#endif
