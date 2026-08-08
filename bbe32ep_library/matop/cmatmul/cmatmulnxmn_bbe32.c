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
#include "cmatmulnxmn_common.h"
#include <string.h>
#define MAX(x,y) ((x)>(y)?(x):(y))
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
#if !(HAVE_MULPC && HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, cmatmulnxmn,( void    *          pScr,
                  complex_fract16 * restrict z, 
            const complex_fract16 * restrict x, 
            const complex_fract16 * restrict y, 
            int N, int M, int L, int Q ))
size_t cmatmulnxmn_getScratchSize(int N, int M, int L) { (void)N; (void)M;  (void)L; return 0; }
#else
/* Block Order, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
     N,M must be multiples of 4
*/
void cmatmulnxmn ( void * pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q )
{
  int MN, MM, _L;

  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN(z, 2*BBE_SIMD_WIDTH );

  NASSERT( N%4==0 && M%4==0 );

  NASSERT( Q>=0 && Q<=16 );
  if (L <= 0 || M <= 0) return;
  if (N <= 0)
  {
    memset(z, 0, 2 * M*M*L*sizeof(int16_t));
    return;
  }

  MN = M*N; MM = M*M;

  if ( ( _L = (L&~7) ) )
  {
    cmatmulnxmn_L8( pScr,
                         z, x, y,
                         N, M, _L, Q );
  }

  if ( ( _L = (L&6) ) )
  {
    cmatmulnxmn_L2( pScr,
                         z + (L&~7)*MM,
                         x + (L&~7)*MN,
                         y + (L&~7)*MN,
                         N, M, _L, Q );
  }

  if ( ( _L = (L&1) ) )
  {
    cmatmulnxmn_tail( pScr,
                         z + (L&~1)*MM,
                         x + (L&~1)*MN,
                         y + (L&~1)*MN,
                         N, M, Q );
  }
} /* cmatmulnxmn() */

/* Return the scratch area size, in bytes. */
size_t cmatmulnxmn_getScratchSize ( int N, int M, int L )
{
  size_t L2_scratch, L8_scratch, sz = 0;
  if (M <= 0 || N <= 0) return 0;
  (void)L;
  if (L >= 8) { L8_scratch = cmatmulnxmn_L8_getScratchSize(N, M); sz = MAX(sz, L8_scratch); }
  L &= 7;
  if (L>0) { L2_scratch = cmatmulnxmn_L2_getScratchSize(N, M); sz = MAX(sz, L2_scratch); }
  return sz;
} /* cmatmulnxmn_getScratchSize() */
#endif
