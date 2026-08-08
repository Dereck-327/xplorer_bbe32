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
    Matrix Hermitian Product
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"

/*-------------------------------------------------------------------------
Matrix Hermitian Product

Description: These functions left multiply each complex input matrix by its
conjugate transpose. The result is a Hermitian (or self-adjoint) matrix. Both
the block order and streaming order are allowed for input/output matrix
sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr      Scratch memory area. To determine the scratch area size required by
          a function <fun>, use the respective helper function 
          <fun>_getScratchSize()
Input:
x[L*Sx]   Complex input matrices
M         Matrix dimension 
N         Matrix dimension (columnar for MxN)
L         Number of matrices 
Q         Position of fractional point in matrix representation, 0..16
Output:
y[L*Sy]   Complex output matrices

Restrictions:
pScr,x,y  Aligned on 32-byte boundary
pScr,x,y  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Streaming Order, 2x2*2x2->2x2, Sx=4, Sy=4
   Restrictions:
     L must be a multiple of 8
*/
void cmatherm2x2s ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L, int Q )
{
  int k;

  xb_vecNx16 y00;
  xb_vecNx16 y01;
  xb_vecNx16 y02;
  xb_vecNx16 y03;
  xb_vecNx16 z00;
  vsaN q = BBE_MOVVSA32(Q);
  xb_vecNx40 z_out0;
  const xb_vecNx16 * restrict px00;
  xb_vecNx16 * restrict pz00;

  px00 = (const xb_vecNx16*)x;
  pz00 = (xb_vecNx16*)y;    
  /* check restrictions */
  NASSERT_ALIGN(x, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, BBE_SIMD_WIDTH);
  NASSERT(L % (BBE_SIMD_WIDTH / 2) == 0);

  if (L <= 0) return;
  __Pragma("ymemory( px00 )");
  __Pragma("loop_count min=1");
  for (k = 0; k<L; k += (BBE_SIMD_WIDTH >> 1))
  {
    /*Load input matrix X */
    BBE_LVNX16_XP(y00, px00, 2 * 2 * L);
    BBE_LVNX16_XP(y01, px00, 2 * 2 * L);
    BBE_LVNX16_XP(y02, px00, 2 * 2 * L);
    BBE_LVNX16_XP(y03, px00, -3 * 2 * 2 * L + 2 * BBE_SIMD_WIDTH);

    z_out0 = BBE_MULNX16J(y00, y00); /* x00*y00 */
    BBE_MULANX16J(z_out0, y02, y02); /* x02*y02 */

    /*Pack and save result */
    z00 = BBE_PACKVNX40(z_out0, q);
    BBE_SVNX16_XP(z00, pz00, 2 * 2 * L); 

    z_out0 = BBE_MULNX16J(y01, y00); /* x00*y01 */
    BBE_MULANX16J(z_out0, y03, y02); /* x02*y03 */

    /*Pack and save result */
    z00 = BBE_PACKVNX40(z_out0, q);
    BBE_SVNX16_XP(z00, pz00, 2 * 2 * L); 

    z_out0 = BBE_CONJNX40C(z_out0); /* x00*y00 */

    /*Pack and save result */
    z00 = BBE_PACKVNX40(z_out0, q);
    BBE_SVNX16_XP(z00, pz00, 2 * 2 * L); 

    z_out0 = BBE_MULNX16J(y01, y01); /* x00*y01 */
    BBE_MULANX16J(z_out0, y03, y03); /* x02*y03 */

    /*Pack and save result */
    z00 = BBE_PACKVNX40(z_out0, q);
    BBE_SVNX16_XP(z00, pz00, -3 * 2 * 2 * L + 2 * BBE_SIMD_WIDTH); 
  }
} /* cmatherm2x2s() */
