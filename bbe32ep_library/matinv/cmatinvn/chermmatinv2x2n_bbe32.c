/* ------------------------------------------------------------------------ */
/* Copyright (c) 2017 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
  NatureDSP_Baseband library. Direct Matrix Inversion
    Direct Matrix Inversion
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "cmatinvn_common.h"

/*-------------------------------------------------------------------------
Direct Matrix Inversion For Hermitian Matrices

Description: perform in-place inversion of 2x2 and 4x4 complex Hermitian 
matrices. 2x2 matrices are inverted by Cramer's rule. For 4x4 matrices we
employ the blockwise inversion algorithm encompassed with a suboptimal 
row/column permutation that gains better conditioning of the block structure.
Due to special properties of Hermitian matrices, the complexity of direct
matrix inversion is much lower than for a generic matrix of the same size.

Data format: block ordered 16-bit fixed-point data

NOTES:
1. In general, accuracy of a matrix inversion algorithm implementation is a
   function of input matrix condition number. Thus it is user's responsibility
   to qualify the reliability of numeric results. Refer to NatureDSP Baseband 
   Library Reference for details. 
2. For blockwise inversion of 4x4 fixed-point matrices, it is reasonable to
   limit the dynamic range of input data by 11..13 significant bits. This
   measure reduces the possibility of an overflow at internal computations.

Parameters:
Temporary:
pScr       Scratch memory area. To determine the scratch area size required by
           a function <fun>, use the respective helper function 
           <fun>_getScratchSize()
Input:
L          Number of matrices
qA         Number of fractional bits for input/output data
Input/Output:
A[L][N*N]  Sequence of L NxN complex input/result matrices

Restrictions:
pScr,A     Aligned on 32-byte boundary

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

#define STEP  (BBE_SIMD_WIDTH/2)

#if !HAVE_HERMMATINV2X2
DISCARD_FUN(void,chermmatinv2x2n,( void * restrict pScr, complex_fract16 * restrict A, int qA, int L ))
#else

void chermmatinv2x2n(void * restrict pScr, complex_fract16 * restrict A, int qA, int L)
{
  const xb_vecNx16 *          A_rd;
        xb_vecNx16 * restrict A_wr;

  xb_vecNx16 a0, a1, a2, a3;
  xb_vecNx16 b0, b1, b2, b3;
  xb_vecNx40 w0, w1, w3;
  xb_vecNx16 t0, t1;
  xb_vecNx16 rd0;
  xb_vecNx16 d_exp;
  xb_vecNx16 c53m2q;

  vboolN vb0;
  vsaN   vsa0, vsa_c23;

  int l;

  NASSERT_ALIGN16( A    );
  NASSERT_ALIGN16( pScr );

  t0 = BBE_MOVVINT16( 23 );

  vsa_c23 = BBE_MOVVSV( t0, 0 );

  c53m2q = BBE_MOVVA16( 53 - 2*qA );

  A_rd = (const xb_vecNx16*)A;
  A_wr = (      xb_vecNx16*)A;

  for ( l=0; l<L/STEP; l++ )
  {
    // Load 16 2x2 Hermitian matrices. Perform formal fixed point conversion:
    // CQ15 <- CQ(qA)/2^(15-qA) + 15 - qA
    BBE_LVNX16_IP( a0, A_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( a1, A_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( a2, A_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( a3, A_rd, +4*BBE_SIMD_WIDTH/2 );

    // Convert data to streaming format by 16x4 -> 4x16 transposition.
    BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_DEINTERLEAVE_2 );
    BBE_DSELNX16I( a3, a2, a3, a2, BBE_DSELI_DEINTERLEAVE_2 );
    BBE_DSELNX16I( a3, a1, a3, a1, BBE_DSELI_DEINTERLEAVE_2 );
    BBE_DSELNX16I( a2, a0, a2, a0, BBE_DSELI_DEINTERLEAVE_2 );

    //
    // Compute the (real) determinant of each matrix. Results are placed to
    // even elements, odd elements hold zero.
    //

    // Diagonal elements are real.
    // Q1.30 <- Q15*Q15
    w0 = BBE_MULNX16( a0, a3 );
    // Off-diagonal elements are adjoint.
    // Q2.30 <- Q1.30 - |CQ15|^2
    BBE_MULSNX16J( w0, a1, a1 );

    //
    // Take absolute value.
    //

    w1 = BBE_MOVWINX40( BBE_MOVWI_ZERO );

    // Remember the sign.
    vb0 = BBE_GTNX40( w1, w0 );

    w0 = BBE_ABSNX40( w0 );

    //
    // Perform 40-bit normalization of determinant values.
    //

    vsa0  = BBE_NSANX40( w0 );
    d_exp = BBE_MOVVVS( vsa0 );

    // Q(2.30+d_exp) <- Q2.30 + d_exp
    w0 = BBE_SLLNX40( w0, vsa0 );

    //
    // Compute reciprocal determinant values (16 even elements)
    //

    // w1: Q(46-d_exp) <- 2^76/Q(30+d_exp)
    // t0: Q(32-2*d_exp) <- 2^92/Q(30+d_exp)^2; 16-bit unsigned
    // t1: Q(14+d_exp) <- Q(30+d_exp) - 16
    BBE_RECIPLUNX40_0( w0,t1,  t0, w0 );

    // Q(46-d_exp) <- Q(46-d_exp) + Q(32-2*d_exp)*Q(14+d_exp)
    BBE_MULUSANX16( w0, t0, t1 );

    w0 = BBE_RNDSADJNX40( w0, vsa_c23 );

    // Q(23-d_exp) <- Q(46-d_exp) - 23 w/ rounding and saturation.
    rd0 = BBE_PACKVNX40( w0, vsa_c23 );

    // Restore the sign.
    BBE_NEGSNX16T( rd0, rd0, vb0 );

    // Duplicate the result in odd positions.
    rd0 = BBE_SHFLNX16I( rd0, BBE_SHFLI_DUPLICATE_1_EVEN );

    //
    // Compute the inverse matrix: invA = 1/d*[A(2,2),-A(1,2);-A(2,1),A(1,1)]
    //

    t0 = BBE_SHFLNX16I( d_exp, BBE_SHFLI_DUPLICATE_1_EVEN );

    t0 = BBE_SUBNX16( c53m2q, t0 );

    vsa0 = BBE_MOVVSV( t0, 0 );

    b0 = a3;
    b1 = BBE_NEGSNX16( a1 );
    b3 = a0;

    // CQ(38-d_exp) <- CQ15*Q(23-d_exp)
    w0 = BBE_MULRNX16( rd0, b0, vsa0 );
    w1 = BBE_MULRNX16( rd0, b1, vsa0 );
    w3 = BBE_MULRNX16( rd0, b3, vsa0 );

    // CQ(qA) <- CQ(38-d_exp)/2^(15-qA) - 38 + d_exp + qA w/ rounding
    b0 = BBE_PACKVNX40( w0, vsa0 );
    b1 = BBE_PACKVNX40( w1, vsa0 );
    b3 = BBE_PACKVNX40( w3, vsa0 );

    // Off-diagonal elements are adjoint.
    b2 = BBE_CONJSNX16C( b1 );

    //
    // Convert data to block format by 4x16 -> 16x4 transposition.
    //

    BBE_DSELNX16I( a1, a0, b1, b0, BBE_DSELI_INTERLEAVE_2 );
    BBE_DSELNX16I( a3, a2, b3, b2, BBE_DSELI_INTERLEAVE_2 );

    // We might use BBE_DSELI_INTERLEAVE4 here, but the schedule appears better
    // with a pair of selects.
    b0 = BBE_SELNX16I( a2, a0, BBE_SELI_INTERLEAVE_4_LO );
    b1 = BBE_SELNX16I( a2, a0, BBE_SELI_INTERLEAVE_4_HI );
    b2 = BBE_SELNX16I( a3, a1, BBE_SELI_INTERLEAVE_4_LO );
    b3 = BBE_SELNX16I( a3, a1, BBE_SELI_INTERLEAVE_4_HI );

    // Save 16 inverted 2x2 matrices.
    BBE_SVNX16_IP( b0, A_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( b1, A_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( b2, A_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( b3, A_wr, +4*BBE_SIMD_WIDTH/2 );
  }

  if (L > 0 && (L&(STEP - 1)))
  {
    valign A_va;

    int n = 16*( L & (STEP-1) );

    A_va = BBE_LAVNX16_PP( A_rd );

    // Read ( L & (STEP-1) ) 2x2 Hermitian matrices.
    // CQ15 <- CQ(qA)/2^(15-qA) + 15 - qA
    BBE_LAVNX16_XP( a0, A_va, A_rd, n ); n -= BBE_SIMD_WIDTH;
    BBE_LAVNX16_XP( a1, A_va, A_rd, n ); n -= BBE_SIMD_WIDTH;
    //BBE_LAVNX16_XP(a2, A_va, A_rd, n);
    if (n > 0) { BBE_LAVNX16_XP(a2, A_va, A_rd, n);} else a2 = BBE_ZERONX16();
    a3 = BBE_MOVVINX16( BBE_MOVVI_ZERO );

    // Convert data to streaming format by 16x4 -> 4x16 transposition.
    BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_DEINTERLEAVE_2 );
    BBE_DSELNX16I( a3, a2, a3, a2, BBE_DSELI_DEINTERLEAVE_2 );
    BBE_DSELNX16I( a3, a1, a3, a1, BBE_DSELI_DEINTERLEAVE_2 );
    BBE_DSELNX16I( a2, a0, a2, a0, BBE_DSELI_DEINTERLEAVE_2 );

    //
    // Compute the (real) determinant of each matrix. Results are placed to
    // even elements, odd elements hold zero.
    //

    // Diagonal elements are real.
    // Q1.30 <- Q15*Q15
    w0 = BBE_MULNX16( a0, a3 );
    // Off-diagonal elements are adjoint.
    // Q2.30 <- Q1.30 - |CQ15|^2
    BBE_MULSNX16J( w0, a1, a1 );

    //
    // Take absolute value.
    //

    w1 = BBE_MOVWINX40( BBE_MOVWI_ZERO );

    // Remember the sign.
    vb0 = BBE_GTNX40( w1, w0 );

    w0 = BBE_ABSNX40( w0 );

    //
    // Perform 40-bit normalization of determinant values.
    //

    vsa0  = BBE_NSANX40( w0 );
    d_exp = BBE_MOVVVS( vsa0 );

    // Q(2.30+d_exp) <- Q2.30 + d_exp
    w0 = BBE_SLLNX40( w0, vsa0 );

    //
    // Compute reciprocal determinant values (16 even elements)
    //

    // w1: Q(46-d_exp) <- 2^76/Q(30+d_exp)
    // t0: Q(32-2*d_exp) <- 2^92/Q(30+d_exp)^2; 16-bit unsigned
    // t1: Q(14+d_exp) <- Q(30+d_exp) - 16
    BBE_RECIPLUNX40_0(w0,  t1, t0, w0 );

    // Q(46-d_exp) <- Q(46-d_exp) + Q(32-2*d_exp)*Q(14+d_exp)
    BBE_MULUSANX16( w0, t0, t1 );

    w0 = BBE_RNDSADJNX40( w0, vsa_c23 );

    // Q(23-d_exp) <- Q(46-d_exp) - 23 w/ rounding and saturation.
    rd0 = BBE_PACKVNX40( w0, vsa_c23 );

    // Restore the sign.
    BBE_NEGSNX16T( rd0, rd0, vb0 );

    // Duplicate the result in odd positions.
    rd0 = BBE_SHFLNX16I( rd0, BBE_SHFLI_DUPLICATE_1_EVEN );

    //
    // Compute the inverse matrix: invA = 1/d*[A(2,2),-A(1,2);-A(2,1),A(1,1)]
    //

    t0 = BBE_SHFLNX16I( d_exp, BBE_SHFLI_DUPLICATE_1_EVEN );

    t0 = BBE_SUBNX16( c53m2q, t0 );

    vsa0 = BBE_MOVVSV( t0, 0 );

    b0 = a3;
    b1 = BBE_NEGSNX16( a1 );
    b3 = a0;

    // CQ(38-d_exp) <- CQ15*Q(23-d_exp)
    w0 = BBE_MULRNX16( rd0, b0, vsa0 );
    w1 = BBE_MULRNX16( rd0, b1, vsa0 );
    w3 = BBE_MULRNX16( rd0, b3, vsa0 );

    // CQ(qA) <- CQ(38-d_exp)/2^(15-qA) - 38 + d_exp + qA w/ rounding
    b0 = BBE_PACKVNX40( w0, vsa0 );
    b1 = BBE_PACKVNX40( w1, vsa0 );
    b3 = BBE_PACKVNX40( w3, vsa0 );

    // Off-diagonal elements are adjoint.
    b2 = BBE_CONJSNX16C( b1 );

    //
    // Convert data to block format by 4x16 -> 16x4 transposition.
    //

    BBE_DSELNX16I( a1, a0, b1, b0, BBE_DSELI_INTERLEAVE_2 );
    BBE_DSELNX16I( a3, a2, b3, b2, BBE_DSELI_INTERLEAVE_2 );

    BBE_DSELNX16I( b1, b0, a2, a0, BBE_DSELI_INTERLEAVE_4 );
    BBE_DSELNX16I( b3, b2, a3, a1, BBE_DSELI_INTERLEAVE_4 );

    //
    // Save ( L & (STEP-1) ) inverted 2x2 matrices.
    //

    n = 16*( L & (STEP-1) );

    A_va = BBE_ZALIGN();

    BBE_SAVNX16_XP( b0, A_va, A_wr, n ); n -= 2*BBE_SIMD_WIDTH;
    BBE_SAVNX16_XP( b1, A_va, A_wr, n ); n -= 2*BBE_SIMD_WIDTH;
    BBE_SAVNX16_XP( b2, A_va, A_wr, n );

    BBE_SAVNX16POS_FP( A_va, A_wr );
  }

} /* chermmatinv2x2n() */
#endif

/* Return the scratch area size, in bytes. */
size_t chermmatinv2x2n_getScratchSize ( int L )
{
  (void)L;
  return (0);
}
