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

#if (HAVE_HERMMATINV4X4)

#define ALIGN_PTR(p)      (void*)( ((uintptr_t)(p)+BBE_SIMD_WIDTH*2-1) & ~(BBE_SIMD_WIDTH*2-1) )
#define ALIGN_SIZE(sz)    ( (sz+BBE_SIMD_WIDTH*2-1) & ~(BBE_SIMD_WIDTH*2-1) )
#define sz_i16            sizeof(int16_t)

#define NORM_S            1
#define NORM_CA           1
#define ROUND_CA          0
#define NORM_AB           1
#define ROUND_AB          0
#define NORM_RECIP_ARG    0 // 0 for PM/EP config, 1 for LP

static const int16_t ALIGN(32) search_tbl[16] =
{
  (int16_t)0x1420,(int16_t)0x2840,(int16_t)0x3c60,(int16_t)0x28c5,(int16_t)0x3ce5,(int16_t)0x3d6a,(int16_t)0x0000,(int16_t)0x0000,
  (int16_t)0x5630,(int16_t)0x6a50,(int16_t)0x7e70,(int16_t)0x6ad5,(int16_t)0x7ef5,(int16_t)0x7f7a,(int16_t)0x0000,(int16_t)0x0000,
};

static const int16_t ALIGN(32) pattern_tbl[6*16] =
{
  (int16_t)0x0000,(int16_t)0x0101,(int16_t)0x0202,(int16_t)0x0303,(int16_t)0x0404,(int16_t)0x0505,(int16_t)0x0606,(int16_t)0x0707,
  (int16_t)0x0808,(int16_t)0x0909,(int16_t)0x0a0a,(int16_t)0x0b0b,(int16_t)0x0c0c,(int16_t)0x0d0d,(int16_t)0x0e0e,(int16_t)0x0f0f,
  (int16_t)0x0000,(int16_t)0x0202,(int16_t)0x0101,(int16_t)0x0303,(int16_t)0x0808,(int16_t)0x0a0a,(int16_t)0x0909,(int16_t)0x0b0b,
  (int16_t)0x0404,(int16_t)0x0606,(int16_t)0x0505,(int16_t)0x0707,(int16_t)0x0c0c,(int16_t)0x0e0e,(int16_t)0x0d0d,(int16_t)0x0f0f,
  (int16_t)0x0000,(int16_t)0x0203,(int16_t)0x0301,(int16_t)0x0102,(int16_t)0x080c,(int16_t)0x0a0f,(int16_t)0x0b0d,(int16_t)0x090e,
  (int16_t)0x0c04,(int16_t)0x0e07,(int16_t)0x0f05,(int16_t)0x0d06,(int16_t)0x0408,(int16_t)0x060b,(int16_t)0x0709,(int16_t)0x050a,
  (int16_t)0x0a05,(int16_t)0x0806,(int16_t)0x0904,(int16_t)0x0b07,(int16_t)0x0209,(int16_t)0x000a,(int16_t)0x0108,(int16_t)0x030b,
  (int16_t)0x0601,(int16_t)0x0402,(int16_t)0x0500,(int16_t)0x0703,(int16_t)0x0e0d,(int16_t)0x0c0e,(int16_t)0x0d0c,(int16_t)0x0f0f,
  (int16_t)0x0a05,(int16_t)0x0807,(int16_t)0x0b04,(int16_t)0x0906,(int16_t)0x020d,(int16_t)0x000f,(int16_t)0x030c,(int16_t)0x010e,
  (int16_t)0x0e01,(int16_t)0x0c03,(int16_t)0x0f00,(int16_t)0x0d02,(int16_t)0x0609,(int16_t)0x040b,(int16_t)0x0708,(int16_t)0x050a,
  (int16_t)0x0a0a,(int16_t)0x0b0b,(int16_t)0x0808,(int16_t)0x0909,(int16_t)0x0e0e,(int16_t)0x0f0f,(int16_t)0x0c0c,(int16_t)0x0d0d,
  (int16_t)0x0202,(int16_t)0x0303,(int16_t)0x0000,(int16_t)0x0101,(int16_t)0x0606,(int16_t)0x0707,(int16_t)0x0404,(int16_t)0x0505,
};

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
void chermmatinv4x4n ( void * restrict pScr, complex_fract16 * restrict A, int qA, int L )
{
  int16_t * restrict a_perm_ix; // Permutation indices for L matrices
  int16_t * restrict a_eX;      // Normalization shift amounts for L matrices
  int16_t * restrict a_qA;      // Fixed point position for L 2x2 subblocks A
  int16_t * restrict a_qT;      // Fixed point position for L 2x2 matrices T
  int16_t * restrict a_qS;      // Fixed point position for L 2x2 matrices S
  int16_t * restrict a_A;       // Upper left 2x2 subblock 
  int16_t * restrict a_B;       // Upper right 2x2 subblock
  int16_t * restrict a_C;       // Lower left 2x2 subblock 
  int16_t * restrict a_D;       // Lower right 2x2 subblock
  int16_t * restrict a_T;       // L 2x2 matrix products C*A^-1, A^-1*B

  int16_t * restrict X = (int16_t * restrict)A;    

  int L_;

  NASSERT_ALIGN( pScr, BBE_SIMD_WIDTH*2 );
  NASSERT_ALIGN( A   , BBE_SIMD_WIDTH*2 );

  if ( L<=0 ) return;

  //
  // Partition the scratch memory area.
  //

  // Number of input matrices rounded up to the next multiple of SIMD_WIDTH/2.
  L_ = L + ( -L & (BBE_SIMD_WIDTH/2-1) );

  {
    void * ptr = pScr;
    
    a_A       = (int16_t*)ALIGN_PTR( ptr ); ptr = a_A  + L_*4*2;
    a_B       = (int16_t*)ALIGN_PTR( ptr ); ptr = a_B  + L_*4*2;
    a_C       = (int16_t*)ALIGN_PTR( ptr ); ptr = a_C  + L_*4*2;
    a_D       = (int16_t*)ALIGN_PTR( ptr ); ptr = a_D  + L_*4*2;
    a_T       = (int16_t*)ALIGN_PTR( ptr ); ptr = a_T  + L_*4*2;
    a_qA      = (int16_t*)ALIGN_PTR( ptr ); ptr = a_qA + L_;
    a_qT      = (int16_t*)ALIGN_PTR( ptr ); ptr = a_qT + L_;
    a_qS      = (int16_t*)ALIGN_PTR( ptr ); ptr = a_qS + L_;
    a_eX      = (int16_t*)ALIGN_PTR( ptr ); ptr = a_eX + L_;
    a_perm_ix = (int16_t*)         ( ptr ); ptr = a_perm_ix + L_;

    NASSERT( (int8_t*)ptr - (int8_t*)pScr <= (int)chermmatinv4x4n_getScratchSize( L ) );
  }

  //----------------------------------------------------------------------------
  // Normalize input 4x4 matrices and search for a suboptimal permutation 
  // for each matrix, so that the 2x2 sub-block at the upper left of a
  // permuted matrix has a large determinant.
  //

  chermmatinv4x4n_searchPermutation( 
        X, a_eX,
        a_perm_ix, search_tbl,
        L );

  cmatinv4x4n_permute_b2bs( 
        a_A, a_B, a_C, a_D,
        X,
        a_perm_ix, pattern_tbl,
        L );

  //----------------------------------------------------------------------------
  // Invert the upper left block A, then compute and invert the Schur
  // complement of A: S = (D - C*A^-1*B)^-1

  // "A" <- A^-1
  // CQ(qA=(15-eA)) <- Q30/( CQ15 + eA )
  chermmatinv4x4n_inv2x2bs( X, a_A, a_qA, L_, 15 );

  // "T" <- C*A^-1
  // CQ(qT=(qA-2+eCA)) <- CQ15*CQ(qA) - ( 17 - eCA ) with/without rounding
  // Exponent is bounded to protect D*2^qT from overflow.
  if ( NORM_CA && !ROUND_CA ) cmatinv4x4n_mul2x2bs_n( a_T,a_qT,a_C,a_A,a_qA,L_,15,10 );
  else NASSERT( !"Option not implemented!" );

  // "D" <- D - C*A^-1*B
  // CQ(qS=(qT-3+eS)) <- {[ CQ15 + qT ] - CQ(qT)*CQ15} - ( 18 - eS ) w/ rounding
  if ( NORM_S ) cmatinv4x4n_mas2x2bs_nr( a_D, a_qS, a_T, a_B, a_qT, L_, 15 );
  else NASSERT( !"Option not implemented!" );

  // "D" <- S = ( D - C*A^-1*B )^-1
  // CQ(30-qS-eSi) <- Q30/( CQ(qS) + eSi )
  // qS <- 30-qS-eSi
  chermmatinv4x4n_inv2x2bs_q( X, a_D, a_qS, L_ );

  //----------------------------------------------------------------------------
  // Compute 2x2 blocks of the inverted matrix.

  // "C" <- -1*S*C*A^-1
  // CQ(qA) <- -1*CQ(qS)*CQ(qT)/2^(15-qA-eX) - ( qS + qT ) + qA
  chermmatinv4x4n_mul2x2bs_srn( a_C, a_B, a_D, a_T, a_eX, a_qS, a_qT, L_, 2*qA-15 );

  // "A" <- A^-1 + A^-1*C'*S*C*A^-1
  // CQ(qA+qT) <- [CQ(qA) - qA + qA + qT]/2^(15-qA-eX) + CQ(qT)*Q(qA)
  // CQ(qA) <- CQ(qA+qT) - qT w/ rounding
  chermmatinv4x4n_mas2x2bs_srj( a_A, a_T, a_C, a_qA, a_qT, a_eX, L_, qA, 15-qA );

  // "D" <- S
  // CQ(qA) <- CQ(qS)/2^(15-qA-eX) - qS + qA
  cmatinv4x4n_rnd2x2bs( a_D, a_eX, a_qS, L_, 15-2*qA );

  //----------------------------------------------------------------------------
  // Combine the inverse matrix from 2x2 sub-blocks and apply the inverse
  // permutation to compensate for initial reordering of the original matrix.

  cmatinv4x4n_permute_bs2b(
        X,
        a_A, a_B, a_C, a_D,
        a_perm_ix, pattern_tbl,
        L );
} /* chermmatinv4x4n() */

/* Return the scratch area size, in bytes. */
size_t chermmatinv4x4n_getScratchSize ( int L )
{
  int L_, sz = 0;

  // Number of input matrices rounded up to the next multiple of SIMD_WIDTH/2.
  L_ = L + ( -L & (BBE_SIMD_WIDTH/2-1) );

  if ( L>0 )
  {
    sz = 5*ALIGN_SIZE( L_*4*2*sz_i16 ) + // a_A, a_B, a_C, a_D, a_T 
         4*ALIGN_SIZE( L_*sz_i16     ) + // a_qA, a_qT, a_qS, a_eX
                     ( L_*sz_i16     );  // a_perm_ix
  }

  return ( (size_t)sz );
}

#else
DISCARD_FUN( void  , chermmatinv4x4n, ( void * restrict pScr, complex_fract16 * restrict A, int qA, int L ) )
size_t chermmatinv4x4n_getScratchSize  ( int L ) { (void)L; return 0;}
#endif /* HAVE_HERMMATINV4X4 */
