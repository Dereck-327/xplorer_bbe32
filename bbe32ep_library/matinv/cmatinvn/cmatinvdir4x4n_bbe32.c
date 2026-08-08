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

#if HAVE_MATINV4X4

#define ALIGN_PTR(p)      (void*)( ((uintptr_t)(p)+BBE_SIMD_WIDTH*2-1) & ~(BBE_SIMD_WIDTH*2-1) )
#define ALIGN_SIZE(sz)    ( (sz+BBE_SIMD_WIDTH*2-1) & ~(BBE_SIMD_WIDTH*2-1) )
#define sz_i16            sizeof(int16_t)

#define NORM_S    1
#define NORM_CA   1
#define ROUND_CA  0
#define NORM_AB   1
#define ROUND_AB  0

static const int16_t ALIGN(32) search_tbl[16*16] =
{
  (int16_t)0x5410,(int16_t)0x6420,(int16_t)0x7430,(int16_t)0x9810,(int16_t)0xa820,(int16_t)0xb830,(int16_t)0xdc10,(int16_t)0xec20,
  (int16_t)0x0000,(int16_t)0x0020,(int16_t)0x0040,(int16_t)0x00c0,(int16_t)0x00e0,(int16_t)0x0100,(int16_t)0x0180,(int16_t)0x01a0,
  (int16_t)0x5410,(int16_t)0x6521,(int16_t)0x7531,(int16_t)0x9810,(int16_t)0xa921,(int16_t)0xb931,(int16_t)0xdc10,(int16_t)0xed21,
  (int16_t)0x0000,(int16_t)0x0060,(int16_t)0x0080,(int16_t)0x00c0,(int16_t)0x0120,(int16_t)0x0140,(int16_t)0x0180,(int16_t)0x01e0,
  (int16_t)0x6420,(int16_t)0x6521,(int16_t)0x7632,(int16_t)0xa820,(int16_t)0xa921,(int16_t)0xba32,(int16_t)0xec20,(int16_t)0xed21,
  (int16_t)0x0020,(int16_t)0x0060,(int16_t)0x00a0,(int16_t)0x00e0,(int16_t)0x0120,(int16_t)0x0160,(int16_t)0x01a0,(int16_t)0x01e0,
  (int16_t)0x7430,(int16_t)0x7531,(int16_t)0x7632,(int16_t)0xb830,(int16_t)0xb931,(int16_t)0xba32,(int16_t)0xfc30,(int16_t)0xfd31,
  (int16_t)0x0040,(int16_t)0x0080,(int16_t)0x00a0,(int16_t)0x0100,(int16_t)0x0140,(int16_t)0x0160,(int16_t)0x01c0,(int16_t)0x0200,
  (int16_t)0x5410,(int16_t)0x6420,(int16_t)0x7430,(int16_t)0x9854,(int16_t)0xa864,(int16_t)0xb874,(int16_t)0xdc54,(int16_t)0xec64,
  (int16_t)0x0000,(int16_t)0x0020,(int16_t)0x0040,(int16_t)0x0240,(int16_t)0x0260,(int16_t)0x0280,(int16_t)0x0300,(int16_t)0x0320,
  (int16_t)0x5410,(int16_t)0x6521,(int16_t)0x7531,(int16_t)0x9854,(int16_t)0xa965,(int16_t)0xb975,(int16_t)0xdc54,(int16_t)0xed65,
  (int16_t)0x0000,(int16_t)0x0060,(int16_t)0x0080,(int16_t)0x0240,(int16_t)0x02a0,(int16_t)0x02c0,(int16_t)0x0300,(int16_t)0x0360,
  (int16_t)0x6420,(int16_t)0x6521,(int16_t)0x7632,(int16_t)0xa864,(int16_t)0xa965,(int16_t)0xba76,(int16_t)0xec64,(int16_t)0xed65,
  (int16_t)0x0020,(int16_t)0x0060,(int16_t)0x00a0,(int16_t)0x0260,(int16_t)0x02a0,(int16_t)0x02e0,(int16_t)0x0320,(int16_t)0x0360,
  (int16_t)0x7430,(int16_t)0x7531,(int16_t)0x7632,(int16_t)0xb874,(int16_t)0xb975,(int16_t)0xba76,(int16_t)0xfc74,(int16_t)0xfd75,
  (int16_t)0x0040,(int16_t)0x0080,(int16_t)0x00a0,(int16_t)0x0280,(int16_t)0x02c0,(int16_t)0x02e0,(int16_t)0x0340,(int16_t)0x0380,
  (int16_t)0x9810,(int16_t)0xa820,(int16_t)0xb830,(int16_t)0x9854,(int16_t)0xa864,(int16_t)0xb874,(int16_t)0xdc98,(int16_t)0xeca8,
  (int16_t)0x00c0,(int16_t)0x00e0,(int16_t)0x0100,(int16_t)0x0240,(int16_t)0x0260,(int16_t)0x0280,(int16_t)0x03c0,(int16_t)0x03e0,
  (int16_t)0x9810,(int16_t)0xa921,(int16_t)0xb931,(int16_t)0x9854,(int16_t)0xa965,(int16_t)0xb975,(int16_t)0xdc98,(int16_t)0xeda9,
  (int16_t)0x00c0,(int16_t)0x0120,(int16_t)0x0140,(int16_t)0x0240,(int16_t)0x02a0,(int16_t)0x02c0,(int16_t)0x03c0,(int16_t)0x0420,
  (int16_t)0xa820,(int16_t)0xa921,(int16_t)0xba32,(int16_t)0xa864,(int16_t)0xa965,(int16_t)0xba76,(int16_t)0xeca8,(int16_t)0xeda9,
  (int16_t)0x00e0,(int16_t)0x0120,(int16_t)0x0160,(int16_t)0x0260,(int16_t)0x02a0,(int16_t)0x02e0,(int16_t)0x03e0,(int16_t)0x0420,
  (int16_t)0xb830,(int16_t)0xb931,(int16_t)0xba32,(int16_t)0xb874,(int16_t)0xb975,(int16_t)0xba76,(int16_t)0xfcb8,(int16_t)0xfdb9,
  (int16_t)0x0100,(int16_t)0x0140,(int16_t)0x0160,(int16_t)0x0280,(int16_t)0x02c0,(int16_t)0x02e0,(int16_t)0x0400,(int16_t)0x0440,
  (int16_t)0xdc10,(int16_t)0xec20,(int16_t)0xfc30,(int16_t)0xdc54,(int16_t)0xec64,(int16_t)0xfc74,(int16_t)0xdc98,(int16_t)0xeca8,
  (int16_t)0x0180,(int16_t)0x01a0,(int16_t)0x01c0,(int16_t)0x0300,(int16_t)0x0320,(int16_t)0x0340,(int16_t)0x03c0,(int16_t)0x03e0,
  (int16_t)0xdc10,(int16_t)0xed21,(int16_t)0xfd31,(int16_t)0xdc54,(int16_t)0xed65,(int16_t)0xfd75,(int16_t)0xdc98,(int16_t)0xeda9,
  (int16_t)0x0180,(int16_t)0x01e0,(int16_t)0x0200,(int16_t)0x0300,(int16_t)0x0360,(int16_t)0x0380,(int16_t)0x03c0,(int16_t)0x0420,
  (int16_t)0xec20,(int16_t)0xed21,(int16_t)0xfe32,(int16_t)0xec64,(int16_t)0xed65,(int16_t)0xfe76,(int16_t)0xeca8,(int16_t)0xeda9,
  (int16_t)0x01a0,(int16_t)0x01e0,(int16_t)0x0220,(int16_t)0x0320,(int16_t)0x0360,(int16_t)0x03a0,(int16_t)0x03e0,(int16_t)0x0420,
  (int16_t)0xfc30,(int16_t)0xfd31,(int16_t)0xfe32,(int16_t)0xfc74,(int16_t)0xfd75,(int16_t)0xfe76,(int16_t)0xfcb8,(int16_t)0xfdb9,
  (int16_t)0x01c0,(int16_t)0x0200,(int16_t)0x0220,(int16_t)0x0340,(int16_t)0x0380,(int16_t)0x03a0,(int16_t)0x0400,(int16_t)0x0440
};

static const int16_t ALIGN(32) pattern_tbl[35*16] =
{
  0x0000,0x0101,0x0202,0x0303,0x0404,0x0505,0x0606,0x0707,0x0808,0x0909,0x0a0a,0x0b0b,0x0c0c,0x0d0d,0x0e0e,0x0f0f,
  0x0000,0x0102,0x0201,0x0303,0x0804,0x0906,0x0a05,0x0b07,0x0408,0x050a,0x0609,0x070b,0x0c0c,0x0d0e,0x0e0d,0x0f0f,
  0x0000,0x0103,0x0201,0x0302,0x0804,0x0907,0x0a05,0x0b06,0x0c08,0x0d0b,0x0e09,0x0f0a,0x040c,0x050f,0x060d,0x070e,
  0x0801,0x0902,0x0a00,0x0b03,0x0005,0x0106,0x0204,0x0307,0x0409,0x050a,0x0608,0x070b,0x0c0d,0x0d0e,0x0e0c,0x0f0f,
  0x0801,0x0903,0x0a00,0x0b02,0x0005,0x0107,0x0204,0x0306,0x0c09,0x0d0b,0x0e08,0x0f0a,0x040d,0x050f,0x060c,0x070e,
  0x0802,0x0903,0x0a00,0x0b01,0x0c06,0x0d07,0x0e04,0x0f05,0x000a,0x010b,0x0208,0x0309,0x040e,0x050f,0x060c,0x070d,
  0x0000,0x0201,0x0102,0x0303,0x0408,0x0609,0x050a,0x070b,0x0804,0x0a05,0x0906,0x0b07,0x0c0c,0x0e0d,0x0d0e,0x0f0f,
  0x0000,0x0202,0x0101,0x0303,0x0808,0x0a0a,0x0909,0x0b0b,0x0404,0x0606,0x0505,0x0707,0x0c0c,0x0e0e,0x0d0d,0x0f0f,
  0x0000,0x0203,0x0101,0x0302,0x0808,0x0a0b,0x0909,0x0b0a,0x0c04,0x0e07,0x0d05,0x0f06,0x040c,0x060f,0x050d,0x070e,
  0x0801,0x0a02,0x0900,0x0b03,0x0009,0x020a,0x0108,0x030b,0x0405,0x0606,0x0504,0x0707,0x0c0d,0x0e0e,0x0d0c,0x0f0f,
  0x0801,0x0a03,0x0900,0x0b02,0x0009,0x020b,0x0108,0x030a,0x0c05,0x0e07,0x0d04,0x0f06,0x040d,0x060f,0x050c,0x070e,
  0x0802,0x0a03,0x0900,0x0b01,0x0c0a,0x0e0b,0x0d08,0x0f09,0x0006,0x0207,0x0104,0x0305,0x040e,0x060f,0x050c,0x070d,
  0x0000,0x0201,0x0302,0x0103,0x040c,0x060d,0x070e,0x050f,0x0804,0x0a05,0x0b06,0x0907,0x0c08,0x0e09,0x0f0a,0x0d0b,
  0x0000,0x0202,0x0301,0x0103,0x080c,0x0a0e,0x0b0d,0x090f,0x0404,0x0606,0x0705,0x0507,0x0c08,0x0e0a,0x0f09,0x0d0b,
  0x0000,0x0203,0x0301,0x0102,0x080c,0x0a0f,0x0b0d,0x090e,0x0c04,0x0e07,0x0f05,0x0d06,0x0408,0x060b,0x0709,0x050a,
  0x0801,0x0a02,0x0b00,0x0903,0x000d,0x020e,0x030c,0x010f,0x0405,0x0606,0x0704,0x0507,0x0c09,0x0e0a,0x0f08,0x0d0b,
  0x0801,0x0a03,0x0b00,0x0902,0x000d,0x020f,0x030c,0x010e,0x0c05,0x0e07,0x0f04,0x0d06,0x0409,0x060b,0x0708,0x050a,
  0x0802,0x0a03,0x0b00,0x0901,0x0c0e,0x0e0f,0x0f0c,0x0d0d,0x0006,0x0207,0x0304,0x0105,0x040a,0x060b,0x0708,0x0509,
  0x0204,0x0005,0x0106,0x0307,0x0608,0x0409,0x050a,0x070b,0x0a00,0x0801,0x0902,0x0b03,0x0e0c,0x0c0d,0x0d0e,0x0f0f,
  0x0204,0x0006,0x0105,0x0307,0x0a08,0x080a,0x0909,0x0b0b,0x0600,0x0402,0x0501,0x0703,0x0e0c,0x0c0e,0x0d0d,0x0f0f,
  0x0204,0x0007,0x0105,0x0306,0x0a08,0x080b,0x0909,0x0b0a,0x0e00,0x0c03,0x0d01,0x0f02,0x060c,0x040f,0x050d,0x070e,
  0x0a05,0x0806,0x0904,0x0b07,0x0209,0x000a,0x0108,0x030b,0x0601,0x0402,0x0500,0x0703,0x0e0d,0x0c0e,0x0d0c,0x0f0f,
  0x0a05,0x0807,0x0904,0x0b06,0x0209,0x000b,0x0108,0x030a,0x0e01,0x0c03,0x0d00,0x0f02,0x060d,0x040f,0x050c,0x070e,
  0x0a06,0x0807,0x0904,0x0b05,0x0e0a,0x0c0b,0x0d08,0x0f09,0x0202,0x0003,0x0100,0x0301,0x060e,0x040f,0x050c,0x070d,
  0x0204,0x0005,0x0306,0x0107,0x060c,0x040d,0x070e,0x050f,0x0a00,0x0801,0x0b02,0x0903,0x0e08,0x0c09,0x0f0a,0x0d0b,
  0x0204,0x0006,0x0305,0x0107,0x0a0c,0x080e,0x0b0d,0x090f,0x0600,0x0402,0x0701,0x0503,0x0e08,0x0c0a,0x0f09,0x0d0b,
  0x0204,0x0007,0x0305,0x0106,0x0a0c,0x080f,0x0b0d,0x090e,0x0e00,0x0c03,0x0f01,0x0d02,0x0608,0x040b,0x0709,0x050a,
  0x0a05,0x0806,0x0b04,0x0907,0x020d,0x000e,0x030c,0x010f,0x0601,0x0402,0x0700,0x0503,0x0e09,0x0c0a,0x0f08,0x0d0b,
  0x0a05,0x0807,0x0b04,0x0906,0x020d,0x000f,0x030c,0x010e,0x0e01,0x0c03,0x0f00,0x0d02,0x0609,0x040b,0x0708,0x050a,
  0x0a06,0x0807,0x0b04,0x0905,0x0e0e,0x0c0f,0x0f0c,0x0d0d,0x0202,0x0003,0x0300,0x0101,0x060a,0x040b,0x0708,0x0509,
  0x0208,0x0309,0x000a,0x010b,0x060c,0x070d,0x040e,0x050f,0x0a00,0x0b01,0x0802,0x0903,0x0e04,0x0f05,0x0c06,0x0d07,
  0x0208,0x030a,0x0009,0x010b,0x0a0c,0x0b0e,0x080d,0x090f,0x0600,0x0702,0x0401,0x0503,0x0e04,0x0f06,0x0c05,0x0d07,
  0x0208,0x030b,0x0009,0x010a,0x0a0c,0x0b0f,0x080d,0x090e,0x0e00,0x0f03,0x0c01,0x0d02,0x0604,0x0707,0x0405,0x0506,
  0x0a09,0x0b0a,0x0808,0x090b,0x020d,0x030e,0x000c,0x010f,0x0601,0x0702,0x0400,0x0503,0x0e05,0x0f06,0x0c04,0x0d07,
  0x0a09,0x0b0b,0x0808,0x090a,0x020d,0x030f,0x000c,0x010e,0x0e01,0x0f03,0x0c00,0x0d02,0x0605,0x0707,0x0404,0x0506
};

/*-------------------------------------------------------------------------
Direct Matrix Inversion For Complex Matrices

Description: perform in-place inversion of 2x2 and 4x4 complex matrices. 
2x2 matrices are inverted by Cramer's rule. For 3x3, 4x4 matrices we employ 
the blockwise inversion algorithm encompassed with a suboptimal row/column
permutation that gains better conditioning of the block structure.

Data format and order options:
  Suffix   Data Order                 Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Notes:
1. In general, accuracy of a matrix inversion algorithm implementation is a
   function of input matrix condition number. Thus it is user's responsibility
   to qualify the reliability of numeric results. Refer to NatureDSP Baseband 
   Library Reference for details. 
2. For blockwise inversion of 4x4 fixed-point matrices, it is reasonable to
   limit the dynamic range of input data by 11..13 significant bits. This
   measure reduces the possibility of an overflow at internal computations.

Parameters:
Temporary:
pScr      Scratch memory area. To determine the scratch area size required by
          a function <fun>, use the respective helper function 
          <fun>_getScratchSize()
Input:
L         Number of matrices
qA        Number of fractional bits for fixed-point input/output data
Input/Output:
A[L][SA]  Sequence of L NxN complex input/result matrices. SA is the number
          of data elements occupied by a single NxN matrix in a block 
          (stream) ordered sequence, see function specifications.

Restrictions:
pScr,A    Aligned on 32-byte boundary

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/
void cmatinv4x4n(void * restrict pScr, complex_fract16 * restrict A, int qA, int L)
{

  int16_t * restrict a_perm_ix; // Permutation indices for L matrices
  int16_t * restrict a_eX;      // Normalization shift amounts for L matrices
  int16_t * restrict a_qA;      // Fixed point position for L 2x2 sub-blocks A
  int16_t * restrict a_qT;      // Fixed point position for L 2x2 matrices T
  int16_t * restrict a_qS;      // Fixed point position for L 2x2 matrices S
  int16_t * restrict a_A;       // Upper left 2x2 sub-block 
  int16_t * restrict a_B;       // Upper right 2x2 sub-block
  int16_t * restrict a_C;       // Lower left 2x2 sub-block 
  int16_t * restrict a_D;       // Lower right 2x2 sub-block
  int16_t * restrict a_T;       // L 2x2 matrix products C*A^-1, A^-1*B

  int16_t * restrict A_ = (int16_t * restrict)A;
  int L_;

  NASSERT_ALIGN32( pScr );
  NASSERT_ALIGN32( A    );

  if (L <= 0) return;

  //
  // Partition the scratch memory area.
  //
  
  // Number of input matrices rounded up to the next multiple of SIMD_WIDTH/2.
  L_ = L + ( -L & (BBE_SIMD_WIDTH/2-1) );

  {
    void * ptr = pScr;
    
    a_A       = (int16_t*)ALIGN_PTR( ptr );    ptr = a_A + L_*4*2;
    a_B       = (int16_t*)ALIGN_PTR( ptr );    ptr = a_B + L_*4*2;
    a_C       = (int16_t*)ALIGN_PTR( ptr );    ptr = a_C + L_*4*2;
    a_D       = (int16_t*)ALIGN_PTR( ptr );    ptr = a_D + L_*4*2;
    a_T       = (int16_t*)ALIGN_PTR( ptr );    ptr = a_T + L_*4*2;
    a_qA      = (int16_t*)ALIGN_PTR( ptr );    ptr = a_qA + L_;
    a_qT      = (int16_t*)ALIGN_PTR( ptr );    ptr = a_qT + L_;
    a_qS      = (int16_t*)ALIGN_PTR( ptr );    ptr = a_qS + L_;
    a_eX      = (int16_t*)ALIGN_PTR( ptr );    ptr = a_eX + L;
    a_perm_ix = (int16_t*)ALIGN_PTR( ptr );    ptr = a_perm_ix + L;

    NASSERT( (int8_t*)ptr - (int8_t*)pScr <= (int)cmatinv4x4n_getScratchSize( L ) );
  }

  //
  // Normalize input 4x4 matrices and search for a suboptimal permutation 
  // for each matrix, so that the 2x2 sub-block at the upper left of a
  // permuted matrix has a large determinant.
  //

  cmatinv4x4n_searchPermutation( A_, a_eX, a_perm_ix, search_tbl, L );

  cmatinv4x4n_permute_b2bs( a_A, a_B, a_C, a_D, A_, a_perm_ix, pattern_tbl, L );

  //----------------------------------------------------------------------------
  // Invert the upper left block A, then compute and invert the Schur
  // complement of A: S = (D - C*A^-1*B)^-1

  // "A" <- A^-1
  // CQ(qA=(15-eA)) <- Q30/( CQ15 + eA )
  cmatinv4x4n_inv2x2bs( A_, a_A, a_qA, L_, 15 );

  // "T" <- C*A^-1
  // CQ(13-eA+eCA) <- CQ15*CQ(15-eA) - ( 17 - eCA ) with/without rounding
  // CQ(qA-2+eCA) <- CQ15*CQ(qA) - ( 17 - eCA ) with/without rounding
  // qT <- 15 + qA - 17 + eCA = qA - 2 + eCA
  // Exponent is bounded to protect D*2^qT from overflow.
  if      ( NORM_CA &&  ROUND_CA ) cmatinv4x4n_mul2x2bs_nr( a_T, a_qT, a_C, a_A, a_qA, L_, 15, 10 );
  else if ( NORM_CA && !ROUND_CA ) cmatinv4x4n_mul2x2bs_n ( a_T, a_qT, a_C, a_A, a_qA, L_, 15, 10 );
  else    NASSERT( !"Option not implemented!" );

  // "D" <- D - C*A^-1*B
  // CQ(qT-3+eS) <- { [ CQ15 + qT ] - CQ(qT)*CQ15 } - ( 18 - eS ) w/ rounding
  // qS <- 15 + qT - 18 + eS = qT - 3 + eS
  if ( NORM_S) cmatinv4x4n_mas2x2bs_nr( a_D, a_qS, a_T, a_B, a_qT, L_, 15 );
  else         NASSERT( !"Option not implemented!" );

  // "D" <- S = ( D - C*A^-1*B )^-1
  // CQ(20+eA-eCA-eS-eSi) <- Q30/( CQ(10-eA+eCA+eS) + eSi )
  // CQ(30-qS) <- Q30/CQ(qS)
  // qS <- 30-qS
  cmatinv4x4n_inv2x2bs_q( A_, a_D, a_qS, L_ );//ref code

  //----------------------------------------------------------------------------
  // Compute 2x2 blocks of the inverted matrix.

  // "C" <- -1*S*C*A^-1
  // CQ(qX) <- -1*CQ(qS)*CQ(qT)/2^(15-qX-eX) - ( qS + qT ) + qX
  cmatinv4x4n_mul2x2bs_srn( a_C, a_D, a_T, a_eX, a_qS, a_qT, L_, 2*qA-15 );

  // "T" <- A^-1*B
  // CQ(qA-2+eAB) <- CQ(qA)*Q15 - ( 17 - eAB )
  // qT <- qA + 15 - 17 + eAB = qA - 2 + eAB
  // Exponent is bounded to protect D*2^qT from overflow.
  if      ( NORM_AB &&  ROUND_AB ) cmatinv4x4n_mul2x2bs_nr( a_T, a_qT, a_A, a_B, a_qA, L_, 15, 10 );
  else if ( NORM_AB && !ROUND_AB ) cmatinv4x4n_mul2x2bs_n ( a_T, a_qT, a_A, a_B, a_qA, L_, 15, 10 );
  else    NASSERT( !"Option not implemented!" );

  // "A" <- A^-1 + A^-1*B*S*C*A^-1
  // CQ(qT+qX) <- [ CQ(qA) + qT + qX - qA ]/2^(15-qX-eX) - CQ(qT)*CQ(qX)
  // CQ(qX) <- CQ(qT+qX) - qT w/ rounding
  cmatinv4x4n_mas2x2bs_sr( a_A, a_T, a_C, a_qA, a_qT, a_eX, L_, qA, 15-qA );

  // "B" <- -1*A^-1*B*S
  // CQ(qX) <- -1*CQ(qS)*CQ(qT)/2^(15-qX-eX) - ( qS + qT ) + qX
  cmatinv4x4n_mul2x2bs_srn( a_B, a_T, a_D, a_eX, a_qS, a_qT, L_, 2*qA-15 );

  // "D" <- S
  // CQ(qX) <- CQ(qS)/2^(15-qX-eX) - qS + qX
  cmatinv4x4n_rnd2x2bs( a_D, a_eX, a_qS, L_, 15-2*qA );


  //----------------------------------------------------------------------------
  // Combine the inverse matrix from 2x2 sub-blocks and apply the inverse
  // permutation to compensate for initial reordering of the original matrix.

  cmatinv4x4n_permute_bs2b( A_, a_A, a_B, a_C, a_D, a_perm_ix, pattern_tbl, L );//desel -6
} /* cmatinv4x4n() */

size_t cmatinv4x4n_getScratchSize ( int L )
{
  int L_, sz = 0;

  // Number of input matrices rounded up to the next multiple of SIMD_WIDTH/2.
  L_ = L + ( -L & (BBE_SIMD_WIDTH/2-1) );

  if ( L>0 )
  {
    sz = 5*ALIGN_SIZE( L_*4*2*sz_i16 ) + // a_A, a_B, a_C, a_D, a_T 
         3*ALIGN_SIZE( L_*sz_i16     ) + // a_qA, a_qT, a_qS
         2*ALIGN_SIZE( L *sz_i16     );  // a_eX, a_perm_ix
  }

  return ( (size_t)sz );
} /* cmatinv4x4n_getScratchSize() */

#else
DISCARD_FUN(void, cmatinv4x4n, ( void * restrict pScr, complex_fract16 * restrict A, int qA, int L ) )
size_t cmatinv4x4n_getScratchSize ( int L ) { (void)L; return 0; }
#endif
