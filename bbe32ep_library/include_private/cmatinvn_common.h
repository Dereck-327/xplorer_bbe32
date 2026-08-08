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
    Common declarations for blockwise inversion of 4x4 complex matrices.
    IntegrIT, 2006-2017
*/

#ifndef CMATINVN_COMMON__
#define CMATINVN_COMMON__

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"

#define HAVE_MATINV2X2 (HAVE_VSAMATH && HAVE_DIV && HAVE_RECIP && 1)
#define HAVE_MATINV4X4 (HAVE_VSAMATH && HAVE_DIV && 1)
#define HAVE_HERMMATINV2X2 (HAVE_VSAMATH && HAVE_RECIP && 1)
#define HAVE_HERMMATINV4X4 (HAVE_VSAMATH && HAVE_DIV && 1)

/*-------------------------------------------------------------------
 Normalize input 4x4 matrices and search for a suboptimal permutation 
 for each matrix, so that the 2x2 subblock at the upper left of a
 permuted matrix has a large determinant.
  Input:
    X[L][16*2]             Input 4x4 complex matrices, CQ(qX)
    cmatinv4x4n_searchPermutation:
        search_tbl[16][2*8]    For each of 16 element locations, enumerates
                               8 variants of the upper left 2x2 subblock
                               selection patterns that relocate the element
                               of interest to that subblock.
    chermmatinv4x4n_searchPermutation:
        search_tbl[8]          Enumerates upper left 2x2 subblock
                               selection patterns for 6 variants of
                               symmetric permutation
  Output:
    X[L][16*2]             Normalized 4x4 complex matrices, CQ(qX+eX)
    eX[L]                  Normalization shift amount for each matrix
    perm_ix[L]             Permutation indices
  Restrictions:
    X  Must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
-------------------------------------------------------------------*/
void cmatinv4x4n_searchPermutation( int16_t * restrict X,
                                    int16_t * restrict eX,
                                    int16_t * restrict perm_ix,
                              const int16_t *          search_tbl,
                                    int L );

/* Simplified version for Hermitian matrices. */
void chermmatinv4x4n_searchPermutation( int16_t * restrict X,
                                        int16_t * restrict eX,
                                        int16_t * restrict perm_ix,
                                  const int16_t *          search_tbl,
                                        int L );

/*-------------------------------------------------------------------
 Apply permutation to each 4x4 complex matrix, then convert the whole
 bulk of matrices to block-streaming format: each group of 8 consecutive
 matrices is converted to streaming format [4*4][8*2], and that groups
 are stored sequentally, i.e. we arrive to format [L/8][4*4][8*2]. Finally,
 2x2 subblocks are separated into 4 arrays A, B, C and D: [L/8][2*2][8*2].
  Input:
    X[L][16*2]              4x4 complex matrices in block format
    perm_ix[L]              Permutation index for each matrix, 0..34
    pattern_tbl[35][16]     35-entry permutation pattern table
  Output:
    A[(L+7)/8][4][8*2]      Upper left 2x2 subblock in block-streaming format
    B[(L+7)/8][4][8*2]      Upper right 2x2 subblock in block-streaming format
    C[(L+7)/8][4][8*2]      Lower left 2x2 subblock in block-streaming format
    D[(L+7)/8][4][8*2]      Lower right 2x2 subblock in block-streaming format
  Note:
    L may be not a multiple of 8, in which case the last block is complemented
    with the respective volume of dummy data.
  Restrictions:
    A,B,C,D,X must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
-------------------------------------------------------------------*/
void cmatinv4x4n_permute_b2bs( int16_t * restrict A,
                               int16_t * restrict B,
                               int16_t * restrict C,
                               int16_t * restrict D,
                               int16_t * restrict X,
                         const int16_t * restrict perm_ix,
                         const int16_t * restrict pattern_tbl,
                               int L );

/*-------------------------------------------------------------------
 For each group of 8 matrices, merge 2x2 subblocks A, B, C and D into
 4x4 matrices, then convert 8 4x4 matrices from streaming format to block
 format. After that, permute each 4x4 matrix according to its permutation
 index.
  Input:
    A[(L+7)/8][4][8*2]      Upper left 2x2 subblock in block-streaming format
    B[(L+7)/8][4][8*2]      Upper right 2x2 subblock in block-streaming format
    C[(L+7)/8][4][8*2]      Lower left 2x2 subblock in block-streaming format
    D[(L+7)/8][4][8*2]      Lower right 2x2 subblock in block-streaming format
    perm_ix[L]              Permutation index for each matrix, 0..34
    pattern_tbl[35][16]     35-entry permutation pattern table
  Output:
    Y[L][16*2]              4x4 complex matrices in block format
  Note:
    L may be not a multiple of 8.
  Restrictions:
    Y,A,B,C,D  Must be aligned on BBE_SIMD_WIDTH*2-byte boundary
-------------------------------------------------------------------*/
void cmatinv4x4n_permute_bs2b( int16_t * restrict Y,
                         const int16_t * restrict A,
                         const int16_t * restrict B,
                         const int16_t * restrict C,
                         const int16_t * restrict D,
                         const int16_t * restrict perm_ix,
                         const int16_t * restrict pattern_tbl,
                               int L );

/*-------------------------------------------------------------------
 In-place matrix inversion for 2x2 complex matrices stored in
 block-streaming format.

 Internal operation for l-th matrix:
   CQ(30-Q-eX[l]) <- (1<<30)/( CQ(Q) + eX[l] )
 where eX[l] is a normalization shift amount computed internally
 for each input matrix.

 Fixed point position for input data Q is actually used only to compute
 the fixed point position for output data Qi[L]!

  Temporary:
    pScr[(L/16)*16*4]  Scratch memory area
  Input:
    X[L/8][4][8*2]     2x2 complex matrices in block-streaming format, CQ(Q)
    Q                  Fixed point position for input matrices
  Output:
    X[L/8][4][8*2]     Inverse 2x2 matrices, CQ(Qi[L])
    Qi[L]              Fixed point position for each inverse matrix
  Restrictions:
    X,Qi,pScr          Must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
    L                  Must be a multiple of 8
-------------------------------------------------------------------*/
void cmatinv4x4n_inv2x2bs( int16_t * restrict pScr,
                           int16_t * restrict X,
                           int16_t * restrict Qi,
                           int L, int Q );

/* Special version for Hermitian matrices. */
void chermmatinv4x4n_inv2x2bs( int16_t * restrict pScr,
                               int16_t * restrict X,
                               int16_t * restrict Qi,
                               int L, int Q );

/*-------------------------------------------------------------------
 In-place matrix inversion for 2x2 complex matrices stored in
 block-streaming format. This function is similar to inv2x2bs(),
 but it allows to specify the fixed point position for each
 input matrix, and it doesn't normalize input matrices.

 Internal operation for the l-th matrix:
   CQ(30-Q[l]-eX[l]) <- (1<<30)/( CQ(Q[l]) + eX[l] )
 where eX[l] is a normalization shift amount computed internally
 for each inverse matrix.

 Fixed point position for input data is actually used only to compute
 the fixed point position for output data!

  Temporary:
    pScr[(L/16)*16*4]  Scratch memory area
  Input:
    X[L/8][4][8*2]     2x2 complex matrices in block-streaming format, CQ(Q)
    Q[L]               Fixed point position for each input matrix
  Output:
    X[L/8][4][8*2]     Inverse 2x2 matrices, CQ(Qi[L])
    Q[L]               Fixed point position for each inverse matrix
  Restrictions:
    X,Q,pScr           Must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
    L                  Must be a multiple of 8
-------------------------------------------------------------------*/
void cmatinv4x4n_inv2x2bs_q( int16_t * restrict pScr,
                             int16_t * restrict X,
                             int16_t * restrict Q,
                             int L );

/* Special version for Hermitian matrices. */
void chermmatinv4x4n_inv2x2bs_q( int16_t * restrict pScr,
                                 int16_t * restrict X,
                                 int16_t * restrict Q,
                                 int L );

/*-------------------------------------------------------------------
 Matrix multiplication for 16-bit complex 2x2 matrices stored in
 block-streaming format. Matrix product is normalized and ROUNDED to
 16 bits.

 Fixed point operation:
   Q(Qx[]+Qy-17+eZ) <- Q(Qx[])*Q(Qy) - 17 + eZ w/ symmetric rounding
 where eZ holds product's normalization shift amount.

  Input:
    X[L/8][4][8*2]  Left-hand multiplicands, L complex 2x2 matrices
    Y[L/8][4][8*2]  Right-hand multiplicands, L complex 2x2 matrices
    Qx[L]           Fixed point position enumerated for left-hand
                    multiplicands
    Qy              Fixed point position for right-hand multiplicands
    Emax            Upper limit for normalization shift amount eZ
  Output:
    Z[L/8][4][8*2]  L complex 2x2 product matrices
    Qz[L]           Fixed point position for product matrices: Qx[L]+Qy-17+eZ
  Notes:
    Qx[L] and Qy may be swapped in the sense that Qx[L] may actually represent
    the fixed point position for Y[] matrices, in which case Qy holds the
    number of fractional bits for X[] matrices.
  Restrictions:
    L                      Multiple of 8
    Z[],Qz[],X[],Qx[],Y[]  Must not overlap
    Z[],X[],Y[]            Must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
    Qz[],Qx[]              Must be aligned on 16-bytes boundary
-------------------------------------------------------------------*/
void cmatinv4x4n_mul2x2bs_nr( int16_t * restrict Z,
                              int16_t * restrict Qz,
                        const int16_t *          X,
                        const int16_t *          Y,
                        const int16_t *          Qx,
                              int L, int Qy, int Emax );

/*-------------------------------------------------------------------
 Matrix multiplication for 16-bit complex 2x2 matrices stored in
 block-streaming format. Matrix product is normalized and TRUNCATED to
 16 bits.

 Fixed point operation:
   Q(Qx[]+Qy-17+eZ) <- Q(Qx[])*Q(Qy) - 17 + eZ w/ saturation
 where eZ holds product's normalization shift amount.

  Input:
    X[L/8][4][8*2]  Left-hand multiplicands, L complex 2x2 matrices
    Y[L/8][4][8*2]  Right-hand multiplicands, L complex 2x2 matrices
    Qx[L]           Fixed point position enumerated for left-hand
                    multiplicands
    Qy              Fixed point position for right-hand multiplicands
    Emax            Upper limit for normalization shift amount eZ
  Output:
    Z[L/8][4][8*2]  L complex 2x2 product matrices
    Qz[L]           Fixed point position for product matrices: Qx[L]+Qy-17+eZ
  Notes:
    Qx[L] and Qy may be swapped in the sense that Qx[L] may actually represent
    the fixed point position for Y[] matrices, in which case Qy holds the
    number of fractional bits for X[] matrices.
  Restrictions:
    L                      Multiple of 8
    Z[],Qz[],X[],Qx[],Y[]  Must not overlap
    Z[],X[],Y[]            Must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
    Qz[],Qx[]              Must be aligned on 16-bytes boundary
-------------------------------------------------------------------*/

void cmatinv4x4n_mul2x2bs_n( int16_t * restrict Z,
                             int16_t * restrict Qz,
                       const int16_t *          X,
                       const int16_t *          Y,
                       const int16_t *          Qx,
                             int L, int Qy, int Emax );

/*-------------------------------------------------------------------
 Matrix multiplication for 16-bit complex 2x2 matrices stored in
 block-streaming format, with scaling, rounding and negation.

   Z[] = -X[]*Y[]*2^(-Qx[]-Qy[]+Qz[]+Ez)

 Fixed point operation:
   Q(Qz[]+Ez) <- Q(Qx[])*Q(Qy[]) - Qx[] - Qy[] + Qz[] + Ez w/ rounding

  Input:
    X[L/8][4][8*2]  Left-hand multiplicands, L complex 2x2 matrices
    Y[L/8][4][8*2]  Right-hand multiplicands, L complex 2x2 matrices
    Qz[L]           Fixed point position for output matrices
    Qx[L]           Fixed point position enumerated for left-hand
                    multiplicands
    Qy[L]           Fixed point position enumerated for right-hand
                    multiplicands
    Ez              Additional scaling shift amount for output matrices
  Output:
    Z[L/8][4][8*2]  L complex 2x2 product matrices
  Notes:
    Qx[L] and Qy[L] may be swapped in the sense that Qx[L] may actually
    represent the fixed point position for Y[] matrices, in which case
    Qy[L] holds the number of fractional bits for X[] matrices.
  Restrictions:
    L                      Multiple of 8
    Z[],X[],Y[],Qx[],Qy[]  Must not overlap
    Z[],X[],Y[]            Must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
-------------------------------------------------------------------*/
void cmatinv4x4n_mul2x2bs_srn( int16_t * restrict Z,
                         const int16_t *          X,
                         const int16_t *          Y,
                         const int16_t *          Qz,
                         const int16_t *          Qx,
                         const int16_t *          Qy,
                               int L, int Ez );

/* Special vatiant for Hermitian inversion: accomplishes product matrices Z[] with
 * conjugate transposes stored in the output array W[]. */
void chermmatinv4x4n_mul2x2bs_srn( int16_t * restrict Z,
                                   int16_t * restrict W,
                             const int16_t *          X,
                             const int16_t *          Y,
                             const int16_t *          Qz,
                             const int16_t *          Qx,
                             const int16_t *          Qy,
                                   int L, int Ez );

/*-------------------------------------------------------------------
 Matrix multiply-and-subtract (MAS) for 16-bit complex 2x2 matrices stored in
 block-streaming format. MAS results are normalized and ROUNDED to 16 bits.

 Fixed point operation:
   1. MAS
       Q(Qzy+Qx[]) <- [ Q(Qzy) + Qx[] ] - Q(Qx[])*Q(Qzy)
   2. Normalization and rounding
       Q(Qzy+Qx[]-18+eZ) <- Q(Qzy+Qx[]) - 18 + eZ w/ symmetric rounding
      where eZ holds accumulator's normalization shift amount.

  Input:
    Z[L/8][4][8*2]  Accumulators, L complex 2x2 matrices, Qzy
    X[L/8][4][8*2]  Left-hand multiplicands, L complex 2x2 matrices, Qx[L]
    Y[L/8][4][8*2]  Right-hand multiplicands, L complex 2x2 matrices, Qzy
    Qx[L]           Fixed point position enumerated for left-hand
                    multiplicands X[]
    Qzy             Fixed point position for INPUT accumulators Z[] and right-
                    hand multiplicands Y[]
  Output:
    Z[L/8][4][8*2]  Accumulators, L complex 2x2 matrices, Qz[L]
    Qz[L]           Fixed point position for accumulators: Qx[L]+Qy-18+eZ
  Notes:
    Qx[L] and Qzy may be swapped in the sense that Qx[L] may actually
    represent the fixed point position for Y[] matrices, in which case
    Qzy holds the fixed point position for X[] matrices and input
    accumulator matrices Z[].
  Restrictions:
    L                      Multiple of 8
    Z[],Qz[],X[],Qx[],Y[]  Must not overlap
    Z[],X[],Y[]            Must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
    Qz[],Qx[]              Must be aligned on 16-bytes boundary
-------------------------------------------------------------------*/
void cmatinv4x4n_mas2x2bs_nr( int16_t * restrict Z,
                              int16_t * restrict Qz,
                        const int16_t *          X,
                        const int16_t *          Y,
                        const int16_t *          Qx,
                        int L, int Qzy );

/*-------------------------------------------------------------------
 Matrix multiply-and-subtract (MAS) for 16-bit complex 2x2 matrices stored in
 block-streaming format. MAS results are SCALED and ROUNDED to 16 bits.

   Z = [ Z*2^(Qx[]+Qxy+Ez[]-rsh-Qz[]) - X*Y ]*2^-Qx[]

 Fixed point operation:
   Q(Qx[]+Qzy+Ez[]-rsh) <- [ Q(Qz[])+Qx[]+Qzy+Ez[]-rsh-Qz[] ] - Q(Qx[])*Q(Qzy)
   Q(Qzy+Ez-rsh) <- Q(Qx[]+Qzy+Ez[]-rsh) - Qx[] w/ rounding

  Input:
    Z[L/8][4][8*2]  Accumulators, L complex 2x2 matrices, Qz[L]
    X[L/8][4][8*2]  Left-hand multiplicands, L complex 2x2 matrices, Qx[L]
    Y[L/8][4][8*2]  Right-hand multiplicands, L complex 2x2 matrices, Qzy
    Qz[L]           Fixed point position enumerated for INPUT accumulators Z[]
    Qx[L]           Fixed point position enumerated for left-hand
                    multiplicands X[]
    Ez[L]           Scaling shift amount enumerated for accumulators Z[]
    Qzy             Fixed point position for OUTPUT accumulators Z[] and
                    right-hand multiplicands Y[]
    rsh             Constant scaling right shift amount for all
                    accumulators Z[]
  Output:
    Z[L/8][4][8*2]  Accumulators, L complex 2x2 matrices, Qzy+Ez
  Notes:
    Qx[L] and Qzy may be swapped in the sense that Qx[L] may actually
    represent the fixed point position for Y[] matrices, in which case
    Qzy holds the fixed point position for X[] matrices and output
    accumulator matrices Z[].
  Restrictions:
    L                            Multiple of 8
    Z[],X[],Y[],Qz[],Qx[],Ez[]   Must not overlap
    Z[],X[],Y[]                  Must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
-------------------------------------------------------------------*/
void cmatinv4x4n_mas2x2bs_sr( int16_t * restrict Z,
                        const int16_t *          X,
                        const int16_t *          Y,
                        const int16_t *          Qz,
                        const int16_t *          Qx,
                        const int16_t *          Ez,
                        int L, int Qzy, int rsh );

/*-------------------------------------------------------------------
 Matrix multiply-and-subtract (MAS) for 16-bit complex 2x2 matrices stored in
 block-streaming format, with left-hand multiplier being conjugate transposed.
 MAS results are SCALED and ROUNDED to 16 bits.

   Z = [ Z*2^(Qx[]+Qzy+Ez[]-rsh-Qz[]) - X'*Y ]*2^-Qx[]

 Fixed point operation:
   Q(Qx[]+Qzy+Ez[]-rsh) <- [ Q(Qz[])+Qx[]+Qzy+Ez[]-rsh-Qz[] ] - Q(Qx[])*Q(Qzy)
   Q(Qzy+Ez-rsh) <- Q(Qx[]+Qzy+Ez[]-rsh) - Qx[] w/ rounding

  Input:
    Z[L/8][4][8*2]  Accumulators, L complex 2x2 matrices, Qz[L]
    X[L/8][4][8*2]  Left-hand multiplicands, L complex 2x2 matrices, Qx[L]
    Y[L/8][4][8*2]  Right-hand multiplicands, L complex 2x2 matrices, Qzy
    Qz[L]           Fixed point position enumerated for INPUT accumulators Z[]
    Qx[L]           Fixed point position enumerated for left-hand
                    multiplicands X[]
    Ez[L]           Scaling shift amount enumerated for accumulators Z[]
    Qzy             Fixed point position for OUTPUT accumulators Z[] and
                    right-hand multiplicands Y[]
    rsh             Constant scaling right shift amount for all
                    accumulators Z[]
  Output:
    Z[L/8][4][8*2]  Accumulators, L complex 2x2 matrices, Qzy+Ez
  Notes:
    Qx[L] and Qzy may be swapped in the sense that Qx[L] may actually
    represent the fixed point position for Y[] matrices, in which case
    Qzy holds the fixed point position for X[] matrices and output
    accumulator matrices Z[].
  Restrictions:
    L                            Multiple of 8
    Z[],X[],Y[],Qz[],Qx[],Ez[]   Must not overlap
    Z[],X[],Y[],Qz[],Qx[],Ez[]   Must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
-------------------------------------------------------------------*/
void chermmatinv4x4n_mas2x2bs_srj( int16_t * Z,
                             const int16_t * X,
                             const int16_t * Y,
                             const int16_t * Qz,
                             const int16_t * Qx,
                             const int16_t * Ez,
                                   int L, int Qzy, int rsh );

/*-------------------------------------------------------------------
 In-place rescaling with rounding for 2x2 complex matrices stored in
 block-streaming format: X[] = X[]*2^(sA[]-sB[]-rsh).
  Input:
    X[L/8][4][8*2]  Input complex 2x2 matrices
    sA[L]           Right shift amount for each matrix
    sB[L]           Left shift amount for each matrix
    rsh             Additional right shift amount for all matrices
  Output:
    X[L/8][4][8*2]  Rescaled and rounded complex 2x2 matrices
  Restrictions:
    L               Multiple of 8
    X[],sA[],sB[]   Must not overlap
    X[]             Must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
-------------------------------------------------------------------*/
void cmatinv4x4n_rnd2x2bs( int16_t * restrict X,
                     const int16_t *          sA,
                     const int16_t *          sB,
                           int L, int rsh );

#endif /* CMATINVN_COMMON__ */
