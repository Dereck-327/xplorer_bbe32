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
    Blockwise inversion of 4x4 Hermitian matrices
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Common declarations for blockwise inversion of 4x4 complex matrices. */
#include "cmatinvn_common.h"

#if HAVE_HERMMATINV4X4
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
/* Special vatiant for Hermitian inversion: accomplishes product matrices Z[] with
 * conjugate transposes stored in the output array W[]. */
void chermmatinv4x4n_mul2x2bs_srn( 
                               int16_t * restrict Z,
                               int16_t * restrict W,
                         const int16_t *          X,
                         const int16_t *          Y,
                         const int16_t *          Qz,
                         const int16_t *          Qx,
                         const int16_t *          Qy,
                         int L, int Ez )
{
        xb_vecNx16 * restrict Z_wr;
        xb_vecNx16 * restrict W_wr;
  const xb_vecNx16 *          X_rd;
  const xb_vecNx16 *          Y_rd;
  const xb_vecNx16 * restrict QZ_rd;
  const xb_vecNx16 * restrict QX_rd;
  const xb_vecNx16 * restrict QY_rd;

  xb_vecNx16 x00, x01, x10, x11;
  xb_vecNx16 y00, y01, y10, y11;
  xb_vecNx16 z00, z01, z10, z11;
  xb_vecNx40 w00, w01, w10, w11;

  xb_vecNx16 qZ, qX, qY, eZ;

  xb_vecNx16 t0, t1;

  vsaN vsa0;

  valign QZ_va, QX_va, QY_va;

  int l;

  NASSERT_ALIGN( Z, (2*BBE_SIMD_WIDTH) );
  NASSERT_ALIGN( W, (2*BBE_SIMD_WIDTH) );
  NASSERT_ALIGN( X, (2*BBE_SIMD_WIDTH) );
  NASSERT_ALIGN( Y, (2*BBE_SIMD_WIDTH) );

  NASSERT( L%(BBE_SIMD_WIDTH/2)==0 );
 
  eZ = BBE_MOVVA16( Ez );

  Z_wr  = (xb_vecNx16*)Z;
  W_wr  = (xb_vecNx16*)W;
  X_rd  = (xb_vecNx16*)X;
  Y_rd  = (xb_vecNx16*)Y;
  QZ_rd = (xb_vecNx16*)Qz;
  QX_rd = (xb_vecNx16*)Qx;
  QY_rd = (xb_vecNx16*)Qy;

  QZ_va = BBE_LAVNX16_PP( QZ_rd );
  QX_va = BBE_LAVNX16_PP( QX_rd );
  QY_va = BBE_LAVNX16_PP( QY_rd );

  for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++ )
  {
    //--------------------------------------------------------------------------
    // Calculate fixed point position for each output matrix.

    BBE_LAVNX16_XP( qZ, QZ_va, QZ_rd, BBE_SIMD_WIDTH );
    BBE_LAVNX16_XP( qX, QX_va, QX_rd, BBE_SIMD_WIDTH );
    BBE_LAVNX16_XP( qY, QY_va, QY_rd, BBE_SIMD_WIDTH );

    t0 = BBE_ADDNX16( qX, qY );
    t1 = BBE_ADDNX16( qZ, eZ );
    t0 = BBE_SUBNX16( t0, t1 );

    t0 = BBE_SHFLNX16I( t0, BBE_SHFLI_DOUBLE_1_LO );

    vsa0 = BBE_MOVVSV( t0, 0 );

    //--------------------------------------------------------------------------
    // Calculate BBE_SIMD_WIDTH/2 2x2*2x2 matrix products, leave results in
    // wide registers.

    BBE_LVNX16_IP( x00, X_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x01, X_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x10, X_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x11, X_rd, 2*BBE_SIMD_WIDTH );
                              
    BBE_LVNX16_IP( y00, Y_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( y01, Y_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( y10, Y_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( y11, Y_rd, 2*BBE_SIMD_WIDTH );

    x00 = BBE_NEGSNX16( x00 );
    x01 = BBE_NEGSNX16( x01 );
    x10 = BBE_NEGSNX16( x10 );
    x11 = BBE_NEGSNX16( x11 );

    // Q(qX+qY) <= Q(qX)*Q(qY)
    w00 = BBE_MULRNX16C( x00, y00, vsa0 );
    w01 = BBE_MULRNX16C( x00, y01, vsa0 );
    w10 = BBE_MULRNX16C( x10, y00, vsa0 );
    w11 = BBE_MULRNX16C( x10, y01, vsa0 );

    BBE_MULANX16C( w00, x01, y10 );
    BBE_MULANX16C( w01, x01, y11 );
    BBE_MULANX16C( w10, x11, y10 );
    BBE_MULANX16C( w11, x11, y11 );

    //--------------------------------------------------------------------------
    // Scale, round and pack BBE_SIMD_WIDTH/2 product matrices and store them
    // to output array.

    // Q(qZ+eZ) <- Q(qX+qY) - qX - qY + qZ + eZ
    z00 = BBE_PACKVNX40( w00, vsa0 );
    z01 = BBE_PACKVNX40( w01, vsa0 );
    z10 = BBE_PACKVNX40( w10, vsa0 );
    z11 = BBE_PACKVNX40( w11, vsa0 );

    BBE_SVNX16_IP( z00, Z_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( z01, Z_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( z10, Z_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( z11, Z_wr, 2*BBE_SIMD_WIDTH );

    /* Store conjugate transpose of product matrices. */
    BBE_SVNX16_IP( BBE_CONJSNX16C(z00), W_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( BBE_CONJSNX16C(z10), W_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( BBE_CONJSNX16C(z01), W_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( BBE_CONJSNX16C(z11), W_wr, 2*BBE_SIMD_WIDTH );
  }

} /* chermmatinv4x4n_mul2x2bs_srn() */

#endif
