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

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Common utility declarations. */
#include "cmatinvn_common.h"

#if HAVE_MATINV4X4

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
                             int L, int Qy, int Emax )
{

        xb_vecNx16 * restrict Z_wr;
        xb_vecNx16 * restrict QZ_wr;
  const xb_vecNx16 *          X_rd;
  const xb_vecNx16 *          Y_rd;
  const xb_vecNx16 *          QX_rd;

  xb_vecNx16 x00, x01, x10, x11;
  xb_vecNx16 y00, y01, y10, y11;
  xb_vecNx16 z00, z01, z10, z11;
  xb_vecNx40 w00, w01, w10, w11;

  xb_vecNx16 t0;
  xb_vecNx16 eZ, qZ, qX;
  valign aQX_rd,aQZ_wr;

  vsaN vsa0, vsa1, vsa2, vsa3, eLim;

  int l;

  NASSERT_ALIGN( Z, (BBE_SIMD_WIDTH*2) );
  NASSERT_ALIGN( X, (BBE_SIMD_WIDTH*2) );
  NASSERT_ALIGN( Y, (BBE_SIMD_WIDTH*2) );

  NASSERT( L%(BBE_SIMD_WIDTH/2)==0 );

  eLim = BBE_MOVVSA32( Emax +7);

  X_rd  = (const xb_vecNx16*)X;
  Y_rd  = (const xb_vecNx16*)Y;
  QX_rd = (const xb_vecNx16*)Qx;
  Z_wr  = (      xb_vecNx16*)Z;
  QZ_wr = (      xb_vecNx16*)Qz;
  aQX_rd=BBE_LA_PP(QX_rd);
  aQZ_wr=BBE_ZALIGN();

  for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++ )
  {
    //--------------------------------------------------------------------------
    // Calculate BBE_SIMD_WIDTH/2 2x2*2x2 matrix products, leave results in wide
    // registers.

    BBE_LVNX16_IP( x00, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x01, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x10, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x11, X_rd, +4*BBE_SIMD_WIDTH/2 );

    BBE_LVNX16_IP( y00, Y_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( y01, Y_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( y10, Y_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( y11, Y_rd, +4*BBE_SIMD_WIDTH/2 );

    // Q(Qx[l]+Qy) <= Q(Qx[l])*Q(Qy)
    w00 = BBE_MULNX16C( x00, y00 );
    w01 = BBE_MULNX16C( x00, y01 );
    w10 = BBE_MULNX16C( x10, y00 );
    w11 = BBE_MULNX16C( x10, y01 );

    BBE_MULANX16C( w00, x01, y10 );
    BBE_MULANX16C( w01, x01, y11 );
    BBE_MULANX16C( w10, x11, y10 );
    BBE_MULANX16C( w11, x11, y11 );

    //--------------------------------------------------------------------------

    {
      vsa0 = BBE_NSANX40C( w00 );
      vsa1 = BBE_NSANX40C( w01 );
      vsa2 = BBE_NSANX40C( w10 );
      vsa3 = BBE_NSANX40C( w11 );

      vsa0 = BBE_MINVSN( vsa0, eLim );
      vsa0 = BBE_MINVSN( vsa0, vsa1 );
      vsa0 = BBE_MINVSN( vsa0, vsa2 );
      vsa0 = BBE_MINVSN( vsa0, vsa3 );
      // 17-eZ <- (17+7) - (eZ+7)
      vsa0 = BBE_SUBSAVSN(24, vsa0);
      eZ   = BBE_MOVVVS( vsa0 );
    }

    // Q(Qx[l]+Qy-17+eZ) <- Q(Qx[l]+Qy) - (17-eZ) w/ saturation
    z00 = BBE_PACKVNX40( w00, vsa0 );
    z01 = BBE_PACKVNX40( w01, vsa0 );
    z10 = BBE_PACKVNX40( w10, vsa0 );
    z11 = BBE_PACKVNX40( w11, vsa0 );

    BBE_SVNX16_IP( z00, Z_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( z01, Z_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( z10, Z_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( z11, Z_wr, +4*BBE_SIMD_WIDTH/2 );

    //--------------------------------------------------------------------------
    // Calculate fixed point position for each output matrix.

    // t0 <- 17-eZ
    t0 = BBE_SELNX16I( eZ, eZ, BBE_SELI_EXTRACT_1_OF_2_OFF_0 );
    BBE_LAVNX16_XP(qX,aQX_rd,QX_rd,BBE_SIMD_WIDTH);

    // qZ <- Qx[l]+Qy
    qZ = BBE_MOVVA16( Qy );
    qZ = BBE_ADDNX16( qZ, qX );
    // qZ <- Qx[l]+Qy-17+eZ
    qZ = BBE_SUBNX16( qZ, t0 );
    BBE_SAVNX16_XP(qZ,aQZ_wr,QZ_wr,BBE_SIMD_WIDTH);
  }

  BBE_SAPOS_FP(aQZ_wr,QZ_wr);

} // cmatinv4x4n_mul2x2bs_n()

#endif
