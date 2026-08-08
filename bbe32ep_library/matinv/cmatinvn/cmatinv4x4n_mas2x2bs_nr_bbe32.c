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
                        int L, int Qzy )
{
  const xb_vecNx16 *          Z_rd;
        xb_vecNx16 * restrict Z_wr;
        xb_vecNx16 * restrict QZ_wr;
  const xb_vecNx16 *          X_rd;
  const xb_vecNx16 *          Y_rd;
  const xb_vecNx16 *          QX_rd;
  valign aQX_rd,aQZ_wr;

  xb_vecNx16 x00, x01, x10, x11;
  xb_vecNx16 y00, y01, y10, y11;
  xb_vecNx16 z00, z01, z10, z11;
  xb_vecNx40 w00, w01, w10, w11;

  xb_vecNx16 qX, eZ, qZ;
  xb_vecNx16 c24;

  vsaN  vsa0;
  int l;

  NASSERT_ALIGN( Z, (2*BBE_SIMD_WIDTH) );
  NASSERT_ALIGN( X, (2*BBE_SIMD_WIDTH) );
  NASSERT_ALIGN( Y, (2*BBE_SIMD_WIDTH) );

  NASSERT( L%(BBE_SIMD_WIDTH/2)==0 );

  c24 = BBE_MOVVINT16( 24 );

  Z_rd  = (const xb_vecNx16*)Z;
  Z_wr  = (      xb_vecNx16*)Z;
  QZ_wr = (      xb_vecNx16*)Qz;
  X_rd  = (const xb_vecNx16*)X;
  Y_rd  = (const xb_vecNx16*)Y;
  QX_rd = (const xb_vecNx16*)Qx;

  aQZ_wr = BBE_ZALIGN();
  aQX_rd = BBE_LA_PP(QX_rd);

  for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++ )
  {
    //--------------------------------------------------------------------------
    // Perform MAS on BBE_SIMD_WIDTH/2 2x2 matrices, with accumulators 
    // pre-scaling.

    // CQ(qX)
    BBE_LVNX16_IP( x00, X_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x01, X_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x10, X_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x11, X_rd, 2*BBE_SIMD_WIDTH );
                              
    // CQ(qZY)                
    BBE_LVNX16_IP( y00, Y_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( y01, Y_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( y10, Y_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( y11, Y_rd, 2*BBE_SIMD_WIDTH );
                              
    // CQ(qZY)                
    BBE_LVNX16_IP( z00, Z_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( z01, Z_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( z10, Z_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( z11, Z_rd, 2*BBE_SIMD_WIDTH );

    BBE_LAVNX16_XP(qX,aQX_rd,QX_rd,BBE_SIMD_WIDTH);
    qX=BBE_SHFLNX16I(qX,BBE_SHFLI_DOUBLE_1_LO);
    vsa0 = BBE_MOVVSV( qX, 0 );
    vsa0 = BBE_SUBSAVSN(15, vsa0);
    // CQ(qZY+qX) <- CQ(qZY) + qX
    w00 = BBE_UNPKNVNX16( z00, vsa0 );
    w01 = BBE_UNPKNVNX16( z01, vsa0 );
    w10 = BBE_UNPKNVNX16( z10, vsa0 );
    w11 = BBE_UNPKNVNX16( z11, vsa0 );

    // CQ(qZY+qX) <- CQ(qZY+qX) - CQ(qX)*CQ(qZY)
    BBE_MULSNX16C( w00, x00, y00 );
    BBE_MULSNX16C( w01, x00, y01 );
    BBE_MULSNX16C( w10, x10, y00 );
    BBE_MULSNX16C( w11, x10, y01 );

    BBE_MULSNX16C( w00, x01, y10 );
    BBE_MULSNX16C( w01, x01, y11 );
    BBE_MULSNX16C( w10, x11, y10 );
    BBE_MULSNX16C( w11, x11, y11 );

    //--------------------------------------------------------------------------
    // Normalize, round and save BBE_SIMD_WIDTH/2 2x2 matrix results.

    {    
      vsaN _24 = BBE_MOVVSA32(24);
      xb_vecNx16 s00,s01,s10,s11;
      vboolN vt;
      s00 = BBE_PACKVNX40(w00,_24);
      s01 = BBE_PACKVNX40(w01,_24);
      s10 = BBE_PACKVNX40(w10,_24);
      s11 = BBE_PACKVNX40(w11,_24);

      BBE_BMAXABSNX16(vt, s00,s00,s01);
      BBE_BMAXABSNX16(vt, s10,s10,s11);
      s00=BBE_MAXNX16(s00,s10);

      // Normalize, round and save 8 matrix products.
      vsa0 = BBE_NSANX16C( s00 );
      // 18-eZ <- 24 - (eZ+6)
      vsa0 = BBE_SUBSAVSN(24,vsa0 );
      eZ =BBE_MOVVVS( vsa0 );
    }

    w00 = BBE_RNDSADJNX40( w00, vsa0 );
    w01 = BBE_RNDSADJNX40( w01, vsa0 );
    w10 = BBE_RNDSADJNX40( w10, vsa0 );
    w11 = BBE_RNDSADJNX40( w11, vsa0 );

    // CQ(qZY+qX-18+eZ) <- CQ(qZY+qX) - (18-eZ)
    z00 = BBE_PACKVNX40( w00, vsa0 );
    z01 = BBE_PACKVNX40( w01, vsa0 );
    z10 = BBE_PACKVNX40( w10, vsa0 );
    z11 = BBE_PACKVNX40( w11, vsa0 );

    BBE_SVNX16_IP( z00, Z_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( z01, Z_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( z10, Z_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( z11, Z_wr, 2*BBE_SIMD_WIDTH );

    //--------------------------------------------------------------------------
    // Calculate fixed point position for each output matrix.

    // qZ <- qZY
    qZ = BBE_MOVVA16( Qzy );
    // qZ <- qZY + qX
    qZ = BBE_ADDNX16( qZ, qX );
    // qZ <- qZY + qZ - 18 + eZ
    qZ = BBE_SUBNX16( qZ, eZ );
    qZ=BBE_SELNX16I(qZ,qZ,BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    BBE_SAVNX16_XP(qZ,aQZ_wr,QZ_wr,BBE_SIMD_WIDTH);
  }

  BBE_SAPOS_FP(aQZ_wr,QZ_wr);

} // cmatinv4x4n_mas2x2bs_nr()

#endif
