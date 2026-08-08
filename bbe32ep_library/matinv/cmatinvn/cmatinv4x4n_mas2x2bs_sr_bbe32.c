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
                              int L, int Qzy, int rsh )
{
#if 0
  const xb_vecNx16 *          Z_rd;
        xb_vecNx16 * restrict Z_wr;
  const xb_vecNx16 *          X_rd;
  const xb_vecNx16 *          Y_rd;
  const xb_vecNx16 *          QZ_rd;
  const xb_vecNx16 *          QX_rd;
  const xb_vecNx16 *          EZ_rd;

  valign QZ_va, QX_va, EZ_va;

  xb_vecNx16 x00, x01, x10, x11;
  xb_vecNx16 y00, y01, y10, y11;
  xb_vecNx16 z00, z01, z10, z11;
  xb_vecNx40 w00, w01, w10, w11;

  xb_vecNx16 csh, qZ, qX, eZ;
  xb_vecNx16 t0;

  vsaN vsa0, vsa1;

  int l;

  NASSERT_ALIGN32( Z );
  NASSERT_ALIGN32( X );
  NASSERT_ALIGN32( Y );

  NASSERT( !(L&7) );

  csh = BBE_MOVVA16( Qzy - rsh ); 

  Z_rd  = (const xb_vecNx16*)Z;
  Z_wr  = (      xb_vecNx16*)Z;
  X_rd  = (const xb_vecNx16*)X;
  Y_rd  = (const xb_vecNx16*)Y;
  QZ_rd = (const xb_vecNx16*)Qz;
  QX_rd = (const xb_vecNx16*)Qx;
  EZ_rd = (const xb_vecNx16*)Ez;

  QZ_va = BBE_LAVNX16_PP( QZ_rd );
  QX_va = BBE_LAVNX16_PP( QX_rd );
  EZ_va = BBE_LAVNX16_PP( EZ_rd );

  for ( l=0; l<L/8; l++ )
  {
    //--------------------------------------------------------------------------
    // Calculate pre-scaling shift amount for 8 complex 2x2 accumulator matrices

    BBE_LAVNX16_XP( qZ, QZ_va, QZ_rd, 8*2 );
    BBE_LAVNX16_XP( qX, QX_va, QX_rd, 8*2 );
    BBE_LAVNX16_XP( eZ, EZ_va, EZ_rd, 8*2 );

    // t0 <- qX + eZ
    t0 = BBE_ADDNX16( qX, eZ );
    // t0 <- qX + eZ - qZ
    t0 = BBE_SUBNX16( t0, qZ );
    // t0 <- qX + qZY + eZ - rsh - qZ
    t0 = BBE_ADDNX16( t0, csh );

    t0 = BBE_SHFLNX16I( t0, BBE_SHFLI_DOUBLE_1_LO );
    qX = BBE_SHFLNX16I( qX, BBE_SHFLI_DOUBLE_1_LO );

    vsa0 = BBE_MOVVSV( t0, 0 );
    vsa1 = BBE_MOVVSV( qX, 0 );

    //--------------------------------------------------------------------------
    // Perform MAS on 8 2x2 matrices, with accumulators pre-scaling.

    // CQ(qX)
    BBE_LVNX16_IP( x00, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x01, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x10, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x11, X_rd, +4*BBE_SIMD_WIDTH/2 );

    // CQ(qZY)
    BBE_LVNX16_IP( y00, Y_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( y01, Y_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( y10, Y_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( y11, Y_rd, +4*BBE_SIMD_WIDTH/2 );

    // CQ(qZ)
    BBE_LVNX16_IP( z00, Z_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( z01, Z_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( z10, Z_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( z11, Z_rd, +4*BBE_SIMD_WIDTH/2 );

    w00 = BBE_UNPKSNX16( z00 );
    w01 = BBE_UNPKSNX16( z01 );
    w10 = BBE_UNPKSNX16( z10 );
    w11 = BBE_UNPKSNX16( z11 );

    // CQ(qX+qZY+eZ-rsh) <- CQ(qZ) + qX + qZY + eZ - rsh - qZ
    w00 = BBE_SLSNX40( w00, vsa0 );
    w01 = BBE_SLSNX40( w01, vsa0 );
    w10 = BBE_SLSNX40( w10, vsa0 );
    w11 = BBE_SLSNX40( w11, vsa0 );

    // CQ(qX+qZY+eZ-rsh) <- CQ(qX+qZY+eZ-rsh) - CQ(qX)*CQ(qZY)
    BBE_MULSNX16C( w00, x00, y00 );
    BBE_MULSNX16C( w01, x00, y01 );
    BBE_MULSNX16C( w10, x10, y00 );
    BBE_MULSNX16C( w11, x10, y01 );

    BBE_MULSNX16C( w00, x01, y10 );
    BBE_MULSNX16C( w01, x01, y11 );
    BBE_MULSNX16C( w10, x11, y10 );
    BBE_MULSNX16C( w11, x11, y11 );

    //--------------------------------------------------------------------------
    // Scale, round and save 8 2x2 matrix results.

    w00 = BBE_RNDSADJNX40( w00, vsa1 );
    w01 = BBE_RNDSADJNX40( w01, vsa1 );
    w10 = BBE_RNDSADJNX40( w10, vsa1 );
    w11 = BBE_RNDSADJNX40( w11, vsa1 );

    // CQ(qZY+eZ-rsh) <- CQ(qX+qZY+eZ-rsh) - qX w/ rounding
    z00 = BBE_PACKVNX40( w00, vsa1 );
    z01 = BBE_PACKVNX40( w01, vsa1 );
    z10 = BBE_PACKVNX40( w10, vsa1 );
    z11 = BBE_PACKVNX40( w11, vsa1 );

    BBE_SVNX16_IP( z00, Z_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( z01, Z_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( z10, Z_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( z11, Z_wr, +4*BBE_SIMD_WIDTH/2 );
  }
#else
  const xb_vecNx16 *          Z_rd;
        xb_vecNx16 * restrict Z_wr;
  const xb_vecNx16 *          X_rd;
  const xb_vecNx16 *          Y_rd;
  const xb_vecNx16 *          QZ_rd;
  const xb_vecNx16 *          QX_rd;
  const xb_vecNx16 *          EZ_rd;

  valign QZ_va, QX_va, EZ_va;

  xb_vecNx16 x00, x01, x10, x11;
  xb_vecNx16 y00, y01, y10, y11;
  xb_vecNx16 z00, z01, z10, z11;
  xb_vecNx40 w00, w01, w10, w11;

  xb_vecNx16 csh, qZ, qX, eZ;
  xb_vecNx16 t0;

  vsaN vsa0, vsa1;

  int l;

  NASSERT_ALIGN( Z, (2*BBE_SIMD_WIDTH) );
  NASSERT_ALIGN( X, (2*BBE_SIMD_WIDTH) );
  NASSERT_ALIGN( Y, (2*BBE_SIMD_WIDTH) );

  NASSERT( L%(BBE_SIMD_WIDTH/2)==0 );

  csh = BBE_MOVVA16( Qzy - rsh ); 

  Z_rd  = (const xb_vecNx16*)Z;
  Z_wr  = (      xb_vecNx16*)Z;
  X_rd  = (const xb_vecNx16*)X;
  Y_rd  = (const xb_vecNx16*)Y;
  QZ_rd = (const xb_vecNx16*)Qz;
  QX_rd = (const xb_vecNx16*)Qx;
  EZ_rd = (const xb_vecNx16*)Ez;

  QZ_va = BBE_LAVNX16_PP( QZ_rd );
  QX_va = BBE_LAVNX16_PP( QX_rd );
  EZ_va = BBE_LAVNX16_PP( EZ_rd );

  for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++ )
  {
    //--------------------------------------------------------------------------
    // Calculate pre-scaling shift amount for BBE_SIMD_WIDTH/2 complex 2x2 
    // accumulator matrices

    BBE_LAVNX16_XP( qZ, QZ_va, QZ_rd, BBE_SIMD_WIDTH );
    BBE_LAVNX16_XP( qX, QX_va, QX_rd, BBE_SIMD_WIDTH );
    BBE_LAVNX16_XP( eZ, EZ_va, EZ_rd, BBE_SIMD_WIDTH );

    // t0 <- qX + eZ
    t0 = BBE_ADDNX16( qX, eZ );
    // t0 <- qX + eZ - qZ
    t0 = BBE_SUBNX16( t0, qZ );
    // t0 <- qX + qZY + eZ - rsh - qZ
    t0 = BBE_ADDNX16( t0, csh );

    t0 = BBE_SHFLNX16I( t0, BBE_SHFLI_DOUBLE_1_LO );
    qX = BBE_SHFLNX16I( qX, BBE_SHFLI_DOUBLE_1_LO );

    vsa0 = BBE_MOVVSV( t0, 0 );
    vsa1 = BBE_MOVVSV( qX, 0 );

    //--------------------------------------------------------------------------
    // Perform MAS on BBE_SIMD_WIDTH/2 2x2 matrices, with accumulators pre-scaling.

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
                              
    // CQ(qZ)                 
    BBE_LVNX16_IP( z00, Z_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( z01, Z_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( z10, Z_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( z11, Z_rd, 2*BBE_SIMD_WIDTH );

    w00 = BBE_UNPKSNX16( z00 );
    w01 = BBE_UNPKSNX16( z01 );
    w10 = BBE_UNPKSNX16( z10 );
    w11 = BBE_UNPKSNX16( z11 );

    // CQ(qX+qZY+eZ-rsh) <- CQ(qZ) + qX + qZY + eZ - rsh - qZ
    w00 = BBE_SLSNX40( w00, vsa0 );
    w01 = BBE_SLSNX40( w01, vsa0 );
    w10 = BBE_SLSNX40( w10, vsa0 );
    w11 = BBE_SLSNX40( w11, vsa0 );

    // CQ(qX+qZY+eZ-rsh) <- CQ(qX+qZY+eZ-rsh) - CQ(qX)*CQ(qZY)
    BBE_MULSNX16C( w00, x00, y00 );
    BBE_MULSNX16C( w01, x00, y01 );
    BBE_MULSNX16C( w10, x10, y00 );
    BBE_MULSNX16C( w11, x10, y01 );

    BBE_MULSNX16C( w00, x01, y10 );
    BBE_MULSNX16C( w01, x01, y11 );
    BBE_MULSNX16C( w10, x11, y10 );
    BBE_MULSNX16C( w11, x11, y11 );

    //--------------------------------------------------------------------------
    // Scale, round and save BBE_SIMD_WIDTH/2 2x2 matrix results.

    w00 = BBE_RNDSADJNX40( w00, vsa1 );
    w01 = BBE_RNDSADJNX40( w01, vsa1 );
    w10 = BBE_RNDSADJNX40( w10, vsa1 );
    w11 = BBE_RNDSADJNX40( w11, vsa1 );

    // CQ(qZY+eZ-rsh) <- CQ(qX+qZY+eZ-rsh) - qX w/ rounding
    z00 = BBE_PACKVNX40( w00, vsa1 );
    z01 = BBE_PACKVNX40( w01, vsa1 );
    z10 = BBE_PACKVNX40( w10, vsa1 );
    z11 = BBE_PACKVNX40( w11, vsa1 );

    BBE_SVNX16_IP( z00, Z_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( z01, Z_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( z10, Z_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( z11, Z_wr, 2*BBE_SIMD_WIDTH );
  }
#endif
} // cmatinv4x4n_mas2x2bs_sr()

#endif
