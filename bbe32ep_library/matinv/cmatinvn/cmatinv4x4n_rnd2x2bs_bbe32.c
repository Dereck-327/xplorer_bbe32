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
                           int L, int rsh )
{
  const xb_vecNx16 *          X_rd;
        xb_vecNx16 * restrict X_wr;
  const xb_vecNx16 *          SA_rd;
  const xb_vecNx16 *          SB_rd;

  valign SA_va, SB_va;

  xb_vecNx16 x0, x1, x2, x3;
  xb_vecNx40 w0, w1, w2, w3;

  xb_vecNx16 t0, t1;
  xb_vecNx16 c1, csh;

  vsaN vsa0;

  int l;

  NASSERT_ALIGN( X, (2*BBE_SIMD_WIDTH) );

  NASSERT( L%(BBE_SIMD_WIDTH/2)==0 );

  csh = BBE_MOVVA16( rsh );

  c1 = BBE_MOVVINX16( BBE_MOVVI_INT16_1 );

  X_rd  = (const xb_vecNx16*)X;
  X_wr  = (      xb_vecNx16*)X;
  SA_rd = (const xb_vecNx16*)sA;
  SB_rd = (const xb_vecNx16*)sB;

  SA_va = BBE_LAVNX16_PP( SA_rd );
  SB_va = BBE_LAVNX16_PP( SB_rd );

  __Pragma( "concurrent" );
  for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++ )
  {
    BBE_LAVNX16_XP( t0, SA_va, SA_rd, BBE_SIMD_WIDTH );
    BBE_LAVNX16_XP( t1, SB_va, SB_rd, BBE_SIMD_WIDTH );

    // t0 <- rsh - sA
    t0 = BBE_SUBNX16( csh, t0 );
    // t0 <- rsh - sA + sB
    t0 = BBE_ADDNX16( t0, t1 );
    t0 = BBE_SHFLNX16I( t0, BBE_SHFLI_DOUBLE_1_LO );

    vsa0 = BBE_MOVVSV( t0, 0 );

    BBE_LVNX16_IP( x0, X_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x1, X_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x2, X_rd, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x3, X_rd, 2*BBE_SIMD_WIDTH );

    w0 = BBE_MULRNX16( x0, c1, vsa0 );
    w1 = BBE_MULRNX16( x1, c1, vsa0 );
    w2 = BBE_MULRNX16( x2, c1, vsa0 );
    w3 = BBE_MULRNX16( x3, c1, vsa0 );

    x0 = BBE_PACKVNX40( w0, vsa0 );
    x1 = BBE_PACKVNX40( w1, vsa0 );
    x2 = BBE_PACKVNX40( w2, vsa0 );
    x3 = BBE_PACKVNX40( w3, vsa0 );

    BBE_SVNX16_IP( x0, X_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( x1, X_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( x2, X_wr, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( x3, X_wr, 2*BBE_SIMD_WIDTH );
  }

} // cmatinv4x4n_rnd2x2bs()

#endif
