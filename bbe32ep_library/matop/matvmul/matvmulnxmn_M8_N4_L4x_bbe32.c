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
  NatureDSP_Baseband library. Matrix operations part.
    Real Matrix-Vector Multiply; Block Order, MxN * Nx1 -> Mx1
    M==8 && N==4 && !(L&3)
    C code optimized for BBE32EP
  IntegrIT, 2006-2013
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#include "matvmulnxmn_common.h"
#if !(HAVE_MULPC && HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, matvmulnxmn_M8_N4_L4x, (void * pScr,
                   int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int N, int M, int L, int Q ))
#else
// M==8 && N==4 && !(L&3)
void matvmulnxmn_M8_N4_L4x( void    *          pScr,
                            int16_t * restrict z, 
                      const int16_t * restrict x, 
                      const int16_t * restrict y, 
                      int N, int M, int L, int Q )
{
  const xb_vecNx16 *          X;
  const xb_vecNx16 *          Y;
        xb_vecNx16 * restrict Z;

  xb_vecNx16 a0;
  xb_vecNx16 x0, x1, x2, x3;
  xb_vecNx16 x4, x5, x6, x7;
  xb_vecNx16 y0, y1, y2, y3;
  xb_vecNx16 t0, t1, t2, t3;
  xb_vecNx16 z0, z1;
  xb_vecNx40 w0, w1, w2, w3;

  vselN sel0;
  vsaN  vsa0;

  int l;

  static const int16_t ALIGN(32) cst[][BBE_SIMD_WIDTH] = {
    { 0,4,8,12,2,6,10,14,16,20,24,28,18,22,26,30 },
  };

  NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);

  NASSERT( !(L&3) );

  NASSERT( M==8 && N==4 );

  NASSERT( Q>=0 && Q<=16 );

  x0 = BBE_LVNX16_I( (const xb_vecNx16*)cst, 0 );

  sel0 = BBE_MOVVSELNX16( x0, 0 );

  vsa0 = BBE_MOVVSA32( Q );

  /*
  * Compute BBE_SIMD_WIDTH/4 matrix-vector products at a time.
  */

  X = (const xb_vecNx16*)x;
  Y = (const xb_vecNx16*)y;
  Z = (      xb_vecNx16*)z;
  __Pragma("loop_count min=1");
  for ( l=0; l<L/(BBE_SIMD_WIDTH/4); l++ )
  {
    BBE_LVNX16_IP( x0, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x1, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x2, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x3, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x4, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x5, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x6, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( x7, X, 2*BBE_SIMD_WIDTH );

    BBE_LVNX16_IP( a0, Y, 2*BBE_SIMD_WIDTH );

    a0 = BBE_CONJSNX16C( a0 );

    y0 = BBE_SHFLNX16I( a0, BBE_SHFLI_REP_0X4 );
    y1 = BBE_SHFLNX16I( a0, BBE_SHFLI_REP_1X4 );
    y2 = BBE_SHFLNX16I( a0, BBE_SHFLI_REP_2X4 );
    y3 = BBE_SHFLNX16I( a0, BBE_SHFLI_REP_3X4 );

    w0 = BBE_MULRNX16PC_0( x0, y0, vsa0 );
    w1 = BBE_MULRNX16PC_0( x2, y1, vsa0 );
    w2 = BBE_MULRNX16PC_0( x4, y2, vsa0 );
    w3 = BBE_MULRNX16PC_0( x6, y3, vsa0 );

    BBE_MULANX16PC_1( w0, x1, y0 );
    BBE_MULANX16PC_1( w1, x3, y1 );
    BBE_MULANX16PC_1( w2, x5, y2 );
    BBE_MULANX16PC_1( w3, x7, y3 );

    t0 = BBE_PACKVNX40( w0, vsa0 );
    t1 = BBE_PACKVNX40( w1, vsa0 );
    t2 = BBE_PACKVNX40( w2, vsa0 );
    t3 = BBE_PACKVNX40( w3, vsa0 );

    z0 = BBE_SELNX16( t1, t0, sel0 );
    z1 = BBE_SELNX16( t3, t2, sel0 );

    BBE_SVNX16_IP( z0, Z, 2*BBE_SIMD_WIDTH );
    BBE_SVNX16_IP( z1, Z, 2*BBE_SIMD_WIDTH );
  }

} // matvmulnxmn_M8_N4_L4x()
#endif
