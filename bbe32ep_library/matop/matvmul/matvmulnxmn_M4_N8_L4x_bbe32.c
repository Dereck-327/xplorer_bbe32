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
    M==4 && N==8 && !(L&3)
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
DISCARD_FUN(void, matvmulnxmn_M4_N8_L4x, (void    *          pScr,
                   int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int N, int M, int L, int Q ))
#else

// M==4 && N==8 && !(L&3)
void matvmulnxmn_M4_N8_L4x( void    *          pScr,
                            int16_t * restrict z, 
                      const int16_t * restrict x, 
                      const int16_t * restrict y, 
                      int N, int M, int L, int Q )
{
  const xb_vecNx16 *          X;
  const xb_vecNx16 *          Y;
        xb_vecNx16 * restrict Z;

  xb_vecNx16 a0, a1, a2, a3, a4, a5, a6, a7;
  xb_vecNx16 b0, b1;
  xb_vecNx16 x0, x1, x2, x3, x4, x5, x6, x7;
  xb_vecNx16 y0, y1, y2, y3, y4, y5, y6, y7;
  xb_vecNx16 t0, t1, z0;
  xb_vecNx40 w0, w1;

  vselN sel0;
  vsaN  vsa0;

  int l;

  static const int16_t ALIGN(32) cst[][BBE_SIMD_WIDTH] = {
    { 0, 4, 8,12,16,20,24,28, 2, 6,10,14,18,22,26,30 },
  };

  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);

  NASSERT( !(L&3) );

  NASSERT( M==4 && N==8 );

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
    BBE_LVNX16_IP( a0, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( a1, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( a2, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( a3, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( a4, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( a5, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( a6, X, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( a7, X, 2*BBE_SIMD_WIDTH );

    BBE_LVNX16_IP( b0, Y, 2*BBE_SIMD_WIDTH );
    BBE_LVNX16_IP( b1, Y, 2*BBE_SIMD_WIDTH );

    b0 = BBE_CONJSNX16C( b0 );
    b1 = BBE_CONJSNX16C( b1 );

    x0 = BBE_SELNX16I( a1, a0, BBE_SELI_EXTRACT_4_OF_8_OFF_0 );
    x1 = BBE_SELNX16I( a1, a0, BBE_SELI_EXTRACT_4_OF_8_OFF_4 );
    x2 = BBE_SELNX16I( a3, a2, BBE_SELI_EXTRACT_4_OF_8_OFF_0 );
    x3 = BBE_SELNX16I( a3, a2, BBE_SELI_EXTRACT_4_OF_8_OFF_4 );
    x4 = BBE_SELNX16I( a5, a4, BBE_SELI_EXTRACT_4_OF_8_OFF_0 );
    x5 = BBE_SELNX16I( a5, a4, BBE_SELI_EXTRACT_4_OF_8_OFF_4 );
    x6 = BBE_SELNX16I( a7, a6, BBE_SELI_EXTRACT_4_OF_8_OFF_0 );
    x7 = BBE_SELNX16I( a7, a6, BBE_SELI_EXTRACT_4_OF_8_OFF_4 );

    y0 = BBE_SHFLNX16I( b0, BBE_SHFLI_REP_0X4 );
    y1 = BBE_SHFLNX16I( b0, BBE_SHFLI_REP_1X4 );
    y2 = BBE_SHFLNX16I( b0, BBE_SHFLI_REP_2X4 );
    y3 = BBE_SHFLNX16I( b0, BBE_SHFLI_REP_3X4 );
    y4 = BBE_SHFLNX16I( b1, BBE_SHFLI_REP_0X4 );
    y5 = BBE_SHFLNX16I( b1, BBE_SHFLI_REP_1X4 );
    y6 = BBE_SHFLNX16I( b1, BBE_SHFLI_REP_2X4 );
    y7 = BBE_SHFLNX16I( b1, BBE_SHFLI_REP_3X4 );

    w0 = BBE_MULRNX16PC_0( x0, y0, vsa0 );
    BBE_MULANX16PC_0(  w0, x1, y1 );
    BBE_MULANX16PC_1(  w0, x4, y4 );
    BBE_MULANX16PC_1(  w0, x5, y5 );

    w1 = BBE_MULRNX16PC_0( x2, y2, vsa0 );
    BBE_MULANX16PC_0(  w1, x3, y3 );
    BBE_MULANX16PC_1(  w1, x6, y6 );
    BBE_MULANX16PC_1(  w1, x7, y7 );

    t0 = BBE_PACKVNX40( w0, vsa0 );
    t1 = BBE_PACKVNX40( w1, vsa0 );

    z0 = BBE_SELNX16( t1, t0, sel0 );

    BBE_SVNX16_IP( z0, Z, 2*BBE_SIMD_WIDTH );
  }

} /* matvmulnxmn_M4_N8_L4x() */

#endif
