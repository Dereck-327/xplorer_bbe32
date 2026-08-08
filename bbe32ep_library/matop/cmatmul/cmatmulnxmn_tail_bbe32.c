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
 * Complex Matrix-Matrix/Matrix-Vector Multiply
 * C code optimized for BBE32
 */

/*  
    Optimized code for matrix multiplication
	Integrit, 2006-2016

    specialized function for:
    L = 1
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matop.h"
#include "cmatmulnxmn_common.h"

#if !HAVE_MULPC
DISCARD_FUN(void, cmatmulnxmn_tail, (void * pScr,
                    complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int N, int M, int Q ))
#else
void cmatmulnxmn_tail(void* pScr, complex_fract16* z, const complex_fract16* x, const complex_fract16* y, int N, int M, int Q)
{
          xb_vecNx16 * restrict Z_scr;
  const xb_vecNx16 *          Z_rd;
        xb_vecNx16 * restrict Z_wr;
        xb_vecNx16 * restrict XY_scr;
  const xb_vecNx16 *          X_rd;
  const xb_vecNx16 *          Y_rd;
        xb_vecNx16 * restrict XY_wr;

  int MN, MM;

  NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);

  NASSERT( N%4==0 && M%4==0 );

  NASSERT( Q>=0 && Q<=16 );

  MN = M*N; MM = M*M;

  /*
  * Partition the scratch memory
  */

  {
    void * ptr = pScr;

    Z_scr  = (xb_vecNx16*)ptr;
    ptr    = (void*)( (uintptr_t)Z_scr + MM*2*4 );
    XY_scr = (xb_vecNx16*)ptr;
    ptr    = (void*)( (uintptr_t)XY_scr + MN*2*4 );
  }

  /*
  * Interleave left- and right-hand matrices x[], y[] so that a vector of
  * BBE_SIMD_WIDTH 16-bit words contains 4 complex elements of x[] and 4
  * complex elements of y[].
  */

  X_rd = (const xb_vecNx16*)x;
  Y_rd = (const xb_vecNx16*)y;

  XY_wr = XY_scr;

  {
    xb_vecNx16 a0, a1, b0, b1;

    int n;

    for ( n=0; n<MN/(BBE_SIMD_WIDTH/2); n++ )
    {
      BBE_LVNX16_IP( a0, X_rd, 4*BBE_SIMD_WIDTH/2 );
      BBE_LVNX16_IP( a1, Y_rd, 4*BBE_SIMD_WIDTH/2 );

      BBE_DSELNX16I( b1, b0, a1, a0, BBE_DSELI_INTERLEAVE_2 );

      BBE_SVNX16_IP( b0, XY_wr, 4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( b1, XY_wr, 4*BBE_SIMD_WIDTH/2 );
    }
  }

  __Pragma( "no_reorder" );

  /*
  * Compute 4x4 sub-blocks of the product matrix z[]. Each complex value of
  * the resulting matrix is padded with zero so that 4 consecutive elements
  * occupy a whole vector of BBE_SIMD_WIDTH 16-bit words.
  */

  {
    xb_vecNx16 x0, x1, x2, x3;
    xb_vecNx16 y0, y1, y2, y3;
    xb_vecNx16 z0, z1, z2, z3;

    xb_vecNx16 x00, x01, x10, x11, x20, x21, x30, x31;

    xb_vecNx40 w0, w1, w2, w3;

    vsaN vsa0;

    int i, j, k;

    vsa0 = BBE_MOVVSA32( Q );

    Z_wr = (xb_vecNx16*)Z_scr;

    for ( i=0; i<M/4; i++ )
    {
      X_rd = (const xb_vecNx16*)( (uintptr_t)XY_scr + 4*i*N*2*4 );

      for ( j=0; j<M/4; j++ )
      {
        Y_rd = (const xb_vecNx16*)( (uintptr_t)XY_scr + 4*j*2*4 );

        w0 = w1 = w2 = w3 = BBE_MOVWA32( (1UL<<Q) >> 1 );

        __Pragma( "loop_count min=1" );
        for ( k=0; k<N/4; k++ )
        {
          BBE_LVNX16_XP( x0, X_rd, N*2*4 );
          BBE_LVNX16_XP( x1, X_rd, N*2*4 );
          BBE_LVNX16_XP( x2, X_rd, N*2*4 );

          BBE_LVNX16_XP( x3, X_rd, -3*N*2*4 + 4*BBE_SIMD_WIDTH/2 );

          x0 = BBE_SELNX16I( x1, x0, BBE_SELI_EXTRACT_2_OF_4_OFF_0 );
          x2 = BBE_SELNX16I( x3, x2, BBE_SELI_EXTRACT_2_OF_4_OFF_0 );

          x00 = BBE_SHFLNX16I( x0, BBE_SHFLI_REP_0X4 );
          x01 = BBE_SHFLNX16I( x0, BBE_SHFLI_REP_1X4 );
          x10 = BBE_SHFLNX16I( x0, BBE_SHFLI_REP_2X4 );
          x11 = BBE_SHFLNX16I( x0, BBE_SHFLI_REP_3X4 );

          x20 = BBE_SHFLNX16I( x2, BBE_SHFLI_REP_0X4 );
          x21 = BBE_SHFLNX16I( x2, BBE_SHFLI_REP_1X4 );
          x30 = BBE_SHFLNX16I( x2, BBE_SHFLI_REP_2X4 );
          x31 = BBE_SHFLNX16I( x2, BBE_SHFLI_REP_3X4 );

          y1 = BBE_LVNX16_X( Y_rd, 1*M*2*4 );
          y2 = BBE_LVNX16_X( Y_rd, 2*M*2*4 );
          y3 = BBE_LVNX16_X( Y_rd, 3*M*2*4 );

          BBE_LVNX16_XP( y0, Y_rd, 4*M*2*4 );

          y0 = BBE_SELNX16I( y1, y0, BBE_SELI_INTERLEAVE_2_ODD );
          y1 = BBE_SELNX16I( y3, y2, BBE_SELI_INTERLEAVE_2_ODD );

          BBE_MULANX16PC_0( w0, x00, y0 );
          BBE_MULANX16PC_0( w0, x01, y1 );

          BBE_MULANX16PC_0( w1, x10, y0 );
          BBE_MULANX16PC_0( w1, x11, y1 );

          BBE_MULANX16PC_0( w2, x20, y0 );
          BBE_MULANX16PC_0( w2, x21, y1 );

          BBE_MULANX16PC_0( w3, x30, y0 );
          BBE_MULANX16PC_0( w3, x31, y1 );
        }

        X_rd = (const xb_vecNx16*)( (uintptr_t)X_rd - (N/4)*4*BBE_SIMD_WIDTH/2 );

        z0 = BBE_PACKVNX40( w0, vsa0 );
        z1 = BBE_PACKVNX40( w1, vsa0 );
        z2 = BBE_PACKVNX40( w2, vsa0 );
        z3 = BBE_PACKVNX40( w3, vsa0 );

        BBE_SVNX16_XP( z0, Z_wr, M*2*4 );
        BBE_SVNX16_XP( z1, Z_wr, M*2*4 );
        BBE_SVNX16_XP( z2, Z_wr, M*2*4 );

        BBE_SVNX16_XP( z3, Z_wr, -3*M*2*4 + 4*BBE_SIMD_WIDTH/2 );
      }

      Z_wr = (xb_vecNx16*)( (uintptr_t)Z_wr + 3*M*2*4 );
    }
  }

  __Pragma( "no_reorder" );

  /*
  * Remove padding zeros and store the product matrix z[] to the output array
  */

  Z_rd = (const xb_vecNx16*)Z_scr;
  Z_wr = (      xb_vecNx16*)z;

  {
    xb_vecNx16 a0, a1;
    xb_vecNx16 b0;

    int n;

    for ( n=0; n<MM/(BBE_SIMD_WIDTH/2); n++ )
    {
      BBE_LVNX16_IP( a0, Z_rd, 4*BBE_SIMD_WIDTH/2 );
      BBE_LVNX16_IP( a1, Z_rd, 4*BBE_SIMD_WIDTH/2 );

      b0 = BBE_SELNX16I( a1, a0, BBE_SELI_EXTRACT_2_OF_4_OFF_0 );

      BBE_SVNX16_IP( b0, Z_wr, 4*BBE_SIMD_WIDTH/2 );
    }
  }
}
#endif


