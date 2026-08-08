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
    Complex Matrix-Vector Multiply; Block Order, MxN * Nx1 -> Mx1
    M<8 && N>=8 && !(L&1)
    C code optimized for BBE32EP
  IntegrIT, 2006-2013
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#include "cmatvmulnxmn_common.h"

#if !(HAVE_MULPC && 1)
DISCARD_FUN(void, cmatvmulnxmn_Mlt8_Ngte8_L2x,( 
                   void    *          pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q ))
#else
void cmatvmulnxmn_Mlt8_Ngte8_L2x( 
                   void    *          pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q )
{
        xb_vecNx16 * restrict Z_wr;
        xb_vecNx16 * restrict X_scr;
  const xb_vecNx16 *          X_rd;
        xb_vecNx16 * restrict X_wr;
        xb_vecNx16 * restrict Y_scr;
  const xb_vecNx16 *          Y_rd;
        xb_vecNx16 * restrict Y_wr;

  int MN, l;
  int _N;

  NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);

  NASSERT( !(L&1) && !(M&3) && !(N&3) );

  NASSERT( M<8 && N>=8 );

  NASSERT( Q>=0 && Q<=16);

  if (L<=0) return;
  MN = 4*N;

  /* Round up to the next multiple of BBE_SIMD_WIDTH/2. */
  _N = ( ( N + (BBE_SIMD_WIDTH/2-1) ) & ~(BBE_SIMD_WIDTH/2-1) );

  /*
  * Partition the scratch memory.
  */

  {
    void * ptr = pScr;

    X_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)X_scr + (BBE_SIMD_WIDTH/8)*MN*4 );
    Y_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)Y_scr + (BBE_SIMD_WIDTH/8)*_N*4 );

  }

  /*
  * Compute BBE_SIMD_WIDTH/8 matrix-vector products at a time.
  */

  Z_wr = (xb_vecNx16*)z;
  __Pragma("loop_count min=1");
  for ( l=0; l<L/(BBE_SIMD_WIDTH/8); l++ )
  {
    /*
    * Convert BBE_SIMD_WIDTH/8 left-hand matrices x[M*N][2] to block-streaming
    * format: x[BBE_SIMD_WIDTH/8][M*N][2] => x[M*N/2][BBE_SIMD_WIDTH/8][2][2].
    */

    X_rd = (const xb_vecNx16*)( (uintptr_t)x + l*(BBE_SIMD_WIDTH/8)*MN*4 );
    X_wr = (      xb_vecNx16*)X_scr;

    {
      xb_vecNx16 a0, a1;

      int n;

      __Pragma( "loop_count min=2, factor=2" );
      for ( n=0; n<MN/(BBE_SIMD_WIDTH/2); n++ )
      {
        a1 = BBE_LVNX16_X( X_rd, 4*MN               );
        BBE_LVNX16_IP( a0, X_rd, 4*BBE_SIMD_WIDTH/2 );

        BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_INTERLEAVE_4 );

        BBE_SVNX16_IP( a0, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a1, X_wr, 4*BBE_SIMD_WIDTH/2 );
      }
    }

    /*
    * Convert BBE_SIMD_WIDTH/8 right-hand vectors y[N][2] to block-streaming
    * format: y[BBE_SIMD_WIDTH/8][N][2] => y[_N/2][BBE_SIMD_WIDTH/8][2][2].
    */

    Y_rd = (const xb_vecNx16*)( (uintptr_t)y + l*(BBE_SIMD_WIDTH/8)*_N*4 );
    Y_wr = (      xb_vecNx16*)Y_scr;

    {
      xb_vecNx16 a0, a1;

      int n;

      for ( n=0; n<_N/(BBE_SIMD_WIDTH/2); n++ )
      {
        a1 = BBE_LVNX16_X( Y_rd, 4*_N );
        BBE_LVNX16_IP( a0, Y_rd, 4*BBE_SIMD_WIDTH/2 );

        BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_INTERLEAVE_4 );

        BBE_SVNX16_IP( a0, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a1, Y_wr, 4*BBE_SIMD_WIDTH/2 );
      }
    }

    __Pragma( "no_reorder" );

    /*
    * Compute BBE_SIMD_WIDTH/8 matrix-vector products; store resulting vectors
    * to the output array.
    */

    X_rd = X_scr;
    Y_rd = Y_scr;

    {
      const xb_vecNx16 * X0_rd;
      const xb_vecNx16 * X1_rd;
      const xb_vecNx16 * X2_rd;
      const xb_vecNx16 * X3_rd;

      xb_vecNx16 x0, x1, x2, x3;
      xb_vecNx16 y0, y1;
      xb_vecNx16 z0, p0;

      xb_vecNx40 w0;

      vsaN vsa0;

      int j;

      vsa0 = BBE_MOVVSA32( Q );

      X0_rd = (const xb_vecNx16*)( (uintptr_t)X_rd + 0*BBE_SIMD_WIDTH/8*N*4 );
      X1_rd = (const xb_vecNx16*)( (uintptr_t)X_rd + 1*BBE_SIMD_WIDTH/8*N*4 );
      X2_rd = (const xb_vecNx16*)( (uintptr_t)X_rd + 2*BBE_SIMD_WIDTH/8*N*4 );
      X3_rd = (const xb_vecNx16*)( (uintptr_t)X_rd + 3*BBE_SIMD_WIDTH/8*N*4 );

      w0 = BBE_MOVWA32( (1UL<<Q) >> 1  );

      __Pragma( "loop_count min=1" );
      for ( j=0; j<N/4; j++ )
      {
        BBE_LVNX16_IP( x0, X0_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x1, X1_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x2, X2_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x3, X3_rd, 4*BBE_SIMD_WIDTH/2 );

        BBE_DSELNX16I( x2, x0, x2, x0, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( x3, x1, x3, x1, BBE_DSELI_INTERLEAVE_4 );

        BBE_LVNX16_IP( p0, Y_rd, 4*BBE_SIMD_WIDTH/2 );

        y0 = BBE_SHFLNX16I( p0, BBE_SHFLI_DOUBLE_4_LO );
        y1 = BBE_SHFLNX16I( p0, BBE_SHFLI_DOUBLE_4_HI );

        BBE_MULANX16PC_0( w0, x0, y0 );
        BBE_MULANX16PC_1( w0, x1, y0 );
        BBE_MULANX16PC_0( w0, x2, y1 );
        BBE_MULANX16PC_1( w0, x3, y1 );
      }

      z0 = BBE_PACKVNX40( w0, vsa0 );

      BBE_SVNX16_IP( z0, Z_wr, 4*BBE_SIMD_WIDTH/2 );
    }
  }

} /* cmatvmulnxmn_Mlt8_Ngte8_L2x() */

/* Return the scratch area size, in bytes. */
size_t cmatvmulnxmn_Mlt8_Ngte8_L2x_getScratchSize(int N, int M)
{
  int _N = ((N + (BBE_SIMD_WIDTH / 2 - 1)) & ~(BBE_SIMD_WIDTH / 2 - 1));
  size_t  sz = 0;
  sz = (4*N + _N ) * 4* sizeof(int16_t);
  return sz;
} /* cmatvmulnxmn_Mlt8_Ngte8_L2x_getScratchSize() */
#endif
