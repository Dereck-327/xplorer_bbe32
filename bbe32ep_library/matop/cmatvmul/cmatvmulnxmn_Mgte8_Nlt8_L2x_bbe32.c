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
    M>=8 && N<8 && !(L&1)
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
DISCARD_FUN(void, cmatvmulnxmn_Mgte8_Nlt8_L2x,( 
                   void    *          pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q ))
#else
void cmatvmulnxmn_Mgte8_Nlt8_L2x( 
                   void    *          pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q )
{
        xb_vecNx16 * restrict Z_scr;
  const xb_vecNx16 *          Z_rd;
        xb_vecNx16 * restrict Z_wr;
        xb_vecNx16 * restrict X_scr;
  const xb_vecNx16 *          X_rd;
        xb_vecNx16 * restrict X_wr;
  const xb_vecNx16 *          Y_rd;

  int MN, l;
  int _M;

  NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);

  NASSERT( !(L&1) && !(M&3) && !(N&3) );

  NASSERT( M>=8 && N<8 );

  NASSERT( Q>=0 && Q<=16 );
  
  if (L<=0) return;
  MN = M*4;

  /* Round up to the next multiple of BBE_SIMD_WIDTH/2. */
  _M = ( ( M + (BBE_SIMD_WIDTH/2-1) ) & ~(BBE_SIMD_WIDTH/2-1) );

  /*
  * Partition the scratch memory.
  */

  {
    void * ptr = pScr;

    X_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)X_scr + (BBE_SIMD_WIDTH/8)*MN*4 );
    Z_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)Z_scr + (BBE_SIMD_WIDTH/8)*_M*4 );

  }

  /*
  * Compute BBE_SIMD_WIDTH/8 matrix-vector products at a time.
  */

  Y_rd = (const xb_vecNx16*)y;
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

    __Pragma( "no_reorder" );

    /*
    * Compute BBE_SIMD_WIDTH/8 matrix-vector products; store resulting vectors
    * to the scratch in block-streaming format: z[M/4][BBE_SIMD_WIDTH/8][4][2].
    */

    X_rd = X_scr;
    Z_wr = Z_scr;

    {
      xb_vecNx16 x0, x1, x2, x3;
      xb_vecNx16 y0, y1;
      xb_vecNx16 z0, p0;

      xb_vecNx40 w0;

      vsaN vsa0;

      int i;

      vsa0 = BBE_MOVVSA32( Q );

      BBE_LVNX16_IP( p0, Y_rd, 4*BBE_SIMD_WIDTH/2 );

      y0 = BBE_SHFLNX16I( p0, BBE_SHFLI_DUPLICATE_4_EVEN );
      y1 = BBE_SHFLNX16I( p0, BBE_SHFLI_DUPLICATE_4_ODD  );

      __Pragma( "loop_count min=1" );
      __Pragma( "ymemory( X_rd )" );
      for ( i=0; i<M/4; i++ )
      {
        w0 = BBE_MOVWA32( (1UL<<Q) >> 1  );

        BBE_LVNX16_IP( x0, X_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x1, X_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x2, X_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x3, X_rd, 4*BBE_SIMD_WIDTH/2 );

        BBE_DSELNX16I( x2, x0, x2, x0, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( x3, x1, x3, x1, BBE_DSELI_INTERLEAVE_4 );

        BBE_MULANX16PC_0( w0, x0, y0 );
        BBE_MULANX16PC_1( w0, x1, y0 );
        BBE_MULANX16PC_0( w0, x2, y1 );
        BBE_MULANX16PC_1( w0, x3, y1 );

        z0 = BBE_PACKVNX40( w0, vsa0 );

        BBE_SVNX16_IP( z0, Z_wr, 4*BBE_SIMD_WIDTH/2 );
      }
    }

    __Pragma( "no_reorder" );

    /*
    * Convert BBE_SIMD_WIDTH/8 resulting vectors to block format: 
    * z[M/4][BBE_SIMD_WIDTH/8][4][2] => Z[BBE_SIMD_WIDTH/8][_M][2].
    */

    Z_rd = (const xb_vecNx16*)Z_scr;
    Z_wr = (      xb_vecNx16*)( (uintptr_t)z + l*BBE_SIMD_WIDTH/8*_M*4 );

    {
      xb_vecNx16 a0, a1, b0, b1;

      int m;

      __Pragma( "loop_count min=1" );
      for ( m=0; m<_M/(BBE_SIMD_WIDTH/2); m++ )
      {
        BBE_LVNX16_IP( a0, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a1, Z_rd, 4*BBE_SIMD_WIDTH/2 );

        b0 = BBE_SELNX16I( a1, a0, BBE_SELI_EXTRACT_LO_HALVES );
        b1 = BBE_SELNX16I( a1, a0, BBE_SELI_EXTRACT_HI_HALVES );

        BBE_SVNX16_X ( b1, Z_wr, 4*_M               );
        BBE_SVNX16_XP( b0, Z_wr, 4*BBE_SIMD_WIDTH/2 );
      }
    }
  }

} /* cmatvmulnxmn_Mgte8_Nlt8_L2x() */

/* Return the scratch area size, in bytes. */
size_t cmatvmulnxmn_Mgte8_Nlt8_L2x_getScratchSize(int N, int M)
{
  int _M = ((M + (BBE_SIMD_WIDTH / 2 - 1)) & ~(BBE_SIMD_WIDTH / 2 - 1));
  size_t  sz = 0;
  sz = (4 * M + _M) * 2 * 2 * sizeof(int16_t);
  return sz;
} /* cmatvmulnxmn_Mgte8_Nlt8_L2x_getScratchSize() */
#endif
