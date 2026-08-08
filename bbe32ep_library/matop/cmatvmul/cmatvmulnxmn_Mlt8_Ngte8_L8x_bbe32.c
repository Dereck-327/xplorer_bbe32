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
/* Cadence, Inc. under Terms and Condition of a Software Lic ense Agreement  */
/* between Cadence, Inc. and IntegrIT, Ltd.                                 */
/* ------------------------------------------------------------------------ */
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
  NatureDSP_Baseband library. Matrix operations part.
    Complex Matrix-Vector Multiply; Block Order, MxN * Nx1 -> Mx1
    M<8 && N>=8 && !(L&7)
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
DISCARD_FUN(void, cmatvmulnxmn_Mlt8_Ngte8_L8x,( 
                   void    *          pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q ))
#else
void cmatvmulnxmn_Mlt8_Ngte8_L8x( 
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
        xb_vecNx16 * restrict Y_scr;
  const xb_vecNx16 *          Y_rd;
        xb_vecNx16 * restrict Y_wr;

  int MN, l;
  int _N;

  /* M<=8 && !(M&3) && M>0 => M==4 */
  const int _M = 4;

  NASSERT_ALIGN( pScr, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN(z, 2*BBE_SIMD_WIDTH );

  NASSERT( !(L&7) && !(M&3) && !(N&3) );

  NASSERT( M<8 && N>=8 );

  NASSERT( Q>=0 && Q<=16 );

  if (L<=0) return;
  MN = M*N;

  /* Round up to the next multiple of BBE_SIMD_WIDTH/2. */
  _N = ( ( N + (BBE_SIMD_WIDTH/2-1) ) & ~(BBE_SIMD_WIDTH/2-1) );

  /*
  * Partition the scratch memory.
  */

  {
    void * ptr = pScr;

    X_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)X_scr + (BBE_SIMD_WIDTH/2)*MN*4 );
    Y_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)Y_scr + (BBE_SIMD_WIDTH/2)*_N*4 );
    Z_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)Z_scr + (BBE_SIMD_WIDTH/2)*_M*4 );
  }

  /*
  * Compute BBE_SIMD_WIDTH/2 matrix-vector products at a time.
  */
  __Pragma("loop_count min=1");
  for ( l=0; l<L/(BBE_SIMD_WIDTH/2); l++ )
  {
    /*
    * Convert BBE_SIMD_WIDTH/2 left-hand matrices x[M*N][2] to streaming
    * format: x[BBE_SIMD_WIDTH/2][M*N][2] => x[M*N][BBE_SIMD_WIDTH/2][2].
    */

    X_rd = (const xb_vecNx16*)( (uintptr_t)x + l*BBE_SIMD_WIDTH/2*MN*4 );
    X_wr = (      xb_vecNx16*)X_scr;

    {
      xb_vecNx16 a0, a1, a2, a3, a4, a5, a6, a7;

      int n;

      __Pragma( "loop_count min=2" );
      for ( n=0; n<MN/(BBE_SIMD_WIDTH/2); n++ )
      {
        BBE_LVNX16_XP( a0, X_rd, MN*4 );
        BBE_LVNX16_XP( a1, X_rd, MN*4 );
        BBE_LVNX16_XP( a2, X_rd, MN*4 );
        BBE_LVNX16_XP( a3, X_rd, MN*4 );
        BBE_LVNX16_XP( a4, X_rd, MN*4 );
        BBE_LVNX16_XP( a5, X_rd, MN*4 );
        BBE_LVNX16_XP( a6, X_rd, MN*4 );

        BBE_LVNX16_XP( a7, X_rd, -7*MN*4 + 4*BBE_SIMD_WIDTH/2 );

        /*
        * Complex-valued 8x8 matrix transpose.
        */

        BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a3, a2, a3, a2, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a5, a4, a5, a4, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a7, a6, a7, a6, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_DSELNX16I( a2, a0, a2, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a6, a4, a6, a4, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a3, a1, a3, a1, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a7, a5, a7, a5, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_DSELNX16I( a4, a0, a4, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a5, a1, a5, a1, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a6, a2, a6, a2, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a7, a3, a7, a3, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_SVNX16_IP( a0, X_wr, +4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a1, X_wr, +4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a2, X_wr, +4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a3, X_wr, +4*BBE_SIMD_WIDTH/2 );

        BBE_SVNX16_IP( a4, X_wr, +4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a5, X_wr, +4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a6, X_wr, +4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a7, X_wr, +4*BBE_SIMD_WIDTH/2 );
      }
    }

    /*
    * Convert BBE_SIMD_WIDTH/2 right-hand vectors y[N][2] to streaming format:
    * y[BBE_SIMD_WIDTH/2][_N][2] => y[N][BBE_SIMD_WIDTH/2][2].
    */

    Y_rd = (const xb_vecNx16*)( (uintptr_t)y + l*BBE_SIMD_WIDTH/2*_N*4 );
    Y_wr = (      xb_vecNx16*)Y_scr;

    {
      xb_vecNx16 a0, a1, a2, a3, a4, a5, a6, a7;

      int n;

      __Pragma( "loop_count min=1" );
      for ( n=0; n<_N/(BBE_SIMD_WIDTH/2); n++ )
      {
        BBE_LVNX16_XP( a0, Y_rd, _N*4 );
        BBE_LVNX16_XP( a1, Y_rd, _N*4 );
        BBE_LVNX16_XP( a2, Y_rd, _N*4 );
        BBE_LVNX16_XP( a3, Y_rd, _N*4 );
        BBE_LVNX16_XP( a4, Y_rd, _N*4 );
        BBE_LVNX16_XP( a5, Y_rd, _N*4 );
        BBE_LVNX16_XP( a6, Y_rd, _N*4 );

        BBE_LVNX16_XP( a7, Y_rd, -7*_N*4 + 4*BBE_SIMD_WIDTH/2 );

        /*
        * Complex-valued 8x8 matrix transpose.
        */

        BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a3, a2, a3, a2, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a5, a4, a5, a4, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a7, a6, a7, a6, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_DSELNX16I( a2, a0, a2, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a6, a4, a6, a4, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a3, a1, a3, a1, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a7, a5, a7, a5, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_DSELNX16I( a4, a0, a4, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a5, a1, a5, a1, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a6, a2, a6, a2, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a7, a3, a7, a3, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_SVNX16_IP( a0, Y_wr, +4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a1, Y_wr, +4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a2, Y_wr, +4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a3, Y_wr, +4*BBE_SIMD_WIDTH/2 );
                        
        BBE_SVNX16_IP( a4, Y_wr, +4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a5, Y_wr, +4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a6, Y_wr, +4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a7, Y_wr, +4*BBE_SIMD_WIDTH/2 );
      }
    }

    __Pragma( "no_reorder" );

    /*
    * Compute BBE_SIMD_WIDTH/2 matrix-vector products; store resulting vectors
    * to the scratch in streaming format: z[M][BBE_SIMD_WIDTH/2][2].
    */

    X_rd = X_scr;
    Y_rd = Y_scr;
    Z_wr = Z_scr;

    {
      const xb_vecNx16 * X0_rd;
      const xb_vecNx16 * X1_rd;
      const xb_vecNx16 * X2_rd;
      const xb_vecNx16 * X3_rd;

      xb_vecNx16 x00, x01, x02, x03;
      xb_vecNx16 x10, x11, x12, x13;
      xb_vecNx16 x20, x21, x22, x23;
      xb_vecNx16 x30, x31, x32, x33;

      xb_vecNx16 y0, y1, y2, y3;
      xb_vecNx40 w0, w1, w2, w3;
      xb_vecNx16 z0, z1, z2, z3;

      vsaN vsa0;

      int j;

      vsa0 = BBE_MOVVSA32( Q );

      X0_rd = (const xb_vecNx16*)( (uintptr_t)X_rd + 0*BBE_SIMD_WIDTH/2*N*4 );
      X1_rd = (const xb_vecNx16*)( (uintptr_t)X_rd + 1*BBE_SIMD_WIDTH/2*N*4 );
      X2_rd = (const xb_vecNx16*)( (uintptr_t)X_rd + 2*BBE_SIMD_WIDTH/2*N*4 );
      X3_rd = (const xb_vecNx16*)( (uintptr_t)X_rd + 3*BBE_SIMD_WIDTH/2*N*4 );

      w0 = w1 = w2 = w3 = BBE_MOVWA32( (1UL<<Q) >> 1 );

      __Pragma( "ymemory( Y_rd )" );
      __Pragma( "loop_count min=1" );
      for ( j=0; j<N/4; j++ )
      {
        BBE_LVNX16_IP( x00, X0_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x01, X0_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x02, X0_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x03, X0_rd, +4*BBE_SIMD_WIDTH/2 );
                             
        BBE_LVNX16_IP( x10, X1_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x11, X1_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x12, X1_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x13, X1_rd, +4*BBE_SIMD_WIDTH/2 );
                             
        BBE_LVNX16_IP( x20, X2_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x21, X2_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x22, X2_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x23, X2_rd, +4*BBE_SIMD_WIDTH/2 );
                             
        BBE_LVNX16_IP( x30, X3_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x31, X3_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x32, X3_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( x33, X3_rd, +4*BBE_SIMD_WIDTH/2 );

        BBE_LVNX16_IP( y0, Y_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( y1, Y_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( y2, Y_rd, +4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( y3, Y_rd, +4*BBE_SIMD_WIDTH/2 );

        BBE_MULANX16C( w0, x00, y0 );
        BBE_MULANX16C( w0, x01, y1 );
        BBE_MULANX16C( w0, x02, y2 );
        BBE_MULANX16C( w0, x03, y3 );

        BBE_MULANX16C( w1, x10, y0 );
        BBE_MULANX16C( w1, x11, y1 );
        BBE_MULANX16C( w1, x12, y2 );
        BBE_MULANX16C( w1, x13, y3 );

        BBE_MULANX16C( w2, x20, y0 );
        BBE_MULANX16C( w2, x21, y1 );
        BBE_MULANX16C( w2, x22, y2 );
        BBE_MULANX16C( w2, x23, y3 );

        BBE_MULANX16C( w3, x30, y0 );
        BBE_MULANX16C( w3, x31, y1 );
        BBE_MULANX16C( w3, x32, y2 );
        BBE_MULANX16C( w3, x33, y3 );
      }

      z0 = BBE_PACKVNX40( w0, vsa0 );
      z1 = BBE_PACKVNX40( w1, vsa0 );
      z2 = BBE_PACKVNX40( w2, vsa0 );
      z3 = BBE_PACKVNX40( w3, vsa0 );

      BBE_SVNX16_IP( z0, Z_wr, +4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( z1, Z_wr, +4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( z2, Z_wr, +4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( z3, Z_wr, +4*BBE_SIMD_WIDTH/2 );
    }

    __Pragma( "no_reorder" );

    /*
    * Convert BBE_SIMD_WIDTH/2 resulting vectors to block format: 
    * z[M][BBE_SIMD_WIDTH/2][2] => Z[BBE_SIMD_WIDTH/2][_M][2].
    */

    Z_rd = (const xb_vecNx16*)Z_scr;
    Z_wr = (      xb_vecNx16*)( (uintptr_t)z + l*BBE_SIMD_WIDTH/2*_M*4 );

    {
      xb_vecNx16 a0, a1, a2, a3;
      xb_vecNx16 b0, b1, b2, b3;

      BBE_LVNX16_IP( a0, Z_rd, 4*BBE_SIMD_WIDTH/2 );
      BBE_LVNX16_IP( a1, Z_rd, 4*BBE_SIMD_WIDTH/2 );
      BBE_LVNX16_IP( a2, Z_rd, 4*BBE_SIMD_WIDTH/2 );
      BBE_LVNX16_IP( a3, Z_rd, 4*BBE_SIMD_WIDTH/2 );

      /*
      * 4x8 => 8x4 (complex valued transpose).
      */

      BBE_DSELNX16I( b1, b0, a1, a0, BBE_DSELI_DEINTERLEAVE_2 );
      BBE_DSELNX16I( b3, b2, a3, a2, BBE_DSELI_DEINTERLEAVE_2 );

      BBE_DSELNX16I( a2, a0, b2, b0, BBE_DSELI_DEINTERLEAVE_2 );
      BBE_DSELNX16I( a3, a1, b3, b1, BBE_DSELI_DEINTERLEAVE_2 );

      BBE_DSELNX16I( b2, b0, a1, a0, BBE_DSELI_DEINTERLEAVE_2 );
      BBE_DSELNX16I( b3, b1, a3, a2, BBE_DSELI_DEINTERLEAVE_2 );

      BBE_SVNX16_IP( b0, Z_wr, 4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( b1, Z_wr, 4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( b2, Z_wr, 4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( b3, Z_wr, 4*BBE_SIMD_WIDTH/2 );
    }
  }

} /* cmatvmulnxmn_Mlt8_Ngte8_L8x() */

/* Return the scratch area size, in bytes. */
size_t cmatvmulnxmn_Mlt8_Ngte8_L8x_getScratchSize(int N, int M)
{
  int _N = ((N + (BBE_SIMD_WIDTH / 2 - 1)) & ~(BBE_SIMD_WIDTH / 2 - 1));
  int _M = ((M + (BBE_SIMD_WIDTH / 2 - 1)) & ~(BBE_SIMD_WIDTH / 2 - 1));
  size_t  sz = 0;
  sz = (M*N + _N + _M ) * 8 * 2 * sizeof(int16_t);
  return sz;
} /* cmatvmulnxmn_Mlt8_Ngte8_L8x_getScratchSize() */
#endif
