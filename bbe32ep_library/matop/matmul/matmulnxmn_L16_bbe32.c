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
  NatureDSP_Baseband library. Matrix Operations
    Real Matrix-Matrix/Matrix-Vector Multiply
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#include "matmulnxmn_common.h"

/* get allocated space per one matrix (real) */
static int getSpaceR(int S)
{
  int m;
  /* compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl ) */
  m=30-XT_NSA(S);
  m=XT_MIN(m,LOG2_BBE_SIMD_WIDTH);
  /* round up to the  next multiple of 32 or lesser degree of 2 */
  S=(((S-1)>>m)+1)<<m;
  return S;
} /* getSpaceR() */

/* Compute L matrix products MxN * NxM -> MxM, where L is a multiple of 16. */
void matmulnxmn_L16( void * pScr,
                       int16_t * restrict z, 
                 const int16_t * restrict x, 
                 const int16_t * restrict y, 
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

  int l;
  int MN, MM;

  NASSERT_ALIGN( pScr,2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);

  NASSERT( N%4==0 && M%4==0 );

  NASSERT( Q>=0 && Q<=16 );

  NASSERT( !( L & (BBE_SIMD_WIDTH/2-1) ) );

  MN = M*N; MM = M*M;

  /*
  * Partition the scratch memory
  */

  {
    void * ptr = pScr;

    Z_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)Z_scr + (BBE_SIMD_WIDTH/2)*MM*4 );
    X_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)X_scr + (BBE_SIMD_WIDTH/2)*MN*4 );
    Y_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)Y_scr + (BBE_SIMD_WIDTH/2)*MN*4 );
  }

  /*
  * Main loop computes BBE_SIMD_WIDTH matrix products at a time.
  */

  for ( l=0; l<L/(BBE_SIMD_WIDTH); l++ )
  {
    /*
    * Convert BBE_SIMD_WIDTH MxN matrices x[] and y[] to streaming format.
    */

    X_rd = (const xb_vecNx16*)( (uintptr_t)x + l*(BBE_SIMD_WIDTH/2)*MN*4 );
    X_wr = (      xb_vecNx16*)X_scr;

    {
      xb_vecNx16 a00, a01, a02, a03, a04, a05, a06, a07;
      xb_vecNx16 a10, a11, a12, a13, a14, a15, a16, a17;
      xb_vecNx16 b00, b01, b02, b03, b04, b05, b06, b07;
      xb_vecNx16 b10, b11, b12, b13, b14, b15, b16, b17;

      int n;

      __Pragma("loop_count min=1")
      for ( n=0; n<MN/(BBE_SIMD_WIDTH); n++ )
      {
        BBE_LVNX16_XP( a00, X_rd, MN*2 ); BBE_LVNX16_XP( a01, X_rd, MN*2 );
        BBE_LVNX16_XP( a02, X_rd, MN*2 ); BBE_LVNX16_XP( a03, X_rd, MN*2 );
        BBE_LVNX16_XP( a04, X_rd, MN*2 ); BBE_LVNX16_XP( a05, X_rd, MN*2 );
        BBE_LVNX16_XP( a06, X_rd, MN*2 ); BBE_LVNX16_XP( a07, X_rd, MN*2 );

        BBE_LVNX16_XP( a10, X_rd, MN*2 ); BBE_LVNX16_XP( a11, X_rd, MN*2 );
        BBE_LVNX16_XP( a12, X_rd, MN*2 ); BBE_LVNX16_XP( a13, X_rd, MN*2 );
        BBE_LVNX16_XP( a14, X_rd, MN*2 ); BBE_LVNX16_XP( a15, X_rd, MN*2 );
        BBE_LVNX16_XP( a16, X_rd, MN*2 );

        BBE_LVNX16_XP( a17, X_rd, -15*MN*2 + 4*BBE_SIMD_WIDTH/2 );

        BBE_DSELNX16I( b01, b00, a01, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b03, b02, a03, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b05, b04, a05, a04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b07, b06, a07, a06, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_DSELNX16I( b11, b10, a11, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b13, b12, a13, a12, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b15, b14, a15, a14, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b17, b16, a17, a16, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_DSELNX16I( a02, a00, b02, b00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a03, a01, b03, b01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a06, a04, b06, b04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a07, a05, b07, b05, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( a12, a10, b12, b10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a11, b13, b11, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a14, b16, b14, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a15, b17, b15, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_DSELNX16I( b04, b00, a04, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b05, b01, a05, a01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b06, b02, a06, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b07, b03, a07, a03, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( b14, b10, a14, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b15, b11, a15, a11, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b16, b12, a16, a12, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b17, b13, a17, a13, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( a10, a00, b10, b00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a11, a01, b11, b01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a12, a02, b12, b02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a03, b13, b03, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( a14, a04, b14, b04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a15, a05, b15, b05, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a06, b16, b06, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a07, b17, b07, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_SVNX16_IP( a00, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a01, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a02, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a03, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a04, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a05, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a06, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a07, X_wr, 4*BBE_SIMD_WIDTH/2 );

        BBE_SVNX16_IP( a10, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a11, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a12, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a13, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a14, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a15, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a16, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a17, X_wr, 4*BBE_SIMD_WIDTH/2 );
      }
    }

    __Pragma("no_reorder")

    Y_rd = (const xb_vecNx16*)( (uintptr_t)y + l*(BBE_SIMD_WIDTH/2)*MN*4 );
    Y_wr = (      xb_vecNx16*)Y_scr;

    {
      xb_vecNx16 a00, a01, a02, a03, a04, a05, a06, a07;
      xb_vecNx16 a10, a11, a12, a13, a14, a15, a16, a17;
      xb_vecNx16 b00, b01, b02, b03, b04, b05, b06, b07;
      xb_vecNx16 b10, b11, b12, b13, b14, b15, b16, b17;

      int n;

      __Pragma("loop_count min=1")
      for ( n=0; n<MN/(BBE_SIMD_WIDTH); n++ )
      {
        BBE_LVNX16_XP( a00, Y_rd, MN*2 ); BBE_LVNX16_XP( a01, Y_rd, MN*2 );
        BBE_LVNX16_XP( a02, Y_rd, MN*2 ); BBE_LVNX16_XP( a03, Y_rd, MN*2 );
        BBE_LVNX16_XP( a04, Y_rd, MN*2 ); BBE_LVNX16_XP( a05, Y_rd, MN*2 );
        BBE_LVNX16_XP( a06, Y_rd, MN*2 ); BBE_LVNX16_XP( a07, Y_rd, MN*2 );
                                                                               
        BBE_LVNX16_XP( a10, Y_rd, MN*2 ); BBE_LVNX16_XP( a11, Y_rd, MN*2 );
        BBE_LVNX16_XP( a12, Y_rd, MN*2 ); BBE_LVNX16_XP( a13, Y_rd, MN*2 );
        BBE_LVNX16_XP( a14, Y_rd, MN*2 ); BBE_LVNX16_XP( a15, Y_rd, MN*2 );
        BBE_LVNX16_XP( a16, Y_rd, MN*2 );

        BBE_LVNX16_XP( a17, Y_rd, -15*MN*2 + 4*BBE_SIMD_WIDTH/2 );

        BBE_DSELNX16I( b01, b00, a01, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b03, b02, a03, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b05, b04, a05, a04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b07, b06, a07, a06, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( b11, b10, a11, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b13, b12, a13, a12, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b15, b14, a15, a14, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b17, b16, a17, a16, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( a02, a00, b02, b00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a03, a01, b03, b01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a06, a04, b06, b04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a07, a05, b07, b05, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( a12, a10, b12, b10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a11, b13, b11, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a14, b16, b14, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a15, b17, b15, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( b04, b00, a04, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b05, b01, a05, a01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b06, b02, a06, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b07, b03, a07, a03, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( b14, b10, a14, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b15, b11, a15, a11, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b16, b12, a16, a12, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b17, b13, a17, a13, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( a10, a00, b10, b00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a11, a01, b11, b01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a12, a02, b12, b02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a03, b13, b03, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( a14, a04, b14, b04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a15, a05, b15, b05, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a06, b16, b06, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a07, b17, b07, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_SVNX16_IP( a00, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a01, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a02, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a03, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a04, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a05, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a06, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a07, Y_wr, 4*BBE_SIMD_WIDTH/2 );

        BBE_SVNX16_IP( a10, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a11, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a12, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a13, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a14, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a15, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a16, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( a17, Y_wr, 4*BBE_SIMD_WIDTH/2 );
      }
    }

    /*
    * Compute BBE_SIMD_WIDTH MxN*NxM matrix products in streaming format.
    */

    {
      xb_vecNx16 x0, x1, y0, y1;

      xb_vecNx16 z00, z01, z10, z11;
      xb_vecNx40 w00, w01, w10, w11;

      vsaN vsa0;

      int i, j, k;

      vsa0 = BBE_MOVVSA32( Q );

      Z_wr = Z_scr;

      for ( i=0; i<M/2; i++ )
      {
        X_rd = (const xb_vecNx16*)( (uintptr_t)X_scr + 2*i*N*4*BBE_SIMD_WIDTH/2 );

        for ( j=0; j<M/2; j++ )
        {
          Y_rd = (const xb_vecNx16*)( (uintptr_t)Y_scr + 2*j*4*BBE_SIMD_WIDTH/2 );

          {
            x1 = BBE_LVNX16_X( X_rd, N*4*BBE_SIMD_WIDTH/2 );
            BBE_LVNX16_IP( x0, X_rd,   4*BBE_SIMD_WIDTH/2 );

            BBE_LVNX16_IP( y0, Y_rd,       4*BBE_SIMD_WIDTH/2 );
            BBE_LVNX16_XP( y1, Y_rd, (M-1)*4*BBE_SIMD_WIDTH/2 );

            w00 = BBE_MULRNX16( x0, y0, vsa0 );
            w01 = BBE_MULRNX16( x0, y1, vsa0 );
            w10 = BBE_MULRNX16( x1, y0, vsa0 );
            w11 = BBE_MULRNX16( x1, y1, vsa0 );
          }

          __Pragma( "ymemory( X_rd )" )
          __Pragma( "ymemory( Y_rd )" )
          __Pragma( "loop_count min=3" )
          for ( k=0; k<N-1; k++ )
          {
            x1 = BBE_LVNX16_X( X_rd, N*4*BBE_SIMD_WIDTH/2 );
            BBE_LVNX16_IP( x0, X_rd,   4*BBE_SIMD_WIDTH/2 );

            y1 = BBE_LVNX16_I( Y_rd,   4*BBE_SIMD_WIDTH/2 );
            BBE_LVNX16_XP( y0, Y_rd, M*4*BBE_SIMD_WIDTH/2 );

            BBE_MULANX16( w00, x0, y0 );
            BBE_MULANX16( w01, x0, y1 );
            BBE_MULANX16( w10, x1, y0 );
            BBE_MULANX16( w11, x1, y1 );
          }

          X_rd = (xb_vecNx16*)( (uintptr_t)X_rd - N*4*BBE_SIMD_WIDTH/2 );

          z00 = BBE_PACKVNX40( w00, vsa0 );
          z01 = BBE_PACKVNX40( w01, vsa0 );
          z10 = BBE_PACKVNX40( w10, vsa0 );
          z11 = BBE_PACKVNX40( w11, vsa0 );

          BBE_SVNX16_X ( z10, Z_wr, M*4*BBE_SIMD_WIDTH/2 );
          BBE_SVNX16_IP( z00, Z_wr,   4*BBE_SIMD_WIDTH/2 );
          BBE_SVNX16_X ( z11, Z_wr, M*4*BBE_SIMD_WIDTH/2 );
          BBE_SVNX16_IP( z01, Z_wr,   4*BBE_SIMD_WIDTH/2 );
        }

        Z_wr = (xb_vecNx16*)( (uintptr_t)Z_wr + M*4*BBE_SIMD_WIDTH/2 );
      }
    }

    /*
    * Convert BBE_SIMD_WIDTH MxM matrices z[] to block format.
    */

    Z_rd = (const xb_vecNx16*)Z_scr;
    Z_wr = (      xb_vecNx16*)( (uintptr_t)z + l*(BBE_SIMD_WIDTH/2)*MM*4 );

    {
      xb_vecNx16 a00, a01, a02, a03, a04, a05, a06, a07;
      xb_vecNx16 a10, a11, a12, a13, a14, a15, a16, a17;
      xb_vecNx16 b00, b01, b02, b03, b04, b05, b06, b07;
      xb_vecNx16 b10, b11, b12, b13, b14, b15, b16, b17;

      int n;

      __Pragma("loop_count min=1")
      for ( n=0; n<MM/(BBE_SIMD_WIDTH); n++ )
      {
        BBE_LVNX16_IP( a00, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a01, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a02, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a03, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a04, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a05, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a06, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a07, Z_rd, 4*BBE_SIMD_WIDTH/2 );

        BBE_LVNX16_IP( a10, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a11, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a12, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a13, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a14, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a15, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a16, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a17, Z_rd, 4*BBE_SIMD_WIDTH/2 );

        BBE_DSELNX16I( b01, b00, a01, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b03, b02, a03, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b05, b04, a05, a04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b07, b06, a07, a06, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( b11, b10, a11, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b13, b12, a13, a12, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b15, b14, a15, a14, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b17, b16, a17, a16, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( a02, a00, b02, b00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a03, a01, b03, b01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a06, a04, b06, b04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a07, a05, b07, b05, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( a12, a10, b12, b10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a11, b13, b11, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a14, b16, b14, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a15, b17, b15, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( b04, b00, a04, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b05, b01, a05, a01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b06, b02, a06, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b07, b03, a07, a03, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( b14, b10, a14, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b15, b11, a15, a11, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b16, b12, a16, a12, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( b17, b13, a17, a13, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( a10, a00, b10, b00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a11, a01, b11, b01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a12, a02, b12, b02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a03, b13, b03, BBE_DSELI_DEINTERLEAVE_1 );
                                                                          
        BBE_DSELNX16I( a14, a04, b14, b04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a15, a05, b15, b05, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a06, b16, b06, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a07, b17, b07, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_SVNX16_XP( a00, Z_wr, MM*2 ); BBE_SVNX16_XP( a01, Z_wr, MM*2 );
        BBE_SVNX16_XP( a02, Z_wr, MM*2 ); BBE_SVNX16_XP( a03, Z_wr, MM*2 );
        BBE_SVNX16_XP( a04, Z_wr, MM*2 ); BBE_SVNX16_XP( a05, Z_wr, MM*2 );
        BBE_SVNX16_XP( a06, Z_wr, MM*2 ); BBE_SVNX16_XP( a07, Z_wr, MM*2 );
                                                                                
        BBE_SVNX16_XP( a10, Z_wr, MM*2 ); BBE_SVNX16_XP( a11, Z_wr, MM*2 );
        BBE_SVNX16_XP( a12, Z_wr, MM*2 ); BBE_SVNX16_XP( a13, Z_wr, MM*2 );
        BBE_SVNX16_XP( a14, Z_wr, MM*2 ); BBE_SVNX16_XP( a15, Z_wr, MM*2 );
        BBE_SVNX16_XP( a16, Z_wr, MM*2 );

        BBE_SVNX16_XP( a17, Z_wr, -15*MM*2 + 4*BBE_SIMD_WIDTH/2 );
      }
    }
  }
} /* matmulnxmn_L16() */


/* Return the scratch area size, in bytes. */
size_t matmulnxmn_L16_getScratchSize ( int N, int M )
{
  int Sx, Sz;
  Sx = getSpaceR(N*M);
  Sz = getSpaceR(M*M);
  return (Sx * 16 * 2 + Sz * 16)*sizeof(int16_t);

} /* matmulnxmn_L16_getScratchSize() */
