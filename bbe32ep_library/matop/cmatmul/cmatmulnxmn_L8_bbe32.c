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
    L - multiple of 8
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matop.h"
#include "cmatmulnxmn_common.h"

/* get allocated space per one matrix (complex) */
static int getSpaceC(int S)
{
  int m;
  /* compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl ) */
  m = 30 - XT_NSA(S);
  m = XT_MIN(m, 3);
  /* round up to the  next multiple of 32 or lesser degree of 2 */
  S = 2 * ((((S - 1) >> m) + 1) << m);
  return S;
}

void cmatmulnxmn_L8(void * pScr, complex_fract16 * restrict z, const complex_fract16 * restrict x, const complex_fract16 * restrict y, int N, int M, int L, int Q)
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

  NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
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
  * Main loop computes BBE_SIMD_WIDTH/2 matrix products at a time.
  */

  for ( l=0; l<L/(BBE_SIMD_WIDTH/2); l++ )
  {
    /*
    * Convert BBE_SIMD_WIDTH/2 MxN matrices x[] and y[] to streaming format.
    */

    X_rd = (const xb_vecNx16*)( (uintptr_t)x + l*(BBE_SIMD_WIDTH/2)*MN*4 );
    X_wr = (      xb_vecNx16*)X_scr;

    {
      xb_vecNx16 a0, a1, a2, a3, a4, a5, a6, a7;
      xb_vecNx16 b0, b1, b2, b3, b4, b5, b6, b7;

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

        BBE_DSELNX16I( b1, b0, a1, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b3, b2, a3, a2, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b5, b4, a5, a4, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b7, b6, a7, a6, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_DSELNX16I( a2, a0, b2, b0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a3, a1, b3, b1, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a6, a4, b6, b4, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a7, a5, b7, b5, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_DSELNX16I( b4, b0, a4, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b5, b1, a5, a1, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b6, b2, a6, a2, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b7, b3, a7, a3, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_SVNX16_IP( b0, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b1, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b2, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b3, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b4, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b5, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b6, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b7, X_wr, 4*BBE_SIMD_WIDTH/2 );
      }
    }

    __Pragma( "no_reorder" );

    Y_rd = (const xb_vecNx16*)( (uintptr_t)y + l*(BBE_SIMD_WIDTH/2)*MN*4 );
    Y_wr = (      xb_vecNx16*)Y_scr;

    {
      xb_vecNx16 a0, a1, a2, a3, a4, a5, a6, a7;
      xb_vecNx16 b0, b1, b2, b3, b4, b5, b6, b7;

      int n;

      __Pragma( "loop_count min=2" );
      for ( n=0; n<MN/(BBE_SIMD_WIDTH/2); n++ )
      {
        BBE_LVNX16_XP( a0, Y_rd, MN*4 );
        BBE_LVNX16_XP( a1, Y_rd, MN*4 );
        BBE_LVNX16_XP( a2, Y_rd, MN*4 );
        BBE_LVNX16_XP( a3, Y_rd, MN*4 );
        BBE_LVNX16_XP( a4, Y_rd, MN*4 );
        BBE_LVNX16_XP( a5, Y_rd, MN*4 );
        BBE_LVNX16_XP( a6, Y_rd, MN*4 );

        BBE_LVNX16_XP( a7, Y_rd, -7*MN*4 + 4*BBE_SIMD_WIDTH/2 );

        BBE_DSELNX16I( b1, b0, a1, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b3, b2, a3, a2, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b5, b4, a5, a4, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b7, b6, a7, a6, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_DSELNX16I( a2, a0, b2, b0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a3, a1, b3, b1, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a6, a4, b6, b4, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a7, a5, b7, b5, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_DSELNX16I( b4, b0, a4, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b5, b1, a5, a1, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b6, b2, a6, a2, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b7, b3, a7, a3, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_SVNX16_IP( b0, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b1, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b2, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b3, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b4, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b5, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b6, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b7, Y_wr, 4*BBE_SIMD_WIDTH/2 );
      }
    }

    /*
    * Compute BBE_SIMD_WIDTH/2 MxN*NxM matrix products in streaming format.
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

            w00 = BBE_MULRNX16C( x0, y0, vsa0 );
            w01 = BBE_MULRNX16C( x0, y1, vsa0 );
            w10 = BBE_MULRNX16C( x1, y0, vsa0 );
            w11 = BBE_MULRNX16C( x1, y1, vsa0 );
          }

          __Pragma( "ymemory( X_rd )" );
          __Pragma( "ymemory( Y_rd )" );
          __Pragma( "loop_count min=3" );
          for ( k=0; k<N-1; k++ )
          {
            x1 = BBE_LVNX16_X( X_rd, N*4*BBE_SIMD_WIDTH/2 );
            BBE_LVNX16_IP( x0, X_rd,   4*BBE_SIMD_WIDTH/2 );

            y1 = BBE_LVNX16_I( Y_rd,   4*BBE_SIMD_WIDTH/2 );
            BBE_LVNX16_XP( y0, Y_rd, M*4*BBE_SIMD_WIDTH/2 );

            BBE_MULANX16C( w00, x0, y0 );
            BBE_MULANX16C( w01, x0, y1 );
            BBE_MULANX16C( w10, x1, y0 );
            BBE_MULANX16C( w11, x1, y1 );
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
    * Convert BBE_SIMD_WIDTH/2 MxM matrices z[] to block format.
    */

    Z_rd = (const xb_vecNx16*)Z_scr;
    Z_wr = (      xb_vecNx16*)( (uintptr_t)z + l*(BBE_SIMD_WIDTH/2)*MM*4 );

    {
      xb_vecNx16 a0, a1, a2, a3, a4, a5, a6, a7;
      xb_vecNx16 b0, b1, b2, b3, b4, b5, b6, b7;

      int n;

      __Pragma( "loop_count min=2" );
      for ( n=0; n<MM/(BBE_SIMD_WIDTH/2); n++ )
      {
        BBE_LVNX16_IP( a0, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a1, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a2, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a3, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a4, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a5, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a6, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a7, Z_rd, 4*BBE_SIMD_WIDTH/2 );

        BBE_DSELNX16I( b1, b0, a1, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b3, b2, a3, a2, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b5, b4, a5, a4, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b7, b6, a7, a6, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_DSELNX16I( a2, a0, b2, b0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a3, a1, b3, b1, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a6, a4, b6, b4, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a7, a5, b7, b5, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_DSELNX16I( b4, b0, a4, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b5, b1, a5, a1, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b6, b2, a6, a2, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( b7, b3, a7, a3, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_SVNX16_XP( b0, Z_wr, MM*4 );
        BBE_SVNX16_XP( b1, Z_wr, MM*4 );
        BBE_SVNX16_XP( b2, Z_wr, MM*4 );
        BBE_SVNX16_XP( b3, Z_wr, MM*4 );
        BBE_SVNX16_XP( b4, Z_wr, MM*4 );
        BBE_SVNX16_XP( b5, Z_wr, MM*4 );
        BBE_SVNX16_XP( b6, Z_wr, MM*4 );

        BBE_SVNX16_XP( b7, Z_wr, -7*MM*4 + 4*BBE_SIMD_WIDTH/2 );
      }
    }
  }
}
/* return scratch size for L8 family functions */
size_t cmatmulnxmn_L8_getScratchSize(int N, int M)
{
    int Sx,Sz;
    Sx = getSpaceC(N*M);
    Sz = getSpaceC(M*M);
    return (Sx*8*2+Sz*8)*sizeof(int16_t);
}

