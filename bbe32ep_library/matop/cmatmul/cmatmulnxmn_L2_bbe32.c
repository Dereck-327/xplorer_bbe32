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
    L - multiple of 2
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matop.h"
#include "cmatmulnxmn_common.h"

/* get allocated space per one matrix (complex) */
static int getSpaceC(int S)
{
    int m;
    /* compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl ) */
    m=30-XT_NSA(S);
    m=XT_MIN(m,3);
    /* round up to the  next multiple of 32 or lesser degree of 2 */
    S=2*((((S-1)>>m)+1)<<m);
    return S;
}

void cmatmulnxmn_L2(void * pScr, complex_fract16 * restrict z, const complex_fract16 * restrict x, const complex_fract16 * restrict y, int N, int M, int L, int Q)
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

  int MN, MM;

  int l;

  NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);

  NASSERT( N%4==0 && M%4==0 );

  NASSERT( Q>=0 && Q<=16 );

  NASSERT( !( L & (BBE_SIMD_WIDTH/8-1) ) );

  MN = M*N; MM = M*M;

  /*
  * Partition the scratch memory
  */

  {
    void * ptr = pScr;

    Z_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)Z_scr + (BBE_SIMD_WIDTH/8)*MM*4 );
    X_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)X_scr + (BBE_SIMD_WIDTH/8)*MN*4 );
    Y_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)Y_scr + (BBE_SIMD_WIDTH/8)*MN*4 );

  }

  /*
  * Main loop computes BBE_SIMD_WIDTH/8 matrix products at a time.
  */

  for ( l=0; l<L/(BBE_SIMD_WIDTH/8); l++ )
  {
    /*
    * Convert BBE_SIMD_WIDTH/8 MxN matrices x[] and y[] to streaming format.
    */

    X_rd = (const xb_vecNx16*)( (uintptr_t)x + l*(BBE_SIMD_WIDTH/8)*MN*4 );
    Y_rd = (const xb_vecNx16*)( (uintptr_t)y + l*(BBE_SIMD_WIDTH/8)*MN*4 );

    X_wr = X_scr;
    Y_wr = Y_scr;

    {
      xb_vecNx16 a00, a01, a10, a11;
      xb_vecNx16 b00, b01, b10, b11;

      int n;

      __Pragma( "loop_count min=2" )
      for ( n=0; n<MN/(BBE_SIMD_WIDTH/2); n++ )
      {
        BBE_LVNX16_XP( a00, X_rd, MN*4                        );
        BBE_LVNX16_XP( a01, X_rd, 4*BBE_SIMD_WIDTH/2 - 1*MN*4 );

        BBE_LVNX16_XP( a10, Y_rd, MN*4                        );
        BBE_LVNX16_XP( a11, Y_rd, 4*BBE_SIMD_WIDTH/2 - 1*MN*4 );

        BBE_DSELNX16I( b01, b00, a01, a00, BBE_DSELI_INTERLEAVE_2 );
        BBE_DSELNX16I( b11, b10, a11, a10, BBE_DSELI_INTERLEAVE_2 );

        BBE_SVNX16_IP( b00, X_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b01, X_wr, 4*BBE_SIMD_WIDTH/2 );

        BBE_SVNX16_IP( b10, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( b11, Y_wr, 4*BBE_SIMD_WIDTH/2 );
      }
    }

    /*
    * Compute BBE_SIMD_WIDTH/8 MxN*NxM matrix products in streaming format.
    */

    {
      xb_vecNx16 a0, a1, a2, a3;
      xb_vecNx16 x0, x1, x2, x3;
      xb_vecNx16 z0, z1, z2, z3;

      xb_vecNx16 y00, y01, y02, y03;
      xb_vecNx16 y10, y11, y12, y13;
      xb_vecNx16 y20, y21, y22, y23;
      xb_vecNx16 y30, y31, y32, y33;

      xb_vecNx40 w0, w1, w2, w3;

      vsaN  vs0;

      int i, j, k;

      vs0 = BBE_MOVVSA32( Q );

      Z_wr = Z_scr;

      for ( i=0; i<M/4; i++ )
      {
        X_rd = (const xb_vecNx16*)( (uintptr_t)X_scr + 4*i*N*(BBE_SIMD_WIDTH/8)*4 );
        X_wr = (      xb_vecNx16*)( (uintptr_t)X_scr + 4*i*N*(BBE_SIMD_WIDTH/8)*4 );

        /*
        * Compute 4x4 sub-block (i,0) of resulting matrices z[]. In the course
        * of multiplication, transpose 4x4 sub-blocks (i,0..N/4) of the left-
        * hand matrices and store them back to scratch memory for reuse in the
        * next loop.
        */

        {
          Y_rd = (const xb_vecNx16*)( (uintptr_t)Y_scr + 4*0*(BBE_SIMD_WIDTH/8)*4 );

          w0 = w1 = w2 = w3 = BBE_MOVWA32( (1UL<<Q) >> 1 );

          __Pragma( "ymemory( X_rd )" );
          __Pragma( "ymemory( Y_rd )" );
          __Pragma( "loop_count min=1" );
          for ( k=0; k<N/4; k++ )
          {
            BBE_LVNX16_XP( a0, X_rd, N*(BBE_SIMD_WIDTH/8)*4 );
            BBE_LVNX16_XP( a1, X_rd, N*(BBE_SIMD_WIDTH/8)*4 );
            BBE_LVNX16_XP( a2, X_rd, N*(BBE_SIMD_WIDTH/8)*4 );

            BBE_LVNX16_XP( a3, X_rd, 4*BBE_SIMD_WIDTH/2 - 3*N*(BBE_SIMD_WIDTH/8)*4 );

            BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_INTERLEAVE_4 );
            BBE_DSELNX16I( a3, a2, a3, a2, BBE_DSELI_INTERLEAVE_4 );

            x0 = BBE_SELNX16I( a2, a0, BBE_SELI_EXTRACT_LO_HALVES );
            x1 = BBE_SELNX16I( a2, a0, BBE_SELI_EXTRACT_HI_HALVES );
            x2 = BBE_SELNX16I( a3, a1, BBE_SELI_EXTRACT_LO_HALVES );
            x3 = BBE_SELNX16I( a3, a1, BBE_SELI_EXTRACT_HI_HALVES );

            BBE_SVNX16_XP( x0, X_wr, N*(BBE_SIMD_WIDTH/8)*4 );
            BBE_SVNX16_XP( x1, X_wr, N*(BBE_SIMD_WIDTH/8)*4 );
            BBE_SVNX16_XP( x2, X_wr, N*(BBE_SIMD_WIDTH/8)*4 );

            BBE_SVNX16_XP( x3, X_wr, 4*BBE_SIMD_WIDTH/2 - 3*N*(BBE_SIMD_WIDTH/8)*4 );

            BBE_LVNX16_XP( a0, Y_rd, M*(BBE_SIMD_WIDTH/8)*4 );
            BBE_LVNX16_XP( a1, Y_rd, M*(BBE_SIMD_WIDTH/8)*4 );
            BBE_LVNX16_XP( a2, Y_rd, M*(BBE_SIMD_WIDTH/8)*4 );
            BBE_LVNX16_XP( a3, Y_rd, M*(BBE_SIMD_WIDTH/8)*4 );

            y00 = BBE_SHFLNX16I( a0, BBE_SHFLI_REP_0X4 );
            y01 = BBE_SHFLNX16I( a0, BBE_SHFLI_REP_1X4 );
            y02 = BBE_SHFLNX16I( a0, BBE_SHFLI_REP_2X4 );
            y03 = BBE_SHFLNX16I( a0, BBE_SHFLI_REP_3X4 );

            y10 = BBE_SHFLNX16I( a1, BBE_SHFLI_REP_0X4 );
            y11 = BBE_SHFLNX16I( a1, BBE_SHFLI_REP_1X4 );
            y12 = BBE_SHFLNX16I( a1, BBE_SHFLI_REP_2X4 );
            y13 = BBE_SHFLNX16I( a1, BBE_SHFLI_REP_3X4 );
                                  
            y20 = BBE_SHFLNX16I( a2, BBE_SHFLI_REP_0X4 );
            y21 = BBE_SHFLNX16I( a2, BBE_SHFLI_REP_1X4 );
            y22 = BBE_SHFLNX16I( a2, BBE_SHFLI_REP_2X4 );
            y23 = BBE_SHFLNX16I( a2, BBE_SHFLI_REP_3X4 );
                                  
            y30 = BBE_SHFLNX16I( a3, BBE_SHFLI_REP_0X4 );
            y31 = BBE_SHFLNX16I( a3, BBE_SHFLI_REP_1X4 );
            y32 = BBE_SHFLNX16I( a3, BBE_SHFLI_REP_2X4 );
            y33 = BBE_SHFLNX16I( a3, BBE_SHFLI_REP_3X4 );

            BBE_MULANX16C( w0, x0, y00 ); BBE_MULANX16C( w0, x1, y10 );
            BBE_MULANX16C( w1, x0, y01 ); BBE_MULANX16C( w1, x1, y11 );
            BBE_MULANX16C( w2, x0, y02 ); BBE_MULANX16C( w2, x1, y12 );
            BBE_MULANX16C( w3, x0, y03 ); BBE_MULANX16C( w3, x1, y13 );

            BBE_MULANX16C( w0, x2, y20 ); BBE_MULANX16C( w0, x3, y30 );
            BBE_MULANX16C( w1, x2, y21 ); BBE_MULANX16C( w1, x3, y31 );
            BBE_MULANX16C( w2, x2, y22 ); BBE_MULANX16C( w2, x3, y32 );
            BBE_MULANX16C( w3, x2, y23 ); BBE_MULANX16C( w3, x3, y33 );
          }

          X_rd = (const xb_vecNx16*)( (uintptr_t)X_rd - N*(BBE_SIMD_WIDTH/8)*4 );

          a0 = BBE_PACKVNX40( w0, vs0 );
          a1 = BBE_PACKVNX40( w1, vs0 );
          a2 = BBE_PACKVNX40( w2, vs0 );
          a3 = BBE_PACKVNX40( w3, vs0 );

          /*
          * 4x4 sub-blocks of resulting matrices z[] appear transposed here.
          */

          BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_INTERLEAVE_4 );
          BBE_DSELNX16I( a3, a2, a3, a2, BBE_DSELI_INTERLEAVE_4 );

          z0 = BBE_SELNX16I( a2, a0, BBE_SELI_EXTRACT_LO_HALVES );
          z1 = BBE_SELNX16I( a2, a0, BBE_SELI_EXTRACT_HI_HALVES );
          z2 = BBE_SELNX16I( a3, a1, BBE_SELI_EXTRACT_LO_HALVES );
          z3 = BBE_SELNX16I( a3, a1, BBE_SELI_EXTRACT_HI_HALVES );

          BBE_SVNX16_XP( z0, Z_wr, M*(BBE_SIMD_WIDTH/8)*4 );
          BBE_SVNX16_XP( z1, Z_wr, M*(BBE_SIMD_WIDTH/8)*4 );
          BBE_SVNX16_XP( z2, Z_wr, M*(BBE_SIMD_WIDTH/8)*4 );

          BBE_SVNX16_XP( z3, Z_wr, 4*BBE_SIMD_WIDTH/2 - 3*M*(BBE_SIMD_WIDTH/8)*4 );
        }

        /*
        * Compute 4x4 sub-blocks (i,1..M/4) of resulting matrices z[].
        */

        for ( j=1; j<M/4; j++ )
        {
          Y_rd = (const xb_vecNx16*)( (uintptr_t)Y_scr + 4*j*(BBE_SIMD_WIDTH/8)*4 );

          w0 = w1 = w2 = w3 = BBE_MOVWA32( (1UL<<Q) >> 1 );

          __Pragma( "ymemory( X_rd )" );
          __Pragma( "ymemory( Y_rd )" );
          __Pragma( "loop_count min=1" );
          for ( k=0; k<N/4; k++ )
          {
            BBE_LVNX16_XP( x0, X_rd, N*(BBE_SIMD_WIDTH/8)*4 );
            BBE_LVNX16_XP( x1, X_rd, N*(BBE_SIMD_WIDTH/8)*4 );
            BBE_LVNX16_XP( x2, X_rd, N*(BBE_SIMD_WIDTH/8)*4 );

            BBE_LVNX16_XP( x3, X_rd, 4*BBE_SIMD_WIDTH/2 - 3*N*(BBE_SIMD_WIDTH/8)*4 );

            BBE_LVNX16_XP( a0, Y_rd, M*(BBE_SIMD_WIDTH/8)*4 );
            BBE_LVNX16_XP( a1, Y_rd, M*(BBE_SIMD_WIDTH/8)*4 );
            BBE_LVNX16_XP( a2, Y_rd, M*(BBE_SIMD_WIDTH/8)*4 );
            BBE_LVNX16_XP( a3, Y_rd, M*(BBE_SIMD_WIDTH/8)*4 );

            y00 = BBE_SHFLNX16I( a0, BBE_SHFLI_REP_0X4 );
            y01 = BBE_SHFLNX16I( a0, BBE_SHFLI_REP_1X4 );
            y02 = BBE_SHFLNX16I( a0, BBE_SHFLI_REP_2X4 );
            y03 = BBE_SHFLNX16I( a0, BBE_SHFLI_REP_3X4 );

            y10 = BBE_SHFLNX16I( a1, BBE_SHFLI_REP_0X4 );
            y11 = BBE_SHFLNX16I( a1, BBE_SHFLI_REP_1X4 );
            y12 = BBE_SHFLNX16I( a1, BBE_SHFLI_REP_2X4 );
            y13 = BBE_SHFLNX16I( a1, BBE_SHFLI_REP_3X4 );
                                  
            y20 = BBE_SHFLNX16I( a2, BBE_SHFLI_REP_0X4 );
            y21 = BBE_SHFLNX16I( a2, BBE_SHFLI_REP_1X4 );
            y22 = BBE_SHFLNX16I( a2, BBE_SHFLI_REP_2X4 );
            y23 = BBE_SHFLNX16I( a2, BBE_SHFLI_REP_3X4 );
                                  
            y30 = BBE_SHFLNX16I( a3, BBE_SHFLI_REP_0X4 );
            y31 = BBE_SHFLNX16I( a3, BBE_SHFLI_REP_1X4 );
            y32 = BBE_SHFLNX16I( a3, BBE_SHFLI_REP_2X4 );
            y33 = BBE_SHFLNX16I( a3, BBE_SHFLI_REP_3X4 );

            BBE_MULANX16C( w0, x0, y00 ); BBE_MULANX16C( w0, x1, y10 );
            BBE_MULANX16C( w1, x0, y01 ); BBE_MULANX16C( w1, x1, y11 );
            BBE_MULANX16C( w2, x0, y02 ); BBE_MULANX16C( w2, x1, y12 );
            BBE_MULANX16C( w3, x0, y03 ); BBE_MULANX16C( w3, x1, y13 );

            BBE_MULANX16C( w0, x2, y20 ); BBE_MULANX16C( w0, x3, y30 );
            BBE_MULANX16C( w1, x2, y21 ); BBE_MULANX16C( w1, x3, y31 );
            BBE_MULANX16C( w2, x2, y22 ); BBE_MULANX16C( w2, x3, y32 );
            BBE_MULANX16C( w3, x2, y23 ); BBE_MULANX16C( w3, x3, y33 );
          }

          X_rd = (const xb_vecNx16*)( (uintptr_t)X_rd - N*(BBE_SIMD_WIDTH/8)*4 );

          a0 = BBE_PACKVNX40( w0, vs0 );
          a1 = BBE_PACKVNX40( w1, vs0 );
          a2 = BBE_PACKVNX40( w2, vs0 );
          a3 = BBE_PACKVNX40( w3, vs0 );

          /*
          * 4x4 sub-blocks of resulting matrices z[] appear transposed here.
          */

          BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_INTERLEAVE_4 );
          BBE_DSELNX16I( a3, a2, a3, a2, BBE_DSELI_INTERLEAVE_4 );

          z0 = BBE_SELNX16I( a2, a0, BBE_SELI_EXTRACT_LO_HALVES );
          z1 = BBE_SELNX16I( a2, a0, BBE_SELI_EXTRACT_HI_HALVES );
          z2 = BBE_SELNX16I( a3, a1, BBE_SELI_EXTRACT_LO_HALVES );
          z3 = BBE_SELNX16I( a3, a1, BBE_SELI_EXTRACT_HI_HALVES );

          BBE_SVNX16_XP( z0, Z_wr, M*(BBE_SIMD_WIDTH/8)*4 );
          BBE_SVNX16_XP( z1, Z_wr, M*(BBE_SIMD_WIDTH/8)*4 );
          BBE_SVNX16_XP( z2, Z_wr, M*(BBE_SIMD_WIDTH/8)*4 );

          BBE_SVNX16_XP( z3, Z_wr, 4*BBE_SIMD_WIDTH/2 - 3*M*(BBE_SIMD_WIDTH/8)*4 );
        }

        Z_wr = (xb_vecNx16*)( (uintptr_t)Z_wr + 3*M*(BBE_SIMD_WIDTH/8)*4 );
      }
    }

    /*
    * Convert BBE_SIMD_WIDTH/8 MxM matrices z[] to block format.
    */

    Z_rd = (const xb_vecNx16*)Z_scr;
    Z_wr = (      xb_vecNx16*)( (uintptr_t)z + l*BBE_SIMD_WIDTH/8*MM*4 );

    {
      xb_vecNx16 a0, a1;
      xb_vecNx16 b0, b1;

      int n;

      for ( n=0; n<MM/(BBE_SIMD_WIDTH/2); n++ )
      {
        BBE_LVNX16_IP( a0, Z_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( a1, Z_rd, 4*BBE_SIMD_WIDTH/2 );

        BBE_DSELNX16I( b1, b0, a1, a0, BBE_DSELI_DEINTERLEAVE_2 );

        BBE_SVNX16_XP( b0, Z_wr, MM*4                        );
        BBE_SVNX16_XP( b1, Z_wr, 4*BBE_SIMD_WIDTH/2 - 1*MM*4 );
      }
    }
  }
}
/* return scratch size for L8 family functions */
size_t cmatmulnxmn_L2_getScratchSize(int N, int M)
{
    int Sx,Sz;
    Sx = getSpaceC(N*M);
    Sz = getSpaceC(M*M);
    return (Sx*2*2+Sz*2)*sizeof(int16_t);
}

