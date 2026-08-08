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
    M>=12 && N>=12 && !(L&15)
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
#if !(HAVE_MULPC && 1)
DISCARD_FUN(void, matvmulnxmn_Mgte12_Ngte12_L16x, (void * pScr,
                   int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int N, int M, int L, int Q ))
size_t matvmulnxmn_Mgte12_Ngte12_L16_getScratchSize(int N, int M) { (void)N; (void)M; return 0; }
#else

/* get allocated space per one matrix (real) */
static int getSpaceR(int S)
{
    int m;
    /* compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl ) */
    m=30-XT_NSA(S);
    m=XT_MIN(m,3);
    /* round up to the  next multiple of 32 or lesser degree of 2 */
    S=(((S-1)>>m)+1)<<m;
    return S;
}

// M>=12 && N>=12 && !(L&15)
void matvmulnxmn_Mgte12_Ngte12_L16x( void    *          pScr,
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

  int _M, _N, MN;
  int off0, off1, off2, off3, off4, off5, off6;
  off0 = (((M&15) > 2 ) || (M%16 == 0)) ? 1:0;
  off1 = (((M&15) > 4 ) || (M%16 == 0)) ? 1:0;
  off2 = (((M&15) > 6 ) || (M%16 == 0)) ? 1:0;
  off3 = (((M&15) > 8 ) || (M%16 == 0)) ? 1:0;
  off4 = (((M&15) > 10) || (M%16 == 0)) ? 1:0;
  off5 = (((M&15) > 12) || (M%16 == 0)) ? 1:0;
  off6 = (((M&15) > 14) || (M%16 == 0)) ? 1:0;

  NASSERT_ALIGN( pScr, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN(z, 2*BBE_SIMD_WIDTH );

  NASSERT( (L>0) && !(L&3) && !(M&3) && !(N&3) );

  NASSERT( M>=12 && N>=12 && !(L&15) );

  NASSERT( Q>=0 && Q<=16 );

  MN = M*N;

  // Round up to the next multiple of BBE_SIMD_WIDTH.
  _M = ( ( M + (BBE_SIMD_WIDTH-1) ) & ~(BBE_SIMD_WIDTH-1) );
  _N = ( ( N + (BBE_SIMD_WIDTH-1) ) & ~(BBE_SIMD_WIDTH-1) );

  //
  // Partition the scratch memory area.
  //

  {
    void * ptr = pScr;

    X_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)X_scr + BBE_SIMD_WIDTH*MN*2 );
    Y_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)Y_scr + BBE_SIMD_WIDTH*_N*2 );
    Z_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)Z_scr + BBE_SIMD_WIDTH*_M*2 );
  }

  //
  // Compute BBE_SIMD_WIDTH matrix-vector products at a time.
  //

  for ( l=0; l<L/BBE_SIMD_WIDTH; l++ )
  {
    //--------------------------------------------------------------------------
    // Convert BBE_SIMD_WIDTH left-hand matrices x[M*N] to streaming format:
    // x[BBE_SIMD_WIDTH][M*N] => x[M*N/4]][BBE_SIMD_WIDTH][4].
    //

    X_rd = (const xb_vecNx16*)( (uintptr_t)x + l*BBE_SIMD_WIDTH*MN*2 );
    X_wr = (      xb_vecNx16*)X_scr;

    {
      int n;

      __Pragma( "loop_count min=1" );
      for ( n=0; n<MN/BBE_SIMD_WIDTH; n++ )
      {
        xb_vecNx16 a00, a01, a02, a03, a10, a11, a12, a13;
        xb_vecNx16 a20, a21, a22, a23, a30, a31, a32, a33;
        xb_vecNx16 b00, b01, b02, b03, b10, b11, b12, b13;
        xb_vecNx16 b20, b21, b22, b23, b30, b31, b32, b33;

        BBE_LVNX16_XP( a00, X_rd, MN*2 );
        BBE_LVNX16_XP( a01, X_rd, MN*2 );
        BBE_LVNX16_XP( a02, X_rd, MN*2 );
        BBE_LVNX16_XP( a03, X_rd, MN*2 );

        BBE_LVNX16_XP( a10, X_rd, MN*2 );
        BBE_LVNX16_XP( a11, X_rd, MN*2 );
        BBE_LVNX16_XP( a12, X_rd, MN*2 );
        BBE_LVNX16_XP( a13, X_rd, MN*2 );

        BBE_LVNX16_XP( a20, X_rd, MN*2 );
        BBE_LVNX16_XP( a21, X_rd, MN*2 );
        BBE_LVNX16_XP( a22, X_rd, MN*2 );
        BBE_LVNX16_XP( a23, X_rd, MN*2 );

        BBE_LVNX16_XP( a30, X_rd, MN*2 );
        BBE_LVNX16_XP( a31, X_rd, MN*2 );
        BBE_LVNX16_XP( a32, X_rd, MN*2 );

        BBE_LVNX16_XP( a33, X_rd, -15*MN*2 + 2*BBE_SIMD_WIDTH );

        //
        // Real-valued ( L64_4 x I4 ):
        //   L64_4 x I4 = ( L16_4 x I16 )*( I4 x L16_8*L16_8 x I4 )
        //

        BBE_DSELNX16I( b01, b00, a02, a00, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( b03, b02, a03, a01, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( a01, a00, b02, b00, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( a03, a02, b03, b01, BBE_DSELI_INTERLEAVE_4 );

        BBE_DSELNX16I( b11, b10, a12, a10, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( b13, b12, a13, a11, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( a11, a10, b12, b10, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( a13, a12, b13, b11, BBE_DSELI_INTERLEAVE_4 );

        BBE_DSELNX16I( b21, b20, a22, a20, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( b23, b22, a23, a21, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( a21, a20, b22, b20, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( a23, a22, b23, b21, BBE_DSELI_INTERLEAVE_4 );

        BBE_DSELNX16I( b31, b30, a32, a30, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( b33, b32, a33, a31, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( a31, a30, b32, b30, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( a33, a32, b33, b31, BBE_DSELI_INTERLEAVE_4 );

        b00 = a00; b01 = a10; b02 = a20; b03 = a30;
        b10 = a01; b11 = a11; b12 = a21; b13 = a31;
        b20 = a02; b21 = a12; b22 = a22; b23 = a32;
        b30 = a03; b31 = a13; b32 = a23; b33 = a33;

        BBE_SVNX16_I( b01, X_wr, 1*2*BBE_SIMD_WIDTH );
        BBE_SVNX16_I( b02, X_wr, 2*2*BBE_SIMD_WIDTH );
        BBE_SVNX16_I( b03, X_wr, 3*2*BBE_SIMD_WIDTH );

        BBE_SVNX16_I( b10, X_wr, 4*2*BBE_SIMD_WIDTH );
        BBE_SVNX16_I( b11, X_wr, 5*2*BBE_SIMD_WIDTH );
        BBE_SVNX16_I( b12, X_wr, 6*2*BBE_SIMD_WIDTH );
        BBE_SVNX16_I( b13, X_wr, 7*2*BBE_SIMD_WIDTH );

        BBE_SVNX16_I( b20, X_wr,  8*2*BBE_SIMD_WIDTH );
        BBE_SVNX16_I( b21, X_wr,  9*2*BBE_SIMD_WIDTH );
        BBE_SVNX16_I( b22, X_wr, 10*2*BBE_SIMD_WIDTH );
        BBE_SVNX16_I( b23, X_wr, 11*2*BBE_SIMD_WIDTH );

        BBE_SVNX16_I( b30, X_wr, 12*2*BBE_SIMD_WIDTH );
        BBE_SVNX16_I( b31, X_wr, 13*2*BBE_SIMD_WIDTH );
        BBE_SVNX16_I( b32, X_wr, 14*2*BBE_SIMD_WIDTH );
        BBE_SVNX16_I( b33, X_wr, 15*2*BBE_SIMD_WIDTH );

        BBE_SVNX16_XP( b00, X_wr, 16*2*BBE_SIMD_WIDTH );
      }
    }

    //--------------------------------------------------------------------------
    // Convert BBE_SIMD_WIDTH right-hand vectors y[N] to streaming format:
    // y[BBE_SIMD_WIDTH][_N] => y[N][BBE_SIMD_WIDTH].
    //

    Y_rd = (const xb_vecNx16*)( (uintptr_t)y + l*BBE_SIMD_WIDTH*_N*2 );
    Y_wr = (      xb_vecNx16*)Y_scr;

    {
      int n;

      __Pragma( "loop_count min=1" );
      for ( n=0; n<_N/BBE_SIMD_WIDTH; n++ )
      {
        xb_vecNx16 a00, a01, a02, a03, a04, a05, a06, a07;
        xb_vecNx16 a10, a11, a12, a13, a14, a15, a16, a17;

        BBE_LVNX16_XP( a00, Y_rd, _N*2 );
        BBE_LVNX16_XP( a01, Y_rd, _N*2 );
        BBE_LVNX16_XP( a02, Y_rd, _N*2 );
        BBE_LVNX16_XP( a03, Y_rd, _N*2 );
        BBE_LVNX16_XP( a04, Y_rd, _N*2 );
        BBE_LVNX16_XP( a05, Y_rd, _N*2 );
        BBE_LVNX16_XP( a06, Y_rd, _N*2 );
        BBE_LVNX16_XP( a07, Y_rd, _N*2 );
                                  
        BBE_LVNX16_XP( a10, Y_rd, _N*2 );
        BBE_LVNX16_XP( a11, Y_rd, _N*2 );
        BBE_LVNX16_XP( a12, Y_rd, _N*2 );
        BBE_LVNX16_XP( a13, Y_rd, _N*2 );
        BBE_LVNX16_XP( a14, Y_rd, _N*2 );
        BBE_LVNX16_XP( a15, Y_rd, _N*2 );
        BBE_LVNX16_XP( a16, Y_rd, _N*2 );

        BBE_LVNX16_XP( a17, Y_rd, -15*_N*2 + 2*BBE_SIMD_WIDTH );

        //
        // Real-valued 16x16 matrix transpose.
        //

        BBE_DSELNX16I( a01, a00, a01, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a03, a02, a03, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a05, a04, a05, a04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a07, a06, a07, a06, BBE_DSELI_DEINTERLEAVE_1 );
                                       
        BBE_DSELNX16I( a11, a10, a11, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a12, a13, a12, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a15, a14, a15, a14, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a16, a17, a16, BBE_DSELI_DEINTERLEAVE_1 );
                                        
        BBE_DSELNX16I( a02, a00, a02, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a06, a04, a06, a04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a03, a01, a03, a01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a07, a05, a07, a05, BBE_DSELI_DEINTERLEAVE_1 );
                                       
        BBE_DSELNX16I( a12, a10, a12, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a14, a16, a14, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a11, a13, a11, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a15, a17, a15, BBE_DSELI_DEINTERLEAVE_1 );
                                        
        BBE_DSELNX16I( a04, a00, a04, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a05, a01, a05, a01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a06, a02, a06, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a07, a03, a07, a03, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_DSELNX16I( a14, a10, a14, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a15, a11, a15, a11, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a12, a16, a12, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a13, a17, a13, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_DSELNX16I( a10, a00, a10, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a11, a01, a11, a01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a12, a02, a12, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a03, a13, a03, BBE_DSELI_DEINTERLEAVE_1 );
                         
        BBE_DSELNX16I( a14, a04, a14, a04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a15, a05, a15, a05, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a06, a16, a06, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a07, a17, a07, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_SVNX16_IP( a00, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a01, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a02, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a03, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a04, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a05, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a06, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a07, Y_wr, +2*BBE_SIMD_WIDTH );

        BBE_SVNX16_IP( a10, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a11, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a12, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a13, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a14, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a15, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a16, Y_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a17, Y_wr, +2*BBE_SIMD_WIDTH );
      }
    }

    __Pragma( "no_reorder" );

    //--------------------------------------------------------------------------
    // Compute BBE_SIMD_WIDTH matrix-vector products; store resulting vectors
    // to the scratch in streaming format: z[M][BBE_SIMD_WIDTH].
    //

    X_rd = X_scr;
    Y_rd = Y_scr;
    Z_wr = Z_scr;

    {
      xb_vecNx16 x00, x01, x02, x03;
      xb_vecNx16 x10, x11, x12, x13;
      xb_vecNx16 x20, x21, x22, x23;
      xb_vecNx16 x30, x31, x32, x33;

      xb_vecNx16 y0, y1, y2, y3;
      xb_vecNx40 w0, w1, w2, w3;
      xb_vecNx16 z0, z1, z2, z3;

      vsaN vsa0;

      int i, j;

      vsa0 = BBE_MOVVSA32( Q );

      for ( i=0; i<M/4; i++ )
      {
        w0 = w1 = w2 = w3 = 0;

        __Pragma( "loop_count min=1" );
        for ( j=0; j<N/4; j++ )
        {
          BBE_LVNX16_IP( x00, X_rd, +1*2*BBE_SIMD_WIDTH                       );
          BBE_LVNX16_IP( x01, X_rd, +1*2*BBE_SIMD_WIDTH                       );
          BBE_LVNX16_IP( x02, X_rd, +1*2*BBE_SIMD_WIDTH                       );
          BBE_LVNX16_XP( x03, X_rd, -3*2*BBE_SIMD_WIDTH + 1*N*BBE_SIMD_WIDTH*2);

          BBE_LVNX16_IP( x10, X_rd, +1*2*BBE_SIMD_WIDTH                       );
          BBE_LVNX16_IP( x11, X_rd, +1*2*BBE_SIMD_WIDTH                       );
          BBE_LVNX16_IP( x12, X_rd, +1*2*BBE_SIMD_WIDTH                       );
          BBE_LVNX16_XP( x13, X_rd, -3*2*BBE_SIMD_WIDTH + 1*N*BBE_SIMD_WIDTH*2);
                                     
          BBE_LVNX16_IP( x20, X_rd, +1*2*BBE_SIMD_WIDTH                       );
          BBE_LVNX16_IP( x21, X_rd, +1*2*BBE_SIMD_WIDTH                       );
          BBE_LVNX16_IP( x22, X_rd, +1*2*BBE_SIMD_WIDTH                       );
          BBE_LVNX16_XP( x23, X_rd, -3*2*BBE_SIMD_WIDTH + 1*N*BBE_SIMD_WIDTH*2);
                                    
          BBE_LVNX16_IP( x30, X_rd, +1*2*BBE_SIMD_WIDTH                       );
          BBE_LVNX16_IP( x31, X_rd, +1*2*BBE_SIMD_WIDTH                       );
          BBE_LVNX16_IP( x32, X_rd, +1*2*BBE_SIMD_WIDTH                       );
          BBE_LVNX16_XP( x33, X_rd, +1*2*BBE_SIMD_WIDTH - 3*N*BBE_SIMD_WIDTH*2);

          BBE_LVNX16_IP( y0, Y_rd, +2*BBE_SIMD_WIDTH );
          BBE_LVNX16_IP( y1, Y_rd, +2*BBE_SIMD_WIDTH );
          BBE_LVNX16_IP( y2, Y_rd, +2*BBE_SIMD_WIDTH );
          BBE_LVNX16_IP( y3, Y_rd, +2*BBE_SIMD_WIDTH );

          //
          // Real-valued I4 x L64_4
          //

          BBE_DSELNX16I( x01, x00, x01, x00, BBE_DSELI_DEINTERLEAVE_1 );
          BBE_DSELNX16I( x03, x02, x03, x02, BBE_DSELI_DEINTERLEAVE_1 );
          BBE_DSELNX16I( x02, x00, x02, x00, BBE_DSELI_DEINTERLEAVE_1 );
          BBE_DSELNX16I( x03, x01, x03, x01, BBE_DSELI_DEINTERLEAVE_1 );

          BBE_DSELNX16I( x11, x10, x11, x10, BBE_DSELI_DEINTERLEAVE_1 );
          BBE_DSELNX16I( x13, x12, x13, x12, BBE_DSELI_DEINTERLEAVE_1 );
          BBE_DSELNX16I( x12, x10, x12, x10, BBE_DSELI_DEINTERLEAVE_1 );
          BBE_DSELNX16I( x13, x11, x13, x11, BBE_DSELI_DEINTERLEAVE_1 );

          BBE_DSELNX16I( x21, x20, x21, x20, BBE_DSELI_DEINTERLEAVE_1 );
          BBE_DSELNX16I( x23, x22, x23, x22, BBE_DSELI_DEINTERLEAVE_1 );
          BBE_DSELNX16I( x22, x20, x22, x20, BBE_DSELI_DEINTERLEAVE_1 );
          BBE_DSELNX16I( x23, x21, x23, x21, BBE_DSELI_DEINTERLEAVE_1 );

          BBE_DSELNX16I( x31, x30, x31, x30, BBE_DSELI_DEINTERLEAVE_1 );
          BBE_DSELNX16I( x33, x32, x33, x32, BBE_DSELI_DEINTERLEAVE_1 );
          BBE_DSELNX16I( x32, x30, x32, x30, BBE_DSELI_DEINTERLEAVE_1 );
          BBE_DSELNX16I( x33, x31, x33, x31, BBE_DSELI_DEINTERLEAVE_1 );

          BBE_MULANX16( w0, x00, y0 );
          BBE_MULANX16( w0, x01, y1 );
          BBE_MULANX16( w0, x02, y2 );
          BBE_MULANX16( w0, x03, y3 );

          BBE_MULANX16( w1, x10, y0 );
          BBE_MULANX16( w1, x11, y1 );
          BBE_MULANX16( w1, x12, y2 );
          BBE_MULANX16( w1, x13, y3 );

          BBE_MULANX16( w2, x20, y0 );
          BBE_MULANX16( w2, x21, y1 );
          BBE_MULANX16( w2, x22, y2 );
          BBE_MULANX16( w2, x23, y3 );

          BBE_MULANX16( w3, x30, y0 );
          BBE_MULANX16( w3, x31, y1 );
          BBE_MULANX16( w3, x32, y2 );
          BBE_MULANX16( w3, x33, y3 );
        }

        X_rd = (const xb_vecNx16*)( (uintptr_t)X_rd + 3*BBE_SIMD_WIDTH*N*2 );
        Y_rd = (const xb_vecNx16*)( (uintptr_t)Y_rd - 1*BBE_SIMD_WIDTH*N*2 );

        w0 = BBE_RNDADJNX40( w0, vsa0 );
        w1 = BBE_RNDADJNX40( w1, vsa0 );
        w2 = BBE_RNDADJNX40( w2, vsa0 );
        w3 = BBE_RNDADJNX40( w3, vsa0 );

        z0 = BBE_PACKVNX40( w0, vsa0 );
        z1 = BBE_PACKVNX40( w1, vsa0 );
        z2 = BBE_PACKVNX40( w2, vsa0 );
        z3 = BBE_PACKVNX40( w3, vsa0 );

        BBE_SVNX16_IP( z0, Z_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( z1, Z_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( z2, Z_wr, +2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( z3, Z_wr, +2*BBE_SIMD_WIDTH );
      }
    }

    __Pragma( "no_reorder" );

    //--------------------------------------------------------------------------
    // Convert BBE_SIMD_WIDTH resulting vectors to block format: 
    // z[M][BBE_SIMD_WIDTH] => Z[BBE_SIMD_WIDTH][_M].
    //

    Z_rd = (const xb_vecNx16*)Z_scr;
    Z_wr = (      xb_vecNx16*)( (uintptr_t)z + l*BBE_SIMD_WIDTH*_M*2 );

    {
      int m;

      //__Pragma( "loop_count min=1" );
      for ( m=1; m<_M/BBE_SIMD_WIDTH; m++ )
      {
        xb_vecNx16 a00, a01, a02, a03, a04, a05, a06, a07;
        xb_vecNx16 a10, a11, a12, a13, a14, a15, a16, a17;

        BBE_LVNX16_IP( a00, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a01, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a02, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a03, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a04, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a05, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a06, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a07, Z_rd, +2*BBE_SIMD_WIDTH );

        BBE_LVNX16_IP( a10, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a11, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a12, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a13, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a14, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a15, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a16, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a17, Z_rd, +2*BBE_SIMD_WIDTH );

        //
        // Real-valued 16x16 matrix transpose.
        //

        BBE_DSELNX16I( a01, a00, a01, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a03, a02, a03, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a05, a04, a05, a04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a07, a06, a07, a06, BBE_DSELI_DEINTERLEAVE_1 );
                                       
        BBE_DSELNX16I( a11, a10, a11, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a12, a13, a12, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a15, a14, a15, a14, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a16, a17, a16, BBE_DSELI_DEINTERLEAVE_1 );
                                        
        BBE_DSELNX16I( a02, a00, a02, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a06, a04, a06, a04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a03, a01, a03, a01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a07, a05, a07, a05, BBE_DSELI_DEINTERLEAVE_1 );
                                       
        BBE_DSELNX16I( a12, a10, a12, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a14, a16, a14, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a11, a13, a11, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a15, a17, a15, BBE_DSELI_DEINTERLEAVE_1 );
                                        
        BBE_DSELNX16I( a04, a00, a04, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a05, a01, a05, a01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a06, a02, a06, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a07, a03, a07, a03, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_DSELNX16I( a14, a10, a14, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a15, a11, a15, a11, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a12, a16, a12, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a13, a17, a13, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_DSELNX16I( a10, a00, a10, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a11, a01, a11, a01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a12, a02, a12, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a03, a13, a03, BBE_DSELI_DEINTERLEAVE_1 );
                         
        BBE_DSELNX16I( a14, a04, a14, a04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a15, a05, a15, a05, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a06, a16, a06, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a07, a17, a07, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_SVNX16_XP( a00, Z_wr, _M*2 );
        BBE_SVNX16_XP( a01, Z_wr, _M*2 );
        BBE_SVNX16_XP( a02, Z_wr, _M*2 );
        BBE_SVNX16_XP( a03, Z_wr, _M*2 );
        BBE_SVNX16_XP( a04, Z_wr, _M*2 );
        BBE_SVNX16_XP( a05, Z_wr, _M*2 );
        BBE_SVNX16_XP( a06, Z_wr, _M*2 );
        BBE_SVNX16_XP( a07, Z_wr, _M*2 );

        BBE_SVNX16_XP( a10, Z_wr, _M*2 );
        BBE_SVNX16_XP( a11, Z_wr, _M*2 );
        BBE_SVNX16_XP( a12, Z_wr, _M*2 );
        BBE_SVNX16_XP( a13, Z_wr, _M*2 );
        BBE_SVNX16_XP( a14, Z_wr, _M*2 );
        BBE_SVNX16_XP( a15, Z_wr, _M*2 );
        BBE_SVNX16_XP( a16, Z_wr, _M*2 );

        BBE_SVNX16_XP( a17, Z_wr, -15*_M*2 + 2*BBE_SIMD_WIDTH );
      }
      {
        xb_vecNx16 a00, a01, a02, a03, a04, a05, a06, a07;
        xb_vecNx16 a10, a11, a12, a13, a14, a15, a16, a17;
        a02 = BBE_LVNX16_X(Z_rd, off0*4*BBE_SIMD_WIDTH );
        a03 = BBE_LVNX16_X(Z_rd, off0*6*BBE_SIMD_WIDTH );
        a04 = BBE_LVNX16_X(Z_rd, off1*8*BBE_SIMD_WIDTH );
        a05 = BBE_LVNX16_X(Z_rd, off1*10*BBE_SIMD_WIDTH );
        a06 = BBE_LVNX16_X(Z_rd, off2*12*BBE_SIMD_WIDTH);
        a07 = BBE_LVNX16_X(Z_rd, off2*14*BBE_SIMD_WIDTH);
        a10 = BBE_LVNX16_X(Z_rd, off3*16*BBE_SIMD_WIDTH);
        a11 = BBE_LVNX16_X(Z_rd, off3*18*BBE_SIMD_WIDTH);
        a12 = BBE_LVNX16_X(Z_rd, off4*20*BBE_SIMD_WIDTH);
        a13 = BBE_LVNX16_X(Z_rd, off4*22*BBE_SIMD_WIDTH);
        a14 = BBE_LVNX16_X(Z_rd, off5*24*BBE_SIMD_WIDTH);
        a15 = BBE_LVNX16_X(Z_rd, off5*26*BBE_SIMD_WIDTH);
        a16 = BBE_LVNX16_X(Z_rd, off6*28*BBE_SIMD_WIDTH);
        a17 = BBE_LVNX16_X(Z_rd, off6*30*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP( a00, Z_rd, +2*BBE_SIMD_WIDTH );
        BBE_LVNX16_XP( a01, Z_rd, +30*BBE_SIMD_WIDTH );

        //
        // Real-valued 16x16 matrix transpose.
        //

        BBE_DSELNX16I( a01, a00, a01, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a03, a02, a03, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a05, a04, a05, a04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a07, a06, a07, a06, BBE_DSELI_DEINTERLEAVE_1 );
                                       
        BBE_DSELNX16I( a11, a10, a11, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a12, a13, a12, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a15, a14, a15, a14, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a16, a17, a16, BBE_DSELI_DEINTERLEAVE_1 );
                                        
        BBE_DSELNX16I( a02, a00, a02, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a06, a04, a06, a04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a03, a01, a03, a01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a07, a05, a07, a05, BBE_DSELI_DEINTERLEAVE_1 );
                                       
        BBE_DSELNX16I( a12, a10, a12, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a14, a16, a14, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a11, a13, a11, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a15, a17, a15, BBE_DSELI_DEINTERLEAVE_1 );
                                        
        BBE_DSELNX16I( a04, a00, a04, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a05, a01, a05, a01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a06, a02, a06, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a07, a03, a07, a03, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_DSELNX16I( a14, a10, a14, a10, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a15, a11, a15, a11, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a12, a16, a12, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a13, a17, a13, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_DSELNX16I( a10, a00, a10, a00, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a11, a01, a11, a01, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a12, a02, a12, a02, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a13, a03, a13, a03, BBE_DSELI_DEINTERLEAVE_1 );
                         
        BBE_DSELNX16I( a14, a04, a14, a04, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a15, a05, a15, a05, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a16, a06, a16, a06, BBE_DSELI_DEINTERLEAVE_1 );
        BBE_DSELNX16I( a17, a07, a17, a07, BBE_DSELI_DEINTERLEAVE_1 );

        BBE_SVNX16_XP( a00, Z_wr, _M*2 );
        BBE_SVNX16_XP( a01, Z_wr, _M*2 );
        BBE_SVNX16_XP( a02, Z_wr, _M*2 );
        BBE_SVNX16_XP( a03, Z_wr, _M*2 );
        BBE_SVNX16_XP( a04, Z_wr, _M*2 );
        BBE_SVNX16_XP( a05, Z_wr, _M*2 );
        BBE_SVNX16_XP( a06, Z_wr, _M*2 );
        BBE_SVNX16_XP( a07, Z_wr, _M*2 );

        BBE_SVNX16_XP( a10, Z_wr, _M*2 );
        BBE_SVNX16_XP( a11, Z_wr, _M*2 );
        BBE_SVNX16_XP( a12, Z_wr, _M*2 );
        BBE_SVNX16_XP( a13, Z_wr, _M*2 );
        BBE_SVNX16_XP( a14, Z_wr, _M*2 );
        BBE_SVNX16_XP( a15, Z_wr, _M*2 );
        BBE_SVNX16_XP( a16, Z_wr, _M*2 );

        BBE_SVNX16_XP( a17, Z_wr, -15*_M*2 + 2*BBE_SIMD_WIDTH );
      }
    }
  }

} // matvmulnxmn_Mgte12_Ngte12_L16x()
size_t matvmulnxmn_Mgte12_Ngte12_L16_getScratchSize(int N, int M)
{
  int Sx, /*_N,*/ _M;
  _M = ((M + (BBE_SIMD_WIDTH - 1)) & ~(BBE_SIMD_WIDTH - 1));
  Sx = getSpaceR(N*M);
  return (Sx * 16 * 2 + _M * 16)*sizeof(int16_t);
} /* matvmulnxmn_Mgte12_Ngte12_L16_getScratchSize() */
#endif
