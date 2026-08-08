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
    M>=12 && N==4 && !(L&3)
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
DISCARD_FUN(void, matvmulnxmn_Mgte12_N4_L4x, (void * pScr,
                   int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int N, int M, int L, int Q ))
size_t matvmulnxmn_Mgte12_N4_getScratchSize(int N, int M) { (void)N; (void)M; return 0; }
#else

/* get allocated space per one matrix (real) */
static int getSpaceR(int S)
{
  int m;
  /* compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl ) */
  m = 30 - XT_NSA(S);
  m = XT_MIN(m, 3);
  /* round up to the  next multiple of 32 or lesser degree of 2 */
  S = (((S - 1) >> m) + 1) << m;
  return S;
}

// M>=12 && N==4 && !(L&3)
void matvmulnxmn_Mgte12_N4_L4x( void    *          pScr,
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
  const xb_vecNx16 *          Y_rd;

  int _M, MN;

  int l;
  int off0, off1, off2;
  off0 = (((M&15) > 4)  || (M%16 == 0)) ? 1:0;
  off1 = (((M&15) > 8)  || (M%16 == 0)) ? 1:0;
  off2 = (((M%16) == 0 ) || (M%16 == 0)) ? 1:0;

  NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);

  NASSERT( !(L&3) && !(M&3) && !(N&3) );

  NASSERT( M>=12 && N==4 );

  NASSERT( Q>=0 && Q<=16 );

  MN = M*4;

  // Round up to the next multiple of BBE_SIMD_WIDTH.
  _M = ( ( M + (BBE_SIMD_WIDTH-1) ) & ~(BBE_SIMD_WIDTH-1) );

  //
  // Partition the scratch memory area.
  //

  {
    void * ptr = pScr;

    X_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)X_scr + BBE_SIMD_WIDTH/4*MN*2 );
    Z_scr = (xb_vecNx16*)ptr;
    ptr   = (void*)( (uintptr_t)Z_scr + BBE_SIMD_WIDTH/4*_M*2 );

  }

  //
  // Compute BBE_SIMD_WIDTH/4 matrix-vector products at a time.
  //

  Y_rd = (xb_vecNx16*)y;

  for ( l=0; l<L/(BBE_SIMD_WIDTH/4); l++ )
  {
    //--------------------------------------------------------------------------
    // Convert BBE_SIMD_WIDTH/4 left-hand matrices x[M*N] to streaming format:
    // x[BBE_SIMD_WIDTH/4][M*N] => x[M*N/2][BBE_SIMD_WIDTH/4][2].
    //

    X_rd = (const xb_vecNx16*)( (uintptr_t)x + l*(BBE_SIMD_WIDTH/4)*MN*2 );
    X_wr = X_scr;

    {
      int n;

      __Pragma( "loop_count min=1" );
      for ( n=0; n<MN/BBE_SIMD_WIDTH; n++ )
      {
        xb_vecNx16 a0, a1, a2, a3;

        BBE_LVNX16_XP( a0, X_rd, MN*2 );
        BBE_LVNX16_XP( a1, X_rd, MN*2 );
        BBE_LVNX16_XP( a2, X_rd, MN*2 );

        BBE_LVNX16_XP( a3, X_rd, -3*MN*2 + 2*BBE_SIMD_WIDTH );

        BBE_DSELNX16I( a2, a0, a2, a0, BBE_DSELI_INTERLEAVE_2 );
        BBE_DSELNX16I( a3, a1, a3, a1, BBE_DSELI_INTERLEAVE_2 );

        BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_INTERLEAVE_2 );
        BBE_DSELNX16I( a3, a2, a3, a2, BBE_DSELI_INTERLEAVE_2 );

        BBE_SVNX16_IP( a0, X_wr, 2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a1, X_wr, 2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a2, X_wr, 2*BBE_SIMD_WIDTH );
        BBE_SVNX16_IP( a3, X_wr, 2*BBE_SIMD_WIDTH );
      }
    }

    __Pragma( "no_reorder" );

    //--------------------------------------------------------------------------
    // Compute BBE_SIMD_WIDTH/4 matrix-vector products; store resulting vectors
    // to the scratch in streaming format: z[M/4][BBE_SIMD_WIDTH/4][4].
    //

    X_rd = X_scr;
    Z_wr = Z_scr;

    {
      xb_vecNx16 x0, x1, x2, x3;
      xb_vecNx16 y0, y1;
      xb_vecNx16 a0, a1;
      xb_vecNx16 z0;
      xb_vecNx40 w0;

      vsaN  vsa0;
      vselN sel0, sel1;

      int i;

      static const int16_t ALIGN(32) cst[][BBE_SIMD_WIDTH] = {
        { 0,16,1,16,4,16,5,16, 8,16, 9,16,12,16,13,16 },
        { 2,16,3,16,6,16,7,16,10,16,11,16,14,16,15,16 }
      };

      vsa0 = BBE_MOVVSA32( Q );

      {
        y0 = BBE_LVNX16_I( (const xb_vecNx16*)cst, 0*2*BBE_SIMD_WIDTH );
        y1 = BBE_LVNX16_I( (const xb_vecNx16*)cst, 1*2*BBE_SIMD_WIDTH );

        sel0 = BBE_MOVVSELNX16( y0, 0 );
        sel1 = BBE_MOVVSELNX16( y1, 0 );

        BBE_LVNX16_IP( a0, Y_rd, 2*BBE_SIMD_WIDTH );

        a1 = 0;

        y0 = BBE_SELNX16( a1, a0, sel0 );
        y1 = BBE_SELNX16( a1, a0, sel1 );
      }

      __Pragma( "loop_count min=1" );
      for ( i=0; i<M/4; i++ )
      {
        BBE_LVNX16_IP( x0, X_rd, 2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( x1, X_rd, 2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( x2, X_rd, 2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( x3, X_rd, 2*BBE_SIMD_WIDTH );

        BBE_DSELNX16I( x1, x0, x1, x0, BBE_DSELI_INTERLEAVE_1 );
        BBE_DSELNX16I( x3, x2, x3, x2, BBE_DSELI_INTERLEAVE_1 );

        w0 = BBE_MULRNX16PC_0( x0, y0, vsa0 );

        BBE_MULANX16PC_0( w0, x1, y1 );
        BBE_MULANX16PC_1( w0, x2, y0 );
        BBE_MULANX16PC_1( w0, x3, y1 );

        z0 = BBE_PACKVNX40( w0, vsa0 );

        BBE_SVNX16_IP( z0, Z_wr, 2*BBE_SIMD_WIDTH );
      }
    }

    __Pragma( "no_reorder" );

    //--------------------------------------------------------------------------
    // Convert BBE_SIMD_WIDTH/4 resulting vectors to block format: 
    // z[M/4][BBE_SIMD_WIDTH/4][4] => Z[BBE_SIMD_WIDTH/4][_M].
    //

    Z_rd = Z_scr;
    Z_wr = (xb_vecNx16*)( (uintptr_t)z + l*(BBE_SIMD_WIDTH/4)*_M*2 );

    {
      int n;

      //__Pragma( "loop_count min=1" );
      for ( n=1; n<_M/BBE_SIMD_WIDTH; n++ )
      {
        xb_vecNx16 a0, a1, a2, a3;

        BBE_LVNX16_IP( a0, Z_rd, 2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a1, Z_rd, 2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a2, Z_rd, 2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a3, Z_rd, 2*BBE_SIMD_WIDTH );

        BBE_DSELNX16I( a2, a0, a2, a0, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( a3, a1, a3, a1, BBE_DSELI_INTERLEAVE_4 );

        BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( a3, a2, a3, a2, BBE_DSELI_INTERLEAVE_4 );

        BBE_SVNX16_XP( a0, Z_wr, _M*2 );
        BBE_SVNX16_XP( a1, Z_wr, _M*2 );
        BBE_SVNX16_XP( a2, Z_wr, _M*2 );

        BBE_SVNX16_XP( a3, Z_wr, -3*_M*2 + BBE_SIMD_WIDTH*2 );
      }
      {
        xb_vecNx16 a0, a1, a2, a3;
        a1 = BBE_LVNX16_X( Z_rd, off0*2*BBE_SIMD_WIDTH );
        a2 = BBE_LVNX16_X( Z_rd, off1*4*BBE_SIMD_WIDTH );
        a3 = BBE_LVNX16_X( Z_rd, off2*6*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a0, Z_rd, 8*BBE_SIMD_WIDTH );

        BBE_DSELNX16I( a2, a0, a2, a0, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( a3, a1, a3, a1, BBE_DSELI_INTERLEAVE_4 );

        BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_INTERLEAVE_4 );
        BBE_DSELNX16I( a3, a2, a3, a2, BBE_DSELI_INTERLEAVE_4 );

        BBE_SVNX16_XP( a0, Z_wr, _M*2 );
        BBE_SVNX16_XP( a1, Z_wr, _M*2 );
        BBE_SVNX16_XP( a2, Z_wr, _M*2 );

        BBE_SVNX16_XP( a3, Z_wr, -3*_M*2 + BBE_SIMD_WIDTH*2 );
      }
    }
  }

} // matvmulnxmn_Mgte12_N4_L4x()

size_t matvmulnxmn_Mgte12_N4_getScratchSize(int N, int M)
{
  int Sx, _M;
  _M = ((M + (BBE_SIMD_WIDTH - 1)) & ~(BBE_SIMD_WIDTH - 1));
  Sx = getSpaceR(N*M);
  return (Sx * 2 * 4 + _M * 4)*sizeof(int16_t);
} /* matvmulnxmn_Mgte12_N4_getScratchSize() */
#endif
