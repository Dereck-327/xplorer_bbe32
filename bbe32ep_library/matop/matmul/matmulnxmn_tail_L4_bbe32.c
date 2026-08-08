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


/* Compute L matrix products MxN * NxM -> MxM, where L equals 1, 2 or 3. */
void matmulnxmn_tail_L4( void * pScr,
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
        xb_vecNx16 zero=0;

  int MN, MM;
  int l;

  NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);

  NASSERT( N%4==0 && M%4==0 );

  NASSERT( Q>=0 && Q<=16 );

  NASSERT( L <= (BBE_SIMD_WIDTH/4-1) );

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



  X_rd = (const xb_vecNx16*)x;
  Y_rd = (const xb_vecNx16*)y;
  X_wr = (      xb_vecNx16*)X_scr;
  Y_wr = (      xb_vecNx16*)Y_scr;
 /*
  * Fill scretch area by zeroes
  */
  {
    int n;
    for ( n=0; n<MN/(BBE_SIMD_WIDTH); n++ )
    {
        BBE_SVNX16_IP( zero, X_wr, 4*BBE_SIMD_WIDTH/2);
        BBE_SVNX16_IP( zero, X_wr, 4*BBE_SIMD_WIDTH/2);
        BBE_SVNX16_IP( zero, X_wr, 4*BBE_SIMD_WIDTH/2);
        BBE_SVNX16_IP( zero, X_wr, 4*BBE_SIMD_WIDTH/2);
            
        BBE_SVNX16_IP( zero, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( zero, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( zero, Y_wr, 4*BBE_SIMD_WIDTH/2 );
        BBE_SVNX16_IP( zero, Y_wr, 4*BBE_SIMD_WIDTH/2 );
    }
  }

 /*
  * Copy input data to the scratch. The rest of function code assumes that
  * we have BBE_SIMD_WIDTH/4 matrix pairs on input.
  */
  {
    xb_vecNx16 x0, y0;

    uint32_t ptrAdj = 4*BBE_SIMD_WIDTH/2 - MN*(BBE_SIMD_WIDTH/8)*4;

    int n;
    X_rd = (const xb_vecNx16*)x;
    Y_rd = (const xb_vecNx16*)y;
    X_wr = (      xb_vecNx16*)X_scr;
    Y_wr = (      xb_vecNx16*)Y_scr;
    for ( l=0; l<L; l++ )
    {
      for ( n=0; n<MN/(BBE_SIMD_WIDTH); n++ )
      {
        BBE_LVNX16_IP( x0, X_rd, 4*BBE_SIMD_WIDTH/2 );
        BBE_LVNX16_IP( y0, Y_rd, 4*BBE_SIMD_WIDTH/2 );

        BBE_SVNX16_XP( x0, X_wr, (BBE_SIMD_WIDTH/8)*4*BBE_SIMD_WIDTH );
        BBE_SVNX16_XP( y0, Y_wr, (BBE_SIMD_WIDTH/8)*4*BBE_SIMD_WIDTH );
      }

      X_wr = (xb_vecNx16*)( (uintptr_t)X_wr + ptrAdj );
      Y_wr = (xb_vecNx16*)( (uintptr_t)Y_wr + ptrAdj );
    }
  }

  /*
  * In-place conversion of BBE_SIMD_WIDTH/4 MxN matrices x[] and y[] to
  * streaming format.
  */

  X_rd = X_wr = X_scr;
  Y_rd = Y_wr = Y_scr;

  {
    xb_vecNx16 a00, a01, a02, a03;
    xb_vecNx16 a10, a11, a12, a13;
    xb_vecNx16 b00, b01, b02, b03;
    xb_vecNx16 b10, b11, b12, b13;

    int n;
    for ( n=0; n<MN/(BBE_SIMD_WIDTH); n++ )
    {
      BBE_LVNX16_IP( a00, X_rd, 4*BBE_SIMD_WIDTH/2);
      BBE_LVNX16_IP( a01, X_rd, 4*BBE_SIMD_WIDTH/2);
      BBE_LVNX16_IP( a02, X_rd, 4*BBE_SIMD_WIDTH/2);
      BBE_LVNX16_IP( a03, X_rd, 4*BBE_SIMD_WIDTH/2);

      BBE_LVNX16_IP( a10, Y_rd, 4*BBE_SIMD_WIDTH/2 );
      BBE_LVNX16_IP( a11, Y_rd, 4*BBE_SIMD_WIDTH/2 );
      BBE_LVNX16_IP( a12, Y_rd, 4*BBE_SIMD_WIDTH/2 );
      BBE_LVNX16_IP( a13, Y_rd, 4*BBE_SIMD_WIDTH/2 );

      BBE_DSELNX16I( b01, b00, a01, a00, BBE_DSELI_DEINTERLEAVE_1 );
      BBE_DSELNX16I( b03, b02, a03, a02, BBE_DSELI_DEINTERLEAVE_1 );
                                        
      BBE_DSELNX16I( a02, a00, b02, b00, BBE_DSELI_DEINTERLEAVE_1 );
      BBE_DSELNX16I( a03, a01, b03, b01, BBE_DSELI_DEINTERLEAVE_1 );
                                        
      BBE_DSELNX16I( b01, b00, a01, a00, BBE_DSELI_DEINTERLEAVE_1 );
      BBE_DSELNX16I( b03, b02, a03, a02, BBE_DSELI_DEINTERLEAVE_1 );
                                                                
      BBE_DSELNX16I( a02, a00, b02, b00, BBE_DSELI_DEINTERLEAVE_1 );
      BBE_DSELNX16I( a03, a01, b03, b01, BBE_DSELI_DEINTERLEAVE_1 );

      BBE_DSELNX16I( b11, b10, a11, a10, BBE_DSELI_DEINTERLEAVE_1 );
      BBE_DSELNX16I( b13, b12, a13, a12, BBE_DSELI_DEINTERLEAVE_1 );
                                                                
      BBE_DSELNX16I( a12, a10, b12, b10, BBE_DSELI_DEINTERLEAVE_1 );
      BBE_DSELNX16I( a13, a11, b13, b11, BBE_DSELI_DEINTERLEAVE_1 );
                                                                
      BBE_DSELNX16I( b11, b10, a11, a10, BBE_DSELI_DEINTERLEAVE_1 );
      BBE_DSELNX16I( b13, b12, a13, a12, BBE_DSELI_DEINTERLEAVE_1 );
                                                                
      BBE_DSELNX16I( a12, a10, b12, b10, BBE_DSELI_DEINTERLEAVE_1 );
      BBE_DSELNX16I( a13, a11, b13, b11, BBE_DSELI_DEINTERLEAVE_1 );

      BBE_SVNX16_IP( a00, X_wr, 4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( a01, X_wr, 4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( a02, X_wr, 4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( a03, X_wr, 4*BBE_SIMD_WIDTH/2 );

      BBE_SVNX16_IP( a10, Y_wr, 4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( a11, Y_wr, 4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( a12, Y_wr, 4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( a13, Y_wr, 4*BBE_SIMD_WIDTH/2 );
    }
  }

  /*
  * Compute BBE_SIMD_WIDTH/4 MxN*NxM matrix products in streaming format.
  */

  {
    xb_vecNx16 a0, a1, a2, a3;
    xb_vecNx16 x0, x1, x2, x3;

    int i, k;

    /*
    * Prior to matrix multiplication loop, transpose each 4x4 sub-block of
    * the left-hand matrices x[].
    */

    for ( i=0; i<M/4; i++ )
    {
      X_rd = (const xb_vecNx16*)( (uintptr_t)X_scr + 4*i*N*(BBE_SIMD_WIDTH/8)*4 );
      X_wr = (      xb_vecNx16*)( (uintptr_t)X_scr + 4*i*N*(BBE_SIMD_WIDTH/8)*4 );

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
      }
    }
  }

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
    vselN vs1, vs2;

    static const int16_t ALIGN(64) seli[][BBE_SIMD_WIDTH] = {
      { 0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3 },
    };

    int i, j, k;

    vs0 = BBE_MOVVSA32( Q );

    a0 = BBE_LVNX16_I( (const xb_vecNx16*)seli, 0 );

    vs1 = BBE_MOVVSELNX16( a0, 0 );
    vs2 = BBE_MOVVSELNX16( a0, 0 );

    Z_wr = Z_scr;

    for ( i=0; i<M/4; i++ )
    {
      X_rd = (const xb_vecNx16*)( (uintptr_t)X_scr + 4*i*N*(BBE_SIMD_WIDTH/8)*4 );

      for ( j=0; j<M/4; j++ )
      {
        Y_rd = (const xb_vecNx16*)( (uintptr_t)Y_scr + 4*j*(BBE_SIMD_WIDTH/8)*4 );

        a0 = 0;

        w0 = w1 = w2 = w3 = BBE_MULRNX16( a0, a0, vs0 );

        __Pragma( "ymemory( X_rd )" )
        __Pragma( "ymemory( Y_rd )" )
        __Pragma( "loop_count min=1" )
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

          BBE_SHFLUNX16( y00, a0, vs1, 4 ); BBE_SHFLUNX16( y10, a1, vs2, 4 );
          BBE_SHFLUNX16( y01, a0, vs1, 4 ); BBE_SHFLUNX16( y11, a1, vs2, 4 );
          BBE_SHFLUNX16( y02, a0, vs1, 4 ); BBE_SHFLUNX16( y12, a1, vs2, 4 );
          BBE_SHFLUNX16( y03, a0, vs1, 4 ); BBE_SHFLUNX16( y13, a1, vs2, 4 );
                                                                         
          BBE_SHFLUNX16( y20, a2, vs1, 4 ); BBE_SHFLUNX16( y30, a3, vs2, 4 );
          BBE_SHFLUNX16( y21, a2, vs1, 4 ); BBE_SHFLUNX16( y31, a3, vs2, 4 );
          BBE_SHFLUNX16( y22, a2, vs1, 4 ); BBE_SHFLUNX16( y32, a3, vs2, 4 );
          BBE_SHFLUNX16( y23, a2, vs1, 4 ); BBE_SHFLUNX16( y33, a3, vs2, 4 );

          BBE_MULANX16( w0, x0, y00 ); BBE_MULANX16( w0, x1, y10 );
          BBE_MULANX16( w1, x0, y01 ); BBE_MULANX16( w1, x1, y11 );
          BBE_MULANX16( w2, x0, y02 ); BBE_MULANX16( w2, x1, y12 );
          BBE_MULANX16( w3, x0, y03 ); BBE_MULANX16( w3, x1, y13 );

          BBE_MULANX16( w0, x2, y20 ); BBE_MULANX16( w0, x3, y30 );
          BBE_MULANX16( w1, x2, y21 ); BBE_MULANX16( w1, x3, y31 );
          BBE_MULANX16( w2, x2, y22 ); BBE_MULANX16( w2, x3, y32 );
          BBE_MULANX16( w3, x2, y23 ); BBE_MULANX16( w3, x3, y33 );
        }

        X_rd = (const xb_vecNx16*)( (uintptr_t)X_rd - N*(BBE_SIMD_WIDTH/8)*4 );

        a0 = BBE_PACKVNX40( w0, vs0 );
        a1 = BBE_PACKVNX40( w1, vs0 );
        a2 = BBE_PACKVNX40( w2, vs0 );
        a3 = BBE_PACKVNX40( w3, vs0 );

        /*
        * 4x4 sub-blocks of BBE_SIMD_WIDTH/8 resulting matrices z[] appear
        * transposed here.
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
  * In-place conversion of BBE_SIMD_WIDTH/4 MxM matrices z[] to block format.
  */

  Z_rd = Z_wr = Z_scr;

  {
    xb_vecNx16 a0, a1, a2, a3;
    xb_vecNx16 b0, b1, b2, b3;

    int n;

    for ( n=0; n<MM/(BBE_SIMD_WIDTH); n++ )
    {
      BBE_LVNX16_IP( a0, Z_rd, 4*BBE_SIMD_WIDTH/2 );
      BBE_LVNX16_IP( a1, Z_rd, 4*BBE_SIMD_WIDTH/2 );
      BBE_LVNX16_IP( a2, Z_rd, 4*BBE_SIMD_WIDTH/2 );
      BBE_LVNX16_IP( a3, Z_rd, 4*BBE_SIMD_WIDTH/2 );

      BBE_DSELNX16I( b1, b0, a1, a0, BBE_DSELI_DEINTERLEAVE_1 );
      BBE_DSELNX16I( b3, b2, a3, a2, BBE_DSELI_DEINTERLEAVE_1 );
                                     
      BBE_DSELNX16I( a2, a0, b2, b0, BBE_DSELI_DEINTERLEAVE_1 );
      BBE_DSELNX16I( a3, a1, b3, b1, BBE_DSELI_DEINTERLEAVE_1 );

      BBE_SVNX16_IP( a0, Z_wr, 4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( a1, Z_wr, 4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( a2, Z_wr, 4*BBE_SIMD_WIDTH/2 );
      BBE_SVNX16_IP( a3, Z_wr, 4*BBE_SIMD_WIDTH/2 );
    }
  }

  /* Copy 1..BBE_SIMD_WIDTH/4-1 product matrices z[] to the output array.*/

  Z_rd = (const xb_vecNx16*)Z_scr;
  Z_wr = (      xb_vecNx16*)z;

  {
    xb_vecNx16 z0;

    uint32_t ptrAdj = 4*BBE_SIMD_WIDTH/2 - MM*(BBE_SIMD_WIDTH/8)*4;

    int n;

    for ( l=0; l<L; l++ )
    {
      for ( n=0; n<MM/(BBE_SIMD_WIDTH); n++ )
      {
        BBE_LVNX16_XP( z0, Z_rd, (BBE_SIMD_WIDTH/8)*4*BBE_SIMD_WIDTH );

        BBE_SVNX16_IP( z0, Z_wr, 4*BBE_SIMD_WIDTH/2 );
      }

      Z_rd = (xb_vecNx16*)( (uintptr_t)Z_rd + ptrAdj );
    }
  }
} /* matmulnxmn_tail_L4() */

/* Return the scratch area size, in bytes. */
size_t matmulnxmn_tail_L4_getScratchSize(int N, int M)
{
  int Sx, Sz;
  Sx = getSpaceR(N*M);
  Sz = getSpaceR(M*M);
  return (Sx * 2 * 4 + Sz * 4)*sizeof(int16_t);
} /* matmulnxmn_tail_L4_getScratchSize() */
