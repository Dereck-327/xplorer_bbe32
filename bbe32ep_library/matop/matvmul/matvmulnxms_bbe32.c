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
/* Common helper macros. */
#include "common.h"


/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#include <string.h>

void matvmulnxms_1(int16_t * restrict z,
            const int16_t * restrict x,
            const int16_t * restrict y,
            int N, int M, int L, int Q);

void matvmulnxms_2(int16_t * restrict z,
            const int16_t * restrict x,
            const int16_t * restrict y,
            int N, int M, int L, int Q);

void matvmulnxms_3(int16_t * restrict z,
            const int16_t * restrict x,
            const int16_t * restrict y,
            int N, int M, int L, int Q);
/*-------------------------------------------------------------------------
Real Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of real matrices or vectors. Both the block order and streaming 
order are allowed for input/output matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand input matrices
y[L*Sy]     Sequence of right-hand input matrices
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Streaming Order, MxN*Nx1->Mx1, Sx=MxN, Sy=N, Sz=M
   Restrictions:
     L must be a multiple of 16
*/
void  matvmulnxms ( int16_t * restrict z,
              const int16_t * restrict x,
              const int16_t * restrict y,
              int N, int M, int L, int Q )
{
/* check restrictions */
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT( L%BBE_SIMD_WIDTH == 0);

  if ((L<=0)||(M<=0)) return;
  if (N <= 0)
  {
    memset(z, 0, M*L*sizeof(int16_t));
    return;
  }
  if (!(M & 1))
  {
    matvmulnxms_1(z, x, y, N, M, L, Q);
  }
  else if (N & 1)
  {
    matvmulnxms_2(z, x, y, N, M, L, Q);
  }
  else
  {
    matvmulnxms_3(z, x, y, N, M, L, Q);
  }
} /* matvmulnxms() */

void matvmulnxms_1(int16_t * restrict z,
             const int16_t * restrict x,
             const int16_t * restrict y,
             int N, int M, int L, int Q)
{
  vsaN  q = BBE_MOVVSA32(Q);
    
  int i, j, k;

  xb_vecNx40 z_out0;
  xb_vecNx40 z_out2;

  xb_vecNx16 x00;
  xb_vecNx16 y00;

  xb_vecNx16 * restrict px00;
  xb_vecNx16 * restrict px01;

  xb_vecNx16 * restrict py00;
  xb_vecNx16 * restrict py01;

  xb_vecNx16 * restrict pz00;

  px00 = (xb_vecNx16*)x;
  py00 = (xb_vecNx16*)y;
  pz00 = (xb_vecNx16*)z;

  /* check restrictions */
  NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
  NASSERT(Q>=0 && Q<=16);
  NASSERT(L%BBE_SIMD_WIDTH==0);

  L>>=1;
  if ( N == 1 )
  {
    for ( i=0; i<M; i+=2 ) 
    {
      px01 = px00;
      py01 = py00;
            
      __Pragma( "ymemory( px01 )" );
      __Pragma( "loop_count min=1" );
      for(j=0; j<2*L; j+=BBE_SIMD_WIDTH) 
      {
        /* load 16 point Y */
        BBE_LVNX16_IP(y00, py01, 2*BBE_SIMD_WIDTH);
        /* load 16 point X */
        x00 = BBE_LVNX16_X(px01, 4*L);
        z_out2 = BBE_MULNX16(x00, y00);
        /* load 16 point X */
        BBE_LVNX16_IP(x00, px01, 2*BBE_SIMD_WIDTH);
        z_out0 = BBE_MULNX16(x00, y00);

        x00 = BBE_PACKVNX40(z_out2, q);
        /* save row 16 result matrix [NxN][L] */
        BBE_SVNX16_X(x00 , pz00, 4*L);

        x00 = BBE_PACKVNX40(z_out0,	q);
        /* save row 16 result matrix [NxN][L] */
        BBE_SVNX16_IP(x00 , pz00, sizeof(*pz00)); 
      }	

      px00 = (xb_vecNx16*)XT_ADDX2(4*L, (int32_t)px00);
      pz00 = (xb_vecNx16*)XT_ADD(4*L, (int32_t)pz00);
    }

    return;
  }
  __Pragma("loop_count min=1");
  for ( i=0; i<M; i+=2 ) 
  {
    px01 = px00;
    py01 = py00;
        
    for(j=0; j<2*L; j+=BBE_SIMD_WIDTH) 
    {
      /* load 16 point Y */
      BBE_LVNX16_XP(y00, py01, 4*L);
      /* load 16 point X */
      x00 = BBE_LVNX16_X(px01, 4*L*N);
      z_out2 = BBE_MULNX16(x00, y00);
      /* load 16 point X */
      BBE_LVNX16_XP(x00, px01, 4*L);
      z_out0 = BBE_MULNX16(x00, y00);

      __Pragma( "ymemory( px01 )" );
      __Pragma( "loop_count min=1" );
      /* loop compute 1 point for 16 result matrix [M][M] */
      for ( k=0; k<N-1; k++ )
      {
        y00 = BBE_LVNX16_I(py01, 0);
        /* load 16 point X */
        x00 = BBE_LVNX16_X(px01, 4*L*N);
        BBE_MULANX16(z_out2 , x00, y00);
        /* load 16 point Y */
        BBE_LVNX16_XP(y00, py01, 4*L);
        /* load 16 point X */
        BBE_LVNX16_XP(x00, px01, 4*L);
        BBE_MULANX16(z_out0 , x00, y00);
      }

      px01 = (xb_vecNx16 *)XT_ADDX2(-2*N*L + BBE_SIMD_WIDTH, (int32_t)px01);
      py01 = (xb_vecNx16 *)XT_ADDX2(-2*N*L + BBE_SIMD_WIDTH, (int32_t)py01);

      x00 = BBE_PACKVNX40(z_out2, q);
      /* save row 16 result matrix [NxN][L] */
      BBE_SVNX16_X(x00 , pz00, 4*L); 

      x00 = BBE_PACKVNX40(z_out0,	q);
      /* save row 16 result matrix [NxN][L] */
      BBE_SVNX16_IP(x00 , pz00, sizeof(*pz00)); 
    }	
    px00 = (xb_vecNx16*)XT_ADDX2(4*N*L, (int32_t)px00);
    pz00 = (xb_vecNx16*)XT_ADD(4*L, (int32_t)pz00);
  }
} /* matvmulnxms_1() */


void matvmulnxms_2(int16_t * restrict z,
             const int16_t * restrict x,
             const int16_t * restrict y,
             int N, int M, int L, int Q)
{
    vsaN  q = BBE_MOVVSA32(Q);
    const int L16 = L>>4;

    int i, j, k;

    xb_vecNx40 z_out0;
    xb_vecNx40 z_out2;

    xb_vecNx16 x00;
    xb_vecNx16 y00;

    xb_vecNx16 * restrict px00;
    xb_vecNx16 * restrict px01;

    xb_vecNx16 * restrict py00;
    xb_vecNx16 * restrict py01;

    xb_vecNx16 * restrict pz00;

    px00 = (xb_vecNx16*)x;
    py00 = (xb_vecNx16*)y;
    pz00 = (xb_vecNx16*)z;

    /* check restrictions */
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT(Q>=0 && Q<=16);
    NASSERT( L%BBE_SIMD_WIDTH==0);

    L>>=1;
    if (N == 1)
    {
      for ( i=0; i<M-1; i+=2 ) 
      {
        __Pragma( "ymemory( px01 )" );
        __Pragma( "loop_count min=1" );
        for(j=0; j<L16; j++) 
        {
          px01 = px00 + j;
          py01 = py00 + j;
          /* load 16 point Y */
          BBE_LVNX16_XP(y00, py01, 4*L);
          /* load 16 point X */
          x00 = BBE_LVNX16_X(px01, 4*L*N);
          z_out2 = BBE_MULNX16(x00, y00);
          /* load 16 point X */
          BBE_LVNX16_XP(x00, px01, 4*L);
          z_out0 = BBE_MULNX16(x00, y00);

          x00 = BBE_PACKVNX40(z_out2, q);
          /* save row 16 result matrix [NxN][L] */
          BBE_SVNX16_X(x00 , pz00, 4*L); 

          x00 = BBE_PACKVNX40(z_out0,	q);
          /* save row 16 result matrix [NxN][L] */
          BBE_SVNX16_IP(x00 , pz00, sizeof(*pz00)); 
        }	

        px00 = (xb_vecNx16*)XT_ADDX2(4*N*L, (int32_t)px00);
        pz00 = (xb_vecNx16*)XT_ADD(4*L, (int32_t)pz00);
      }

      __Pragma( "loop_count min=1" );
      for(j=0; j<L16; j++) 
      {
        px01 = px00 + j;
        py01 = py00 + j;
        /* load 16 point X */
        BBE_LVNX16_XP(x00, px01, 2*L);
        /* load 16 point Y */
        BBE_LVNX16_XP(y00, py01, 2*L);
        z_out0 = BBE_MULNX16(x00, y00);

        x00 = BBE_PACKVNX40(z_out0, q);
        /* save row 16 result matrix [NxN][L] */
        BBE_SVNX16_IP(x00 , pz00, sizeof(*pz00)); 
      }	

      return;
    }

    for ( i=0; i<M-1; i+=2 ) 
    {
      px01 = px00;
      py01 = py00;

      for(j=0; j<L16; j++) 
      {
        /* load 16 point Y */
        BBE_LVNX16_XP(y00, py01, 4*L);
        /* load 16 point X */
        x00 = BBE_LVNX16_X(px01, 4*L*N);
        z_out2 = BBE_MULNX16(x00, y00);
        /* load 16 point X */
        BBE_LVNX16_XP(x00, px01, 4*L);
        z_out0 = BBE_MULNX16(x00, y00);

        __Pragma( "ymemory( px01 )" );
        __Pragma( "loop_count min=1" );
        /* loop compute 1 point for 16 result matrix [M][M] */
        for ( k=0; k<N-1; k++ )
        {
          y00 = BBE_LVNX16_I(py01, 0);
          /* load 16 point X */
          x00 = BBE_LVNX16_X(px01, 4*L*N);
          BBE_MULANX16(z_out2 , x00, y00);
          /* load 16 point Y */
          BBE_LVNX16_XP(y00, py01, 4*L);
          /* load 16 point X */
          BBE_LVNX16_XP(x00, px01, 4*L);
          BBE_MULANX16(z_out0 , x00, y00);
        }

        px01 = (xb_vecNx16 *)XT_ADDX2(-2*N*L + BBE_SIMD_WIDTH, (int32_t)px01);
        py01 = (xb_vecNx16 *)XT_ADDX2(-2*N*L + BBE_SIMD_WIDTH, (int32_t)py01);

        x00 = BBE_PACKVNX40(z_out2, q);
        /* save row 16 result matrix [NxN][L] */
        BBE_SVNX16_X(x00 , pz00, 4*L); 

        x00 = BBE_PACKVNX40(z_out0,	q);
        /* save row 16 result matrix [NxN][L] */
        BBE_SVNX16_IP(x00 , pz00, sizeof(*pz00)); 
      }	

      px00 = (xb_vecNx16*)XT_ADDX2(4*N*L, (int32_t)px00);
      pz00 = (xb_vecNx16*)XT_ADD(4*L, (int32_t)pz00);
    }

    px01 = px00;
    py01 = py00;

    __Pragma( "loop_count min=1" );
    for(j=0; j<L16; j++) 
    {
      /* load 16 point X */
      BBE_LVNX16_XP(x00, px01, 4*L);
      /* load 16 point Y */
      BBE_LVNX16_XP(y00, py01, 4*L);
      z_out0 = BBE_MULNX16(x00, y00);

      __Pragma( "ymemory( px01 )" );
      __Pragma( "loop_count min=1" );
      /* loop compute 1 point for 16 result matrix [M][M] */
      for ( k=0; k<N-1; k+=2 )
      {
        /* load 16 point X */
        BBE_LVNX16_XP(x00, px01, 4*L);
        /* load 16 point Y */
        BBE_LVNX16_XP(y00, py01, 4*L);
        BBE_MULANX16(z_out0 , x00, y00);

        /* load 16 point X */
        BBE_LVNX16_XP(x00, px01, 4*L);
        /* load 16 point Y */
        BBE_LVNX16_XP(y00, py01, 4*L);
        BBE_MULANX16(z_out0 , x00, y00);
      }

      px01 = (xb_vecNx16 *)XT_ADDX2(-2*N*L + BBE_SIMD_WIDTH, (int32_t)px01);
      py01 = (xb_vecNx16 *)XT_ADDX2(-2*N*L + BBE_SIMD_WIDTH, (int32_t)py01);

      x00 = BBE_PACKVNX40(z_out0, q);
      /* save row 16 result matrix [NxN][L] */
      BBE_SVNX16_IP(x00 , pz00, sizeof(*pz00)); 
    }	
} /* matvmulnxms_2() */


void matvmulnxms_3(int16_t * restrict z,
             const int16_t * restrict x,
             const int16_t * restrict y,
             int N, int M, int L, int Q)
{
    vsaN  q = BBE_MOVVSA32(Q);
    const int L16 = L>>4;

    int i, j, k;

    xb_vecNx40 z_out0;
    xb_vecNx40 z_out2;

    xb_vecNx16 x00;
    xb_vecNx16 y00;

    xb_vecNx16 * restrict px00;
    xb_vecNx16 * restrict px01;

    xb_vecNx16 * restrict py00;
    xb_vecNx16 * restrict py01;

    xb_vecNx16 * restrict pz00;

    px00 = (xb_vecNx16*)x;
    py00 = (xb_vecNx16*)y;
    pz00 = (xb_vecNx16*)z;

    /* check restrictions */
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT(Q>=0 && Q<=16);
    NASSERT(L%BBE_SIMD_WIDTH==0);

    L>>=1;
    for ( i=0; i<M-1; i+=2 ) 
    {
      px01 = px00;
      py01 = py00;

      for(j=0; j<L16; j++) 
      {
        /* load 16 point Y */
        BBE_LVNX16_XP(y00, py01, 4*L);
        /* load 16 point X */
        x00 = BBE_LVNX16_X(px01, 4*L*N);
        z_out2 = BBE_MULNX16(x00, y00);
        /* load 16 point X */
        BBE_LVNX16_XP(x00, px01, 4*L);
        z_out0 = BBE_MULNX16(x00, y00);

        __Pragma( "ymemory( px01 )" );
        __Pragma( "loop_count min=1" );
        /* loop compute 1 point for 16 result matrix [M][M] */
        for ( k=0; k<N-1; k++ )
        {
          y00 = BBE_LVNX16_I(py01, 0);
          /* load 16 point X */
          x00 = BBE_LVNX16_X(px01, 4*L*N);	
          BBE_MULANX16(z_out2 , x00, y00);
          /* load 16 point Y */
          BBE_LVNX16_XP(y00, py01, 4*L);
          /* load 16 point X */
          BBE_LVNX16_XP(x00, px01, 4*L);
          BBE_MULANX16(z_out0 , x00, y00);
        }

        px01 = (xb_vecNx16 *)XT_ADDX2(-2*N*L + BBE_SIMD_WIDTH, (int32_t)px01);
        py01 = (xb_vecNx16 *)XT_ADDX2(-2*N*L + BBE_SIMD_WIDTH, (int32_t)py01);

        x00 = BBE_PACKVNX40(z_out2, q);
        /* save row 16 result matrix [NxN][L] */
        BBE_SVNX16_X(x00 , pz00, 4*L); 

        x00 = BBE_PACKVNX40(z_out0,	q);
        /* save row 16 result matrix [NxN][L] */
        BBE_SVNX16_IP(x00 , pz00, sizeof(*pz00)); 
      }	
      px00 = (xb_vecNx16*)XT_ADDX2(4*N*L, (int32_t)px00);
      pz00 = (xb_vecNx16*)XT_ADD(4*L, (int32_t)pz00);
    }

    __Pragma( "loop_count min=1" );
    for(j=0; j<L16; j++) 
    {
      px01 = px00 + j;
      py01 = py00 + j;
      /* load 16 point X */
      BBE_LVNX16_XP(x00, px01, 4*L);
      /* load 16 point Y */
      BBE_LVNX16_XP(y00, py01, 4*L);
      z_out0 = BBE_MULNX16(x00, y00);

      __Pragma( "ymemory( px01 )" );
      for ( k=0; k<N-2; k+=2 )
      {
        /* load 16 point X */
        BBE_LVNX16_XP(x00, px01, 4*L);
        /* load 16 point Y */
        BBE_LVNX16_XP(y00, py01, 4*L);
        BBE_MULANX16(z_out0 , x00, y00);
        /* load 16 point X */
        BBE_LVNX16_XP(x00, px01, 4*L);
        /* load 16 point Y */
        BBE_LVNX16_XP(y00, py01, 4*L);
        BBE_MULANX16(z_out0 , x00, y00);
      }
      /* load 16 point X */
      BBE_LVNX16_XP(x00, px01, 4*L);
      /* load 16 point Y */
      BBE_LVNX16_XP(y00, py01, 4*L);
      BBE_MULANX16(z_out0 , x00, y00);

      x00 = BBE_PACKVNX40(z_out0, q);
      /* save row 16 result matrix [NxN][L] */
      BBE_SVNX16_IP(x00 , pz00, sizeof(*pz00)); 
    }	
} /* matvmulnxms_3() */
