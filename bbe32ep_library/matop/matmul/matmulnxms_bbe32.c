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

/* Streaming Order, MxN * NxM ->MxM, Sx=MxN,Sy=NxM,Sz=MxM
   Restrictions:
   Q=0..15, L is a multiple of 8
*/
static void matmulnxms_1(int16_t * restrict z, 
                   const int16_t * restrict x, 
                   const int16_t * restrict y, 
                   int N, int M, int L, int Q);

/* Streaming Order, MxN * NxM ->MxM, Sx=MxN,Sy=NxM,Sz=MxM
   Restrictions:
   Q=0..15, L is a multiple of 8
*/
static void matmulnxms_2(int16_t * restrict z, 
                   const int16_t * restrict x, 
                   const int16_t * restrict y, 
                   int N, int M, int L, int Q);

/* Streaming Order, MxN * NxM ->MxM, Sx=MxN,Sy=NxM,Sz=MxM
   Restrictions:
   Q=0..15, L is a multiple of 8
*/
static void matmulnxms_3(int16_t * restrict z, 
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

/* Streaming Order, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
     L must be a multiple of 16
*/
void matmulnxms ( int16_t * restrict z,
            const int16_t * restrict x,
            const int16_t * restrict y,
            int N, int M, int L, int Q )
{
  NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
  NASSERT(Q>=0 && Q<=16);
  NASSERT(L%BBE_SIMD_WIDTH==0);

  if (L<=0) return;
  if (N <= 0)
  {
    memset(z, 0, M*M*L*sizeof(int16_t));
    return;
  }
  if(!(N&1) && !(M&1) )
  {
      matmulnxms_3(z, x, y, N, M, L, Q);
  }
  else if ( M&1 )
  {
      matmulnxms_1(z, x, y, N, M, L, Q);
  }
  else
  {
      matmulnxms_2(z, x, y, N, M, L, Q);
  }
} /* matmulnxms() */

void matmulnxms_1(int16_t * restrict z,
            const int16_t * restrict x,
            const int16_t * restrict y,
            int N, int M, int L, int Q)
{
  int i, j, k;

  int offset_y = 0;
  int offset_x = 0;
  const uint32_t mod_y = (BBE_SIMD_WIDTH * 2) + ((L*M * 2) << 16);
  const uint32_t mod_x = (BBE_SIMD_WIDTH * 2) + ((L * 2) << 16);
  vsaN q = BBE_MOVVSA32(Q);

  xb_vecNx16 x00, x01;
  xb_vecNx16 y00, y10;
  xb_vecNx40 z_out0, z_out1;

  xb_vecNx16 * restrict px00;
  xb_vecNx16 * restrict py00;
  xb_vecNx16 * restrict px01;
  xb_vecNx16 * restrict py01;
  xb_vecNx16 * restrict pz00;

  px00 = (xb_vecNx16*)x;
  py00 = (xb_vecNx16*)y;
  pz00 = (xb_vecNx16*)z;

  /* check restrictions */
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT( L%BBE_SIMD_WIDTH == 0);
  NASSERT(Q >= 0 && Q <= 16);

  pz00 = (xb_vecNx16*)(z);
  for (i = 0; i<M; i++)
  {
    for (j = 0; j<(M - 1)*L; j += 32)
    {
      z_out0 = 0;
      z_out1 = 0;

      py00 = (xb_vecNx16*)(offset_y + (uintptr_t)y);
      px00 = (xb_vecNx16*)(offset_x + (uintptr_t)x);

      offset_y = BBE_ADDMOD16U(offset_y, mod_y);
      offset_x = BBE_ADDMOD16U(offset_x, mod_x);

      py01 = (xb_vecNx16*)(offset_y + (uintptr_t)y);
      px01 = (xb_vecNx16*)(offset_x + (uintptr_t)x);

      offset_y = BBE_ADDMOD16U(offset_y, mod_y);
      offset_x = BBE_ADDMOD16U(offset_x, mod_x);

      __Pragma("ymemory( px00 )");
      __Pragma("ymemory( px01 )");
      __Pragma("loop_count min=1");
      for (k = 0; k<N; k++)
      {
        /* load 16 point X */
        BBE_LVNX16_XP(x01, px01, L * 2);
        /* load 16 point Y */
        BBE_LVNX16_XP(y10, py01, L*M * 2);
        BBE_MULANX16(z_out1, x01, y10);
        /* load 16 point X */
        BBE_LVNX16_XP(x00, px00, L * 2);
        /* load 16 point Y */
        BBE_LVNX16_XP(y00, py00, L*M * 2);
        BBE_MULANX16(z_out0, x00, y00);
      }
      /* save row 16 result matrix [NxN][L] */
      x00 = BBE_PACKVNX40(z_out0, q);
      BBE_SVNX16_IP(x00, pz00, BBE_SIMD_WIDTH * 2);                            
      /* save row 16 result matrix [NxN][L] */
      x01 = BBE_PACKVNX40(z_out1, q);
      BBE_SVNX16_IP(x01, pz00, BBE_SIMD_WIDTH * 2);                            
    }

    for (j = 0; j<L; j += 16)
    {
      z_out0 = 0;

      py00 = (xb_vecNx16*)(offset_y + (uintptr_t)y);
      px00 = (xb_vecNx16*)(offset_x + (uintptr_t)x);

      for (k = 0; k<N; k++)
      {
        /* load 16 point X */
        BBE_LVNX16_XP(x00, px00, L * 2);
        /* load 16 point Y */
        BBE_LVNX16_XP(y00, py00, L*M * 2);
        BBE_MULANX16(z_out0, x00, y00);
      }

      x00 = BBE_PACKVNX40(z_out0, q);
      /* save row 16 result matrix [NxN][L] */
      BBE_SVNX16_IP(x00, pz00, BBE_SIMD_WIDTH * 2);                            

      offset_y = BBE_ADDMOD16U(offset_y, mod_y);
      offset_x = BBE_ADDMOD16U(offset_x, mod_x);
    }

    x += N*L;
  }
} /* matmulnxms_1() */

void matmulnxms_2(int16_t * restrict z,
            const int16_t * restrict x,
            const int16_t * restrict y,
            int N, int M, int L, int Q)
{
  const int M_2 = M / 2;
  int i, j, k;

  int offset_y = 0;
  int offset_x = 0;
  const uint32_t mod_y_2 = (BBE_SIMD_WIDTH * 2) + ((L*M_2 * 2) << 16);
  const uint32_t mod_x = (BBE_SIMD_WIDTH * 2) + ((L * 2) << 16);

  vsaN q = BBE_MOVVSA32(Q);

  xb_vecNx16 x00;
  xb_vecNx16 x10;
  xb_vecNx16 y00;
  xb_vecNx16 y01;
  xb_vecNx40 z_out0;
  xb_vecNx40 z_out10, z_out01, z_out11;

  xb_vecNx16 * restrict px00;
  xb_vecNx16 * restrict py00;
  xb_vecNx16 * restrict pz00;

  px00 = (xb_vecNx16*)x;
  py00 = (xb_vecNx16*)y;
  pz00 = (xb_vecNx16*)z;

  /* check restrictions */
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(L%BBE_SIMD_WIDTH == 0);
  NASSERT(Q >= 0 && Q <= 16);

  pz00 = (xb_vecNx16*)z;
  px00 = (xb_vecNx16*)x;

  if (N == 1)
  {
    for (i = 0; i<M; i += 2)
    {
      __Pragma("ymemory( px00 )");
      for (j = 0; j<M*L; j += 32)
      {
        py00 = (xb_vecNx16*)(offset_y + (uintptr_t)y);
        px00 = (xb_vecNx16*)(offset_x + (uintptr_t)x);

        x10 = BBE_LVNX16_X(px00, M*L*N);
        y01 = BBE_LVNX16_X(py00, M*L);

        BBE_LVNX16_XP(x00, px00, L * 2);
        BBE_LVNX16_XP(y00, py00, L*M * 2);

        z_out0 = BBE_MULNX16(x00, y00);
        z_out10 = BBE_MULNX16(x10, y00);

        z_out01 = BBE_MULNX16(x00, y01);
        z_out11 = BBE_MULNX16(x10, y01);

        x00 = BBE_PACKVNX40(z_out10, q);
        /* save row 16 result matrix [NxN][L] */
        BBE_SVNX16_X(x00, pz00, M*L*M);                           

        x00 = BBE_PACKVNX40(z_out01, q);
        /* save row 16 result matrix [NxN][L] */
        BBE_SVNX16_X(x00, pz00, M*L);                            

        x00 = BBE_PACKVNX40(z_out11, q);
        /* save row 16 result matrix [NxN][L] */
        BBE_SVNX16_X(x00, pz00, M*L*(M + 1));                            

        x00 = BBE_PACKVNX40(z_out0, q);
        /* save row 16 result matrix [NxN][L] */
        BBE_SVNX16_IP(x00, pz00, BBE_SIMD_WIDTH * 2);                           

        offset_y = BBE_ADDMOD16U(offset_y, mod_y_2);
        offset_x = BBE_ADDMOD16U(offset_x, mod_x);
      }
      pz00 = (xb_vecNx16 *)XT_ADD(M*L, (int32_t)pz00);
      x = (int16_t *)XT_ADDX2(N*L, (int32_t)x);
    }

    return;
  }

  for (i = 0; i<M; i += 2)
  {
    for (j = 0; j<M*L; j += 32)
    {
      py00 = (xb_vecNx16*)(offset_y + (uintptr_t)y);
      px00 = (xb_vecNx16*)(offset_x + (uintptr_t)x);

      x10 = BBE_LVNX16_X(px00, M*L*N);
      y01 = BBE_LVNX16_X(py00, M*L);

      BBE_LVNX16_XP(x00, px00, L * 2);
      BBE_LVNX16_XP(y00, py00, L*M * 2);

      z_out0 = BBE_MULNX16(x00, y00);
      z_out10 = BBE_MULNX16(x10, y00);

      z_out01 = BBE_MULNX16(x00, y01);
      z_out11 = BBE_MULNX16(x10, y01);

      __Pragma("ymemory( px00 )");
      __Pragma("ymemory( py00 )");
      __Pragma("loop_count min=1");
      for (k = 0; k<N - 1; k++)
      {
        x10 = BBE_LVNX16_X(px00, M*L*N);
        y01 = BBE_LVNX16_X(py00, M*L);

        BBE_LVNX16_XP(x00, px00, L * 2);
        BBE_LVNX16_XP(y00, py00, L*M * 2);

        BBE_MULANX16(z_out0, x00, y00);
        BBE_MULANX16(z_out10, x10, y00);

        BBE_MULANX16(z_out01, x00, y01);
        BBE_MULANX16(z_out11, x10, y01);
      }

      x00 = BBE_PACKVNX40(z_out10, q);
      /* save row 16 result matrix [NxN][L] */
      BBE_SVNX16_X(x00, pz00, M*L*M);                           

      x00 = BBE_PACKVNX40(z_out01, q);
      /* save row 16 result matrix [NxN][L] */
      BBE_SVNX16_X(x00, pz00, M*L);                            

      x00 = BBE_PACKVNX40(z_out11, q);
      /* save row 16 result matrix [NxN][L] */
      BBE_SVNX16_X(x00, pz00, M*L*(M + 1));                            

      x00 = BBE_PACKVNX40(z_out0, q);
      /* save row 16 result matrix [NxN][L] */
      BBE_SVNX16_IP(x00, pz00, BBE_SIMD_WIDTH * 2);                            

      offset_y = BBE_ADDMOD16U(offset_y, mod_y_2);
      offset_x = BBE_ADDMOD16U(offset_x, mod_x);
    }
    pz00 = (xb_vecNx16 *)XT_ADD(M*L, (int32_t)pz00);
    x = (int16_t *)XT_ADDX2(N*L, (int32_t)x);
  }
} /* matmulnxms_2() */


void matmulnxms_3(int16_t * restrict z,
            const int16_t * restrict x,
            const int16_t * restrict y,
            int N, int M, int L, int Q)
{
  int i, j, k;

  int offset_x = 0;

  const uint32_t mod_x = (BBE_SIMD_WIDTH * 2) + ((L * 2) << 16);

  vsaN q = BBE_MOVVSA32(Q);

  xb_vecNx16 x00, x10, y00, y01;
  xb_vecNx40 z00, z10, z01, z11;

  const xb_vecNx16 *          px;
  const xb_vecNx16 *          py;
  xb_vecNx16 * restrict pz;

  /* check restrictions */
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(L%BBE_SIMD_WIDTH == 0);
  NASSERT(Q >= 0 && Q <= 16);

  py = (xb_vecNx16*)y;
  pz = (xb_vecNx16*)z;

  for (i = 0; i<M; i += 2)
  {
    __Pragma("loop_count min=1");
    for (j = 0; j<M*L; j += 2 * BBE_SIMD_WIDTH)
    {
      px = (xb_vecNx16*)(offset_x + (uintptr_t)x);

      x10 = BBE_LVNX16_X(px, M*L*N);
      y01 = BBE_LVNX16_X(py, M*L);

      BBE_LVNX16_XP(x00, px, L * 2);
      BBE_LVNX16_XP(y00, py, L*M * 2);

      z00 = BBE_MULNX16(x00, y00);
      z10 = BBE_MULNX16(x10, y00);

      z01 = BBE_MULNX16(x00, y01);
      z11 = BBE_MULNX16(x10, y01);

      __Pragma("loop_count min=1");
      for (k = 1; k<N; k++)
      {
        x10 = BBE_LVNX16_X(px, M*L*N);
        y01 = BBE_LVNX16_X(py, M*L);

        BBE_LVNX16_XP(x00, px, L * 2);
        BBE_LVNX16_XP(y00, py, L*M * 2);

        BBE_MULANX16(z00, x00, y00);
        BBE_MULANX16(z10, x10, y00);

        BBE_MULANX16(z01, x00, y01);
        BBE_MULANX16(z11, x10, y01);
      }

      x00 = BBE_PACKVNX40(z10, q);
      BBE_SVNX16_X(x00, pz, M*L*M);

      x00 = BBE_PACKVNX40(z01, q);
      BBE_SVNX16_X(x00, pz, M*L);

      x00 = BBE_PACKVNX40(z11, q);
      BBE_SVNX16_X(x00, pz, M*L*(M + 1));

      x00 = BBE_PACKVNX40(z00, q);
      BBE_SVNX16_IP(x00, pz, BBE_SIMD_WIDTH * 2);

      py = (xb_vecNx16*)((uintptr_t)py - N*M*L * 2 + 2 * BBE_SIMD_WIDTH);

      offset_x = BBE_ADDMOD16U(offset_x, mod_x);
    }

    pz = (xb_vecNx16 *)XT_ADD((int32_t)pz, M*L);
    py = (xb_vecNx16 *)XT_SUB((int32_t)py, M*L);

    x = (int16_t *)XT_ADDX2(N*L, (int32_t)x);
  }
} /* matmulnxms_3() */
