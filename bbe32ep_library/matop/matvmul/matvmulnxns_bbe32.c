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

/* Streaming Order, NxN*Nx1->Nx1, Sx=NxN, Sy=N, Sz=N
   Restrictions:
     L must be a multiple of 16
*/
void matvmulnxns ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int N, int L, int Q )
{
  const int M = N;

  int i, j, l;
  vsaN  q = BBE_MOVVSA32(Q);
  xb_vecNx16 * restrict px;
  xb_vecNx16 * restrict px0;
  xb_vecNx16 * restrict py;
  xb_vecNx16 * restrict pz0 = (xb_vecNx16 *)z;
  xb_vecNx16 x0, y0, z0;
  xb_vecNx40 acc0;

  /* check restrictions */
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  NASSERT(L%BBE_SIMD_WIDTH == 0);

  if (L <= 0) return;
  if ((M & 3) == 0)
  {
    xb_vecNx16 x1, x2, x3;
    xb_vecNx16 z1, z2, z3;
    xb_vecNx40 acc1, acc2, acc3;
    for (i = 0; i<M / 4; i++)
    {
      px0 = (xb_vecNx16*)(x + i*N*L);
      py = (xb_vecNx16*)(y);

      for (l = 0; l<L; l += BBE_SIMD_WIDTH)
      {
        x2 = BBE_LVNX16_X(px0, M / 2 * N*L*sizeof(int16_t));
        x1 = BBE_LVNX16_X(px0, M / 4 * N*L*sizeof(int16_t));
        x3 = BBE_LVNX16_X(px0, 3 * M / 4 * N*L*sizeof(int16_t));
        BBE_LVNX16_XP(x0, px0, L*sizeof(int16_t));
        BBE_LVNX16_XP(y0, py, L*sizeof(int16_t));

        acc0 = BBE_MULRNX16(x0, y0, q);
        acc2 = BBE_MULRNX16(x2, y0, q);
        acc1 = BBE_MULRNX16(x1, y0, q);
        acc3 = BBE_MULRNX16(x3, y0, q);

        __Pragma("ymemory( px0 )");
        __Pragma("loop_count min=1");
        for (j = 0; j<N - 2; j += 2)
        {
          x2 = BBE_LVNX16_X(px0, M / 2 * N*L*sizeof(int16_t));
          x1 = BBE_LVNX16_X(px0, M / 4 * N*L*sizeof(int16_t));
          x3 = BBE_LVNX16_X(px0, 3 * M / 4 * N*L*sizeof(int16_t));

          BBE_LVNX16_XP(x0, px0, L*sizeof(int16_t));
          BBE_LVNX16_XP(y0, py, L*sizeof(int16_t));

          BBE_MULANX16(acc0, x0, y0);
          BBE_MULANX16(acc2, x2, y0);

          BBE_MULANX16(acc1, x1, y0);
          BBE_MULANX16(acc3, x3, y0);

          x2 = BBE_LVNX16_X(px0, M / 2 * N*L*sizeof(int16_t));
          x1 = BBE_LVNX16_X(px0, M / 4 * N*L*sizeof(int16_t));
          x3 = BBE_LVNX16_X(px0, 3 * M / 4 * N*L*sizeof(int16_t));

          BBE_LVNX16_XP(x0, px0, L*sizeof(int16_t));
          BBE_LVNX16_XP(y0, py, L*sizeof(int16_t));

          BBE_MULANX16(acc0, x0, y0);
          BBE_MULANX16(acc2, x2, y0);

          BBE_MULANX16(acc1, x1, y0);
          BBE_MULANX16(acc3, x3, y0);
        }

        x2 = BBE_LVNX16_X(px0, M / 2 * N*L*sizeof(int16_t));
        x1 = BBE_LVNX16_X(px0, M / 4 * N*L*sizeof(int16_t));
        x3 = BBE_LVNX16_X(px0, 3 * M / 4 * N*L*sizeof(int16_t));

        BBE_LVNX16_XP(x0, px0, (BBE_SIMD_WIDTH - (N - 1)*L)*sizeof(int16_t));
        BBE_LVNX16_XP(y0, py, (BBE_SIMD_WIDTH - (N - 1)*L)*sizeof(int16_t));

        BBE_MULANX16(acc0, x0, y0);
        BBE_MULANX16(acc2, x2, y0);
        BBE_MULANX16(acc1, x1, y0);
        BBE_MULANX16(acc3, x3, y0);

        z1 = BBE_PACKVNX40(acc1, q);
        z3 = BBE_PACKVNX40(acc3, q);
        z2 = BBE_PACKVNX40(acc2, q);
        z0 = BBE_PACKVNX40(acc0, q);

        BBE_SVNX16_X(z1, pz0, M / 4 * L*sizeof(int16_t));
        BBE_SVNX16_X(z3, pz0, 3 * M / 4 * L*sizeof(int16_t));

        BBE_SVNX16_X(z2, pz0, M / 2 * L*sizeof(int16_t));
        BBE_SVNX16_IP(z0, pz0, sizeof(*pz0));
      }
    }
  }
  else
  {
    for (i = 0; i<M; i++)
    {
      for (l = 0; l<L; l += BBE_SIMD_WIDTH)
      {
        acc0 = BBE_ZERONX40();
        px = (xb_vecNx16*)(x + l + i*N*L);
        py = (xb_vecNx16*)(y + l);

        __Pragma("loop_count min=1");
        for (j = 0; j<N; j++)
        {
          BBE_LVNX16_XP(x0, px, L*sizeof(int16_t));
          BBE_LVNX16_XP(y0, py, L*sizeof(int16_t));
          BBE_MULANX16(acc0, x0, y0);
        }

        z0 = BBE_PACKVNX40(acc0, q);
        BBE_SVNX16_IP(z0, pz0, sizeof(*pz0));
      }
    }
  }
} /* matvmulnxns() */
