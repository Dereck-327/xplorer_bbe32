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
    Matrix Hermitian Product
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"
#define ITER(a1,a0,y1,y0,rep,py)              \
{                                             \
    xb_vecNx16 x0;                            \
    BBE_LVNX16_IP(y0, py,2*BBE_SIMD_WIDTH);   \
    BBE_LVNX16_IP(y1, py,2*BBE_SIMD_WIDTH);   \
    x0=BBE_SELNX16(y1,y0,rep);                \
    BBE_MULANX16J(a0,y0,x0);                  \
    BBE_MULANX16J(a1,y1,x0);                  \
}

/*-------------------------------------------------------------------------
Matrix Hermitian Product

Description: These functions left multiply each complex input matrix by its
conjugate transpose. The result is a Hermitian (or self-adjoint) matrix. Both
the block order and streaming order are allowed for input/output matrix
sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr      Scratch memory area. To determine the scratch area size required by
          a function <fun>, use the respective helper function 
          <fun>_getScratchSize()
Input:
x[L*Sx]   Complex input matrices
M         Matrix dimension 
N         Matrix dimension (columnar for MxN)
L         Number of matrices 
Q         Position of fractional point in matrix representation, 0..16
Output:
y[L*Sy]   Complex output matrices

Restrictions:
pScr,x,y  Aligned on 32-byte boundary
pScr,x,y  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, 16x16*16x16->16x16, Sx=256, Sy=256
   Restrictions:
     None
*/
void cmatherm16x16n ( complex_fract16 * restrict y, 
                const complex_fract16 * restrict x, 
                int L, int Q )
{
  int k, i, ystride; 
	vsaN q;

	const xb_vecNx16 * restrict py = (const xb_vecNx16 *)x;  
	xb_vecNx16 * restrict pz = (xb_vecNx16 *)y; 
	xb_vecNx16  x0, y0,y1,x1,y2,y3,x2,y4,y5,x3,y6,y7, z0,_rep; 
	xb_vecNx40 a0,a1; 
  vselN rep0,rep1;

  NASSERT_ALIGN(x, BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, BBE_SIMD_WIDTH);
  if (L <= 0) return;

	q = BBE_MOVVSA32(Q);
  _rep=BBE_MOVVA16C(1<<16); /* pair 1,0 */
  rep1=rep0=BBE_MOVVSELNX16(_rep,0);
  /* prologue */
  BBE_LVNX16_IP(y0, py,   2*BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(y1, py,   2*BBE_SIMD_WIDTH); 
  x0=BBE_SELNX16(y1,y0,rep1);
  BBE_LVNX16_IP(y2, py,   2*BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(y3, py,   2*BBE_SIMD_WIDTH); 
  x1=BBE_SELNX16(y3,y2,rep1);
  BBE_LVNX16_IP(y4, py,   2*BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(y5, py,   2*BBE_SIMD_WIDTH); 
  x2=BBE_SELNX16(y5,y4,rep1);
  BBE_LVNX16_IP(y6, py,   2*BBE_SIMD_WIDTH);
  BBE_LVNX16_IP(y7, py,   2*BBE_SIMD_WIDTH); 
  BBE_SELUNX16(x3,y7,y6,rep1,2);
  __Pragma( "loop_count min=15" ); 
  for(i=k=0; k<L*16-1; k++)
  {
    i=BBE_ADDMOD16U(i,0x100001);
    ystride=-31*2*BBE_SIMD_WIDTH;
    XT_MOVEQZ(ystride,2*BBE_SIMD_WIDTH,i);
    a0 = BBE_MULRNX16J(y0,x0,  q); 
    a1 = BBE_MULRNX16J(y1,x0,  q); 
    BBE_MULANX16J(a0,y2,x1); 
    BBE_MULANX16J(a1,y3,x1); 
    BBE_MULANX16J(a0,y4,x2); 
    BBE_MULANX16J(a1,y5,x2); 
    BBE_MULANX16J(a0,y6,x3); 
    BBE_MULANX16J(a1,y7,x3); 
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    BBE_LVNX16_IP(y0, py,    2*BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(y1, py,    ystride); 
    BBE_SELUNX16(x0,y1,y0,rep0,2);
    BBE_MULANX16J(a0,y0,x0); 
    BBE_MULANX16J(a1,y1,x0); 

    BBE_LVNX16_IP(y0, py,   2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y1, py,   2*BBE_SIMD_WIDTH); 
    x0=BBE_SELNX16(y1,y0,rep1);
    BBE_LVNX16_IP(y2, py,   2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y3, py,   2*BBE_SIMD_WIDTH); 
    x1=BBE_SELNX16(y3,y2,rep1);
    BBE_LVNX16_IP(y4, py,   2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y5, py,   2*BBE_SIMD_WIDTH); 
    x2=BBE_SELNX16(y5,y4,rep1);
    BBE_LVNX16_IP(y6, py,   2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y7, py,   2*BBE_SIMD_WIDTH); 
    BBE_SELUNX16(x3,y7,y6,rep1,2);
    z0 = BBE_PACKVNX40(a0, q); BBE_SVNX16_IP(z0, pz, 2*BBE_SIMD_WIDTH);
    z0 = BBE_PACKVNX40(a1, q); BBE_SVNX16_IP(z0, pz, 2*BBE_SIMD_WIDTH);
  }
  /* epilogue */
  {
    a0 = BBE_MULRNX16J(y0,x0,  q); 
    a1 = BBE_MULRNX16J(y1,x0,  q); 
    BBE_MULANX16J(a0,y2,x1); 
    BBE_MULANX16J(a1,y3,x1); 
    BBE_MULANX16J(a0,y4,x2); 
    BBE_MULANX16J(a1,y5,x2); 
    BBE_MULANX16J(a0,y6,x3); 
    BBE_MULANX16J(a1,y7,x3); 
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    ITER(a1,a0,y1,y0,rep0,py)
    BBE_LVNX16_IP(y0, py,    2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(y1, py,    2*BBE_SIMD_WIDTH); 
    x0=BBE_SELNX16(y1,y0,rep0);
    BBE_MULANX16J(a0,y0,x0); 
    BBE_MULANX16J(a1,y1,x0); 

    z0 = BBE_PACKVNX40(a0, q); BBE_SVNX16_IP(z0, pz, 2*BBE_SIMD_WIDTH);
    z0 = BBE_PACKVNX40(a1, q); BBE_SVNX16_IP(z0, pz, 2*BBE_SIMD_WIDTH);
  }
} /* cmatherm16x16n() */
