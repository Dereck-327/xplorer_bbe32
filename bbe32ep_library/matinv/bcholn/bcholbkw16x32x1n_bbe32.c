/* ------------------------------------------------------------------------ */
/* Copyright (c) 2017 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
  NatureDSP_Baseband library. Banded Cholesky backward recursion for pseudo-inversion API (complex data)
    These functions make backward recursion stage of pseudo-inversion. They use
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "bcholn_common.h"

/*-------------------------------------------------------------------------
These functions make backward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices and results of forward recursion. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. Specifically, matrix sizes SR,SD,SY,SX are selected as usual for 
complex block ordered matrix sequencies, i.e. total size is rounded up to 
the closest bigger multiple of BBE_SIMD_WIDTH/2==8 elements. 
SR=size(W*N)
SD=size(N)
SY=size(N*P)
SX=size(N*P)

Input:
W             Band width
N             Matrix dimension (number of columns in matrices R)
P             Number of columns in right-side matrices B
L             Number of matrices
Rt[L][SR]     Cholesky upper triangular matrices R represented in the compact
              form (saved only elements on the main diagonal and above in 
              such a way that diagoanal elements are in the last raw)
D[L][SD]      Sequence of L reciprocals of main diagonal R represented in the  
              block floating point (mantissa and exponent). N' is computed as 
              for complex block ordered matrices of size N
Yt[L][SY]     Results of forward recursion stage. SY is computed as for complex 
              block ordered matrices of size N*P
qA,qX,qY      Fixed point representation of matrices A(or R which is the same), 
              x and y
Output:
Xt[L][SX]     Decision matrix x. SX is computed as for complex block ordered 
              matrices of size N*P

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. For banded matrices, W must be 4, 8, 12 or 16, N must be a positive multiple of 4
3. P>=1
---------------------------------------------------------------------------*/
#if !HAVE_BCHOLN
DISCARD_FUN(void,bcholbkw16x32x1n ,(
                  complex_fract16* restrict Xt, 
            const complex_fract16* restrict Rt, 
            const complex_fract16* restrict D, 
            const complex_fract16* restrict Yt, 
            int qA, int qY, int qX,
            int L))
#else
void bcholbkw16x32x1n (complex_fract16* restrict Xt, const complex_fract16* restrict Rt, const complex_fract16* restrict D, const complex_fract16* restrict Yt, int qA, int qY, int qX,int L)
{
    if (L <= 0) return;
    bcholbkw16xnx1n ((int16_t*)Xt, (const int16_t*)Rt,(const int16_t*)D, (const int16_t*)Yt,qX-qY+qA, 32,L);
} /* bcholbkw16x32x1n() */
#endif
