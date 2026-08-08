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
/*          Copyright (C) 2009-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
  NatureDSP_Baseband Library API
  Matrix Decomposition and Inversion Functions
  QR decomposition, floating point, complex data, stream format
  C code optimized for BBE32EP with VFPU
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include "cqrsf_common.h"

#if HAVE_VFPU
/*-------------------------------------------------------------------------
[c]qr_calc_qbMxNxPsf

These functions apply Householder reflections to L MxP matrices B in the
course of solving a set of complex-valued linear problems A*X=B through
the QR decomposition of matrices A: A*X=B, A=Q*R => Q*R*X=B => R*X=Q'*B.
Instead of direct multiplication of each matrix B by conjugate transpose
of the corresponding matrix Q, we use a set of vectors V to perform a
sequence of Householder  reflections (see QR decomposition routines
[c]qr_build_rMxNsf.

Data format: IEEE-754 Std single precision floating-point

Data transform is performed in-place.

NOTE:
Data layout for matrices is selected as for other matrices written in a stream 
order. 

Input:
B[M*P]L]                Matrices B (L matrices of size MxP)
V[((2*M-N+1)*N/2+N)*L]  L sets of Householder vectors
Output:
B[M*P][L]               Matrices Q'*B (L matrices of size MxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 4 for complex data and 
   8 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,P,L)
4. M must be greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
---------------------------------------------------------------------------*/
void cqr_calc_qb2x2x1sf  (void *pScr, complex_float *_B, const complex_float *_V , int L)
{
    float32_t *B=(float32_t *)_B;
    float32_t *V=(float32_t *)_V;
    NASSERT_ALIGN(pScr,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(B,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(V,(2*BBE_SIMD_WIDTH));
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/4)==0);
    cqrsfUpdateB2(B+2*0*L,V,L);   V+=2*(2-0)*L;
    cqrsfUpdateB1(B+2*1*L,V,L);   V+=2*(2-1)*L;
    cqrsfRotateBconj1(B,V,2,L);
}
size_t cqr_calc_qb2x2x1sf_getScratchSize(int M, int P, int L) 
{ 
    (void)M,(void)P,(void)L;
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/4)==0);
    return 0; 
}
#else
DISCARD_FUN(void,cqr_calc_qb2x2x1sf, (void *pScr, complex_float *B, const complex_float *V , int L))
size_t cqr_calc_qb2x2x1sf_getScratchSize(int M, int P, int L) 
{ 
    (void)M,(void)P,(void)L;
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/4)==0);
    return 0; 
}
#endif
