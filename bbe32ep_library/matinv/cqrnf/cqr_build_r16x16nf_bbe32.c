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
  QR decomposition, floating point, complex data, block format
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "cqrnf_common.h"
#include "common.h"

#if HAVE_VFPU

/*-------------------------------------------------------------------------
Make QR decomposition for block ordered matrices.
Matrix sizes SA,SV are selected as usual for block ordered matrix 
sequencies of corresponding type, i.e. total size is rounded up to the 
closest bigger multiple of 
- BBE_SIMD_WIDTH/2==8 elements for float32_t
- BBE_SIMD_WIDTH/4==4 elements for complex_float
or, if it is less, to the closest bigger 
multiple of degree of 2. 
SA=size(M*N)
SV=size(((2*M-N+1)*N/2+N)*L)
SD=size(N)
Scratch size in bytes is defined by functions xxxx_getScratchSize()

Input:
 M, N         Dimensional parameters
 L            Number of matrices
Input/output:
 A[L][SA]     On input it is the sequence of L matrices A. 
              At the end of the process, matrices R replace input
              matrices A. In a case of non-square matrices (N!=M),
              only N*N elements of each output matrix are valid.
Output:
 V[SV]        Sequence of L Housholder rotation vectors 
 D[L*SD]      Reciprocals of main diagonal in a special format

Restrictions:
1. A, V, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. Number of columns for input matrices A must not exceed the number
   of rows: N <= M.
---------------------------------------------------------------------------*/
void  cqr_build_r16x16nf(void *pScr,complex_float* _A,complex_float* _V,complex_float* _D,int L)
{
    static const tCqrnIteration it[]=
    {
        {cqrnfTakeColumn16, cqrnfHousholder16, cqrnfUpdateR16},
        {cqrnfTakeColumn16, cqrnfHousholder16, cqrnfUpdateR16},
        {cqrnfTakeColumn16, cqrnfHousholder16, cqrnfUpdateR16},
        {cqrnfTakeColumn16, cqrnfHousholder16, cqrnfUpdateR16},
        {cqrnfTakeColumn12, cqrnfHousholder16, cqrnfUpdateR12},
        {cqrnfTakeColumn12, cqrnfHousholder16, cqrnfUpdateR12},
        {cqrnfTakeColumn12, cqrnfHousholder16, cqrnfUpdateR12},
        {cqrnfTakeColumn12, cqrnfHousholder16, cqrnfUpdateR12},
        {cqrnfTakeColumn8, cqrnfHousholder8, cqrnfUpdateR8},
        {cqrnfTakeColumn8, cqrnfHousholder8, cqrnfUpdateR8},
        {cqrnfTakeColumn8, cqrnfHousholder8, cqrnfUpdateR8},
        {cqrnfTakeColumn8, cqrnfHousholder8, cqrnfUpdateR8},
        {cqrnfTakeColumn4, cqrnfHousholder4, cqrnfUpdateR4},
        {cqrnfTakeColumn4, cqrnfHousholder4, cqrnfUpdateR4},
        {cqrnfTakeColumn4, cqrnfHousholder4, cqrnfUpdateR4}
    };
    float32_t* A=(float32_t*)_A;
    float32_t* V=(float32_t*)_V;
    float32_t* D=(float32_t*)_D;
    size_t Hsize;
    int m;
    int SA=(16*16<<1);
    int SD=(16   <<1);
    float32_t* Z;
    float32_t* Zt;
    float32_t* Fi;
    float32_t* pV;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    if (L<=0) return;
    Hsize=cqrnfHousholder_getScratchSize(16,L);
    Hsize=(Hsize+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    Z=(float32_t*)(((uintptr_t)pScr)+Hsize);
    Zt=Z+2*16*L;
    Fi=V+(2*16-16+1)*16*L;
    pV=V;
    for(m=0; m<16-1; m++)
    {
        it[m].takeColumn (Z,Zt,A+2*m*(16+1),16-m,16,SA,L);
        it[m].housholder (pScr,pV,Fi+2*m*L,D+2*m,Z,Zt,16-m,SD,L);
        it[m].updateR (Z,A+2*m*(16+1),pV,SA,16-m,16-m,16,L);
        pV+=2*(16-m)*L;
    }
    cqrnfTakeColumn1(Z,A+2*m*(16+1),SA,L);
    cqrnfHousholder1(pV,Fi+2*m*L,D+2*m,Z,SD,L);
    cqrnfUpdateR1   (A+2*m*(16+1),pV,SA,L);

    // final rotation Fi'*R
    cqrnfrotateR16(A,Fi,L);
}
size_t  cqr_build_r16x16nf_getScratchSize(int M, int N,int L)
{
    size_t Zsize,Hsize;
    NASSERT(N==16 && M==16);
    NASSERT(L>0);
    L=XT_MAX(L,0);
    Hsize=cqrnfHousholder_getScratchSize(M,L);
    Hsize=(Hsize+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    Zsize=M*L*sizeof(complex_float);
    Zsize=(Zsize+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    return 2*Zsize+Hsize;
}
#else
DISCARD_FUN(void, cqr_build_r16x16nf, (void *pScr,complex_float* _A,complex_float* _V,complex_float* _D,int L))
size_t  cqr_build_r16x16nf_getScratchSize(int M, int N, int L)
{
    NASSERT(N==16 && M==16);
    NASSERT(L>0);
    L=XT_MAX(L,0);
    return 0;
}
#endif
