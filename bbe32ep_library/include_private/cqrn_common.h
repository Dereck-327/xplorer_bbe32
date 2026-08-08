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
    common definitions/routines for cqrn part
    IntegrIT, 2006-2017
*/
#ifndef CQRN_COMMON_H__
#define CQRN_COMMON_H__

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"

#define HAVE_CQRN (HAVE_RSQRT && HAVE_VSAMATH && HAVE_DIV)

#define EXTRA_SHIFT 0

/*-------------------------------------------------------
    rotate R[L][SA] by diagonal matrix Fi'[L][SV]
    Input:
    Fi[L][SV] diagonal rotation matrix (N elements per matrix)
    N         number of columns in R
    extraShift extra shift right
    Input/output:
    R[L][SA]  sequence of upper-triangle matrices of size MxN. 
              Note: we may rotate only NxN elements because 
              lower (M-N)xN elements in upper trinagle matrix 
              are zeroed!
-------------------------------------------------------*/
void cqrnRotateR(int16_t* restrict R,const int16_t* restrict Fi,int N,int SA,int L,int extraShift);
/* optimized case for specific conditions */
/* M=N == 8  */
void cqrnRotateR8(int16_t* restrict R,const int16_t* restrict Fi,int L,int extraShift);
/* M=N == 16 */
void cqrnRotateR16(int16_t* restrict R,const int16_t* restrict Fi,int L,int extraShift);

/*-------------------------------------------------------
    partial update of R matrix
    Fi[L][SV] diagonal rotation matrix (only 0-th element filled)
    v[L][SV]  Housholder vector (M elements filled)
    SV,SD     strides for V/Fi and D
    M,N       matrix size
    L         number of matrices
    Input/output:
    R[L][SA]    L matrices (MxN columns updated with stride N0)
    Temporary:
    Z[2*N*L]
-------------------------------------------------------*/
void cqrnUpdateR(int16_t* restrict Z,
                 int16_t* restrict R,
                 const int16_t* restrict v,
                 int SA, int M,int N, int N0, int L);
/* optimized cases for specific conditions */
/* N<=8 */
void cqrnUpdateR8(int16_t* restrict Z,
                 int16_t* restrict R,
                 const int16_t* restrict v,
                 int SA, int M,int N, int N0, int L);
/* N<=16 */
void cqrnUpdateR16(int16_t* restrict Z,
                 int16_t* restrict R,
                 const int16_t* restrict v,
                 int SA, int M,int N, int N0, int L);
/* N<=32 */
void cqrnUpdateR32(int16_t* restrict Z,
                 int16_t* restrict R,
                 const int16_t* restrict v,
                 int SA, int M,int N, int N0, int L);

/*-------------------------------------------------------
    partial update of R matrix
    Fi[L][SV] diagonal rotation matrix (only 0-th element filled)
    v[L][SV]  Housholder vector (M elements filled)
    SV,SD     strides for V/Fi and D
    M,N       matrix size
    L         number of matrices
    Input/output:
    R[L][SA]    L matrices (MxN columns updated with stride N0)
    special case: N=M==1
-------------------------------------------------------*/
void cqrnUpdateR1(int16_t* restrict R,
                        const int16_t* restrict v,
                        int SA, int L);

/*-------------------------------------------------------
    find Housholder vectors (V and Fi), diagonal element D
    Input:
    x[L*M]    L input columnar vectors of length M
    SV,SD     strides for V/Fi and D
    M         vector length
    Output:
    D[L][SD]  reciprocals of main diagonal (only 0-th element filled)
    Fi[L]     diagonal rotation matrix
    V[M][L]   Housholder vectors (M elements each)

    special case: M==1
-------------------------------------------------------*/
void cqrnHousholder1      (int16_t* restrict v,
                           int16_t* restrict Fi,
                           int16_t *restrict D,
                           const int16_t* restrict x, 
                           int SD, int L);

/*-------------------------------------------------------
    find Housholder vectors (V and Fi), diagonal element D
    Input:
    x[L*M]    L input columnar vectors of length M
    SV,SD     strides for V/Fi and D
    M         vector length
    Output:
    D[L][SD]  reciprocals of main diagonal (only 0-th element filled)
    Fi[L]     diagonal rotation matrix
    V[M][L]   Housholder vectors (M elements each)
    temporary:
    pScr[]    scratch, defined by cqrnHousholder_getScratchSz()
-------------------------------------------------------*/
void cqrnHousholder(void* pScr,
                    int16_t* restrict v,
                    int16_t* restrict Fi,
                    int16_t *restrict D,
                    const int16_t* restrict x, 
                    int M, int SD, int L);
/* specific variants for some conditions */
/* M<=8 */
void cqrnHousholder8(void* pScr,
                    int16_t* restrict v,
                    int16_t* restrict Fi,
                    int16_t *restrict D,
                    const int16_t* restrict x, 
                    int M, int SD, int L);
/* M<=16 */
void cqrnHousholder16(void* pScr,
                    int16_t* restrict v,
                    int16_t* restrict Fi,
                    int16_t *restrict D,
                    const int16_t* restrict x, 
                    int M, int SD, int L);
/* M<=32 */
void cqrnHousholder32(void* pScr,
                    int16_t* restrict v,
                    int16_t* restrict Fi,
                    int16_t *restrict D,
                    const int16_t* restrict x, 
                    int M, int SD, int L);





size_t cqrnHousholder_getScratchSz(int M,int L);

typedef void (*fnTakeColumn)(int16_t* restrict x,const int16_t* restrict A,int M,int N,int SA,int L);
typedef void (*fnHousholder)(void* pScr,
                    int16_t* restrict v,
                    int16_t* restrict Fi,
                    int16_t *restrict D,
                    const int16_t* restrict x, 
                    int M, int SD, int L);
typedef void (*fnUpdateR)(int16_t* restrict Z,
                 int16_t* restrict R,
                 const int16_t* restrict v,
                 int SA, int M,int N, int N0, int L);

typedef struct
{
    fnTakeColumn takeColumn;
    fnHousholder housholder;
    fnUpdateR    updateR;
}
tCqrn_buildr_Iteration;


/*-------------------------------------------------------
    take column from sequence of block ordered matrices
    and put it to the linear array
    Input:
    A[L][SA]    matrices
    output:
    x[L*M]      contingious array
-------------------------------------------------------*/
void cqrnTakeColumn(int16_t* restrict x,const int16_t* restrict A,int M,int N,int SA,int L);
/* specific optimized cases */
/* N,M==1 */
void cqrnTakeColumn1(int16_t* restrict x,const int16_t* restrict A,int SA,int L);
/* M=1...4 */
void cqrnTakeColumn4(int16_t* restrict x,const int16_t* restrict A,int M,int N,int SA,int L);
/* M=5...8 */
void cqrnTakeColumn8(int16_t* restrict x,const int16_t* restrict A,int M,int N,int SA,int L);
/* M=9...12 */
void cqrnTakeColumn12(int16_t* restrict x,const int16_t* restrict A,int M,int N,int SA,int L);
/* M=13...16 */
void cqrnTakeColumn16(int16_t* restrict x,const int16_t* restrict A,int M,int N,int SA,int L);
/* M=17...32 */
void cqrnTakeColumn32(int16_t* restrict x,const int16_t* restrict A,int M,int N,int SA,int L);


/*------------------------------------------------
    rotate B[L][SB] by diagonal matrix Fi'[L][SV]
    introduces additional shift right !
    Input:
    Fi[L][SV]
    Input/output:
    B[L][SB]
    Temporary:
    pScr    - size in bytes N*2*BBE_SIMD_WIDTH
------------------------------------------------*/
/*
    rotate B[L][SB] by diagonal matrix Fi'[L][SV]
    introduces additional shift right !
*/
void cqrnRotateB(int16_t* B,const int16_t* Fi,int N,int P,int SB,int L);
/* specific cases */
/* P=1 */
void cqrnRotateB1(void* pScr, int16_t* B,const int16_t* Fi,int N,int SB,int L);
/* P=1, N==8 */
void cqrnRotateB8(void* pScr, int16_t* B,const int16_t* Fi,int L);
/* P=1, N==16 */
void cqrnRotateB16(void* pScr, int16_t* B,const int16_t* Fi,int L);
/* P=1, N==32 */
void cqrnRotateB32(void* pScr, int16_t* B,const int16_t* Fi,int L);

/*-------------------------------------------------------
    backward recursion: P==1

   Input:
    M, N, P      dimensional parameters
    L            number of matrices
    qABX         qA-qB+qX
   Input/output:
    x[L][SB][2]  at the input it is the sequence of L updated right parts Z=Q'B.
                 They will be replaced with MMSE solution vectors X (only N*P 
                 elements are used)
   Input:
    R[L][SA][2]  Upper triangle matrices R (only N*N 
                 elements of each matrix are used)
    D[L][SD][2]  reciprocal of main diagonal (mantissa, exponent) 
                 in the special format
-------------------------------------------------------*/
void cqrnBkwnx1(
            int16_t* restrict x, 
            const int16_t* restrict R,
            const int16_t* restrict D,
            int qABX,
            int M,int N,int L);

/*-------------------------------------------------------------------------
  Scale B[SB] by 1 bit left with rounding
  Input/output:
  B[SB]        input/output matrix
---------------------------------------------------------------------------*/
void cqrnScaleB(int16_t* B,int SB);

#endif
