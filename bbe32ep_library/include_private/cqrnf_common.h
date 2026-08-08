/* ------------------------------------------------------------------------ */
/* Copyright (c) 2017 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs (“Cadence    */
/* Libraries”) are the copyrighted works of Cadence Design Systems Inc.     */
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
  NatureDSP_Baseband Library API
  Matrix Decomposition and Inversion Functions
  QR decomposition, floating point, complex data, block format
  common functions
  IntegrIT, 2006-2017
*/
#ifndef CQRNF_COMMON_H__
#define CQRNF_COMMON_H__

typedef struct
{
    void (*takeColumn)(float32_t* x,float32_t* xt,const float32_t* A,int M,int N,int SA,int L);
    void (*housholder)(void* pScr,
                             float32_t* v,
                             float32_t* Fi,
                             float32_t* D,
                       const float32_t* x, 
                       const float32_t* xt, 
                       int M, int SD, int L);
    void (*updateR)   (       float32_t* Z,
                              float32_t* R,
                        const float32_t* v,
                        int SA, int M,int N, int N0, int L);
}
tCqrnIteration;

/*------------------------------------------------
    rotate complex B[L][SB] by complex diagonal matrix Fi'[L][SV]
    Input:
    B[L][SB][2]
    Fi'[L][SV][2]
    Output:
    B[L][SB][2]
    Temporary:
    pScr - 2N vectors

    Note, Fi might be non-aligned!
------------------------------------------------*/
void cqrnfRotateB1(void *pScr,float32_t* B,const float32_t* Fi,int N,int SB,int L);

/*-------------------------------------------------------
    take column from sequence of block ordered matrices
    and put it to the linear array
    M==1
    Input:
    A[L][SA]    matrices
    output:
    x[L*M]      contingious array
-------------------------------------------------------*/
void cqrnfTakeColumn1(float32_t* restrict x,const float32_t* restrict A,int SA,int L);
/*-------------------------------------------------------
    take column from sequence of block ordered matrices
    and put it to the linear array
    Input:
    A[L][SA]    matrices
    output:
    x[L*M]      contingious array
-------------------------------------------------------*/
// M<=4     
void cqrnfTakeColumn4 (float32_t* restrict x,float32_t* restrict xt, const float32_t* restrict A,int M,int N,int SA,int L);
// M=5..8   
void cqrnfTakeColumn8 (float32_t* restrict x,float32_t* restrict xt, const float32_t* restrict A,int M,int N,int SA,int L);
// M=9...12 
void cqrnfTakeColumn12(float32_t* restrict x,float32_t* restrict xt, const float32_t* restrict A,int M,int N,int SA,int L);
// M=13...16
void cqrnfTakeColumn16(float32_t* restrict x,float32_t* restrict xt, const float32_t* restrict A,int M,int N,int SA,int L);
// generic  
void cqrnfTakeColumn  (float32_t* restrict x,float32_t* restrict xt, const float32_t* restrict A,int M,int N,int SA,int L);

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
void cqrnfUpdateR(            float32_t* restrict Z,
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int SA, int M,int N, int N0, int L);
// N<=4
void cqrnfUpdateR4(           float32_t* restrict Z,
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int SA, int M,int N, int N0, int L);
// N=5..8
void cqrnfUpdateR8(           float32_t* restrict Z,
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int SA, int M,int N, int N0, int L);
// N=9...12
void cqrnfUpdateR12(           float32_t* restrict Z,
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int SA, int M,int N, int N0, int L);
// N=13...16
void cqrnfUpdateR16(           float32_t* restrict Z,
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int SA, int M,int N, int N0, int L);

/*-------------------------------------------------------
    partial update of R matrix
    Fi[L][SV][2] diagonal rotation matrix (only 0-th element filled)
    v[L][SV][2]  Housholder vector (M elements filled)
    SV,SD        strides for V/Fi and D
    M,N          matrix size
    L            number of matrices
    Input/output:
    R[L][SA]    L matrices (MxN columns updated with stride N0)
    special case: N=M==1
-------------------------------------------------------*/
void cqrnfUpdateR1(       float32_t* restrict R,
                    const float32_t* restrict v,
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
void cqrnfHousholder1(       float32_t* restrict v,
                             float32_t* restrict Fi,
                             float32_t* restrict D,
                       const float32_t* restrict x, 
                       int SD, int L);


/*-------------------------------------------------------
    find Housholder vectors (V and Fi), diagonal element D
    Input:
    x[L*M]    L input columnar vectors of length M
    xt[L*M]   transposed x
    SV,SD     strides for V/Fi and D
    M         vector length
    Output:
    D[L][SD]  reciprocals of main diagonal (only 0-th element filled)
    Fi[L]     diagonal rotation matrix
    V[M][L]   Housholder vectors (M elements each)
    temporary:
    pScr[]    scratch, defined by cqrnHousholder_getScratchSz()
-------------------------------------------------------*/
void cqrnfHousholder(void* pScr,
                                 float32_t*  v,
                                 float32_t*  Fi,
                                 float32_t*  D,
                           const float32_t*  x, 
                           const float32_t*  xt, 
                           int M, int SD, int L);
//M<=4
void cqrnfHousholder4(void* pScr,
                                 float32_t*  v,
                                 float32_t*  Fi,
                                 float32_t*  D,
                           const float32_t*  x, 
                           const float32_t*  xt, 
                           int M, int SD, int L);
//M<=8
void cqrnfHousholder8(void* pScr,
                                 float32_t*  v,
                                 float32_t*  Fi,
                                 float32_t*  D,
                           const float32_t*  x, 
                           const float32_t*  xt, 
                           int M, int SD, int L);
//M<=16
void cqrnfHousholder16(void* pScr,
                                 float32_t*  v,
                                 float32_t*  Fi,
                                 float32_t*  D,
                           const float32_t*  x, 
                           const float32_t*  xt, 
                           int M, int SD, int L);
/*-------------------------------------------------------
    return scratch size of cqrnfHousholderxxx
    Input:
    M,L  matrix size
-------------------------------------------------------*/
size_t cqrnfHousholder_getScratchSize(int M, int L);

/*-------------------------------------------------------
    rotate R[L][SA] by diagonal matrix Fi'[L][SV]
    Input:
    Fi[L][SV] diagonal rotation matrix (N elements per matrix)
    N         number of columns in R
    Input/output:
    R[L][SA]  sequence of upper-triangle matrices of size MxN. 
              Note: we may rotate only NxN elements because 
              lower (M-N)xN elements in upper trinagle matrix 
              are zeroed!
-------------------------------------------------------*/
void cqrnfrotateR  (float32_t* restrict R,const float32_t* restrict Fi,int N,int SA,int L);
//N=8
void cqrnfrotateR8 (float32_t* restrict R,const float32_t* restrict Fi,int L);
//N=16
void cqrnfrotateR16(float32_t* restrict R,const float32_t* restrict Fi,int L);

#endif
