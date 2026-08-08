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
  NatureDSP_Baseband Library API
  Matrix Decomposition and Inversion Functions
  QR decomposition, floating point, complex data, stream format
  common data and definitions
  IntegrIT, 2006-2017
*/
#ifndef CQRSF_COMMON_H__
#define CQRSF_COMMON_H__

/*-----------------------------------------------------------
    find Householder vectors (V and Fi), diagonal element D
    Input:
    x[]       pointer to the diagonal element of original matrix A given 
              in streaming order
    NL,SV,SD  strides for x, V/Fi and D
    M         vector length
    Output:
    D[L][SD]  reciprocals of main diagonal (only 0-th element filled)
    Fi[L]     diagonal rotation matrix
    V[M][L]   Householder vectors (M elements each)
-----------------------------------------------------------*/
void cqrfsHousholder(void *pScr,
                            float32_t* restrict v,
                            float32_t* restrict Fi,
                            float32_t *restrict D,
                    const float32_t* restrict x, 
                    int SV, int M, int N, int L);
// M==1
void cqrfsHousholder1(    float32_t* restrict v,
                                 float32_t* restrict Fi,
                                 float32_t *restrict D,
                           const float32_t* restrict x, 
                           int SV, int M, int N, int L);
// M==2
void cqrfsHousholder2(void *pScr,
                            float32_t* restrict v,
                            float32_t* restrict Fi,
                            float32_t *restrict D,
                    const float32_t* restrict x, 
                    int SV, int M, int N, int L);
// M==3
void cqrfsHousholder3(void *pScr,
                            float32_t* restrict v,
                            float32_t* restrict Fi,
                            float32_t *restrict D,
                    const float32_t* restrict x, 
                    int SV, int M, int N, int L);
// M==4
void cqrfsHousholder4(void *pScr,
                            float32_t* restrict v,
                            float32_t* restrict Fi,
                            float32_t *restrict D,
                    const float32_t* restrict x, 
                    int SV, int M, int N, int L);
size_t cqrfsHousholder_getScratchSize(int M, int N, int L);

/*-----------------------------------------------------------
    partial update of R matrix - update begins from 
    v         pointer to m-th Householder vector 
    K         current number of column in the matrix R
    M,N       matrix size
    L         number of matrices
    Input/output:
    R         pointer to the beginning of column  
-----------------------------------------------------------*/
void cqrfsUpdateR(
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int K,int M,int N, int L);

// M==1
void cqrfsUpdateR1(
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int K,int M,int N, int L);
// M==2
void cqrfsUpdateR2(
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int K,int M,int N, int L);
// M==3
void cqrfsUpdateR3(
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int K,int M,int N, int L);
// M==4
void cqrfsUpdateR4(
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int K,int M,int N, int L);

/*-----------------------------------------------------------
    rotate R[L][MxN] by diagonal matrix Fi'[L][N]
    Input:
    Fi[L][N]  diagonal rotation matrix (N elements per matrix)
    N         number of columns in R
    Input/output:
    R[L][MxN] sequence of upper-triangle matrices of size MxN. 
              Note: we may rotate only NxN elements because 
              lower (M-N)xN elements in upper triangle matrix 
              remain unchanged
-----------------------------------------------------------*/
void cqrfsRotateR(float32_t* restrict R,const float32_t* restrict Fi,int M, int N,int L);
//N==2
void cqrfsRotateR2(float32_t* restrict R,const float32_t* restrict Fi,int M, int N,int L);
//N==3
void cqrfsRotateR3(float32_t* restrict R,const float32_t* restrict Fi,int M, int N,int L);
//N==4
void cqrfsRotateR4(float32_t* restrict R,const float32_t* restrict Fi,int M, int N,int L);

/*
    rotate B[L][NxP] by diagonal matrix Fi'[L][N]
    P==1
*/
void cqrsfRotateBconj1(float32_t* B,const float32_t* Fi,int N,int L);

/*---------------------------------------------------------
    update matrix B[L][MxP] by housholder vectors V[L][SV]

    Input:
    M,N,P,L     dimensions
    V[L][SV]    Householder vectors
    Input/output:
    B[L][MxP]   B matrices MxP
  
---------------------------------------------------------*/
void cqrsfUpdateB(float32_t* B,const float32_t* V,int M,int P,int L);
// M=1,P=1
void cqrsfUpdateB1(float32_t* B,const float32_t* V,int L);
// M=2,P=1
void cqrsfUpdateB2(float32_t* B,const float32_t* V,int L);
// M=3,P=1
void cqrsfUpdateB3(float32_t* B,const float32_t* V,int L);
// M=4,P=1
void cqrsfUpdateB4(float32_t* B,const float32_t* V,int L);

#endif
