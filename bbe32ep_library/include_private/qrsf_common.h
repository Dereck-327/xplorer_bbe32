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
  QR decomposition, floating point, real data, stream format
  common data and definitions
  IntegrIT, 2006-2017
*/
#ifndef QRSF_COMMON_H__
#define QRSF_COMMON_H__

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
void qrfsHousholder(void *pScr,
                            float32_t* restrict v,
                            float32_t* restrict Fi,
                            float32_t *restrict D,
                    const float32_t* restrict x, 
                    int SV, int M, int N, int L);
// M==1
void qrfsHousholder1(    float32_t* restrict v,
                                 float32_t* restrict Fi,
                                 float32_t *restrict D,
                           const float32_t* restrict x, 
                           int SV, int M, int N, int L);
// M==2
void qrfsHousholder2(void *pScr,
                            float32_t* restrict v,
                            float32_t* restrict Fi,
                            float32_t *restrict D,
                    const float32_t* restrict x, 
                    int SV, int M, int N, int L);
// M==3
void qrfsHousholder3(void *pScr,
                            float32_t* restrict v,
                            float32_t* restrict Fi,
                            float32_t *restrict D,
                    const float32_t* restrict x, 
                    int SV, int M, int N, int L);
// M==4
void qrfsHousholder4(void *pScr,
                            float32_t* restrict v,
                            float32_t* restrict Fi,
                            float32_t *restrict D,
                    const float32_t* restrict x, 
                    int SV, int M, int N, int L);
/*-----------------------------------------------------------
  return scratch size for qrfsHousholderxxx
  Input:
  M,N,L  matrix sizes
-----------------------------------------------------------*/
size_t qrfsHousholder_getScratchSize(int M, int N, int L);

/*-----------------------------------------------------------
    partial update of R matrix - update begins from 
    v         pointer to m-th Householder vector 
    K         current number of column in the matrix R
    M,N       matrix size
    L         number of matrices
    Input/output:
    R         pointer to the beginning of column  
-----------------------------------------------------------*/
void qrfsUpdateR(
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int K,int M,int N, int L);

// M==1
void qrfsUpdateR1(
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int K,int M,int N, int L);
// M==2
void qrfsUpdateR2(
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int K,int M,int N, int L);
// M==3
void qrfsUpdateR3(
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int K,int M,int N, int L);
// M==4
void qrfsUpdateR4(
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
void qrfsRotateR(float32_t* restrict R,const float32_t* restrict Fi,int M, int N,int L);
//N==2
void qrfsRotateR2(float32_t* restrict R,const float32_t* restrict Fi,int M, int N,int L);
//N==3
void qrfsRotateR3(float32_t* restrict R,const float32_t* restrict Fi,int M, int N,int L);
//N==4
void qrfsRotateR4(float32_t* restrict R,const float32_t* restrict Fi,int M, int N,int L);

/*-----------------------------------------------------------
    rotate B[L][NxP] by diagonal matrix Fi'[L][N]
    Input:
    B
    Fi
    Output:
    B
    P==1
-----------------------------------------------------------*/
void qrsfRotateBconj1(float32_t* B,const float32_t* Fi,int N,int L);

/*---------------------------------------------------------
    update matrix B[L][MxP] by housholder vectors V[L][SV]

    Input:
    M,N,P,L     dimensions
    V[L][SV]    Householder vectors
    Input/output:
    B[L][MxP]   B matrices MxP
  
---------------------------------------------------------*/
void qrsfUpdateB(float32_t* B,const float32_t* V,int M,int P,int L);
// M=1,P=1
void qrsfUpdateB1(float32_t* B,const float32_t* V,int L);
// M=2,P=1
void qrsfUpdateB2(float32_t* B,const float32_t* V,int L);
// M=3,P=1
void qrsfUpdateB3(float32_t* B,const float32_t* V,int L);
// M=4,P=1
void qrsfUpdateB4(float32_t* B,const float32_t* V,int L);

#endif
