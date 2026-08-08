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
    Internal functions/definitions for block Cholesky (complex floating point)
    IntegrIT, 2006-2017
*/
#ifndef CHOLNF_COMMON_H__
#define CHOLNF_COMMON_H__
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"

// bcholn_xxx require BBE_RSQRTLUNX40_0 and VSA arithmetic (BBE_SUBSR1SAVSN, etc.)
#define HAVE_CHOLNF (HAVE_VFPU)

/*---------------------------------------------------
   compute n-th column of A'*A for all L matrices

   Input:
   A[L][SA]     L matrices MxN
   sigma2[L]    regularization term
   n            number of column
   Output:
   Z[L][n+1][2] results
---------------------------------------------------*/
/* n=0...3   */ void cholnfConv0_3  (float32_t* Z,const float32_t* A,const float32_t* sigma2,int n,int M,int N,int L,int SA);
/* n=4...7   */ void cholnfConv4_7  (float32_t* Z,const float32_t* A,const float32_t* sigma2,int n,int M,int N,int L,int SA);
/* n=8...11  */ void cholnfConv8_11 (float32_t* Z,const float32_t* A,const float32_t* sigma2,int n,int M,int N,int L,int SA);
/* n=12...15 */ void cholnfConv12_15(float32_t* Z,const float32_t* A,const float32_t* sigma2,int n,int M,int N,int L,int SA);

/*-------------------------------------------------
   update N-th diagonal element
   Input:
   Z[L][N+1][2]  convolutions in N-th column
   Input/output:
   y             pointer to the begining of column 
                 in matrix R[L][SR] (N+1 elements 
                 is written)
   Output:
   D[L][SD]      reciprocals of main diagonal 
                 (pointer to the N-th element
-------------------------------------------------*/
void cholnfDiagUpd(float32_t* y,float32_t* D,const float32_t* Z,int N,int L,int SR,int SD);
/* the same as above but with limitations on N */
/* N<4  */void cholnfDiagUpd4 (float32_t* y,float32_t* D,const float32_t* Z,int N,int L,int SR,int SD);
/* N<8  */void cholnfDiagUpd8 (float32_t* y,float32_t* D,const float32_t* Z,int N,int L,int SR,int SD);
/* N<12 */void cholnfDiagUpd12(float32_t* y,float32_t* D,const float32_t* Z,int N,int L,int SR,int SD);
/* N<16 */void cholnfDiagUpd16(float32_t* y,float32_t* D,const float32_t* Z,int N,int L,int SR,int SD);

/* one iteration for Cholesky */
typedef void (*fnConv)(float32_t* Z,const float32_t* A,const float32_t *sigma2,int n,int M,int N,int L,int SA);
typedef void (*fnFwdrec)(float32_t* y,const float32_t* R,const float32_t* D,const float32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);
typedef void (*fnDiagUpd)(float32_t* y,float32_t* D,const float32_t* Z,int N,int L,int SR,int SD);
typedef struct
{
    fnConv    conv;
    fnFwdrec  fwdrec;
    fnDiagUpd diagUpd;
} tCholnfIteration;

/* --------------------------------------------------
   make forward recursion to update n column elements
   Input:
   Z[L][SZ]  convolutions in N-th column
   D[L][SD]  reciprocals of main diagonal
   Output:
   y[L][SY]  result of recursion (N elements filled)
--------------------------------------------------*/
void cholnfFwdrec(float32_t* y,const float32_t* R,const float32_t* D,const float32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);
/* the same as above but with limitations on N */
/* N== 0 */void cholnfFwdrec0 (float32_t* y,const float32_t* R,const float32_t* D,const float32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);
/* N<= 4 */void cholnfFwdrec4 (float32_t* y,const float32_t* R,const float32_t* D,const float32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);
/* N<= 8 */void cholnfFwdrec8 (float32_t* y,const float32_t* R,const float32_t* D,const float32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);
/* N<=12 */void cholnfFwdrec12(float32_t* y,const float32_t* R,const float32_t* D,const float32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);
/* N<=16 */void cholnfFwdrec16(float32_t* y,const float32_t* R,const float32_t* D,const float32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);

/*--------------------------------------------------
    backward recursion: P==1
    Input:
    Rt
    D
    y
    N,L  size of matrices
    Output:
    X
--------------------------------------------------*/
void cholnfBkwnx1( float32_t* restrict x, 
             const float32_t* restrict Rt,
             const float32_t* restrict D,
             const float32_t* restrict y, 
                  int N,int L);

/*----------------------------------------------------------------------------------
   reversing R matrices for easier readings by rows (diagonal elements are omitted):
   original R    transformed R
   0 1 3 6 a     d 8 c 4 7 b 1 3 6 a
     2 4 7 b
       5 8 c
         9 d
           e

   Input:
   R[L][SR]        L input matrices
   Rt[L*N*(N-1)]   stream of L trasposed matrices
----------------------------------------------------------------------------------*/
void cholnfTransformR(float32_t* Rt,const float32_t* R,int N,int L);

#endif
