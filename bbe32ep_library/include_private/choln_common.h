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
    Internal functions/definitions for block Cholesky
    IntegrIT, 2006-2017
*/
#ifndef CHOLN_COMMON_H__
#define CHOLN_COMMON_H__
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"

// bcholn_xxx require BBE_RSQRTLUNX40_0 and VSA arithmetic (BBE_SUBSR1SAVSN, etc.)
#define HAVE_CHOLN (HAVE_RSQRT && HAVE_VSAMATH)

/*---------------------------------------------------
   compute n-th column of A'*A for all L matrices

   Input:
   A[L][SA]     L matrices MxN
   n            number of column
   Output:
   Z[L][n+1][2] results
---------------------------------------------------*/
// n=0...7
void cholnConv0_7(int32_t* Z,const int16_t* A,int n,int M,int N,int L,int SA);
// n=8...15
void cholnConv8_15(int32_t* Z,const int16_t* A,int n,int M,int N,int L,int SA);
// n=16...23
void cholnConv16_23(int32_t* Z,const int16_t* A,int n,int M,int N,int L,int SA);
// n=24...31
void cholnConv24_31(int32_t* Z,const int16_t* A,int n,int M,int N,int L,int SA);

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
void cholnDiagUpd(int16_t* y,int16_t* D,const int32_t* Z,int N,int L,int SR,int SD);
/* the same as above but with limitations on N */
/* N<8 */
void cholnDiagUpd8(int16_t* y,int16_t* D,const int32_t* Z,int N,int L,int SR,int SD);
/* N<16 */
void cholnDiagUpd16(int16_t* y,int16_t* D,const int32_t* Z,int N,int L,int SR,int SD);
/* N<24 */
void cholnDiagUpd24(int16_t* y,int16_t* D,const int32_t* Z,int N,int L,int SR,int SD);
/* N<32 */
void cholnDiagUpd32(int16_t* y,int16_t* D,const int32_t* Z,int N,int L,int SR,int SD);

/* one iteration for Cholesky */
typedef void (*fnConv)(int32_t* Z,const int16_t* A,int n,int M,int N,int L,int SA);
typedef void (*fnFwdrec)(int16_t* y,const int16_t* R,const int16_t* D,const int32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);
typedef void (*fnDiagUpd)(int16_t* y,int16_t* D,const int32_t* Z,int N,int L,int SR,int SD);
typedef struct
{
    fnConv    conv;
    fnFwdrec  fwdrec;
    fnDiagUpd diagUpd;
} tCholnIteration;

/* --------------------------------------------------
   make forward recursion to update n column elements
   Input:
   Z[L][SZ]  convolutions in N-th column
   D[L][SD]  reciprocals of main diagonal
   Output:
   y[L][SY]  result of recursion (N elements filled)
--------------------------------------------------*/
void cholnFwdrec(int16_t* y,const int16_t* R,const int16_t* D,const int32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);
/* the same as above but with limitations on N */
/* N==0 */
void cholnFwdrec0(int16_t* y,const int16_t* R,const int16_t* D,const int32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);
/* N<=8 */
void cholnFwdrec8(int16_t* y,const int16_t* R,const int16_t* D,const int32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);
/* N<=16 */
void cholnFwdrec16(int16_t* y,const int16_t* R,const int16_t* D,const int32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);
/* N<=24 */
void cholnFwdrec24(int16_t* y,const int16_t* R,const int16_t* D,const int32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);
/* N<=32 */
void cholnFwdrec32(int16_t* y,const int16_t* R,const int16_t* D,const int32_t* Z,int N,int L,int SR,int SD,int SY,int SZ);

/*
    backward recursion: P==1
*/
void cholnBkwnx1( int16_t* restrict x, 
            const int16_t* restrict Rt,
            const int16_t* restrict D,
            const int16_t* restrict y, 
                  int qXYA,
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
void cholnTransformR(int16_t* Rt,const int16_t* R,int N,int L);

#endif
