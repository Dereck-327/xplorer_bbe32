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
 * NatureDSP_Baseband Library API
 * Matrix Decomposition and Inversion Functions

    LU decomposition for complex matrices (block ordered)
    IntegrIT, 2006-2017
*/
#ifndef CLUNF_COMMON_H__
#define CLUNF_COMMON_H__
#include "NatureDSP_types.h"

/*--------------------------------------------------
update rows below i-th in matrix NxN (8x8, 16x16 only)

Input:
norm[L]   - normalization value for diagonal Aji, j=i+1...N-1
N           matrix size (8 or 16)
L           number of matrices
Input/output:
A[(N-i)*N]  matrix - pointer to the i-th row
Returns 
none:
--------------------------------------------------*/
void clunf_update(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
// specialized versions for different N-i
void clunf_update2(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update3(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update4(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update5(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update6(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update7(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update8(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update9(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update10(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update11(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update12(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update13(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update14(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update15(complex_float* A, const complex_float* restrict norm, int i, int N, int L);
void clunf_update16(complex_float* A, const complex_float* restrict norm, int i, int N, int L);

#endif
