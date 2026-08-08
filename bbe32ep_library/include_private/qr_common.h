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

#ifndef _QR_COMMON_H_
#define _QR_COMMON_H_

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"

#define HAVE_QR (XCHAL_HAVE_BBEN_VECDIVIDE && XCHAL_HAVE_BBEN_RSQRT)

/* divide half of vector (even elements only) */
#define QUO16X32(res,dvdr,dvsr) {       \
    BBE_DIVNX32S_5STEP0_0(dvdr,dvsr);   \
    BBE_DIVNX16S_4STEP_0(dvsr);         \
    BBE_DIVNX16S_4STEP_0(dvsr);         \
    res=BBE_DIVNX16S_3STEPN_0(dvsr);    \
}

/* divide full vector  */
#define QUO32X32(res,dvdr,dvsr) {       \
    BBE_DIVNX32S_5STEP0_0(dvdr,dvsr);   \
    BBE_DIVNX32S_5STEP0_1(dvdr,dvsr);   \
    BBE_DIVNX16S_4STEP_0(dvsr);         \
    BBE_DIVNX16S_4STEP_1(dvsr);         \
    BBE_DIVNX16S_4STEP_0(dvsr);         \
    BBE_DIVNX16S_4STEP_1(dvsr);         \
    res=BBE_DIVNX16S_3STEPN_0(dvsr);    \
    res=BBE_DIVNX16S_3STEPN_1(dvsr);    \
}

#define SIZE_OF_V(M,N) (((M)==(N)) ? (((M)-1)*((M)+2)/2) : ((M)*(N)-(((N)-1)*(N))/2))

#ifdef __cplusplus
extern "C" {
#endif
    void qr_calc_qb8x8xps(int16_t *B, const int16_t *V, int P, int L);


    void qrHouseholder2(
        const int16_t * restrict A, /* input matrix         */
        int16_t * restrict  v,     /* output streaming order       */
        const int M,
        const int N,
        const int L
        );
    void qrHouseholder3(
        const int16_t * restrict A, /* input matrix         */
        int16_t * restrict  v,     /* output streaming order       */
        const int M,
        const int N,
        const int L
        );
    void qrHouseholder4(
        const int16_t * restrict A, /* input matrix         */
        int16_t * restrict  v,     /* output streaming order       */
        const int M,
        const int N,
        const int L
        );
    void qrHouseholder5(
        const int16_t * restrict A, /* input matrix         */
        int16_t * restrict  v,     /* output streaming order       */
        const int M,
        const int N,
        const int L
        );
    void qrHouseholder6(
        const int16_t * restrict A, /* input matrix         */
        int16_t * restrict  v,     /* output streaming order       */
        const int M,
        const int N,
        const int L
        );
    void qrHouseholder7(
        const int16_t * restrict A, /* input matrix         */
        int16_t * restrict  v,     /* output streaming order       */
        const int M,
        const int N,
        const int L
        );
    void qrHouseholder8(
        const int16_t * restrict A, /* input matrix         */
        int16_t * restrict  v,     /* output streaming order       */
        const int M,
        const int N,
        const int L
        );

    void qrUpdateR2(int16_t* restrict R,
        const int16_t * restrict V,     /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void qrUpdateR3(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void qrUpdateR4(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void qrUpdateR5(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void qrUpdateR6(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void qrUpdateR7(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void qrUpdateR8(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
#ifdef __cplusplus
}
#endif
#endif // _QR_COMMON_H_
