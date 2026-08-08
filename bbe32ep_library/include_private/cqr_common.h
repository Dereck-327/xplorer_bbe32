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
#ifndef CQR_COMMON_H__
#define CQR_COMMON_H__

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"


#define MAX_M 16 
#define SIZE_OF_V(M,N) (((M)==(N)) ? (((M)-1)*((M)+2)/2) : ((M)*(N)+(((N)-1)*(N))/2))
#define SIZE_OF_FI(M)  (M)

#define CLIP_TO_ZERO(X) (((X)<0)?0:X)

/* divide half of vector (even elements only) */
#define QUO8X32(res,dvdr,dvsr) {       \
    BBE_DIVNX32S_5STEP0_0(dvdr,dvsr);   \
    BBE_DIVNX16S_4STEP_0(dvsr);         \
    BBE_DIVNX16S_4STEP_0(dvsr);         \
    res=BBE_DIVNX16S_3STEPN_0(dvsr);    \
}

#ifdef __cplusplus
extern "C" {
#endif
    int cqr_calc_qb5x5xps(int16_t *B, const int16_t *V, int P, int L);
    int cqr_calc_qb6x6xps(int16_t *B, const int16_t *V, int P, int L);
    int cqr_calc_qb7x7xps(int16_t *B, const int16_t *V, int P, int L);


    void cqrComputeLastRot(int16_t* pFi, const int16_t* R, int M, int L);
    void cqrHouseholder2(
        const int16_t *   A, /* input matrix               */
        int16_t *    v,     /* output (housholder vector) */
        int16_t*    Fi_out,  /* output phase rotator       */
        const int M,         /* number of rows             */
        const int N,         /* number of columns          */
        const int L          /* number of matrices         */
        );
    void cqrHouseholder3(
        const int16_t *   A, /* input matrix               */
        int16_t *    v,     /* output (housholder vector) */
        int16_t*    Fi_out,  /* output phase rotator       */
        const int M,         /* number of rows             */
        const int N,         /* number of columns          */
        const int L          /* number of matrices         */
        );
    void cqrHouseholder4(
        const int16_t *   A, /* input matrix               */
        int16_t *    v,     /* output (housholder vector) */
        int16_t*    Fi_out,  /* output phase rotator       */
        const int M,         /* number of rows             */
        const int N,         /* number of columns          */
        const int L          /* number of matrices         */
        );
    void cqrHouseholder5(
        void* pScr,          /* scratch (4*L bytes)        */
        const int16_t *   A, /* input matrix               */
        int16_t *    v,     /* output (housholder vector) */
        int16_t*    Fi_out,  /* output phase rotator       */
        const int M,         /* number of rows             */
        const int N,         /* number of columns          */
        const int L          /* number of matrices         */
        );
    void cqrHouseholder6(
        void* pScr,          /* scratch (4*L bytes)        */
        const int16_t *   A, /* input matrix               */
        int16_t *    v,     /* output (housholder vector) */
        int16_t*    Fi_out,  /* output phase rotator       */
        const int M,         /* number of rows             */
        const int N,         /* number of columns          */
        const int L          /* number of matrices         */
        );
    void cqrHouseholder7(
        void* pScr,          /* scratch (4*L bytes)        */
        const int16_t *   A, /* input matrix               */
        int16_t *    v,     /* output (housholder vector) */
        int16_t*    Fi_out,  /* output phase rotator       */
        const int M,         /* number of rows             */
        const int N,         /* number of columns          */
        const int L          /* number of matrices         */
        );
    void cqrHouseholder8(
        void* pScr,          /* scratch (4*L bytes)        */
        const int16_t *   A, /* input matrix               */
        int16_t *    v,     /* output (housholder vector) */
        int16_t*    Fi_out,  /* output phase rotator       */
        const int M,         /* number of rows             */
        const int N,         /* number of columns          */
        const int L          /* number of matrices         */
        );
    void cqrUpdateR2(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void cqrUpdateR3(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void cqrUpdateR4(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void cqrUpdateR5(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void cqrUpdateR6(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void cqrUpdateR7(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void cqrUpdateR8(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void cqrUpdateR9(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void cqrUpdateR10(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void cqrUpdateR11(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void cqrUpdateR12(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
    void cqrUpdateR13(int16_t* restrict R,
        const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
        int m, int M, int N, int L);
#ifdef __cplusplus
}
#endif
#endif  //CQR_COMMON_H__
