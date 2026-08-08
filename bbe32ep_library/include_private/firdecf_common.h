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
	NatureDSP_Baseband library. FIR filters part.
    Decimating Block Complex FIR Filter with Real Coefficients: common definitions
    Optimized code for BBE32EP
	IntegrIT, 2006-2017
*/
#ifndef FIRDECF_COMMON_H__
#define FIRDECF_COMMON_H__

#include "NatureDSP_types.h"
#include "common.h"

typedef void(*proc_fxdxn_t)(complex_float * restrict y, const complex_float * restrict x, const float32_t * restrict coef, float32_t * restrict delayLine, int M, int N, int D);

typedef struct
{
    int coefNum;
    int delLength;
}
tFilterLayout;

typedef struct
{
    void (*fnalloc)(tFilterLayout * pFltr, int M, int D);
    void (*fninit )(float32_t * restrict coef, const float32_t * restrict h, int M, int D);
}
tFirDecAlloc;

typedef struct
{
    const tFirDecAlloc * initAllocFxdxn;
    proc_fxdxn_t         procFxdxn;
}
tFirFxdxns;

extern const tFirDecAlloc firdecf_alloc_dx     ;
extern const tFirDecAlloc firdecf_alloc_d2_4_mx;
extern const tFirDecAlloc firdecf_alloc_d3_mx  ;
extern const tFirDecAlloc firdecf_alloc_dx_mx  ;

extern const tFirFxdxns firdecf_2d_2m ;
extern const tFirFxdxns firdecf_2d_4m ;
extern const tFirFxdxns firdecf_2d_8m ;
extern const tFirFxdxns firdecf_2d_16m;
extern const tFirFxdxns firdecf_2d_xm ;
extern const tFirFxdxns firdecf_3d_2m ;
extern const tFirFxdxns firdecf_3d_4m ;
extern const tFirFxdxns firdecf_3d_8m ;
extern const tFirFxdxns firdecf_3d_16m;
extern const tFirFxdxns firdecf_3d_xm ;
extern const tFirFxdxns firdecf_4d_2m ;
extern const tFirFxdxns firdecf_4d_4m ;
extern const tFirFxdxns firdecf_4d_8m ;
extern const tFirFxdxns firdecf_4d_16m;
extern const tFirFxdxns firdecf_4d_xm ;
extern const tFirFxdxns firdecf_xd_xm ;


#ifdef __cplusplus
extern "C"
{
#endif
void firdecf_proc_D2_M2( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D2_M4( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D2_M8( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D2_M16( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D2_MX( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D3_M2( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D3_M4( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D3_M8( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D3_M16( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D3_MX( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D4_M2( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D4_M4( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D4_M8( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D4_M16( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_D4_MX( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
void firdecf_proc_DX_MX( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D );
#ifdef __cplusplus
}
#endif

#endif //FIRDECF_COMMON_H__
