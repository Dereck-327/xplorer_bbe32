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
#ifndef FIRDEC_COMMON_H__
#define FIRDEC_COMMON_H__

#include "NatureDSP_types.h"
#include "common.h"

typedef void (*proc_fxdxn_t)(int16_t * restrict y, const int16_t * restrict x, const int16_t * restrict coef, int16_t * restrict delayLine, int M, int N, int D);

typedef struct
{
  int coefNum;
  int delLength;
}
tFilterLayout;

typedef struct
{
    void (*fnalloc)(tFilterLayout* pFltr, int M, int D);
    void (*fninit )(int16_t* restrict coef, const int16_t* restrict h, int M, int D);
}
tFirDecAlloc;

typedef struct
{
  const tFirDecAlloc *   initAllocFxdxn;
  proc_fxdxn_t           procFxdxn;
}
tFirFxdxns;

extern const tFirDecAlloc firdec_alloc_gen   ;
extern const tFirDecAlloc firdec_alloc_d2_m2 ;
extern const tFirDecAlloc firdec_alloc_d2_mx ;
extern const tFirDecAlloc firdec_alloc_d3_m16;
extern const tFirDecAlloc firdec_alloc_d3_m32;
extern const tFirDecAlloc firdec_alloc_d3_mx ;
extern const tFirDecAlloc firdec_alloc_d3_m8 ;
extern const tFirDecAlloc firdec_alloc_d3_m4 ;
extern const tFirDecAlloc firdec_alloc_dx    ;

extern const tFirFxdxns firdec_2d_2_8n    ;
extern const tFirFxdxns firdec_2d_4_8n    ;
extern const tFirFxdxns firdec_2d_8_8n    ;
extern const tFirFxdxns firdec_2d_16_8n   ;
extern const tFirFxdxns firdec_2d_32_8n   ;
extern const tFirFxdxns firdec_2d_x_8n    ;
extern const tFirFxdxns firdec_3d_2_8n    ;
extern const tFirFxdxns firdec_3d_4_8n    ;
extern const tFirFxdxns firdec_3d_8_8n    ;
extern const tFirFxdxns firdec_3d_16_8n   ;
extern const tFirFxdxns firdec_3d_32_8n   ;
extern const tFirFxdxns firdec_3d_x_8n    ;
extern const tFirFxdxns firdec_4d_2_8n    ;
extern const tFirFxdxns firdec_4d_4_8n    ;
extern const tFirFxdxns firdec_4d_8_8n    ;
extern const tFirFxdxns firdec_4d_16_8n   ;
extern const tFirFxdxns firdec_4d_32_8n   ;
extern const tFirFxdxns firdec_4d_xeven_8n;
extern const tFirFxdxns firdec_4d_xodd_8n ;
extern const tFirFxdxns firdec_dx_8n      ;

#ifdef __cplusplus
extern "C"
{
#endif
    void firdec_proc_D2_M2(int16_t * restrict y, const int16_t * restrict x, const int16_t * restrict coef, int16_t * restrict delayLine, int M, int N, int D);
    void firdec_proc_D3_M2(int16_t * restrict y, const int16_t * restrict x, const int16_t * restrict coef, int16_t * restrict delayLine, int M, int N, int D);
    void firdec_proc_D4_M2(int16_t * restrict y, const int16_t * restrict x, const int16_t * restrict coef, int16_t * restrict delayLine, int M, int N, int D);
#ifdef __cplusplus
}
#endif

#endif
