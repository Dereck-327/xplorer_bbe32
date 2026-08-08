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

#ifndef FIRINTERPF_COMMON_H__
#define FIRINTERPF_COMMON_H__

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
tFirInitAlloc_int;

typedef struct
{
    const tFirInitAlloc_int * initAllocFxdxn;
    proc_fxdxn_t              procFxdxn;
}
tFirFxdxns;

extern const tFirInitAlloc_int firinterpf_dx        ;
extern const tFirInitAlloc_int firinterpf_d_2_3_4_mx;
extern const tFirInitAlloc_int firinterpf_d_6_12    ;
extern const tFirInitAlloc_int firinterpf_d_6_12_mx ;
extern const tFirInitAlloc_int firinterpf_dx_mx     ;

extern const tFirFxdxns interpf_2d_2m  ;
extern const tFirFxdxns interpf_2d_4m  ;
extern const tFirFxdxns interpf_2d_8m  ;
extern const tFirFxdxns interpf_2d_16m ;
extern const tFirFxdxns interpf_2d_xm  ;
extern const tFirFxdxns interpf_3d_2m  ;
extern const tFirFxdxns interpf_3d_4m  ;
extern const tFirFxdxns interpf_3d_8m  ;
extern const tFirFxdxns interpf_3d_16m ;
extern const tFirFxdxns interpf_3d_xm  ;
extern const tFirFxdxns interpf_4d_2m  ;
extern const tFirFxdxns interpf_4d_4m  ;
extern const tFirFxdxns interpf_4d_8m  ;
extern const tFirFxdxns interpf_4d_16m ;
extern const tFirFxdxns interpf_4d_xm  ;
extern const tFirFxdxns interpf_6d_2m  ;
extern const tFirFxdxns interpf_6d_4m  ;
extern const tFirFxdxns interpf_6d_8m  ;
extern const tFirFxdxns interpf_6d_16m ;
extern const tFirFxdxns interpf_6d_xm  ;
extern const tFirFxdxns interpf_12d_2m ;
extern const tFirFxdxns interpf_12d_4m ;
extern const tFirFxdxns interpf_12d_8m ;
extern const tFirFxdxns interpf_12d_16m;
extern const tFirFxdxns interpf_12d_xm ;
extern const tFirFxdxns interpf_xd_xm  ;

/* interleaving functions */
void firinterpf_12d_intlv(complex_float * y, int N);
void firinterpf_6d_intlv (complex_float * y, int N);

#ifdef __cplusplus
extern "C"
{
#endif
void firinterpf_proc_D2_M2( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D2_M4( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D2_M8( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D2_M16( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D2_MX( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D3_M2( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D3_M4( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D3_M8( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D3_M16( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D3_MX( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D4_M2( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D4_M4( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D4_M8( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D4_M16( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D4_MX( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D6_M2( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D6_M4( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D6_M8( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D6_M16( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D6_MX( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D12_M2( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D12_M4( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D12_M8( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D12_M16( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_D12_MX( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );
void firinterpf_proc_DX_MX( complex_float * restrict y,const complex_float * restrict x,const float32_t     * restrict coef,float32_t     * restrict delayLine,int M, int N, int D );

#ifdef __cplusplus
}
#endif
#endif //FIRINTERPF_COMMON_H__


