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
#ifndef BKFIRF_COMMON_H__
#define BKFIRF_COMMON_H__
#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C"
{
#endif
        void bkfirf_process_8m(float32_t * restrict delay, float32_t * restrict y, const float32_t * restrict x, const float32_t * restrict h, int M, int N);
        void bkfirf_process_16(float32_t * restrict delay, float32_t * restrict y, const float32_t * restrict x, const float32_t * restrict h, int M, int N);
        void bkfirf_process_8 (float32_t * restrict delay, float32_t * restrict y, const float32_t * restrict x, const float32_t * restrict h, int M, int N);
        void bkfirf_process_4 (float32_t * restrict delay, float32_t * restrict y, const float32_t * restrict x, const float32_t * restrict h, int M, int N);
        void bkfirf_process_2 (float32_t * restrict delay, float32_t * restrict y, const float32_t * restrict x, const float32_t * restrict h, int M, int N);
#ifdef __cplusplus
}
#endif

typedef void (*proc_fxn_t)    (float32_t * restrict delay, float32_t * restrict y, const float32_t * restrict x, const float32_t * restrict h, int M, int N);

#endif //BKFIRF_COMMON_H__
