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
#ifndef CXFIR_COMMON_H__
#define CXFIR_COMMON_H__
#include "NatureDSP_types.h"
#include "common.h"

        void cxfir_process_2 (complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict coef, int16_t * restrict delayLine, int M, int N);
        void cxfir_process_4 (complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict coef, int16_t * restrict delayLine, int M, int N);
        void cxfir_process_8 (complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict coef, int16_t * restrict delayLine, int M, int N);
        void cxfir_process_16(complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict coef, int16_t * restrict delayLine, int M, int N);
        void cxfir_process_8m(complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict coef, int16_t * restrict delayLine, int M, int N);
typedef void (*proc_fxn_t)   (complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict coef, int16_t * restrict delayLine, int M, int N);

#endif //CXFIR_COMMON_H__
