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
#ifndef SRCFIR_COMMON_H__
#define SRCFIR_COMMON_H__
#include "NatureDSP_types.h"
#include "common.h"

#ifdef __cplusplus
extern "C"
{
#endif
        void srcfir_process_33 (void* _srcfir, int16_t * restrict delay, complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict h, int M, int N);
        void srcfir_process_32 (void* _srcfir, int16_t * restrict delay, complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict h, int M, int N);
        void srcfir_process_16 (void* _srcfir, int16_t * restrict delay, complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict h, int M, int N);
        void srcfir_process_17 (void* _srcfir, int16_t * restrict delay, complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict h, int M, int N);
        void srcfir_process_16m(void* _srcfir, int16_t * restrict delay, complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict h, int M, int N);
typedef void (*proc_fxn_t)     (void* _srcfir, int16_t * restrict delay, complex_fract16 * restrict y, const complex_fract16 * restrict x, const int16_t * restrict h, int M, int N);
#ifdef __cplusplus
}
#endif

/* Filter instance structure. */
typedef struct
{
  uint32_t        magic;     // Instance pointer validation number
  int             M;         // Number of filter coefficients
  const int16_t * coef;      // (M+1)/2 filter coefficients
  int16_t *       delayLine; // Delay line for complex samples
  int16_t *       p0;        // Beginning of the delay line 
  int16_t *       p1;        // Ending of the delay line
  proc_fxn_t      process;
} srcfir_t, *srcfir_ptr_t;

#endif //SRCFIR_COMMON_H__
