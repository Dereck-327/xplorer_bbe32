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
#ifndef BQRIIRS_COMMON__
#define BQRIIRS_COMMON__

#include "NatureDSP_types.h"

/* Biquad real block IIR processing function, Fast fixed-point implementation. */
void bqriirs_sp_proc( int16_t * restrict r,
                      int16_t * restrict sect, // 2*M
                const int16_t *          x,
                const int16_t *          coef,
                      int16_t            gain,
                      int N, int L, int M );

/* Biquad real block IIR processing function, Low Noise fixed-point implementation. */
void bqriirs_dp_proc( int16_t * restrict r,
                      int16_t * restrict sect, // 6*M
                const int16_t *          x,
                      int16_t * restrict scr,  // 4*N
                const int16_t *          coef,
                      int16_t            gain,
                      int N, int L, int M );
#endif // BQRIIRS_COMMON__
