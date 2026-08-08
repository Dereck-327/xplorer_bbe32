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
    common data code for matrix multiplication
	Integrit, 2006-2016
*/
#ifndef RCMATMULNXMN_COMMON__
#define RCMATMULNXMN_COMMON__
#include "NatureDSP_types.h"
#include "NatureDSP_Math.h"
#include "common.h"
#ifdef __cplusplus
extern "C" {
#endif
/*    specialized functions*/
void rcmatmulnxmn_N4(void * pScr, complex_fract16 * restrict z, const int16_t * restrict x, const complex_fract16 * restrict y, int M, int L, int Q);
void rcmatmulnxmn_M4(void * pScr, complex_fract16 * restrict z, const int16_t * restrict x, const complex_fract16 * restrict y, int M, int L, int Q);
void rcmatmulnxmn_gen(void * pScr, complex_fract16 * restrict z, const int16_t * restrict x, const complex_fract16 * restrict y, int N, int M, int L, int Q);

/* return scratch size */
size_t rcmatmulnxmn_M4_getScratchSize(int N);
size_t rcmatmulnxmn_N4_getScratchSize(int M);
size_t rcmatmulnxmn_gen_getScratchSize(int N, int M);
#ifdef __cplusplus
}
#endif

#endif
