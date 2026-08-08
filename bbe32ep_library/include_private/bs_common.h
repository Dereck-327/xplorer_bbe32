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
#ifndef BS_COMMON_H__
#define BS_COMMON_H__

/*
	NatureDSP_Baseband library. matrix operations
	packed/streaming conversion
	Integrit, 2006-2016
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matop.h"
#include "common.h"

#if 1
#define INTLV(Y1,Y0,X1,X0) BBE_DSELNX16I(Y1,Y0,X1,X0,BBE_DSELI_DEINTERLEAVE_2)
#define DEINTLV1(Y1,Y0,X1,X0) BBE_DSELNX16I(Y1,Y0,X1,X0,BBE_DSELI_INTERLEAVE_2)
#define DEINTLV2(Y1,Y0,X1,X0) BBE_DSELNX16I(Y1,Y0,X1,X0,BBE_DSELI_INTERLEAVE_4)
#define DEINTLV3(Y1,Y0,X1,X0) {                            \
        Y0=BBE_SELNX16I(X1,X0,BBE_SELI_EXTRACT_LO_HALVES); \
        Y1=BBE_SELNX16I(X1,X0,BBE_SELI_EXTRACT_HI_HALVES);}
#else
#define INTLV(Y1,Y0,X1,X0) {                              \
    Y0=BBE_SELNX16I(X1,X0,BBE_SELI_EXTRACT_2_OF_4_OFF_0); \
    Y1=BBE_SELNX16I(X1,X0,BBE_SELI_EXTRACT_2_OF_4_OFF_2); }

#define DEINTLV1(Y1,Y0,X1,X0) {                          \
        Y0=BBE_SELNX16I(X1,X0,BBE_SELI_INTERLEAVE_2_LO); \
        Y1=BBE_SELNX16I(X1,X0,BBE_SELI_INTERLEAVE_2_HI);}

#define DEINTLV2(Y1,Y0,X1,X0) {                          \
        Y0=BBE_SELNX16I(X1,X0,BBE_SELI_INTERLEAVE_4_LO); \
        Y1=BBE_SELNX16I(X1,X0,BBE_SELI_INTERLEAVE_4_HI);}

#define DEINTLV3(Y1,Y0,X1,X0) {                          \
        Y0=BBE_SELNX16I(X1,X0,BBE_SELI_EXTRACT_LO_HALVES); \
        Y1=BBE_SELNX16I(X1,X0,BBE_SELI_EXTRACT_HI_HALVES);}
#endif

// internal functions

// convert 3x1
void cbs3x1(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L);
// convert 5x1,6x1,7x1
void cbs5_7x1(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L);
// converts 9x1...15x1
void cbs9_15x1(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L);
// mxn function that works with MN not a multiple of 16
void cbsmxn_odd(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L);
void cbsmxn_even(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L);
// convert 3x1
void csb3x1(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L);
// convert 5x1,6x1,7x1
void csb5_7x1(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L);
/* convert 9x1...15x1 */
void csb9_15x1(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L);
// for MN not a multiple of 16
void csbmxn_odd(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L);
void csbmxn_even(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L);
// convert 3x1
void rbs3x1(int16_t* restrict y, const int16_t* restrict x, int MN,int L);
// convert 5x1,6x1,7x1
void rbs5_7x1(int16_t* restrict y, const int16_t* restrict x, int MN, int L);
// converts 9x1...15x1
void rbs9_15x1(int16_t* restrict yy, const int16_t* restrict x, int MN, int L);
// converts 17x1...31x1
void rbs17_31x1(int16_t* restrict yy, const int16_t* restrict x, int MN, int L);

// convert 3x1
void rsb3x1(int16_t* restrict y, const int16_t* restrict x, int MN, int L);
// convert 5x1,6x1,7x1
void rsb5_7x1(int16_t* restrict y, const int16_t* restrict x, int MN, int L);
/* convert 9x1...15x1 */
void rsb9_15x1(int16_t* restrict y, const int16_t* restrict x, int MN, int L);
// converts 17x1...31x1
void rsb17_31x1(int16_t* restrict yy, const int16_t* restrict x, int MN, int L);

#endif  //BS_COMMON_H__
