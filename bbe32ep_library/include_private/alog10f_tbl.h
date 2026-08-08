/* ------------------------------------------------------------------------ */
/* IntegrIT, Ltd.   www.integrIT.com, info@integrIT.com                     */
/*                                                                          */
/* DSP Library                                                              */
/*                                                                          */
/* This library contains copyrighted materials, trade secrets and other     */
/* proprietary information of IntegrIT, Ltd. This software is licensed for  */
/* use with Cadence processor cores only and must not be used for any other */
/* processors and platforms. The license to use these sources was given to  */
/* Cadence, Inc. under Terms and Condition of a Software License Agreement  */
/* between Cadence, Inc. and IntegrIT, Ltd.                                 */
/* ------------------------------------------------------------------------ */
/*          Copyright (C) 2014-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */

/*
    tables for alog10f(x) approximation
*/
#ifndef ALOG10FTBL_H__
#define ALOG10FTBL_H__

/* Portable data types. */
#include "NatureDSP_types.h"
#include "common.h"

externC const union ufloat32uint32 alog10fminmax[2];  /* minimum and maximum arguments of alog10f() input */
externC const int32_t invlog10_2_Q29; /* 1/log10(2), Q29 */
externC const union ufloat32uint32 ALIGN(32) log2_10[2];
#endif /* ALOG10FTBL_H__ */
