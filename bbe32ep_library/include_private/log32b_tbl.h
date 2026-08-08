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
    tables for log(x) approximation, Q31
*/
#ifndef LOG32B_THL_H__
#define LOG32B_THL_H__
#include "NatureDSP_types.h"
#include "common.h"

/* A=round(pow2(2.^((1-(2:17))/16),31)); 
   First Cody coefficient is exact 1Q31 and can be omitted
*/
externC const int32_t ALIGN(32) Alog2_Q31_tbl[];
externC const int32_t ALIGN(32) Plog2_Q31_tbl[];

#endif
