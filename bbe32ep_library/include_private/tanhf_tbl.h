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
    tables for tanhf() approximation
*/
#ifndef TANHF_TBL_H__
#define TANHF_TBL_H__
#include "NatureDSP_types.h"
#include "common.h"

#define TANHF_ALG 0 /* 0 - 2 ULP, 1 - 1 ULP */

externC const union ufloat32uint32 ALIGN(32) polytanhf_tbl[];

externC const union ufloat32uint32 halfln3 ; /* log(3)/2 - tanh(log(3)/2)==0.5 */

#endif
