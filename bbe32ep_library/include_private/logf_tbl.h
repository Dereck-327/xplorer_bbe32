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
    tables for logf(x) approximation
*/
#ifndef LOGF_TBL_H__
#define LOGF_TBL_H__
#include "NatureDSP_types.h"
#include "logf_tbl.h"
#include "common.h"

#define LOGF_ALG 0  /* 0 - 2 ULP algorithm, 1 - 1 ULP algorithm */

#if LOGF_ALG==0
/* 
   polynomial coefficients for ln(x)/(1-x) 
   derived by MATLAB code:
   x=(sqrt(0.5):pow2(1,-16):sqrt(2));
   z=1-x;
   y=log(x)./z;
   p=polyfit(z,y,8);

*/
externC const union ufloat32uint32 ALIGN(32) logf_tbl[];
#elif LOGF_ALG==1
externC const int32_t ALIGN(32) logf_tbl_q31[];
#else
#error wrong LOGF_ALG
#endif
externC const union ufloat32uint32 ln2;

#endif /* LOG2F_TBL_H__ */
