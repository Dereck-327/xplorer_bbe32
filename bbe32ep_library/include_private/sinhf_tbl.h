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
	tables for sinhf(x) approximation
*/
#ifndef SINHF_TBL_H__
#define SINHF_TBL_H__

#include "NatureDSP_types.h"
#include "common.h"

#define SINHF_ALG 0 /* 0 - 2 ULP, 1- 1 ULP */

/* Maximum magnitude of input value such that sinhf does not overflow. */
externC const union ufloat32uint32 sinhf_maxarg;

/* polynomial coefficients for sinh(x) in range +/-1. 
  
   s=pow2(2,-16); x=(s:s:1); x=[-x(end:-1:1) x]; y=sinh(x)./x;
   p=polyfit(x,y,6);
   p=p(1:2:end); p(end)=[];  
*/
externC const union ufloat32uint32 ALIGN(32) polysinhf_tbl[];

#endif 
