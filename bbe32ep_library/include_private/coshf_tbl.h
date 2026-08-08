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
	tables for coshf(x) approximation
*/
#ifndef COSHF_THBL_H__
#define COSHF_THBL_H__
#include "NatureDSP_types.h"
#include "common.h"

/* cosh algorithm options: 
 *   0 - direct polynomial approximation for small argument: |x|<=0.25
 *   1 - computation via exp(x) estimation, for any argument. */
#define COSHF_ALG   1

/* Maximum magnitude of input value such that coshf does not overflow. */
externC const union ufloat32uint32 coshf_maxarg;

#if COSHF_ALG==0
/* polynomial coefficients for cosh(x) in range +/-0.25 

   x=(-0.25:pow2(2,-16):0.25); y=cosh(x);
   polyfit(x,y,4)
   p=p(1:2:end); p(end)=[];
*/
externC const union ufloat32uint32 ALIGN(32) polycoshf_tbl[];
#endif

#endif
