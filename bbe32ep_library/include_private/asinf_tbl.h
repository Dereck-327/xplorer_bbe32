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
    tables for asinf(x) approximation
*/
#ifndef ASINFTBL_H__
#define ASINFTBL_H__

/* Portable data types. */
#include "NatureDSP_types.h"
#include "common.h"

#define ASINF_ALG  0 /* 0 - 2 ULP code, 1 - 1 ULP code */

externC const union  ufloat32uint32 ALIGN(32) asinftbl[5];

#endif /* ASINFTBL_H__ */
