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
    tables for atanf(x) approximation
*/
#ifndef ATANFTBL_H__
#define ATANFTBL_H__

/* Portable data types. */
#include "NatureDSP_types.h"
#include "common.h"

externC const union ufloat32uint32 atanftbl1[8];
externC const union ufloat32uint32 atanftbl2[8];
externC const union ufloat32uint32 atanftbl1a[8]; 
externC const union ufloat32uint32 atanftbl2a[8];

externC const union ufloat32uint32 atanftbl_10ord[10];

#endif /* ATANFTBL_H__ */
