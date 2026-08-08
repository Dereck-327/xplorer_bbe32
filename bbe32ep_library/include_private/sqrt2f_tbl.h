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
    sqrt(2) and sqrt(0.5) constants for single precision routines
*/
#ifndef SQRT2F_H__
#define SQRT2F_H__

/* Portable data types. */
#include "NatureDSP_types.h"
#include "common.h"

externC const union ufloat32uint32 sqrt2f  ; /* sqrt(2)   */
externC const union ufloat32uint32 sqrt0_5f; /* sqrt(0.5) */

#endif /* SQRT2F_H__ */
