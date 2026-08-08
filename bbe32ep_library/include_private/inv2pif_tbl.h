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
    1/(2*pi) constant
*/
#ifndef INV2PIF_H__
#define INV2PIF_H__

#include "NatureDSP_types.h"
#include "common.h"

externC const int64_t inv2pif_Q53; /* 1/(2pi) in Q53 */

externC const union ufloat32uint32 inv2pif; /* 2/pi */
externC const union ufloat32uint32 inv4pif; /* 4/pi */

#endif /* INV2PIF_H__ */
