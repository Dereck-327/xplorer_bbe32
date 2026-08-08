/* ------------------------------------------------------------------------ */
/* Copyright (c) 2017 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs ('Cadence    */
/* Libraries') are the copyrighted works of Cadence Design Systems Inc.     */
/* Cadence IP is licensed for use with Cadence processor cores only and     */
/* must not be used for any other processors and platforms. Your use of the */
/* Cadence Libraries is subject to the terms of the license agreement you   */
/* have entered into with Cadence Design Systems, or a sublicense granted   */
/* to you by a direct Cadence licensee.                                     */
/* ------------------------------------------------------------------------ */
/*  IntegrIT, Ltd.   www.integrIT.com, info@integrIT.com                    */
/*                                                                          */
/* NatureDSP_Baseband Library                                               */
/*                                                                          */
/* This library contains copyrighted materials, trade secrets and other     */
/* proprietary information of IntegrIT, Ltd. This software is licensed for  */
/* use with Cadence processor cores only and must not be used for any other */
/* processors and platforms. The license to use these sources was given to  */
/* Cadence, Inc. under Terms and Condition of a Software License Agreement  */
/* between Cadence, Inc. and IntegrIT, Ltd.                                 */
/* ------------------------------------------------------------------------ */
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */

#ifndef __NATUREDSP_BASEBAND_H
#define __NATUREDSP_BASEBAND_H
   
#include "NatureDSP_types.h"
/* FIR Filters and Related Functions */
#include "NatureDSP_Baseband_fir.h"
/* IIR Filters */
#include "NatureDSP_Baseband_iir.h"
/* FFT Routines */
#include "NatureDSP_Baseband_fft.h"
/* Matrix Operations */
#include "NatureDSP_Baseband_matop.h"
/* Matrix Decomposition and Inversion Functions */
#include "NatureDSP_Baseband_matinv.h"
/* Vector Operations */
#include "NatureDSP_Baseband_vector.h"
/* Math Functions */
#include "NatureDSP_Baseband_math.h"
/* Complex Math Functions */
#include "NatureDSP_Baseband_complex.h"
/* Communications */
#include "NatureDSP_Baseband_comm.h"
/* Fitting and Interpolation Routines */
#include "NatureDSP_Baseband_fit.h"
/* library identification */
#include "NatureDSP_Baseband_id.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
  Aligned data definition
===========================================================================*/
/* aligned arrays */
#if (defined COMPILER_XTENSA)||(defined COMPILER_GNU)
  #define ALIGNED_ARRAY(n,name) \
    typedef int16_t tInt16a_##name[n] __attribute__ ((aligned(32)));
  #define ALIGNED_ARRAY32(n,name) \
    typedef int32_t tInt32a_##name[n] __attribute__ ((aligned(32)));
#elif defined COMPILER_MSVC
  #define ALIGNED_ARRAY(n,name) \
    typedef __declspec(align(64)) int16_t tInt16a_##name[n];
  #define ALIGNED_ARRAY32(n,name) \
    typedef __declspec(align(64)) int32_t tInt32a_##name[n];
#endif

#ifdef __cplusplus
}
#endif

#endif /* __NATUREDSP_BASEBAND_H */
