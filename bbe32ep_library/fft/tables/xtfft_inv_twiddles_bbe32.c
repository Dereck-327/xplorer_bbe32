/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
/*          Copyright (C) 2009-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_fft.h"
#include "common.h"
#include "fft_tw.h"

ALIGN(32) const short __xtfft_inv_twiddles[32] = {
      0,   28378,  /* radix 3 inverse  (0,  0.8660254037844385) */
  -9630,   15582,  /* radix 5 forward ( -0.29389262614623656, 0.47552825814757676) */
  18318,       0,  /* radix 5 forward/inverse (0.55901699437494751, 0) */
   9630,   15582,  /* radix 5 inverse (  0.29389262614623656,  0.47552825814757676) */
      0,   32767,  /* radix 8 inverse (0, 1) */
  23170,   23170,  /* radix 8 inverse (0.70710678118654757,  0.70710678118654757) */
      0,       0,
      0,       0,
      0,       0,
      0,       0,
      0,       0,
      0,       0,
      0,       0,
      0,       0,
      0,       0,
      0,       0
};
