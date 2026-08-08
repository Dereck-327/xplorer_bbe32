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
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */

/*
    NatureDSP_Baseband library. IIR part
    Lattice complex block IIR w/ real coefficients, floating point
    C code optimized for BBE32
    Integrit, 2006-2016
*/

/* M=7...8 */

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Filter instance and processing functions. */
#include "latcf_common.h"

#if HAVE_VFPU

/*-------------------------------------------------------------------------
  Pass a block of complex input samples to the filter and return the 
  response samples.
  Input:
    _latcf  Filter handle
    N       Input/output signal chunk size, in complex samples
    x[N]    Complex input samples.
  Output:
    r[N]    Complex output samples.
-------------------------------------------------------------------------*/
void latcf_process8( latcf_t       *          latc,
                     complex_float * restrict r,
               const complex_float *          x,
                     int                      N )
{
#ifndef complex_float_
  #ifdef COMPILER_XTENSA
    #define complex_float_ complex_float
  #else
    #define complex_float_ xtcomplexfloat
  #endif
#endif

  const complex_float_ * restrict px;
        complex_float_ * restrict pr;
        xb_vecN_2xf32 * restrict pdelay;
  const xb_vecN_2xf32 * restrict pcoef;
  xb_vecN_4xcf32 t;
  xb_vecN_2xf32 xin, rout, gain, zero;
  xb_vecN_2xf32 r0123, r456;
  xb_vecN_2xf32 t00, t01, t10, t11;
  xb_vecN_2xf32 d0123, d4567;
  xb_vecN_2xf32 cf0123, cf4567;
  int i;

  NASSERT((latc->M==7) || (latc->M==8));

  zero = BBE_ZERON_2XF32();
  gain = latc->gain;
  gain = BBE_REPN_2XF32(gain, 0);
  pcoef  = (const xb_vecN_2xf32 *)latc->coef;
  pdelay = (      xb_vecN_2xf32 *)latc->delLine;
  px = (const complex_float_ *)x;
  pr = (      complex_float_ *)r;

  // Load coefficients
  t00 = BBE_LVN_2XF32_I(pcoef, 0);
  BBE_DSELN_2XF32I(cf4567, cf0123, t00, t00, BBE_DSELI_INTERLEAVE_2);
  // Load delay line elements
  d0123 = BBE_LVN_2XF32_I(pdelay, 0);
  d4567 = BBE_LVN_2XF32_I(pdelay, 4*sizeof(complex_float));
  
  __Pragma("loop_count min=1");
  for ( i=0; i<N; i++ )
  {
    // Load and scale the input sample.
    BBE_LSN_4XCF32_IP(t, px, sizeof(complex_float));
    t = BBE_REPN_4XCF32(t, 0);
    xin = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(t));
    xin = BBE_MULN_2XF32(xin, gain);

    t01 = BBE_MULN_2XF32(d4567, cf4567);
    t00 = BBE_MULN_2XF32(d0123, cf0123);

    t11 = BBE_SELN_2XF32I(zero, t01, BBE_SELI_ROTATE_RIGHT_4);
    t10 = BBE_SELN_2XF32I(t01 , t00, BBE_SELI_ROTATE_RIGHT_4);
    t01 = BBE_ADDN_2XF32(t01, t11);
    t00 = BBE_ADDN_2XF32(t00, t10);

    t11 = BBE_SELN_2XF32I(zero, t01, BBE_SELI_ROTATE_RIGHT_8);
    t10 = BBE_SELN_2XF32I(t01 , t00, BBE_SELI_ROTATE_RIGHT_8);
    t01 = BBE_ADDN_2XF32(t01, t11);
    t00 = BBE_ADDN_2XF32(t00, t10);

    r456  = BBE_SUBN_2XF32(xin  , t01);
    r0123 = BBE_SUBN_2XF32(xin  , t00);
    r0123 = BBE_SUBN_2XF32(r0123, t01);
    
    // Update the delay line elements.
    BBE_MULAN_2XF32(d4567, r456 , cf4567);
    BBE_MULAN_2XF32(d0123, r0123, cf0123);
    d4567 = BBE_SELN_2XF32I(d4567, d0123, BBE_SELI_ROTATE_LEFT_4);
    d0123 = BBE_SELN_2XF32I(d0123, r0123, BBE_SELI_PACK_4);
    // Store the output sample.
    rout = r0123;
    t = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(rout));
    BBE_SSN_4XCF32_IP(t, pr, sizeof(complex_float));
  }

  // Store updated delay line elements
  BBE_SVN_2XF32_I(d0123, pdelay, 0);
  BBE_SVN_2XF32_I(d4567, pdelay, 4*sizeof(complex_float));
} /* latcf_process8() */

#endif /* if HAVE_VFPU */
