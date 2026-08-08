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

/* M=5...6 */

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
void latcf_process6( latcf_t       *          latc,
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
        xb_vecN_4xcf32 * restrict pdelay;
  const xb_vecN_2xf32  * restrict pcoef;
  xb_vecN_4xcf32 t;
  xb_vecN_2xf32 xin, rout, gain;
  xb_vecN_2xf32 d0, d1, d2, d3, d4, d5,
                cf0, cf1, cf2, cf3, cf4, cf5,
                r0, r1, r2, r3, r4;
  int i;

  NASSERT((latc->M==5) || (latc->M==6));

  pcoef  = (const xb_vecN_2xf32 *)latc->coef;
  pdelay = (xb_vecN_4xcf32 *)latc->delLine;
  px = (const complex_float_ *)x;
  pr = (      complex_float_ *)r;
  gain = latc->gain;
  cf0 = BBE_LVN_2XF32_I(pcoef, 0);
  cf5 = BBE_REPN_2XF32(cf0, 5);
  cf4 = BBE_REPN_2XF32(cf0, 4);
  cf3 = BBE_REPN_2XF32(cf0, 3);
  cf2 = BBE_REPN_2XF32(cf0, 2);
  cf1 = BBE_REPN_2XF32(cf0, 1);
  cf0 = BBE_REPN_2XF32(cf0, 0);
  t = BBE_LVN_4XCF32_I(pdelay, BBE_SIMD_WIDTH/4*sizeof(complex_float));
  d5 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(BBE_REPN_4XCF32(t, 1)));
  d4 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(BBE_REPN_4XCF32(t, 0)));
  t = BBE_LVN_4XCF32_I(pdelay, 0);
  d3 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(BBE_REPN_4XCF32(t, 3)));
  d2 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(BBE_REPN_4XCF32(t, 2)));
  d1 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(BBE_REPN_4XCF32(t, 1)));
  d0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(BBE_REPN_4XCF32(t, 0)));
  
  __Pragma("loop_count min=1");
  for ( i=0; i<N; i++ )
  {
    BBE_LSN_4XCF32_IP(t, px, sizeof(complex_float));
    xin = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(t));
    // Scale the input sample.
    r4 = BBE_MULN_2XF32(xin, gain);
    // Update the delay line elements.
    BBE_MULSN_2XF32(r4, d5, cf5);
    d5 = d4;
    BBE_MULSN_2XF32(r4, d4, cf4);
    BBE_MULAN_2XF32(d5, r4, cf4);
    d4 = d3;
    r3 = r4;
    BBE_MULSN_2XF32(r3, d3, cf3);
    BBE_MULAN_2XF32(d4, r3, cf3);
    d3 = d2;
    r2 = r3;
    BBE_MULSN_2XF32(r2, d2, cf2);
    BBE_MULAN_2XF32(d3, r2, cf2);
    d2 = d1;
    r1 = r2;
    BBE_MULSN_2XF32(r1, d1, cf1);
    BBE_MULAN_2XF32(d2, r1, cf1);
    d1 = d0;
    r0 = r1;
    BBE_MULSN_2XF32(r0, d0, cf0);
    BBE_MULAN_2XF32(d1, r0, cf0);
    d0 = r0;
    
    // Store the output sample.
    rout = r0;
    t = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(rout));
    BBE_SSN_4XCF32_IP(t, pr, sizeof(complex_float));
  }
  
  d5 = BBE_SELN_2XF32I(d5, d4, BBE_SELI_PACK_4);
  d3 = BBE_SELN_2XF32I(d3, d2, BBE_SELI_PACK_4);
  d3 = BBE_SELN_2XF32I(d3, d1, BBE_SELI_PACK_4);
  d3 = BBE_SELN_2XF32I(d3, d0, BBE_SELI_PACK_4);
  t = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(d5));
  BBE_SVN_4XCF32_I(t, pdelay, BBE_SIMD_WIDTH/4*sizeof(complex_float));
  t = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(d3));
  BBE_SVN_4XCF32_I(t, pdelay, 0);
} /* latcf_process6() */

#endif /* if HAVE_VFPU */
