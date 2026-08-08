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

/* M=3...4 */

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
void latcf_process4( latcf_t       *          latc,
                     complex_float * restrict r,
               const complex_float *          x,
                     int                      N )
{
#ifndef complex_float_
  #ifdef COMPILER_XTENSA
    #define float32_t_ float32_t
    #define complex_float_ complex_float
  #else
    #define float32_t_     xtfloat
    #define complex_float_ xtcomplexfloat
  #endif
#endif

  const complex_float_ * restrict px;
        complex_float_ * restrict pr;
        xb_vecN_4xcf32 * restrict pdelay;
  const float32_t_     * restrict pcoef;
  xb_vecN_4xcf32 t;
  xb_vecN_2xf32 xin, rout, gain;
  xb_vecN_2xf32 d0, d1, d2, d3, cf0, cf1, cf2, cf3, r0, r1, r2;
  int i;

  NASSERT((latc->M==3) || (latc->M==4));

  pcoef  = (const xtfloat *)latc->coef;
  pdelay = (xb_vecN_4xcf32 *)latc->delLine;
  px = (const complex_float_ *)x;
  pr = (      complex_float_ *)r;
  gain = latc->gain;
  cf3 = BBE_LSN_2XF32_I(pcoef, 3*sizeof(float32_t));
  cf2 = BBE_LSN_2XF32_I(pcoef, 2*sizeof(float32_t));
  cf1 = BBE_LSN_2XF32_I(pcoef, 1*sizeof(float32_t));
  cf0 = BBE_LSN_2XF32_I(pcoef, 0*sizeof(float32_t));
  cf3 = BBE_REPN_2XF32(cf3, 0);
  cf2 = BBE_REPN_2XF32(cf2, 0);
  cf1 = BBE_REPN_2XF32(cf1, 0);
  cf0 = BBE_REPN_2XF32(cf0, 0);
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
    r2 = BBE_MULN_2XF32(xin, gain);
    // Update the delay line elements.
    BBE_MULSN_2XF32(r2, d3, cf3);
    d3 = d2;  
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
  
  d3 = BBE_SELN_2XF32I(d3, d2, BBE_SELI_PACK_4);
  d3 = BBE_SELN_2XF32I(d3, d1, BBE_SELI_PACK_4);
  d3 = BBE_SELN_2XF32I(d3, d0, BBE_SELI_PACK_4);
  t = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(d3));
  BBE_SVN_4XCF32_I(t, pdelay, 0);
} /* latcf_process4() */

#endif /* if HAVE_VFPU */
