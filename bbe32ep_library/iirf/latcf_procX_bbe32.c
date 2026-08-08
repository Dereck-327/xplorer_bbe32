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

/* M>8 */

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
void latcf_processX( latcf_t       *          latc,
                     complex_float * restrict r,
               const complex_float *          x,
                     int                      N )
#if 0
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
  const complex_float_ * restrict d_ld;
        complex_float_ * restrict d_st;
  const float32_t_     * restrict coef;
  xb_vecN_4xcf32 tcf;
  xb_vecN_2xf32 gain, accr, xin;
  xb_vecN_2xf32 cf, d0, d1;
  int     n, m, M;

  M = latc->M;
  px   = (const complex_float_ *)x;
  pr   = (      complex_float_ *)r;
  gain = latc->gain;
  
  for ( n=0; n<N; n++ )
  {
    d_ld = (const complex_float_ *)latc->delLine;
    d_st = (      complex_float_ *)latc->delLine;
    coef = (const float32_t_     *)latc->coef + M-1;
    // Scale the input sample.
    BBE_LSN_4XCF32_IP(tcf, px, sizeof(complex_float));
    xin = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(tcf));
    accr = BBE_MULN_2XF32(xin, gain);
    
    BBE_LSN_4XCF32_IP(tcf, d_ld, sizeof(complex_float));
    d1 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(tcf));
    cf = *coef--;

    BBE_MULSN_2XF32(accr, d1, cf);

    __Pragma("loop_count min=4");
    for ( m=M-2; m>=0; m-- )
    {
      BBE_LSN_4XCF32_IP(tcf, d_ld, sizeof(complex_float));
      d0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(tcf));
      cf = *coef--;
      
      BBE_MULSN_2XF32(accr, d0, cf);
      // Update the (m+1)-th delay line element.
      d1 = d0;
      BBE_MULAN_2XF32(d1, accr, cf);
      tcf = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(d1));
      BBE_SSN_4XCF32_IP(tcf, d_st, sizeof(complex_float));
    }
    
    // Update the first delay line element with the resulting sample
    tcf = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(accr));
    BBE_SSN_4XCF32_IP(tcf, d_st, sizeof(complex_float));
    // Make the output sample.
    tcf = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(accr));
    BBE_SSN_4XCF32_IP(tcf, pr, sizeof(complex_float));
  }
} /* latcf_processX() */
#else
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
  const complex_float_ * restrict d_ld;
        complex_float_ * restrict d_st;
  const float32_t_     * restrict coef;
  xb_vecN_4xcf32 tcplx;
  xb_vecN_2xf32 gain, accr0, accr1, xin0, xin1;
  xb_vecN_2xf32 cf0, cf1, cf2, d0, d1, d2;
  int     n, m, M;

  NASSERT((latc->M)>8);

  M = latc->M;
  px   = (const complex_float_ *)x;
  pr   = (      complex_float_ *)r;
  gain = latc->gain;
  
  for ( n=0; n<(N>>1); n++ )
  {
    d_ld = (const complex_float_ *)latc->delLine;
    d_st = (      complex_float_ *)latc->delLine;
    coef = (const float32_t_     *)latc->coef + M-1;
    // Load and scale 2 input samples.
    BBE_LSN_4XCF32_IP(tcplx, px, sizeof(complex_float));  xin0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(tcplx));
    BBE_LSN_4XCF32_IP(tcplx, px, sizeof(complex_float));  xin1 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(tcplx));

    accr0 = BBE_MULN_2XF32(xin0, gain);
    accr1 = BBE_MULN_2XF32(xin1, gain);
    
    BBE_LSN_4XCF32_IP(tcplx, d_ld, sizeof(complex_float));  d2 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(tcplx));
    BBE_LSN_4XCF32_IP(tcplx, d_ld, sizeof(complex_float));  d1 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(tcplx));
    BBE_LSN_2XF32_IP(cf2, coef, -(int)sizeof(float32_t));   cf2 = BBE_REPN_2XF32(cf2, 0);
    BBE_LSN_2XF32_IP(cf1, coef, -(int)sizeof(float32_t));   cf1 = BBE_REPN_2XF32(cf1, 0);

    BBE_MULSN_2XF32(accr0, d2, cf2);
    BBE_MULSN_2XF32(accr0, d1, cf1);

    d2 = d1;
    BBE_MULAN_2XF32(d2, accr0, cf1);

    BBE_MULSN_2XF32(accr1, d2, cf2);
    
    __Pragma("loop_count min=4");
    for ( m=M-3; m>=0; m-- )
    {
      BBE_LSN_4XCF32_IP(tcplx, d_ld, sizeof(complex_float));  d0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(tcplx));
      BBE_LSN_2XF32_IP(cf0, coef, -(int)sizeof(float32_t));   cf0 = BBE_REPN_2XF32(cf0, 0);
      
      BBE_MULSN_2XF32(accr0, d0, cf0);
      d1 = d0;
      BBE_MULAN_2XF32(d1, accr0, cf0);
      
      BBE_MULSN_2XF32(accr1, d1, cf1);
      d2 = d1;
      BBE_MULAN_2XF32(d2, accr1, cf1);
      
      tcplx = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(d2));  BBE_SSN_4XCF32_IP(tcplx, d_st, sizeof(complex_float));

      cf1 = cf0;
    }
    
    d1 = accr0;
    BBE_MULSN_2XF32(accr1, d1, cf1);
    // Update the first delay line element with the resulting sample
    BBE_MULAN_2XF32(d1, accr1, cf1);
    d0 = accr1;
    tcplx = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(d1));  BBE_SSN_4XCF32_IP(tcplx, d_st, sizeof(complex_float));
    tcplx = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(d0));  BBE_SSN_4XCF32_IP(tcplx, d_st, sizeof(complex_float));

    // Make the output sample.
    tcplx = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(accr0));  BBE_SSN_4XCF32_IP(tcplx, pr, sizeof(complex_float));
    tcplx = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(accr1));  BBE_SSN_4XCF32_IP(tcplx, pr, sizeof(complex_float));
  }
  if (N&1)
  {
    d_ld = (const complex_float_ *)latc->delLine;
    d_st = (      complex_float_ *)latc->delLine;
    coef = (const float32_t_     *)latc->coef + M-1;
    // Scale the input sample.
    BBE_LSN_4XCF32_IP(tcplx, px, sizeof(complex_float));
    xin0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(tcplx));
    accr0 = BBE_MULN_2XF32(xin0, gain);
    
    BBE_LSN_4XCF32_IP(tcplx, d_ld, sizeof(complex_float));
    d1 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(tcplx));
    BBE_LSN_2XF32_IP(cf1, coef, -(int)sizeof(float32_t));
    cf1 = BBE_REPN_2XF32(cf1, 0);

    BBE_MULSN_2XF32(accr0, d1, cf1);

    for ( m=M-2; m>=0; m-- )
    {
      BBE_LSN_4XCF32_IP(tcplx, d_ld, sizeof(complex_float));
      d0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(tcplx));
      BBE_LSN_2XF32_IP(cf0, coef, -(int)sizeof(float32_t));
      cf0 = BBE_REPN_2XF32(cf0, 0);
      
      BBE_MULSN_2XF32(accr0, d0, cf0);
      // Update the (m+1)-th delay line element.
      d1 = d0;
      BBE_MULAN_2XF32(d1, accr0, cf0);
      tcplx = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(d1));
      BBE_SSN_4XCF32_IP(tcplx, d_st, sizeof(complex_float));
    }
    
    // Update the first delay line element with the resulting sample
    tcplx = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(accr0));
    BBE_SSN_4XCF32_IP(tcplx, d_st, sizeof(complex_float));
    // Make the output sample.
    tcplx = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(accr0));
    BBE_SSN_4XCF32_IP(tcplx, pr, sizeof(complex_float));
  }
} /* latcf_processX() */
#endif

#endif /* if HAVE_VFPU */
