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

/* M=1...2 */

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
void latcf_process2( latcf_t       *          latc,
                     complex_float * restrict r,
               const complex_float *          x,
                     int                      N )
{
#if 1

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
        complex_float_ * restrict pdelay;
  const float32_t_     * restrict pcoef;
  xb_vecN_4xcf32 t;
  xb_vecN_2xf32 xin, rout, gain, d0, d1, cf0, cf1, r0;
  int i;

  NASSERT((latc->M==1) || (latc->M==2));

  pcoef = (const xtfloat *)latc->coef;
  pdelay = (complex_float_ *)latc->delLine;
  px = (const complex_float_ *)x;
  pr = (      complex_float_ *)r;
  gain = latc->gain;
  cf0 = BBE_LSN_2XF32_I(pcoef, 0);
  cf1 = BBE_LSN_2XF32_I(pcoef, sizeof(float32_t));
  cf0 = BBE_REPN_2XF32(cf0, 0);
  cf1 = BBE_REPN_2XF32(cf1, 0);
  t = pdelay[0]; d0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(t));
  t = pdelay[1]; d1 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(t));
  
  for ( i=0; i<N; i++ )
  {
    BBE_LSN_4XCF32_IP(t, px, sizeof(complex_float));
    xin = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(t));
    // Scale the input sample.
    r0 = BBE_MULN_2XF32(xin, gain);
    BBE_MULSN_2XF32(r0, d0, cf0);
    BBE_MULSN_2XF32(r0, d1, cf1);
    
    // Update the delay line elements.
    d1 = d0;
    BBE_MULAN_2XF32(d1, r0, cf0);
    d0 = r0;
    // Store the output sample.
    rout = r0;
    t = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(rout));
    BBE_SSN_4XCF32_IP(t, pr, sizeof(complex_float));
  }
  
  t = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(d0)); pdelay[0] = t;
  t = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(d1)); pdelay[1] = t;
#else
    float32_t d0,d1,c0,c1;
          float32_t * restrict delLine;
    const float32_t * restrict coef;

    float32_t t0;
    float32_t scale;

    int n;
    
    NASSERT((latr->M==1) || (latr->M==2));

    delLine = latr->delLine;
    coef    = latr->coef;
    scale   = latr->gain;
    d0=delLine[0];
    d1=delLine[1];
    c0=coef[0];
    c1=coef[1];

    for ( n=0; n<N; n++ )
    {
        t0 = x[n]*scale;
        t0 -= d0 * c0;
        t0 -= d1 * c1;
        d1 = d0 + t0 * c0;
        d0 = t0;
        r[n] = t0;
    }
    delLine[0]=d0;
    delLine[1]=d1;
#endif
} /* latcf_process2() */

#endif /* if HAVE_VFPU */
