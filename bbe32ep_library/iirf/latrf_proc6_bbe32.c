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
    Lattice real block IIR, floating point
    C code optimized for BBE32
    Integrit, 2006-2016
*/

/* M=5...6 */

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Filter instance and processing functions. */
#include "latrf_common.h"

#if HAVE_VFPU

/*-------------------------------------------------------------------------
   Pass a block of input samples to the filter 
   Input:
     _latr  Filter handle
     N      Input/output signal chunk size, in samples
     x[N]   Input samples
   Output:
     r[N]   Output samples
   Restrictions:
     None 
-------------------------------------------------------------------------*/
void latrf_process6( latrf_t   *          latr,
                     float32_t * restrict r,
               const float32_t *          x,
                     int                  N )
{
#if 1
  const float32_t     * restrict px;
        float32_t     * restrict pr;
  const xb_vecN_2xf32 * restrict pcoef;
        xb_vecN_2xf32 * restrict pdelay;
  xb_vecN_2xf32 xin, rout, gain;
  xb_vecN_2xf32 r0, r1, r2, r3, r4;
  xb_vecN_2xf32 d0, d1, d2, d3, d4, d5;
  xb_vecN_2xf32 cf0, cf1, cf2, cf3, cf4, cf5;
  xb_vecN_2xf32 D, CF;
  int i;

  NASSERT((latr->M==5) || (latr->M==6));

  gain = latr->gain;
  pcoef = (const xb_vecN_2xf32 *)latr->coef;
  pdelay = (     xb_vecN_2xf32 *)latr->delLine;
  px = x;
  pr = r;
  CF = BBE_LVN_2XF32_I(pcoef, 0);
  cf5 = BBE_REPN_2XF32(CF, 5);
  cf4 = BBE_REPN_2XF32(CF, 4);
  cf3 = BBE_REPN_2XF32(CF, 3);
  cf2 = BBE_REPN_2XF32(CF, 2);
  cf1 = BBE_REPN_2XF32(CF, 1);
  cf0 = BBE_REPN_2XF32(CF, 0);
  D  = BBE_LVN_2XF32_I(pdelay, 0);
  d5 = BBE_REPN_2XF32(D, 5);
  d4 = BBE_REPN_2XF32(D, 4);
  d3 = BBE_REPN_2XF32(D, 3);
  d2 = BBE_REPN_2XF32(D, 2);
  d1 = BBE_REPN_2XF32(D, 1);
  d0 = BBE_REPN_2XF32(D, 0);
  
  __Pragma("loop_count min=1");
  for ( i=0; i<N; i++ )
  {
    BBE_LSN_2XF32_IP(xin, px, sizeof(float32_t));
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
    BBE_SSN_2XF32_IP(rout, pr, sizeof(float32_t));
  }

  D = BBE_SELN_2XF32I(D, d5, BBE_SELI_PACK_2);
  D = BBE_SELN_2XF32I(D, d4, BBE_SELI_PACK_2);
  D = BBE_SELN_2XF32I(D, d3, BBE_SELI_PACK_2);
  D = BBE_SELN_2XF32I(D, d2, BBE_SELI_PACK_2);
  D = BBE_SELN_2XF32I(D, d1, BBE_SELI_PACK_2);
  D = BBE_SELN_2XF32I(D, d0, BBE_SELI_PACK_2);
  BBE_SVN_2XF32_I(D, pdelay, 0);
#elif 1
          float32_t * restrict  delLine;
    const float32_t * restrict  coef;

    float32_t t0;
    float32_t scale;

    int n;
    
    NASSERT((latr->M==5) || (latr->M==6));
    delLine = latr->delLine;
    coef    = latr->coef;
    scale   = latr->gain;
    __Pragma("loop_count min=1");
    for (n=0; n<N; n++)
    {
        t0 = x[n]*scale;
        t0 -=         delLine[5]*coef[5];
        t0 -=         delLine[4]       * coef[4];
        delLine[4+1]= delLine[4] +  t0 * coef[4];
        t0 -=         delLine[3]       * coef[3];
        delLine[3+1]= delLine[3] +  t0 * coef[3];
        t0 -=         delLine[2]       * coef[2];
        delLine[2+1]= delLine[2] +  t0 * coef[2];
        t0 -=         delLine[1]       * coef[1];
        delLine[1+1]= delLine[1] +  t0 * coef[1];
        t0 -=         delLine[0]       * coef[0];
        delLine[0+1]= delLine[0] +  t0 * coef[0];
        delLine[0] = t0;
        r[n] = t0;
    }
#else
  const float32_t * restrict px;
        float32_t * restrict pr;
  const xb_vecN_2xf32 * restrict pcoef;
        xb_vecN_2xf32 * restrict pdelay;
  xb_vecN_2xf32 xin, rout, gain, zero;
  xb_vecN_2xf32 r01234;
  xb_vecN_2xf32 d012345, d0, d1, d2, d3, d4, d5;
  xb_vecN_2xf32 cf012345, cf0, cf1, cf2, cf3, cf4, cf5;
  int i;

  NASSERT((latr->M==5) || (latr->M==6));

  zero = BBE_ZERON_2XF32();
  gain = latr->gain;
  gain = BBE_REPN_2XF32(gain, 0);
  pcoef  = (const xb_vecN_2xf32 *)latr->coef;
  pdelay = (      xb_vecN_2xf32 *)latr->delLine;
  px = x;
  pr = r;

  // Load coefficients
  cf012345 = BBE_LVN_2XF32_I(pcoef, 0);
  cf0 = BBE_REPN_2XF32(cf012345, 0);
  cf1 = BBE_REPN_2XF32(cf012345, 1);
  cf2 = BBE_REPN_2XF32(cf012345, 2);
  cf3 = BBE_REPN_2XF32(cf012345, 3);
  cf4 = BBE_REPN_2XF32(cf012345, 4);
  cf5 = BBE_REPN_2XF32(cf012345, 5);
  cf0 = BBE_SELN_2XF32I(zero, cf0, BBE_SELI_ROTATE_LEFT_2);
  cf1 = BBE_SELN_2XF32I(zero, cf1, BBE_SELI_ROTATE_LEFT_4);
  cf2 = BBE_SELN_2XF32I(zero, cf2, BBE_SELI_ROTATE_LEFT_6);
  cf3 = BBE_SELN_2XF32I(zero, cf3, BBE_SELI_ROTATE_LEFT_8);
  cf4 = BBE_SELN_2XF32I(zero, cf4, BBE_SELI_ROTATE_LEFT_10);
  cf5 = BBE_SELN_2XF32I(zero, cf5, BBE_SELI_ROTATE_LEFT_10);
  // Load delay line elements
  d012345 = BBE_LVN_2XF32_I(pdelay, 0);
  d0 = BBE_REPN_2XF32(d012345, 0);
  d1 = BBE_REPN_2XF32(d012345, 1);
  d2 = BBE_REPN_2XF32(d012345, 2);
  d3 = BBE_REPN_2XF32(d012345, 3);
  d4 = BBE_REPN_2XF32(d012345, 4);
  d5 = BBE_REPN_2XF32(d012345, 5);
  
  __Pragma("loop_count min=1");
  for ( i=0; i<N; i++ )
  {
    BBE_LSN_2XF32_IP(xin, px, sizeof(float32_t));
    xin = BBE_REPN_2XF32(xin, 0);
    // Scale the input sample.
    r01234 = BBE_MULN_2XF32(xin, gain);
    BBE_MULSN_2XF32(r01234, d0, cf0);
    BBE_MULSN_2XF32(r01234, d1, cf1);
    BBE_MULSN_2XF32(r01234, d2, cf2);
    BBE_MULSN_2XF32(r01234, d3, cf3);
    BBE_MULSN_2XF32(r01234, d4, cf4);
    BBE_MULSN_2XF32(r01234, d5, cf5);
    
    // Update the delay line elements.
    BBE_MULAN_2XF32(d012345, r01234, cf012345);
    d0 = r01234;
    d1 = BBE_REPN_2XF32(d012345, 0);
    d2 = BBE_REPN_2XF32(d012345, 1);
    d3 = BBE_REPN_2XF32(d012345, 2);
    d4 = BBE_REPN_2XF32(d012345, 3);
    d5 = BBE_REPN_2XF32(d012345, 4);
    d012345 = BBE_SELN_2XF32I(d012345, r01234, BBE_SELI_PACK_2);
    // Store the output sample.
    rout = r01234;
    BBE_SSN_2XF32_IP(rout, pr, sizeof(float32_t));
  }

  // Store updated delay line elements
  BBE_SVN_2XF32_I(d012345, pdelay, 0);
#endif
} /* latrf_process6() */

#endif /* if HAVE_VFPU */
