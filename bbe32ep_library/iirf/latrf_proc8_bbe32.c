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

/* M=7...8 */

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
void latrf_process8( latrf_t   *          latr,
                     float32_t * restrict r,
               const float32_t *          x,
                     int                  N )
{
#if 0
        float32_t * restrict  delLine;
  const float32_t * restrict  coef;

  float32_t t0,t1,t2,t3;
  float32_t scale;

  int n;

  NASSERT((latr->M==7)||(latr->M==8));
  delLine = latr->delLine;
  coef    = latr->coef;
  scale   = latr->gain;

  for ( n=0; n<(N&~3); n+=4 )
  {
    // Scale the input sample.
    t0 = x[n+0]*scale;
    t1 = x[n+1]*scale;
    t2 = x[n+2]*scale;
    t3 = x[n+3]*scale;

    t0 -=          delLine[8-1]       * coef[8-1];
    t0 -=          delLine[8-2]       * coef[8-2];
    delLine[8-1] = delLine[8-2] +  t0 * coef[8-2];
    t0 -=          delLine[8-3]       * coef[8-3];
    delLine[8-2] = delLine[8-3] +  t0 * coef[8-3];
    t0 -=          delLine[8-4]       * coef[8-4];
    delLine[8-3] = delLine[8-4] +  t0 * coef[8-4];
    t1 -=          delLine[8-1]       * coef[8-1];
    t1 -=          delLine[8-2]       * coef[8-2];
    delLine[8-1] = delLine[8-2] +  t1 * coef[8-2];
    t1 -=          delLine[8-3]       * coef[8-3];
    delLine[8-2] = delLine[8-3] +  t1 * coef[8-3];
    t2 -=          delLine[8-1]       * coef[8-1];
    t2 -=          delLine[8-2]       * coef[8-2];
    delLine[8-1] = delLine[8-2] +  t2 * coef[8-2];
    t3 -=          delLine[8-1]       * coef[8-1];

    t0 -=          delLine[3+0]       * coef[3+0];
    delLine[3+1] = delLine[3+0] +  t0 * coef[3+0];
    t1 -=          delLine[3+1]       * coef[3+1];
    delLine[3+2] = delLine[3+1] +  t1 * coef[3+1];
    t2 -=          delLine[3+2]       * coef[3+2];
    delLine[3+3] = delLine[3+2] +  t2 * coef[3+2];
    t3 -=          delLine[3+3]       * coef[3+3];
    delLine[3+4] = delLine[3+3] +  t3 * coef[3+3];
    t0 -=          delLine[2+0]       * coef[2+0];
    delLine[2+1] = delLine[2+0] +  t0 * coef[2+0];
    t1 -=          delLine[2+1]       * coef[2+1];
    delLine[2+2] = delLine[2+1] +  t1 * coef[2+1];
    t2 -=          delLine[2+2]       * coef[2+2];
    delLine[2+3] = delLine[2+2] +  t2 * coef[2+2];
    t3 -=          delLine[2+3]       * coef[2+3];
    delLine[2+4] = delLine[2+3] +  t3 * coef[2+3];
    t0 -=          delLine[1+0]       * coef[1+0];
    delLine[1+1] = delLine[1+0] +  t0 * coef[1+0];
    t1 -=          delLine[1+1]       * coef[1+1];
    delLine[1+2] = delLine[1+1] +  t1 * coef[1+1];
    t2 -=          delLine[1+2]       * coef[1+2];
    delLine[1+3] = delLine[1+2] +  t2 * coef[1+2];
    t3 -=          delLine[1+3]       * coef[1+3];
    delLine[1+4] = delLine[1+3] +  t3 * coef[1+3];
    t0 -=          delLine[0+0]       * coef[0+0];
    delLine[0+1] = delLine[0+0] +  t0 * coef[0+0];
    t1 -=          delLine[0+1]       * coef[0+1];
    delLine[0+2] = delLine[0+1] +  t1 * coef[0+1];
    t2 -=          delLine[0+2]       * coef[0+2];
    delLine[0+3] = delLine[0+2] +  t2 * coef[0+2];
    t3 -=          delLine[0+3]       * coef[0+3];
    delLine[0+4] = delLine[0+3] +  t3 * coef[0+3];
    r[n+0] = t0;
    t1 -=        t0       * coef[0];
    r[n+1] = t1;
    t0 =         t0 +  t1 * coef[0];
    t2 -=        t0       * coef[1];
    t0 =         t0 +  t2 * coef[1];
    t2 -=        t1       * coef[0];
    r[n+2] = t2;
    t1 =         t1 +  t2 * coef[0];
    t3 -=        t0       * coef[2];
    delLine[3] = t0 +  t3 * coef[2];
    t3 -=        t1       * coef[1];
    delLine[2] = t1 +  t3 * coef[1];
    t3 -=        t2       * coef[0];
    r[n+3] = t3;
    delLine[1] = t2 +  t3 * coef[0];
    delLine[0] = t3;

  }
  if (N&3)
  {
      for ( ; n<N; n++ )
      {
        t0 = x[n]*scale;
        t0 -= delLine[8-1]*coef[8-1];
        t0 -=         delLine[6]       * coef[6];
        delLine[6+1]= delLine[6] +  t0 * coef[6];
        t0 -=         delLine[5]       * coef[5];
        delLine[5+1]= delLine[5] +  t0 * coef[5];
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
  }
#else
  const float32_t * restrict px;
        float32_t * restrict pr;
  const xb_vecN_2xf32 * restrict pcoef;
        xb_vecN_2xf32 * restrict pdelay;
  xb_vecN_2xf32 xin, rout, gain, zero;
  xb_vecN_2xf32 r0123456;
  xb_vecN_2xf32 t0, t1, t2, t3;
  xb_vecN_2xf32 d01234567;
  xb_vecN_2xf32 cf01234567;
  int i;

  NASSERT((latr->M==7) || (latr->M==8));

  zero = BBE_ZERON_2XF32();
  gain = latr->gain;
  gain = BBE_REPN_2XF32(gain, 0);
  pcoef  = (const xb_vecN_2xf32 *)latr->coef;
  pdelay = (      xb_vecN_2xf32 *)latr->delLine;
  px = x;
  pr = r;

  // Load coefficients
  cf01234567 = BBE_LVN_2XF32_I(pcoef, 0);
  // Load delay line elements
  d01234567 = BBE_LVN_2XF32_I(pdelay, 0);
  
  __Pragma("loop_count min=1");
  for ( i=0; i<N; i++ )
  {
    // Load and scale the input sample.
    BBE_LSN_2XF32_IP(xin, px, sizeof(float32_t));
    xin = BBE_REPN_2XF32(xin, 0);
    xin = BBE_MULN_2XF32(xin, gain);

    t0 = BBE_MULN_2XF32(d01234567, cf01234567);
    t1 = BBE_SELN_2XF32I(zero, t0, BBE_SELI_ROTATE_RIGHT_2);
    t0 = BBE_ADDN_2XF32(t0, t1);
    t1 = BBE_SELN_2XF32I(zero, t0, BBE_SELI_ROTATE_RIGHT_4);
    t2 = BBE_SELN_2XF32I(zero, t0, BBE_SELI_ROTATE_RIGHT_8);
    t3 = BBE_SELN_2XF32I(zero, t0, BBE_SELI_ROTATE_RIGHT_12);
    t2 = BBE_ADDN_2XF32(t2, t3);

    r0123456 = BBE_SUBN_2XF32(xin     , t0);
    r0123456 = BBE_SUBN_2XF32(r0123456, t1);
    r0123456 = BBE_SUBN_2XF32(r0123456, t2);
    
    // Update the delay line elements.
    BBE_MULAN_2XF32(d01234567, r0123456, cf01234567);
    d01234567 = BBE_SELN_2XF32I(d01234567, r0123456, BBE_SELI_PACK_2);
    // Store the output sample.
    rout = r0123456;
    BBE_SSN_2XF32_IP(rout, pr, sizeof(float32_t));
  }

  // Store updated delay line elements
  BBE_SVN_2XF32_I(d01234567, pdelay, 0);
#endif
} /* latrf_process8() */

#endif /* if HAVE_VFPU */
