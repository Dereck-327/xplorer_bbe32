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
    Lattice real block IIR, floating point, streaming version
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* M>8 */

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Filter instance and processing functions. */
#include "latrsf_common.h"

#if HAVE_VFPU

/*---------------------------------------------------------------------------
  Pass a block of input samples to the filter 
  Input:
    latrs   Filter handle
    N       Input/output signal chunk size, in samples (per stream)
    x[N*L]  Input samples
  Output:
    r[N*L]  Output samples

  Restrictions:
    x,r     Must be aligned on 32-byte boundary
    L       Must be a multiple of 8
---------------------------------------------------------------------------*/
void latrsf_processX( latrsf_t  *          latrs,
                      float32_t * restrict r,
                const float32_t *          x,
                      int                  N )
{
  const xb_vecN_2xf32 * restrict px;
        xb_vecN_2xf32 * restrict pr;
  const xb_vecN_2xf32 * restrict d_ld;
        xb_vecN_2xf32 * restrict d_st;
  const float32_t     * restrict coef;
  xb_vecN_2xf32 gain, accr0, accr1, xin0, xin1;
  xb_vecN_2xf32 cf0, cf1, cf2, d0, d1, d2;
  int     n, m, l, M, L;

  NASSERT((latrs->M)>8);
  
  L = latrs->L;
  M = latrs->M;
  gain = latrs->gain;
  
  for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++ )
  {
    px = (const xb_vecN_2xf32 *)(x + l*BBE_SIMD_WIDTH/2);
    pr = (      xb_vecN_2xf32 *)(r + l*BBE_SIMD_WIDTH/2);
    for ( n=0; n<(N>>1); n++ )
    {
      d_ld = (const xb_vecN_2xf32 *)(latrs->delLine + l*BBE_SIMD_WIDTH/2*M);
      d_st = (      xb_vecN_2xf32 *)(latrs->delLine + l*BBE_SIMD_WIDTH/2*M);
      coef = (const float32_t     *)(latrs->coef + M-1);
      // Load and scale 2 input samples.
      BBE_LVN_2XF32_XP(xin0, px, L*sizeof(float32_t));
      BBE_LVN_2XF32_XP(xin1, px, L*sizeof(float32_t));

      accr0 = BBE_MULN_2XF32(xin0, gain);
      accr1 = BBE_MULN_2XF32(xin1, gain);
    
      BBE_LVN_2XF32_IP(d2, d_ld, BBE_SIMD_WIDTH/2*sizeof(float32_t));
      BBE_LVN_2XF32_IP(d1, d_ld, BBE_SIMD_WIDTH/2*sizeof(float32_t));
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
        BBE_LVN_2XF32_IP(d0, d_ld, BBE_SIMD_WIDTH/2*sizeof(float32_t));
        BBE_LSN_2XF32_IP(cf0, coef, -(int)sizeof(float32_t));   cf0 = BBE_REPN_2XF32(cf0, 0);
      
        BBE_MULSN_2XF32(accr0, d0, cf0);
        d1 = d0;
        BBE_MULAN_2XF32(d1, accr0, cf0);
      
        BBE_MULSN_2XF32(accr1, d1, cf1);
        d2 = d1;
        BBE_MULAN_2XF32(d2, accr1, cf1);
      
        BBE_SVN_2XF32_IP(d2, d_st, BBE_SIMD_WIDTH/2*sizeof(float32_t));

        cf1 = cf0;
      }
    
      d1 = accr0;
      BBE_MULSN_2XF32(accr1, d1, cf1);
      // Update the first delay line element with the resulting sample
      BBE_MULAN_2XF32(d1, accr1, cf1);
      d0 = accr1;
      BBE_SVN_2XF32_IP(d1, d_st, BBE_SIMD_WIDTH/2*sizeof(float32_t));
      BBE_SVN_2XF32_IP(d0, d_st, BBE_SIMD_WIDTH/2*sizeof(float32_t));

      // Make the output sample.
      BBE_SVN_2XF32_XP(accr0, pr, L*sizeof(float32_t));
      BBE_SVN_2XF32_XP(accr1, pr, L*sizeof(float32_t));
    }
    if (N&1)
    {
      d_ld = (const xb_vecN_2xf32 *)(latrs->delLine + l*BBE_SIMD_WIDTH/2*M);
      d_st = (      xb_vecN_2xf32 *)(latrs->delLine + l*BBE_SIMD_WIDTH/2*M);
      coef = (const float32_t     *)(latrs->coef + M-1);
      // Scale the input sample.
      BBE_LVN_2XF32_XP(xin0, px, L*sizeof(float32_t));
      accr0 = BBE_MULN_2XF32(xin0, gain);
    
      BBE_LVN_2XF32_IP(d1, d_ld, BBE_SIMD_WIDTH/2*sizeof(float32_t));
      BBE_LSN_2XF32_IP(cf1, coef, -(int)sizeof(float32_t));
      cf1 = BBE_REPN_2XF32(cf1, 0);

      BBE_MULSN_2XF32(accr0, d1, cf1);

      for ( m=M-2; m>=0; m-- )
      {
        BBE_LVN_2XF32_IP(d0, d_ld, BBE_SIMD_WIDTH/2*sizeof(float32_t));
        BBE_LSN_2XF32_IP(cf0, coef, -(int)sizeof(float32_t));
        cf0 = BBE_REPN_2XF32(cf0, 0);
      
        BBE_MULSN_2XF32(accr0, d0, cf0);
        // Update the (m+1)-th delay line element.
        d1 = d0;
        BBE_MULAN_2XF32(d1, accr0, cf0);
        BBE_SVN_2XF32_IP(d1, d_st, BBE_SIMD_WIDTH/2*sizeof(float32_t));
      }
    
      // Update the first delay line element with the resulting sample
      BBE_SVN_2XF32_IP(accr0, d_st, BBE_SIMD_WIDTH/2*sizeof(float32_t));
      // Make the output sample.
      BBE_SVN_2XF32_XP(accr0, pr, L*sizeof(float32_t));
    }
  }
} /* latrsf_processX() */

#endif /* if HAVE_VFPU */
