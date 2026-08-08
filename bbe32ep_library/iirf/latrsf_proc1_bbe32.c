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

/* M=1 */

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
void latrsf_process1( latrsf_t  *          latrs,
                      float32_t * restrict r,
                const float32_t *          x,
                      int                  N )
{
    const xb_vecN_2xf32 * restrict px;
          xb_vecN_2xf32 * restrict pr;
          xb_vecN_2xf32 * restrict delLine;
    const float32_t     * restrict coef;
    xb_vecN_2xf32 Xin,Rout,DL,CF;

    xb_vecN_2xf32 ACCR;
    xb_vecN_2xf32 Gain;

    int n, l, L;
    
    NASSERT(latrs->M==1);
    
    L = latrs->L;
    px = (const xb_vecN_2xf32 *)x;
    pr = (      xb_vecN_2xf32 *)r;
    delLine = ( xb_vecN_2xf32 *)latrs->delLine;
    coef    = latrs->coef;

    Gain    = latrs->gain;
    CF      = *coef;

    for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++ )
    {
        DL = BBE_LVN_2XF32_I(delLine, 0);

        __Pragma("loop_count min=1");
        for ( n=0; n<N; n++ )
        {
            BBE_LVN_2XF32_XP(Xin, px, L*sizeof(float32_t));
            ACCR = BBE_MULN_2XF32(Xin, Gain);
            BBE_MULSN_2XF32(ACCR, DL, CF);
            Rout = ACCR;
            DL = ACCR;
            BBE_SVN_2XF32_XP(Rout, pr, L*sizeof(float32_t));
        }

        px = (const xb_vecN_2xf32 *)((float32_t *)px - L*N + BBE_SIMD_WIDTH/2);
        pr = (      xb_vecN_2xf32 *)((float32_t *)pr - L*N + BBE_SIMD_WIDTH/2);

        BBE_SVN_2XF32_IP(DL, delLine, BBE_SIMD_WIDTH/2*sizeof(float32_t));
    }
} /* latrsf_process1() */

#endif /* if HAVE_VFPU */
