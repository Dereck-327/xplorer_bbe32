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

/*
  NatureDSP_Baseband library. IIR part
    Lattice real block IIR, floating point, streaming version
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* M=8 */

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
void latrsf_process8( latrsf_t  *          latrs,
                      float32_t * restrict r,
                const float32_t *          x,
                      int                  N )
{
    const xb_vecN_2xf32 * restrict px;
          xb_vecN_2xf32 * restrict pr;
          xb_vecN_2xf32 * restrict delLine;
    const float32_t     * restrict coef;
    xb_vecN_2xf32 Xin,
                  DL0,DL1,DL2,DL3,DL4,DL5,DL6,DL7,
                  CF0,CF1,CF2,CF3,CF4,CF5,CF6,CF7;

    xb_vecN_2xf32 ACCR;
    xb_vecN_2xf32 Gain;

    int n, l, L;
    
    NASSERT(latrs->M==8);
    
    L = latrs->L;
    px = (const xb_vecN_2xf32 *)x;
    pr = (      xb_vecN_2xf32 *)r;
    delLine = ( xb_vecN_2xf32 *)latrs->delLine;
    coef    = latrs->coef;

    Gain    = latrs->gain;
    CF0     = coef[0];
    CF1     = coef[1];
    CF2     = coef[2];
    CF3     = coef[3];
    CF4     = coef[4];
    CF5     = coef[5];
    CF6     = coef[6];
    CF7     = coef[7];

    for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++ )
    {
        DL7 = BBE_LVN_2XF32_I(delLine, 0*BBE_SIMD_WIDTH/2*sizeof(float32_t));
        DL6 = BBE_LVN_2XF32_I(delLine, 1*BBE_SIMD_WIDTH/2*sizeof(float32_t));
        DL5 = BBE_LVN_2XF32_I(delLine, 2*BBE_SIMD_WIDTH/2*sizeof(float32_t));
        DL4 = BBE_LVN_2XF32_I(delLine, 3*BBE_SIMD_WIDTH/2*sizeof(float32_t));
        DL3 = BBE_LVN_2XF32_I(delLine, 4*BBE_SIMD_WIDTH/2*sizeof(float32_t));
        DL2 = BBE_LVN_2XF32_I(delLine, 5*BBE_SIMD_WIDTH/2*sizeof(float32_t));
        DL1 = BBE_LVN_2XF32_I(delLine, 6*BBE_SIMD_WIDTH/2*sizeof(float32_t));
        DL0 = BBE_LVN_2XF32_I(delLine, 7*BBE_SIMD_WIDTH/2*sizeof(float32_t));

        __Pragma("loop_count min=1");
        for ( n=0; n<N; n++ )
        {
            BBE_LVN_2XF32_XP(Xin, px, L*sizeof(float32_t));
            ACCR = BBE_MULN_2XF32(Xin, Gain);
            BBE_MULSN_2XF32(ACCR, DL7, CF7);

            DL7 = DL6;
            BBE_MULSN_2XF32(ACCR, DL6, CF6);
            BBE_MULAN_2XF32(DL7, ACCR, CF6);

            DL6 = DL5;
            BBE_MULSN_2XF32(ACCR, DL5, CF5);
            BBE_MULAN_2XF32(DL6, ACCR, CF5);

            DL5 = DL4;
            BBE_MULSN_2XF32(ACCR, DL4, CF4);
            BBE_MULAN_2XF32(DL5, ACCR, CF4);

            DL4 = DL3;
            BBE_MULSN_2XF32(ACCR, DL3, CF3);
            BBE_MULAN_2XF32(DL4, ACCR, CF3);

            DL3 = DL2;
            BBE_MULSN_2XF32(ACCR, DL2, CF2);
            BBE_MULAN_2XF32(DL3, ACCR, CF2);

            DL2 = DL1;
            BBE_MULSN_2XF32(ACCR, DL1, CF1);
            BBE_MULAN_2XF32(DL2, ACCR, CF1);

            DL1 = DL0;
            BBE_MULSN_2XF32(ACCR, DL0, CF0);
            BBE_MULAN_2XF32(DL1, ACCR, CF0);
            
            DL0 = ACCR;
            BBE_SVN_2XF32_XP(ACCR, pr, L*sizeof(float32_t));
        }

        px = (const xb_vecN_2xf32 *)((float32_t *)px - L*N + BBE_SIMD_WIDTH/2);
        pr = (      xb_vecN_2xf32 *)((float32_t *)pr - L*N + BBE_SIMD_WIDTH/2);
        
        BBE_SVN_2XF32_IP(DL7, delLine, BBE_SIMD_WIDTH/2*sizeof(float32_t));
        BBE_SVN_2XF32_IP(DL6, delLine, BBE_SIMD_WIDTH/2*sizeof(float32_t));
        BBE_SVN_2XF32_IP(DL5, delLine, BBE_SIMD_WIDTH/2*sizeof(float32_t));
        BBE_SVN_2XF32_IP(DL4, delLine, BBE_SIMD_WIDTH/2*sizeof(float32_t));
        BBE_SVN_2XF32_IP(DL3, delLine, BBE_SIMD_WIDTH/2*sizeof(float32_t));
        BBE_SVN_2XF32_IP(DL2, delLine, BBE_SIMD_WIDTH/2*sizeof(float32_t));
        BBE_SVN_2XF32_IP(DL1, delLine, BBE_SIMD_WIDTH/2*sizeof(float32_t));
        BBE_SVN_2XF32_IP(DL0, delLine, BBE_SIMD_WIDTH/2*sizeof(float32_t));
    }
} /* latrsf_process8() */

#endif /* if HAVE_VFPU */
