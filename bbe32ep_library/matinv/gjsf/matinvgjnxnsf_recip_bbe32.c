/* ------------------------------------------------------------------------ */
/* Copyright (c) 2017 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
    Real Matrix Gauss-Jordan inversion, floating point real data, block 
    format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "matinvgjnxnsf_common.h"

#if HAVE_VFPU
/*---------------------------------------------
reciprocal of main diagonal
Input/output:
x[L]    - diagonal element (streaming format)

Returns:
none
---------------------------------------------*/
void matinvgjnxnsf_recip(float32_t * x,int L)
{
    int l;
    const xb_vecN_2xf32 * restrict pZrd0;
    const xb_vecN_2xf32 * restrict pZrd1;
          xb_vecN_2xf32 * restrict pZwr;
    pZrd0  =(const xb_vecN_2xf32*)(x);
    pZrd1  =(const xb_vecN_2xf32*)(x);
    pZwr   =(      xb_vecN_2xf32*)(x);
    __Pragma("loop_count min=1")
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++) 
    {
        xb_vecN_2xf32 z4,zkk0,zkk1;
        BBE_LVN_2XF32_IP(zkk0,pZrd0,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(zkk1,pZrd1,2*BBE_SIMD_WIDTH);
        {
            xb_vecN_2xf32 vt1, vt2, vt3;
            vt1=BBE_RECIP0N_2XF32 (zkk0);
            vt2=BBE_CONSTN_2XF32 (1);
            BBE_MULSN_2XF32 (vt2, zkk0, vt1);
            BBE_MULANN_2XF32(vt1, vt1, vt2);
            vt3=BBE_CONSTN_2XF32 (1);
            BBE_MULSN_2XF32 (vt3, zkk1, vt1);
            BBE_MULANN_2XF32 (vt1, vt1, vt3);
            z4=vt1;
        }

        BBE_SVN_2XF32_IP(z4,pZwr,2*BBE_SIMD_WIDTH);
    }
}
#endif
