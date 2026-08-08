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
    Complex Matrix Gauss-Jordan inversion, floating point complex data, block 
    format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "cmatinvgjnxnsf_common.h"

#if HAVE_VFPU

/*---------------------------------------------
reciprocal of main diagonal
Input/output:
x[L]    - diagonal element (streaming format)

Returns:
none
---------------------------------------------*/
void cmatinvgjnxnsf_recip(complex_float * x,int L)
#if 0
{
    int l;
    const xb_vecN_4xcf32 * restrict pZrd;
          xb_vecN_4xcf32 * restrict pZwr;
    pZrd   =(const xb_vecN_4xcf32*)(x);
    pZwr   =(      xb_vecN_4xcf32*)(x);
    __Pragma("loop_count min=1")
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_4xcf32 z4;
        BBE_LVN_4XCF32_IP(z4,pZrd,2*BBE_SIMD_WIDTH);
        z4=BBE_RECIPN_4XCF32(z4);
        BBE_SVN_4XCF32_IP(z4,pZwr,2*BBE_SIMD_WIDTH);
    }
}
#else
{
    int l;
    const xb_vecN_2xf32 * restrict pZrd1;
    const xb_vecN_2xf32 * restrict pZrd0;
          xb_vecN_2xf32 * restrict pZwr;
    NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    NASSERT((L%(BBE_SIMD_WIDTH/4))==0);
    pZrd0   =(const xb_vecN_2xf32*)(x);
    pZrd1   =(const xb_vecN_2xf32*)(x);
    pZwr   =(      xb_vecN_2xf32*)(x);
    __Pragma("loop_count min=1")
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++) 
    {
        xb_vecN_2xf32 c,t8,t1,t2,t5;
        BBE_LVN_2XF32_IP(c,pZrd0,2*BBE_SIMD_WIDTH);
        t8=BBE_MULN_2XF32(c, c);
        BBE_MULMASN_2XF32 (t8, c, c, 0, 3);
        t1=BBE_RECIP0N_2XF32 (t8);
        t2=BBE_CONSTN_2XF32 (1);
        BBE_MULSN_2XF32 (t2, t8, t1);
        BBE_MULANN_2XF32(t1, t1, t2);
        t5=BBE_CONSTN_2XF32 (1);
        BBE_MULSN_2XF32 (t5, t8, t1);
        BBE_MULANN_2XF32(t1, t1, t5);
        BBE_LVN_2XF32_IP(c,pZrd1,2*BBE_SIMD_WIDTH);
        c=BBE_MULMN_2XF32(t1, c, 2, 12);
        BBE_SVN_2XF32_IP(c,pZwr,2*BBE_SIMD_WIDTH);
    }
}
#endif

#endif
