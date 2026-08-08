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
    NatureDSP_Baseband library. Direct Matrix Inversion
    Direct inversion of 4x4 floating point matrices 
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
#include "matinv4x4_common.h"

#define __OPTIMIZED__ 1

#if HAVE_VFPU

/*-------------------------------------------------
    inplace stream-to-block/block-to-stream conversion:
    Input:
    X[16*L] - L matrices 4x4 in stream order
    Output:
    X       - L matrices written in strange format
              elements of 1-st matrix is written to the 
              place of x00,x01 of original input etc.
              so, 8 matrices replaces 
-------------------------------------------------*/
void cmatinv4x4f_csb4x4_inplace(complex_float *X,int L,eLayout layout)
#if !__OPTIMIZED__
{
    int p,n,l;
    complex_float x[16][4];
    if (layout==e4x4_stream)
    {
        for (l=0; l<L; l+=4,X+=4)
        {
            for (p=0; p<4 ; p++)
            for (n=0; n<16; n++) x[n][p]=X[p+L*n];

            for (n=0; n<16; n+=4) 
            for (p=0; p<4 ; p++)
            {
                X[p+L*(n+0)]= x[p   ][n>>2];
                X[p+L*(n+1)]= x[p+ 4][n>>2];
                X[p+L*(n+2)]= x[p+ 8][n>>2];
                X[p+L*(n+3)]= x[p+12][n>>2];
            }
        }
    }
    else
    {
        NASSERT(layout==e4x4_block);
        for (l=0; l<L; l+=BBE_SIMD_WIDTH/4,X+=16*BBE_SIMD_WIDTH/4)
        {
            for (p=0; p<BBE_SIMD_WIDTH/4; p++)
            {
                for (n=0; n<16; n++) x[n][p]=X[p+n*BBE_SIMD_WIDTH/4];
            }
            for (p=0; p<BBE_SIMD_WIDTH/4; p++)
            {
                for (n=0; n<16; n++) X[p*16+n]=x[n][p];
            }
        }
    }
}
#else
{
    const xb_vecN_4xcf32 * restrict px;
        xb_vecN_4xcf32 * restrict py;
    xb_vecN_4xcf32 x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15;
    int l,estride,mstride;

    NASSERT_ALIGN(X, 2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(X, 2*BBE_SIMD_WIDTH);
    NASSERT((L & (BBE_SIMD_WIDTH/4-1)) == 0);
    if (L<=0) return;

    px = (const xb_vecN_4xcf32 *)(X);
    py = (      xb_vecN_4xcf32 *)(X);
    if(layout==e4x4_stream)
    {
        estride=(L << 3); mstride=2*BBE_SIMD_WIDTH; 
    }
    else
    {
        NASSERT(layout==e4x4_block);
        estride=2*BBE_SIMD_WIDTH; mstride=16*2*BBE_SIMD_WIDTH; 
    }

  __Pragma("loop_count min=1");
    for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++ )
    {
        BBE_LVN_4XCF32_XP(x0 , px, estride);
        BBE_LVN_4XCF32_XP(x1 , px, estride);
        BBE_LVN_4XCF32_XP(x2 , px, estride);
        BBE_LVN_4XCF32_XP(x3 , px, estride);
        BBE_LVN_4XCF32_XP(x4 , px, estride);
        BBE_LVN_4XCF32_XP(x5 , px, estride);
        BBE_LVN_4XCF32_XP(x6 , px, estride);
        BBE_LVN_4XCF32_XP(x7 , px, estride);
        BBE_LVN_4XCF32_XP(x8 , px, estride);
        BBE_LVN_4XCF32_XP(x9 , px, estride);
        BBE_LVN_4XCF32_XP(x10, px, estride);
        BBE_LVN_4XCF32_XP(x11, px, estride);
        BBE_LVN_4XCF32_XP(x12, px, estride);
        BBE_LVN_4XCF32_XP(x13, px, estride);
        BBE_LVN_4XCF32_XP(x14, px, estride);
        BBE_LVN_4XCF32_XP(x15, px, -15*estride+mstride);

        BBE_DSELN_4XCF32I(x2, x0, x2, x0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(x3, x1, x3, x1, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(x1, x0, x1, x0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(x3, x2, x3, x2, BBE_DSELI_INTERLEAVE_4);

        BBE_DSELN_4XCF32I(x6, x4, x6, x4, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(x7, x5, x7, x5, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(x5, x4, x5, x4, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(x7, x6, x7, x6, BBE_DSELI_INTERLEAVE_4);

        BBE_DSELN_4XCF32I(x10, x8 , x10, x8 , BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(x11, x9 , x11, x9 , BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(x9 , x8 , x9 , x8 , BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(x11, x10, x11, x10, BBE_DSELI_INTERLEAVE_4);

        BBE_DSELN_4XCF32I(x14, x12, x14, x12, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(x15, x13, x15, x13, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(x13, x12, x13, x12, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(x15, x14, x15, x14, BBE_DSELI_INTERLEAVE_4);

        BBE_SVN_4XCF32_XP(x0 , py,  4*estride);
        BBE_SVN_4XCF32_XP(x1 , py,  4*estride);
        BBE_SVN_4XCF32_XP(x2 , py,  4*estride);
        BBE_SVN_4XCF32_XP(x3 , py,-11*estride);
        BBE_SVN_4XCF32_XP(x4 , py,  4*estride);
        BBE_SVN_4XCF32_XP(x5 , py,  4*estride);
        BBE_SVN_4XCF32_XP(x6 , py,  4*estride);
        BBE_SVN_4XCF32_XP(x7 , py,-11*estride);
        BBE_SVN_4XCF32_XP(x8 , py,  4*estride);
        BBE_SVN_4XCF32_XP(x9 , py,  4*estride);
        BBE_SVN_4XCF32_XP(x10, py,  4*estride);
        BBE_SVN_4XCF32_XP(x11, py,-11*estride);
        BBE_SVN_4XCF32_XP(x12, py,  4*estride);
        BBE_SVN_4XCF32_XP(x13, py,  4*estride);
        BBE_SVN_4XCF32_XP(x14, py,  4*estride);
        BBE_SVN_4XCF32_XP(x15, py,-15*estride+mstride);
  }
}
#endif

#endif
