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
    Direct inversion of 3x3 floating point matrices 
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
#include "matinv3x3_common.h"
#if HAVE_VFPU
#define __OPTIMIZED__ 1
/*-------------------------------------------------
    inplace stream-to-block/block-to-stream conversion:
    Input:
    X[9*L] - L matrices 3x3 in stream order 
             (last element is not processed!)
    Output:
    X       - L matrices written in strange format
              elements of 1-st matrix is written to the 
              place of x00,x01 of original input etc.
              so, 8 matrices replaces 8 original elements
              last element is left unchanged
-------------------------------------------------*/
void matinv3x3f_rsb3x3_inplace(float32_t *X,int L)
#if !__OPTIMIZED__
{
    int p,n,l;
    float32_t x[8][8];
    for (l=0; l<L; l+=8,X+=8)
    {
        for (p=0; p<8 ; p++)
        for (n=0; n<8; n++) x[n][p]=X[p+L*n];

        for (n=0; n<8; n++) 
        for (p=0; p<8 ; p++)
        {
            X[p+L*n]= x[p][n];
        }
    }
}
#else
{
    xb_vecN_2xf32* restrict w = (xb_vecN_2xf32*)X;
    const xb_vecN_2xf32* restrict z = (const xb_vecN_2xf32*)X;
    int l;
    xb_vecN_2xf32 X0, X1, X2, X3, X4, X5, X6, X7;
    xb_vecN_2xf32 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
    NASSERT_ALIGN(X, 2 * BBE_SIMD_WIDTH);
    NASSERT((L&(BBE_SIMD_WIDTH / 2 - 1)) == 0);
    for (l = 0; l<(L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        X1=BBE_LVN_2XF32_X ( z, 1*L*sizeof(float32_t));
        X2=BBE_LVN_2XF32_X ( z, 2*L*sizeof(float32_t));
        X3=BBE_LVN_2XF32_X ( z, 3*L*sizeof(float32_t));
        X4=BBE_LVN_2XF32_X ( z, 4*L*sizeof(float32_t));
        X5=BBE_LVN_2XF32_X ( z, 5*L*sizeof(float32_t));
        X6=BBE_LVN_2XF32_X ( z, 6*L*sizeof(float32_t));
        X7=BBE_LVN_2XF32_X ( z, 7*L*sizeof(float32_t));
        BBE_LVN_2XF32_IP(X0, z, 2*BBE_SIMD_WIDTH);

        BBE_DSELN_2XF32I(Y6, Y2, X6, X2,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(Y7, Y3, X7, X3,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(Y4, Y0, X4, X0,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(Y5, Y1, X5, X1,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(X6, X4, Y6, Y4,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(X7, X5, Y7, Y5,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(Y5, Y4, X5, X4,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(Y7, Y6, X7, X6,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(X2, X0, Y2, Y0,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(X3, X1, Y3, Y1,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(Y1, Y0, X1, X0,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(Y3, Y2, X3, X2,BBE_DSELI_INTERLEAVE_2);

        BBE_SVN_2XF32_X (Y1, w, 1*L*sizeof(float32_t));
        BBE_SVN_2XF32_X (Y2, w, 2*L*sizeof(float32_t));
        BBE_SVN_2XF32_X (Y3, w, 3*L*sizeof(float32_t));
        BBE_SVN_2XF32_X (Y4, w, 4*L*sizeof(float32_t));
        BBE_SVN_2XF32_X (Y5, w, 5*L*sizeof(float32_t));
        BBE_SVN_2XF32_X (Y6, w, 6*L*sizeof(float32_t));
        BBE_SVN_2XF32_X (Y7, w, 7*L*sizeof(float32_t));
        BBE_SVN_2XF32_IP(Y0, w, 2*BBE_SIMD_WIDTH);
    }
}
#endif

#endif
