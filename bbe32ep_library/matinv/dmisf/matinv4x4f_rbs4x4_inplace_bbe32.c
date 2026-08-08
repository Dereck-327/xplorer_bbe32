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
#include "bs_common.h"
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
void matinv4x4f_rbs4x4_inplace(float32_t *X,int L,eLayout layout)
#if !__OPTIMIZED__
{
    int p,n,l;
    float32_t x[16][BBE_SIMD_WIDTH/2];
    if (layout==e4x4_stream)
    {
        for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
        {
            for (n=0; n<16; n+=2) 
            for (p=0; p<BBE_SIMD_WIDTH/2; p++)
            {
                x[p  ][n>>1]= X[l+p+L*(n+0)];
                x[p+8][n>>1]= X[l+p+L*(n+1)];
            }
            for (p=0; p<BBE_SIMD_WIDTH/2; p++)
            for (n=0; n<16; n++) X[l+p+L*n]=x[n][p];
        }
    }
    else
    {
        NASSERT(layout==e4x4_block);
        for (l=0; l<L; l+=BBE_SIMD_WIDTH/2,X+=16*BBE_SIMD_WIDTH/2)
        {
            for (p=0; p<BBE_SIMD_WIDTH/2; p++)
            {
                for (n=0; n<16; n++) x[n][p]=X[p*16+n];
            }

            for (p=0; p<BBE_SIMD_WIDTH/2; p++)
            {
                for (n=0; n<16; n++) X[p+n*BBE_SIMD_WIDTH/2]=x[n][p];
            }
        }
    }
}
#else
{
    const xb_vecNx16* restrict pX = (const xb_vecNx16*)X;
          xb_vecNx16* restrict pY = (xb_vecNx16*)X;
    xb_vecNx16 x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15;
    xb_vecNx16 y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15;
    int l,estride,mstride;
    NASSERT_ALIGN(X, 2 * BBE_SIMD_WIDTH);
    NASSERT((L&(BBE_SIMD_WIDTH / 2 - 1)) == 0);
    if(layout==e4x4_stream)
    {
        estride=(L << 2); mstride=2*BBE_SIMD_WIDTH; 
    }
    else
    {
        NASSERT(layout==e4x4_block);
        estride=2*BBE_SIMD_WIDTH; mstride=16*2*BBE_SIMD_WIDTH; 
    }
    __Pragma("loop_count min=1")
    for (l = 0; l<(L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        BBE_LVNX16_XP(x0 , pX, 2*estride);
        BBE_LVNX16_XP(x1 , pX, 2*estride); 
        BBE_LVNX16_XP(x2 , pX, 2*estride); 
        BBE_LVNX16_XP(x3 , pX, 2*estride); 
        BBE_LVNX16_XP(x4 , pX, 2*estride); 
        BBE_LVNX16_XP(x5 , pX, 2*estride); 
        BBE_LVNX16_XP(x6 , pX, 2*estride); 
        BBE_LVNX16_XP(x7 , pX, -13*estride); 

        BBE_LVNX16_XP(x8 , pX, 2*estride); 
        BBE_LVNX16_XP(x9 , pX, 2*estride); 
        BBE_LVNX16_XP(x10, pX, 2*estride); 
        BBE_LVNX16_XP(x11, pX, 2*estride); 
        BBE_LVNX16_XP(x12, pX, 2*estride); 
        BBE_LVNX16_XP(x13, pX, 2*estride); 
        BBE_LVNX16_XP(x14, pX, 2*estride); 
        BBE_LVNX16_XP(x15, pX, -15*estride+mstride);

        BBE_DSELNX16I(y1 ,y0 ,x1 ,x0 ,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(y3 ,y2 ,x3 ,x2 ,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(y5 ,y4 ,x5 ,x4 ,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(y7 ,y6 ,x7 ,x6 ,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(y9 ,y8 ,x9 ,x8 ,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(y11,y10,x11,x10,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(y13,y12,x13,x12,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(y15,y14,x15,x14,BBE_DSELI_INTERLEAVE_2);

        BBE_DSELNX16I(x1 ,x0 ,y2 ,y0 ,BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(x3 ,x2 ,y3 ,y1 ,BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(x5 ,x4 ,y6 ,y4 ,BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(x7 ,x6 ,y7 ,y5 ,BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(x9 ,x8 ,y10,y8 ,BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(x11,x10,y11,y9 ,BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(x13,x12,y14,y12,BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(x15,x14,y15,y13,BBE_DSELI_INTERLEAVE_4);

        DEINTLV3(y1 , y0 , x4 , x0 );
        DEINTLV3(y3 , y2 , x5 , x1 );
        DEINTLV3(y5 , y4 , x6 , x2 );
        DEINTLV3(y7 , y6 , x7 , x3 );
        DEINTLV3(y9 , y8 , x12, x8 );
        DEINTLV3(y11, y10, x13, x9 );
        DEINTLV3(y13, y12, x14, x10);
        DEINTLV3(y15, y14, x15, x11);

        BBE_SVNX16_XP(y0 , pY, estride);
        BBE_SVNX16_XP(y1 , pY, estride);
        BBE_SVNX16_XP(y2 , pY, estride);
        BBE_SVNX16_XP(y3 , pY, estride);
        BBE_SVNX16_XP(y4 , pY, estride);
        BBE_SVNX16_XP(y5 , pY, estride);
        BBE_SVNX16_XP(y6 , pY, estride);
        BBE_SVNX16_XP(y7 , pY, estride);
        BBE_SVNX16_XP(y8 , pY, estride);
        BBE_SVNX16_XP(y9 , pY, estride);
        BBE_SVNX16_XP(y10, pY, estride);
        BBE_SVNX16_XP(y11, pY, estride);
        BBE_SVNX16_XP(y12, pY, estride);
        BBE_SVNX16_XP(y13, pY, estride);
        BBE_SVNX16_XP(y14, pY, estride);
        BBE_SVNX16_XP(y15, pY, -15*estride+mstride);
    }
}
#endif

#endif
