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
  NatureDSP_Baseband library. Matrix Operations
    Streaming to Packed Conversion for Real and Complex Matrices
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Math.h"
#include "NatureDSP_Baseband_matop.h"
#include "common.h"

#define BBE_MOVVW(hi,lo,w) { hi=BBE_MOVVWH(w);lo=BBE_MOVVWL(w);}

/*   MN=9...15, Sx=MN, Sy=16 */
void rsb9_15x1(int16_t* restrict y, const int16_t* restrict x, int MN, int L)
{
    const xb_vecNx16* restrict X=(const xb_vecNx16*)x;
          xb_vecNx16* restrict Y=(xb_vecNx16* )y;
    xb_vecNx16 X0,X1,X2,X3,X4,X5,X6,X7,X8,X9,XA,XB,XC,XD,XE,XF;
    xb_vecNx16 Y0,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,YA,YB,YC,YD,YE,YF;
    xb_vecNx40 W0,W1; 
    int k,_2L,off9,off10,off11,off12,off13,off14;
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT((L&(BBE_SIMD_WIDTH - 1)) == 0);
    NASSERT(MN>=9 && MN<=15);
    _2L=(L<<1);
    off9 =MN>9  ?       _2L:0;
    off10=MN>10 ? off9 +_2L:0;
    off11=MN>11 ? off10+_2L:0;
    off12=MN>12 ? off11+_2L:0;
    off13=MN>13 ? off12+_2L:0;
    off14=MN>14 ? off13+_2L:0;
    if (L <= 0) return;
    for (k = 0; k<(L >> LOG2_BBE_SIMD_WIDTH); k++)
    {
        BBE_LVNX16_XP(Y0,X, _2L);
        BBE_LVNX16_XP(Y1,X, _2L);
        BBE_DSELNX16I(X8,X0,Y1,Y0,BBE_DSELI_INTERLEAVE_1);
        W0=BBE_MOVWV(X8, X0);
        BBE_LVNX16_XP(Y2,X, _2L);
        BBE_LVNX16_XP(Y3,X, _2L);
        BBE_DSELNX16I(X9,X1,Y3,Y2,BBE_DSELI_INTERLEAVE_1);
        W1=BBE_MOVWV(X9, X1);
        BBE_LVNX16_XP(Y4,X, _2L);
        BBE_LVNX16_XP(Y5,X, _2L);
        BBE_DSELNX16I(XA,X2,Y5,Y4,BBE_DSELI_INTERLEAVE_1);
        BBE_LVNX16_XP(Y6,X, _2L);
        BBE_LVNX16_XP(Y7,X, _2L);
        BBE_DSELNX16I(XB,X3,Y7,Y6,BBE_DSELI_INTERLEAVE_1);
        Y8=BBE_LVNX16_I (X, 0);
        Y9=BBE_LVNX16_X (X, off9);
        BBE_DSELNX16I(XC,X4,Y9,Y8,BBE_DSELI_INTERLEAVE_1);
        YA=BBE_LVNX16_X (X, off10);
        YB=BBE_LVNX16_X (X, off11);
        BBE_DSELNX16I(XD,X5,YB,YA,BBE_DSELI_INTERLEAVE_1);
        YC=BBE_LVNX16_X (X, off12);
        YD=BBE_LVNX16_X (X, off13);
        BBE_DSELNX16I(XE,X6,YD,YC,BBE_DSELI_INTERLEAVE_1);
        YE=BBE_LVNX16_X (X, off14);
        BBE_DSELNX16I(XF,X7,YF,YE,BBE_DSELI_INTERLEAVE_1);
        X=X+(2*BBE_SIMD_WIDTH-8*_2L)/sizeof(*X);

        BBE_MOVVW(X8, X0, W0);  
        BBE_DSELNX16I(Y4,Y0,X4,X0,BBE_DSELI_INTERLEAVE_2);
        BBE_MOVVW(X9, X1, W1);  
        BBE_DSELNX16I(Y5,Y1,X5,X1,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(Y6,Y2,X6,X2,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(Y7,Y3,X7,X3,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(YC,Y8,XC,X8,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(YD,Y9,XD,X9,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(YE,YA,XE,XA,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(YF,YB,XF,XB,BBE_DSELI_INTERLEAVE_2);

        BBE_DSELNX16I(X2,X0,Y2,Y0,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(X3,X1,Y3,Y1,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(X6,X4,Y6,Y4,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(X7,X5,Y7,Y5,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(XA,X8,YA,Y8,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(XB,X9,YB,Y9,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(XE,XC,YE,YC,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(XF,XD,YF,YD,BBE_DSELI_INTERLEAVE_2);

        BBE_DSELNX16I(Y1,Y0,X1,X0,BBE_DSELI_INTERLEAVE_2);
        BBE_SVNX16_IP(Y0,Y, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(Y1,Y, 2*BBE_SIMD_WIDTH);
        BBE_DSELNX16I(Y3,Y2,X3,X2,BBE_DSELI_INTERLEAVE_2);
        BBE_SVNX16_IP(Y2,Y, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(Y3,Y, 2*BBE_SIMD_WIDTH);
        BBE_DSELNX16I(Y5,Y4,X5,X4,BBE_DSELI_INTERLEAVE_2);
        BBE_SVNX16_IP(Y4,Y, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(Y5,Y, 2*BBE_SIMD_WIDTH);
        BBE_DSELNX16I(Y7,Y6,X7,X6,BBE_DSELI_INTERLEAVE_2);
        BBE_SVNX16_IP(Y6,Y, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(Y7,Y, 2*BBE_SIMD_WIDTH);
        BBE_DSELNX16I(Y9,Y8,X9,X8,BBE_DSELI_INTERLEAVE_2);
        BBE_SVNX16_IP(Y8,Y, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(Y9,Y, 2*BBE_SIMD_WIDTH);
        BBE_DSELNX16I(YB,YA,XB,XA,BBE_DSELI_INTERLEAVE_2);
        BBE_SVNX16_IP(YA,Y, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(YB,Y, 2*BBE_SIMD_WIDTH);
        BBE_DSELNX16I(YD,YC,XD,XC,BBE_DSELI_INTERLEAVE_2);
        BBE_SVNX16_IP(YC,Y, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(YD,Y, 2*BBE_SIMD_WIDTH);
        BBE_DSELNX16I(YF,YE,XF,XE,BBE_DSELI_INTERLEAVE_2);
        BBE_SVNX16_IP(YE,Y, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(YF,Y, 2*BBE_SIMD_WIDTH);
    }
}
