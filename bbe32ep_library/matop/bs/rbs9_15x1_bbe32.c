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
	NatureDSP_Baseband library. Communication part
    packed/streaming conversion
	Integrit, 2006-2013
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Math.h"
#include "common.h"

#define BBE_MOVVW(hi,lo,w) { hi=BBE_MOVVWH(w);lo=BBE_MOVVWL(w);}

/*   MN=9..15, Sx=16, Sy=MN */
void rbs9_15x1(int16_t* restrict y, const int16_t* restrict x, int MN, int L)
{
  const xb_vecNx16* restrict X=(const xb_vecNx16*)x;
        xb_vecNx16* restrict Y=(xb_vecNx16* )y;
  xb_vecNx16 X0,X1,X2,X3,X4,X5,X6,X7,X8,X9,XA,XB,XC,XD,XE,XF;
  xb_vecNx16 Y0,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,YA,YB,YC,YD,YE,YF;
  xb_vecNx40 W0; 
  int k,_2L,off9,off10,off11,off12,off13,off14;
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT((L&(BBE_SIMD_WIDTH - 1)) == 0);
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
    BBE_LVNX16_IP(Y0,X, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y1,X, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y2,X, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y3,X, 2*BBE_SIMD_WIDTH);
    BBE_DSELNX16I(X1,X0,Y1,Y0,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X3,X2,Y3,Y2,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y2,Y0,X2,X0,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y3,Y1,X3,X1,BBE_DSELI_DEINTERLEAVE_2);

    BBE_LVNX16_IP(Y4,X, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y5,X, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y6,X, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y7,X, 2*BBE_SIMD_WIDTH);
    BBE_DSELNX16I(X5,X4,Y5,Y4,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X7,X6,Y7,Y6,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y6,Y4,X6,X4,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(Y7,Y5,X7,X5,BBE_DSELI_DEINTERLEAVE_2);

    BBE_LVNX16_IP(Y8,X, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(Y9,X, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(YA,X, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(YB,X, 2*BBE_SIMD_WIDTH);
    BBE_DSELNX16I(X9,X8,Y9,Y8,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(XB,XA,YB,YA,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YA,Y8,XA,X8,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YB,Y9,XB,X9,BBE_DSELI_DEINTERLEAVE_2);

    BBE_LVNX16_IP(YC,X, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(YD,X, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(YE,X, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(YF,X, 2*BBE_SIMD_WIDTH);
    BBE_DSELNX16I(XD,XC,YD,YC,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(XF,XE,YF,YE,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YE,YC,XE,XC,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(YF,YD,XF,XD,BBE_DSELI_DEINTERLEAVE_2);

    BBE_DSELNX16I(X4,X0,Y4,Y0,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X5,X1,Y5,Y1,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X6,X2,Y6,Y2,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(X7,X3,Y7,Y3,BBE_DSELI_DEINTERLEAVE_2);

    BBE_DSELNX16I(XC,X8,YC,Y8,BBE_DSELI_DEINTERLEAVE_2);
    W0=BBE_MOVWV(XC, X4);
    BBE_DSELNX16I(XD,X9,YD,Y9,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(XE,XA,YE,YA,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELNX16I(XF,XB,YF,YB,BBE_DSELI_DEINTERLEAVE_2);

    BBE_DSELNX16I(Y1,Y0,X8,X0,BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y0,Y, _2L);
    BBE_SVNX16_XP(Y1,Y, _2L);
    BBE_DSELNX16I(Y3,Y2,X9,X1,BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y2,Y, _2L);
    BBE_SVNX16_XP(Y3,Y, _2L);
    BBE_DSELNX16I(Y5,Y4,XA,X2,BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y4,Y, _2L);
    BBE_SVNX16_XP(Y5,Y, _2L);
    BBE_DSELNX16I(Y7,Y6,XB,X3,BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_XP(Y6,Y, _2L);
    BBE_SVNX16_XP(Y7,Y, _2L);
    BBE_DSELNX16I(YF,YE,XF,X7,BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X (YE,Y, off14);
    BBE_DSELNX16I(YD,YC,XE,X6,BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X (YC,Y, off12);
    BBE_SVNX16_X (YD,Y, off13);
    BBE_DSELNX16I(YB,YA,XD,X5,BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X (YA,Y, off10);
    BBE_SVNX16_X (YB,Y, off11);
    BBE_MOVVW(XC, X4, W0);  
    BBE_DSELNX16I(Y9,Y8,XC,X4,BBE_DSELI_DEINTERLEAVE_1);
    BBE_SVNX16_X (Y9,Y, off9);
    BBE_SVNX16_XP(Y8,Y, 2*BBE_SIMD_WIDTH-8*_2L);
  }
}
