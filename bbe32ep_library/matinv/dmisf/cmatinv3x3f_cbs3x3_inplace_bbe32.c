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
/* NatureDSP_Baseband library API. */
#include "matinv3x3_common.h"
#define __OPTIMIZED__ 1

#if HAVE_VFPU

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
void cmatinv3x3f_cbs3x3_inplace(complex_float *X,int L, eLayout layout)
#if !__OPTIMIZED__
{
    int p,n,l;
    if (layout==e3x3_stream)
    {
        complex_float x[8][4];
        for (l=0; l<L; l+=4,X+=4)
        {
            for (n=0; n<8; n+=2) 
            for (p=0; p<4 ; p++)
            {
                x[p  ][n>>1]= X[p+L*(n+0)];
                x[p+4][n>>1]= X[p+L*(n+1)];
            }
            for (p=0; p<4; p++)
            for (n=0; n<8; n++) X[p+L*n]=x[n][p];
        }
    }
    else
    {
        for (l=0; l<L; l+=BBE_SIMD_WIDTH/4)
        {
            complex_float x[BBE_SIMD_WIDTH/4][12];
            memcpy(x,X,(BBE_SIMD_WIDTH/4)*12*sizeof(complex_float));
            for (p=0; p<BBE_SIMD_WIDTH/4; p++) 
            {
                X[p+ 0*(BBE_SIMD_WIDTH/4)]=x[p][ 0]; X[p+ 1*(BBE_SIMD_WIDTH/4)]=x[p][ 1]; X[p+ 2*(BBE_SIMD_WIDTH/4)]=x[p][ 2]; 
                X[p+ 3*(BBE_SIMD_WIDTH/4)]=x[p][ 3]; X[p+ 4*(BBE_SIMD_WIDTH/4)]=x[p][ 4]; X[p+ 5*(BBE_SIMD_WIDTH/4)]=x[p][ 5]; 
                X[p+ 6*(BBE_SIMD_WIDTH/4)]=x[p][ 6]; X[p+ 7*(BBE_SIMD_WIDTH/4)]=x[p][ 7]; X[p+ 8*(BBE_SIMD_WIDTH/4)]=x[p][ 8]; 

            }
            X+=12*(BBE_SIMD_WIDTH/4);
        }
    }
}
#else
{
  const xb_vecN_4xcf32 * restrict pX;
        xb_vecN_4xcf32 * restrict pY;
  int l;
  xb_vecN_4xcf32 X0, X1, X2, X3;
  xb_vecN_4xcf32 X4, X5, X6, X7;
  xb_vecN_4xcf32 X8, X9, X10, X11;

  NASSERT_ALIGN(X, 2*BBE_SIMD_WIDTH);
  NASSERT((L & (BBE_SIMD_WIDTH/4-1)) == 0);
  if (L<=0) return;


  if (layout==e3x3_stream)
  {
      pX = (const xb_vecN_4xcf32 *)(X+7*L);
      pY = (      xb_vecN_4xcf32 *)(X);
      /* Convert by 8 values */
      for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++ )
      {
            X0 = BBE_LVN_4XCF32_X(pX,-7*L*sizeof(complex_float));
            X4 = BBE_LVN_4XCF32_X(pX,-6*L*sizeof(complex_float));
            X1 = BBE_LVN_4XCF32_X(pX,-5*L*sizeof(complex_float));
            X5 = BBE_LVN_4XCF32_X(pX,-4*L*sizeof(complex_float));
            X2 = BBE_LVN_4XCF32_X(pX,-3*L*sizeof(complex_float));
            X6 = BBE_LVN_4XCF32_X(pX,-2*L*sizeof(complex_float));
            X3 = BBE_LVN_4XCF32_X(pX,-1*L*sizeof(complex_float));
            BBE_LVN_4XCF32_IP(X7, pX,2*BBE_SIMD_WIDTH);

            BBE_DSELN_4XCF32I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
            BBE_DSELN_4XCF32I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
            BBE_DSELN_4XCF32I(X1, X0, X1, X0, BBE_DSELI_INTERLEAVE_4);
            BBE_DSELN_4XCF32I(X3, X2, X3, X2, BBE_DSELI_INTERLEAVE_4);

            BBE_DSELN_4XCF32I(X6, X4, X6, X4, BBE_DSELI_INTERLEAVE_4);
            BBE_DSELN_4XCF32I(X7, X5, X7, X5, BBE_DSELI_INTERLEAVE_4);
            BBE_DSELN_4XCF32I(X5, X4, X5, X4, BBE_DSELI_INTERLEAVE_4);
            BBE_DSELN_4XCF32I(X7, X6, X7, X6, BBE_DSELI_INTERLEAVE_4);

            BBE_SVN_4XCF32_XP(X0,pY,L*sizeof(complex_float));
            BBE_SVN_4XCF32_XP(X1,pY,L*sizeof(complex_float));
            BBE_SVN_4XCF32_XP(X2,pY,L*sizeof(complex_float));
            BBE_SVN_4XCF32_XP(X3,pY,L*sizeof(complex_float));
            BBE_SVN_4XCF32_XP(X4,pY,L*sizeof(complex_float));
            BBE_SVN_4XCF32_XP(X5,pY,L*sizeof(complex_float));
            BBE_SVN_4XCF32_XP(X6,pY,L*sizeof(complex_float));
            BBE_SVN_4XCF32_XP(X7,pY,-7*L*sizeof(complex_float)+2*BBE_SIMD_WIDTH);
        }
  }
  else
  {
      NASSERT(layout==e3x3_block);
    pX = (const xb_vecN_4xcf32 *)(X);
    pY = (      xb_vecN_4xcf32 *)(X);
    __Pragma("loop_count min=1");
    for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++ )
    {
        X1 =BBE_LVN_4XCF32_I ( pX,  3*2*BBE_SIMD_WIDTH); //3
        X2 =BBE_LVN_4XCF32_I ( pX,  6*2*BBE_SIMD_WIDTH); //6
        X3 =BBE_LVN_4XCF32_I ( pX,  9*2*BBE_SIMD_WIDTH); //9
        X4 =BBE_LVN_4XCF32_I ( pX,  1*2*BBE_SIMD_WIDTH); //1
        X5 =BBE_LVN_4XCF32_I ( pX,  4*2*BBE_SIMD_WIDTH); //4
        X6 =BBE_LVN_4XCF32_I ( pX,  7*2*BBE_SIMD_WIDTH); //7
        X7 =BBE_LVN_4XCF32_I ( pX, 10*2*BBE_SIMD_WIDTH); //10
        X8 =BBE_LVN_4XCF32_I ( pX,  2*2*BBE_SIMD_WIDTH); //2
        X9 =BBE_LVN_4XCF32_I ( pX,  5*2*BBE_SIMD_WIDTH); //5
        X10=BBE_LVN_4XCF32_I ( pX,  8*2*BBE_SIMD_WIDTH); //8
        X11=BBE_LVN_4XCF32_I ( pX, 11*2*BBE_SIMD_WIDTH); //11
        BBE_LVN_4XCF32_IP( X0, pX, 12*2*BBE_SIMD_WIDTH); //0

        BBE_DSELN_4XCF32I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(X1, X0, X1, X0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(X3, X2, X3, X2, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(X6, X4, X6, X4, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(X7, X5, X7, X5, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(X5, X4, X5, X4, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_4XCF32I(X7, X6, X7, X6, BBE_DSELI_INTERLEAVE_4);

        X8 = BBE_SELN_4XCF32I(BBE_SELN_4XCF32I(X11, X10, BBE_SELI_PACK_4), BBE_SELN_4XCF32I(X9, X8, BBE_SELI_PACK_4), BBE_SELI_PACK_8);

        BBE_SVN_4XCF32_XP(X0, pY, 2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_XP(X1, pY, 2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_XP(X2, pY, 2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_XP(X3, pY, 2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_XP(X4, pY, 2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_XP(X5, pY, 2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_XP(X6, pY, 2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_IP(X7, pY, 2*BBE_SIMD_WIDTH);
        BBE_SVN_4XCF32_IP(X8, pY, 4*2*BBE_SIMD_WIDTH);
    }
  }
}
#endif

#endif
