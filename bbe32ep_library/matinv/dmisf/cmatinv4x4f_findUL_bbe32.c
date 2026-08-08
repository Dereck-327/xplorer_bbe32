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
#include "matinv4x4Tbl.h"
#include "matinv4x4_common.h"

#define __OPTIMIZED__ 1

#if HAVE_VFPU

#if !__OPTIMIZED__
#include <complex.h>
#include <math.h>
#endif

/*------------------------------------------
find position of upper left corner
Look for an element of maximum absolute value.
Input:
X[16*L]  -L matrices - interleaved in intermediate format
Output:
pos[L]   - index of position of element with 
maximum absolute value
------------------------------------------*/
void cmatinv4x4f_findUL(int16_t *pos,const complex_float *X, int L, eLayout layout)
#if !__OPTIMIZED__
{
    int l;
    int        n,maxIx;
    float32_t  w,maxW;
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    NASSERT(layout==e4x4_stream ||layout==e4x4_block);
    for ( l=0; l<L; l++ )
    {
        complex_float x[16];
        if (layout==e4x4_stream)
        {
            for (n=0; n<4 ; n++)
            {
               x[n+ 0]= X[(l&~3)+n+L*((l&3)*4+0)];
               x[n+ 4]= X[(l&~3)+n+L*((l&3)*4+1)];
               x[n+ 8]= X[(l&~3)+n+L*((l&3)*4+2)];
               x[n+12]= X[(l&~3)+n+L*((l&3)*4+3)];
            }
        }
        else
        {
            for (n=0; n<16; n++) x[n]=X[16*l+n];
        }
        for ( maxW=0.f, maxIx=0, n=0; n<16; n++ )
        {
            complex_float t=x[n];
            w=fabsf(crealf(t)*crealf(t)+cimagf(t)*cimagf(t));
            if ( maxW < w) { maxW = w, maxIx = n; }
        }
        pos[l]=maxIx;
    }
}
#else
{
    xb_vecN_2xc16 ix0,ix8;
    short * restrict pPos=(short *)pos;
    const xb_vecN_2xf32* restrict pX=(const xb_vecN_2xf32*)X;
    int l,estride,mstride;
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);

    ix0 = BBE_SRLIN_2XC16(BBE_MOVN_2XC16_FROMNX16(BBE_SEQNX16()),1);
    ix8 = BBE_SRLIN_2XC16(BBE_MOVN_2XC16_FROMNX16(BBE_ADDNX16(BBE_SEQNX16(),BBE_SIMD_WIDTH)),1);
    if(layout==e4x4_block)
    {
        estride=2*BBE_SIMD_WIDTH;  mstride=16*2*BBE_SIMD_WIDTH; 
    }
    else
    {
        NASSERT(layout==e4x4_stream);
        estride=L<<3;  mstride=2*BBE_SIMD_WIDTH; 
    }
    for ( l=0; l<L; l++ )
    {
        xb_vecN_2xf32 a,x0,x1,x2,x3;
        xb_vecN_2xc16 maxix;
        vboolN_2 bmax;
        int addx;
        int dummy;
        xtfloat MAXW;
        vselN vmax;

        addx = 4*estride;
        XT_MOVEQZ(addx, -12*estride+mstride, (l&3)^3);
        x1=BBE_LVN_2XF32_X (pX,1*estride);
        x2=BBE_LVN_2XF32_X (pX,2*estride);
        x3=BBE_LVN_2XF32_X (pX,3*estride);
        BBE_LVN_2XF32_XP(x0,pX,addx);

        a=BBE_MULMN_2XF32(x0,x0,0,0); BBE_MULMASN_2XF32(a,x0,x0,0,15); x0=a;
        a=BBE_MULMN_2XF32(x1,x1,0,0); BBE_MULMASN_2XF32(a,x1,x1,0,15); x1=a;
        a=BBE_MULMN_2XF32(x2,x2,0,0); BBE_MULMASN_2XF32(a,x2,x2,0,15); x2=a;
        a=BBE_MULMN_2XF32(x3,x3,0,0); BBE_MULMASN_2XF32(a,x3,x3,0,15); x3=a;
        x0=BBE_SELN_2XF32I(x1,x0,BBE_SELI_16B_EXTRACT_2_OF_4_OFF_0);
        x1=BBE_SELN_2XF32I(x3,x2,BBE_SELI_16B_EXTRACT_2_OF_4_OFF_0);
        maxix=BBE_MOVN_2XC16T(ix8,ix0,BBE_OGTN_2XF32(x1,x0)); 
        x0=BBE_MAXNUMN_2XF32(x0,x1);
        BBE_RBMAXNUMN_2XF32(bmax,MAXW,x0);
        (void)MAXW;
        BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));(void)dummy;
        BBE_SSNX16_IP(BBE_SHFLNX16(BBE_MOVNX16_FROMN_2XC16(maxix),vmax),pPos,sizeof(int16_t));
    }
}
#endif

#endif
