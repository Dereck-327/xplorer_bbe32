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
#include "matinv3x3Tbl.h"
#include "matinv3x3_common.h"
#if HAVE_VFPU
#define __OPTIMIZED__ 1

#if !__OPTIMIZED__
#include <math.h>
#endif
/*-----------------------------------------------
    Search for permutation in block/stream ordered 
    arrays
    Input:
    X[]  - 3x3 block/stream ordered data
    Output:
    permIx[L] - permutation index
-----------------------------------------------*/
void matinv3x3f_search(int16_t *permIx,const float32_t *X, int L, eLayout layout)
#if !__OPTIMIZED__
{
    int l;
    const int16_t *perm=matinv3x3sf_searchTbl;
    for ( l=0; l<L; l++ )
    {
        float32_t x[9],a[8],b[8],c[8],d[8];
        int p, n;
        float32_t  w,maxW;
        int      maxIx;

        maxIx=permIx[l];
       if (layout==e3x3_stream)
       {
           for (n=0; n<8 ; n++)
           {
               x[n  ]= X[(l&~7)+n+L*(l&7)];
           }
           x[8]=X[l+8*L];
       }
       else
       {
           NASSERT(layout==e3x3_block);
           for (n=0; n<9 ; n++) x[n]=X[16*l+n];
       }

        for (n=0; n<8; n++)
        {
            a[n]=x[perm[16*0+n*2]>>1];
            b[n]=x[perm[16*1+n*2]>>1];
            c[n]=x[perm[16*2+n*2]>>1];
            d[n]=x[perm[16*3+n*2]>>1];
        }

        for ( maxW=0.f, p=0, n=0; n<8; n++ )
        {
            w=fabsf(a[n]*d[n]-b[n]*c[n]);
            if ( maxW < w ) { maxW = w, p = n; }
        }
        permIx[l]=p;
    }
}
#else
{
    xb_vecNx16 permA,permB,permC,permD;
    vselN_2 selA,selB,selC,selD;
    short * restrict pPermWr;
    const xb_vecN_2xf32* restrict pXrd;
    const xtfloat      * restrict pXrd8;
    const xb_vecNx16   * restrict perm;
    int l,addx,estride,mstride,estride8;
    pXrd=(const xb_vecN_2xf32*)X;
    if (layout==e3x3_block)
    {
        estride = 4*BBE_SIMD_WIDTH;
        mstride =16*2*BBE_SIMD_WIDTH;
        pXrd8=( xtfloat*)(X+8);
        estride8=4*BBE_SIMD_WIDTH;
    }
    else
    {
        NASSERT(layout==e3x3_stream);
        estride = L*sizeof(float32_t);
        mstride = 2*BBE_SIMD_WIDTH;
        pXrd8=(     xtfloat*)(X+8*L);
        estride8=sizeof(float32_t);
    }
    pPermWr=(short*)permIx;
    perm=(const xb_vecNx16*)matinv3x3sf_searchTbl;

    permA=BBE_LVNX16_I(perm,0*2*BBE_SIMD_WIDTH);
    permB=BBE_LVNX16_I(perm,1*2*BBE_SIMD_WIDTH);
    permC=BBE_LVNX16_I(perm,2*2*BBE_SIMD_WIDTH);
    permD=BBE_LVNX16_I(perm,3*2*BBE_SIMD_WIDTH);
    selA=BBE_MOVVSELN_2NX16(permA,0);
    selB=BBE_MOVVSELN_2NX16(permB,0);
    selC=BBE_MOVVSELN_2NX16(permC,0);
    selD=BBE_MOVVSELN_2NX16(permD,0);

    for ( l=0; l<L; l++ )
    {
        xb_vecN_2xf32 X0,X1,A,B,C,D,W;
        xtfloat MAXW;
        vboolN_2 bmax;
        vselN vmax;
        int dummy;
        xb_vecNx16 maxix;

        addx = estride;
        XT_MOVEQZ(addx, -7*estride+mstride, (l&7)^7);

        BBE_LVN_2XF32_XP(X0,pXrd,addx);
        BBE_LSN_2XF32_XP(X1,pXrd8,estride8);
        maxix=BBE_LVNX16_I(perm,4*2*BBE_SIMD_WIDTH);
        A=BBE_SELN_2XF32(X1,X0,selA);
        B=BBE_SELN_2XF32(X1,X0,selB);
        C=BBE_SELN_2XF32(X1,X0,selC);
        D=BBE_SELN_2XF32(X1,X0,selD);
        W=BBE_MULN_2XF32(A,D);
        BBE_MULSN_2XF32(W,B,C);
        W=BBE_ABSN_2XF32(W);
        BBE_RBMAXNUMN_2XF32(bmax,MAXW,W);
        (void)MAXW;
        BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));(void)dummy;
        BBE_SSNX16_IP(BBE_SHFLNX16(maxix,vmax),pPermWr,sizeof(int16_t));
    }
}
#endif


#endif
