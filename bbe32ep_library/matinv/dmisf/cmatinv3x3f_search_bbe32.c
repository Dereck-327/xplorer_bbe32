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

#define __OPTIMIZED__ 1
#if HAVE_VFPU

#if !__OPTIMIZED__
#include <math.h>
#include <complex.h>

static complex_float subc(complex_float x,complex_float y)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    z.s.re=crealf(x)-crealf(y);
    z.s.im=cimagf(x)-cimagf(y);
    return z.u;
}

static complex_float mulc(complex_float x,complex_float y)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    z.s.re=crealf(x)*crealf(y) - cimagf(x)*cimagf(y);
    z.s.im=crealf(x)*cimagf(y) + cimagf(x)*crealf(y);
    return z.u;
}
#endif
/*-----------------------------------------------
    Search for permutation in block/stream ordered 
    arrays
    Input:
    X[]  - 3x3 block/stream ordered data
    Output:
    permIx[L] - permutation index
-----------------------------------------------*/
void cmatinv3x3f_search(int16_t *permIx,const complex_float *X, int L, eLayout layout)
#if !__OPTIMIZED__
{
    int l;
    const int16_t *perm=matinv3x3sf_searchTbl;
    NASSERT (layout==e3x3_stream || layout==e3x3_block);

    for ( l=0; l<L; l++ )
    {
        complex_float x[9],a[8],b[8],c[8],d[8];
        int p, n;
        float32_t  w,maxW;
        int      maxIx;

        maxIx=permIx[l];
        if (layout==e3x3_stream)
        {
            for (n=0; n<4 ; n++)
            {
                x[n  ]= X[(l&~3)+n+L*((l&3)*2+0)];
                x[n+4]= X[(l&~3)+n+L*((l&3)*2+1)];
            }
            x[8]=X[8*L+l];
        }
        else
        {
            for (n=0; n<9 ; n++) x[n]=X[12*l+n];
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
            complex_float t;
            t=subc(mulc(a[n],d[n]),mulc(b[n],c[n]));
            w=fabsf(crealf(t)*crealf(t)+cimagf(t)*cimagf(t));
            if ( maxW < w ) { maxW = w, p = n; }
        }
        permIx[l]=p;
    }
}
#else
{
    int estride,estride2,estride8,mstride;
    xb_vecNx16 permA,permB,permC,permD;
    vselN_2 selA,selB,selC,selD;
    short * restrict pPermWr;
    const xb_vecN_2xf32* restrict pXrd;
    const xtfloat      * restrict pXrd8;
    const xb_vecNx16   * restrict perm;
    int l,addx;
    pXrd=(const xb_vecN_2xf32*)X;

    if (layout==e3x3_stream)
    {
        estride =  L*sizeof(complex_float);
        estride2=2*L*sizeof(complex_float);
        mstride = 2*BBE_SIMD_WIDTH;
        pXrd8=(     xtfloat*)(X+8*L);
        estride8=sizeof(complex_float);
    }
    else
    {
        NASSERT(layout==e3x3_block);
        estride =   2*BBE_SIMD_WIDTH;
        estride2=   6*BBE_SIMD_WIDTH;
        mstride =12*2*BBE_SIMD_WIDTH;
        pXrd8=( xtfloat*)(X+8);
        estride8=6*BBE_SIMD_WIDTH;
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
        xb_vecN_2xf32 reX0,imX0,reX1,imX1,W;
        xb_vecN_2xf32 imA,reA,imB,reB,imC,reC,imD,reD;
        xtfloat MAXW;
        vboolN_2 bmax;
        vselN vmax;
        int dummy;
        xb_vecNx16 maxix;

        addx = estride2;
        XT_MOVEQZ(addx, -3*estride2+mstride, (l&3)^3);

        imX0=BBE_LVN_2XF32_X (pXrd,estride);
        BBE_LVN_2XF32_XP(reX0,pXrd,addx);
        imX1=BBE_LSN_2XF32_I (pXrd8,sizeof(float32_t));
        BBE_LSN_2XF32_XP(reX1,pXrd8,estride8);
        BBE_DSELN_2XF32I(imX0,reX0,imX0,reX0,BBE_DSELI_DEINTERLEAVE_2);

        reA=BBE_SELN_2XF32(reX1,reX0,selA);
        reB=BBE_SELN_2XF32(reX1,reX0,selB);
        reC=BBE_SELN_2XF32(reX1,reX0,selC);
        reD=BBE_SELN_2XF32(reX1,reX0,selD);

        imA=BBE_SELN_2XF32(imX1,imX0,selA);
        imB=BBE_SELN_2XF32(imX1,imX0,selB);
        imC=BBE_SELN_2XF32(imX1,imX0,selC);
        imD=BBE_SELN_2XF32(imX1,imX0,selD);

        //Are*Dre-Aim*Dim
        reX0=BBE_MULN_2XF32(reA,reD);
        imX0=BBE_MULN_2XF32(reA,imD);
        BBE_MULSN_2XF32(reX0,imA,imD);
        BBE_MULAN_2XF32(imX0,imA,reD);

        reX1=BBE_MULN_2XF32(reB,reC);
        imX1=BBE_MULN_2XF32(reB,imC);
        BBE_MULSN_2XF32(reX1,imB,imC);
        BBE_MULAN_2XF32(imX1,imB,reC);

        reX0=BBE_SUBN_2XF32(reX0,reX1);
        imX0=BBE_SUBN_2XF32(imX0,imX1);

        W=BBE_MULN_2XF32(reX0,reX0);
        BBE_MULAN_2XF32(W,imX0,imX0);

        maxix=BBE_LVNX16_I(perm,4*2*BBE_SIMD_WIDTH);
        BBE_RBMAXNUMN_2XF32(bmax,MAXW,W);
        (void)MAXW;
        BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));(void)dummy;
        BBE_SSNX16_IP(BBE_SHFLNX16(maxix,vmax),pPermWr,sizeof(int16_t));
    }
}
#endif


#endif
