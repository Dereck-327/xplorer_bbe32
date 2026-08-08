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
static complex_float addc(complex_float x,complex_float y)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    z.s.re=crealf(x)+crealf(y);
    z.s.im=cimagf(x)+cimagf(y);
    return z.u;
}
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
#define TRUNCATED_SEARCH 1
/*------------------------------------------
vectorized search algorithm for BBE32EP
(truncated variant - look for 8 possible
combinations of possible 9 for given upper 
left position
Input:
X[16*L]   - L matrices
permIx[L] - upper left positions
Output:
permIx[L] - permutation index from 
matinv4x4_fwd_perm_tbl[]
------------------------------------------*/
void cmatinv4x4f_truncatedSearch(int16_t* permIx, const complex_float *X, int L, eLayout layout)
#if !__OPTIMIZED__
{
    int l;
    for ( l=0; l<L; l++ )
    {
        complex_float x[16],a[8],b[8],c[8],d[8];
        const int16_t *perm;
        int p, n;
        float32_t  w,maxW;
        int      maxIx;

        maxIx=permIx[l];
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
           NASSERT(layout==e4x4_block);
           for (n=0; n<16; n++) x[n]=X[l*16+n];
       }

        perm=&matinv4x4sf_searchTbl[16*5*maxIx];
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
        permIx[l]=perm[16*4+p*2];
    }
}
#else
{
    short* restrict pPermWr;
    const xb_vecN_2xf32* restrict pXrd;
    const xb_vecNx16   * restrict pTbl;
    int l,addx,estride,mstride;
    pXrd=(const xb_vecN_2xf32*)X;
    pPermWr=(short*)permIx;
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
        xb_vecN_2xf32 imX0,reX0,imX1,reX1,A,B,C,D,W;
        xb_vecN_2xf32 imA,reA,imB,reB,imC,reC,imD,reD;

        xb_vecNx16 permA,permB,permC,permD;
        xtfloat MAXW;
        vboolN_2 bmax;
        vselN vmax;
        int dummy,maxIx;
        xb_vecNx16 maxix;

        maxIx=permIx[0];permIx++;
        addx = 4*estride;
        XT_MOVEQZ(addx, -12*estride+mstride, (l&3)^3);
        B=BBE_LVN_2XF32_X (pXrd,1*estride);
        C=BBE_LVN_2XF32_X (pXrd,2*estride);
        D=BBE_LVN_2XF32_X (pXrd,3*estride);
        BBE_LVN_2XF32_XP(A,pXrd,addx);
        BBE_DSELN_2XF32I(imX0,reX0,B,A,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(imX1,reX1,D,C,BBE_DSELI_DEINTERLEAVE_2);

        pTbl=(const xb_vecNx16*)&matinv4x4sf_searchTbl[16*5*maxIx];
        permA=BBE_LVNX16_I(pTbl,0*2*BBE_SIMD_WIDTH);
        permB=BBE_LVNX16_I(pTbl,1*2*BBE_SIMD_WIDTH);
        permC=BBE_LVNX16_I(pTbl,2*2*BBE_SIMD_WIDTH);
        permD=BBE_LVNX16_I(pTbl,3*2*BBE_SIMD_WIDTH);
        maxix=BBE_LVNX16_I(pTbl,4*2*BBE_SIMD_WIDTH);
        reA=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permA,0));
        reB=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permB,0));
        reC=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permC,0));
        reD=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permD,0));

        imA=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permA,0));
        imB=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permB,0));
        imC=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permC,0));
        imD=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permD,0));

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

        BBE_RBMAXNUMN_2XF32(bmax,MAXW,W);
        (void)MAXW;
        BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));(void)dummy;
        maxix=BBE_SELNX16(maxix,maxix,vmax);
        maxix=maxix;
        BBE_SSNX16_IP(maxix,pPermWr,sizeof(int16_t));
    }
}
#endif

#endif
