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
#if HAVE_VFPU

#define __OPTIMIZED__ 1

#define TRUNCATED_SEARCH 1

#if !__OPTIMIZED__
#include <math.h>
#endif

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
void matinv4x4f_truncatedSearch(int16_t* permIx, const float32_t *X, int L, eLayout layout)
#if !__OPTIMIZED__
{
    int l;
    for ( l=0; l<L; l++ )
    {
        float32_t x[16],a[8],b[8],c[8],d[8];
        const int16_t *perm;
        int p, n;
        float32_t  w,maxW;
        int      maxIx;

        maxIx=permIx[l];
        if (layout==e4x4_stream)
        {
           for (n=0; n<8 ; n++)
           {
               x[n  ]= X[(l&~7)+n+L*((l&7)*2+0)];
               x[n+8]= X[(l&~7)+n+L*((l&7)*2+1)];
           }
        }
        else
        {
            NASSERT(layout==e4x4_block);
            for (n=0; n<16; n++) x[n]=X[16*l+n];
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
            w=fabsf(a[n]*d[n]-b[n]*c[n]);
            if ( maxW < w ) { maxW = w, p = n; }
        }
        permIx[l]=perm[16*4+p*2];
    }
}
#else
{
    short * restrict pPermWr;
    const xb_vecN_2xf32* restrict pXrd;
    const xb_vecNx16   * restrict perm;
    int l,addx,inc0,inc1,inc2;
    pXrd=(const xb_vecN_2xf32*)X;
    pPermWr=(short*)permIx;
    if(layout==e4x4_block)
    {
        inc0 =inc1=4*BBE_SIMD_WIDTH; inc2=2*BBE_SIMD_WIDTH; 
    }
    else
    {
        NASSERT(layout==e4x4_stream);
        inc0 =L*2*sizeof(float32_t); inc1=(-14*L+8)*sizeof(float32_t); inc2=L*sizeof(float32_t); 
    }
    for ( l=0; l<L; l++ )
    {
        xb_vecN_2xf32 X0,X1,A,B,C,D,W;
        xb_vecNx16 permA,permB,permC,permD;
        xtfloat MAXW;
        vboolN_2 bmax;
        vselN vmax;
        int dummy,maxIx;
        xb_vecNx16 maxix;

        addx = inc0;
        XT_MOVEQZ(addx, inc1, (l&7)^7);

        X1=BBE_LVN_2XF32_X(pXrd,inc2);
        BBE_LVN_2XF32_XP(X0,pXrd,addx);
        maxIx=permIx[0];
        perm=(const xb_vecNx16*)&matinv4x4sf_searchTbl[16*5*maxIx];
        permA=BBE_LVNX16_I(perm,0*2*BBE_SIMD_WIDTH);
        permB=BBE_LVNX16_I(perm,1*2*BBE_SIMD_WIDTH);
        permC=BBE_LVNX16_I(perm,2*2*BBE_SIMD_WIDTH);
        permD=BBE_LVNX16_I(perm,3*2*BBE_SIMD_WIDTH);
        maxix=BBE_LVNX16_I(perm,4*2*BBE_SIMD_WIDTH);
        A=BBE_SELN_2XF32(X1,X0,BBE_MOVVSELN_2NX16(permA,0));
        B=BBE_SELN_2XF32(X1,X0,BBE_MOVVSELN_2NX16(permB,0));
        C=BBE_SELN_2XF32(X1,X0,BBE_MOVVSELN_2NX16(permC,0));
        D=BBE_SELN_2XF32(X1,X0,BBE_MOVVSELN_2NX16(permD,0));
        W=BBE_MULN_2XF32(A,D);
        BBE_MULSN_2XF32(W,B,C);
        W=BBE_ABSN_2XF32(W);
        BBE_RBMAXNUMN_2XF32(bmax,MAXW,W);
        (void)MAXW;
        BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));(void)dummy;
        maxix=BBE_SELNX16(maxix,maxix,vmax);
        BBE_SSNX16_IP(maxix,pPermWr,sizeof(int16_t));
        permIx++;
    }
}
#endif
#endif
