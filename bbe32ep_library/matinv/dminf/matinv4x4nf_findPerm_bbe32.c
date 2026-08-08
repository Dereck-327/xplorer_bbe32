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

#if !__OPTIMIZED__ 
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
static void findUL(int16_t *pos,const float *X, int L)
#if !__OPTIMIZED__ 
{
    int l;
    int        n,maxIx;
    float32_t  w,maxW;
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);

    for ( l=0; l<L; l++ )
    {
        float32_t x[16];
        for (n=0; n<16; n++)
        {
           x[n+ 0]= X[l*16+n];
        }
        for ( maxW=0.f, maxIx=0, n=0; n<16; n++ )
        {
            float32_t t=x[n];
            w=fabsf(t);
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
    int l;
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);

    ix0 = BBE_SRLIN_2XC16(BBE_MOVN_2XC16_FROMNX16(BBE_SEQNX16()),1);
    ix8 = BBE_SRLIN_2XC16(BBE_MOVN_2XC16_FROMNX16(BBE_ADDNX16(BBE_SEQNX16(),BBE_SIMD_WIDTH)),1);

    for ( l=0; l<L; l++ )
    {
        xb_vecN_2xf32 x0,x1;
        xb_vecN_2xc16 maxix;
        vboolN_2 bmax;
        int dummy;
        xtfloat MAXW;
        vselN vmax;

        BBE_LVN_2XF32_IP(x0,pX,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(x1,pX,2*BBE_SIMD_WIDTH);
        x0=BBE_ABSN_2XF32(x0);
        x1=BBE_ABSN_2XF32(x1);
        maxix=BBE_MOVN_2XC16T(ix8,ix0,BBE_OGTN_2XF32(x1,x0)); 
        x0=BBE_MAXNUMN_2XF32(x0,x1);
        BBE_RBMAXNUMN_2XF32(bmax,MAXW,x0);
        (void)MAXW;
        BBE_SQZN(vmax, dummy, BBE_MOVN_FROMN_2(bmax));(void)dummy;
        BBE_SSNX16_IP(BBE_SHFLNX16(BBE_MOVNX16_FROMN_2XC16(maxix),vmax),pPos,sizeof(int16_t));
    }
}
#endif
/*-------------------------------------------------------------------------
    find permutation for L matrices written in block/stream order
    Input:
    X[L*16]  input matrices
    L        number of matrices
    Output:
    permIx[L] permutation indices
-------------------------------------------------------------------------*/
void matinv4x4nf_findPerm(int16_t *permIx,const float32_t *X, int L)
{
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    findUL(permIx,X,L);
    matinv4x4f_truncatedSearch(permIx, X, L, e4x4_block);
}
#endif
