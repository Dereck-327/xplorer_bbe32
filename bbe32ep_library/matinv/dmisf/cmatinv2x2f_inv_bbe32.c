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
    Direct inversion of 2x2 complex floating point matrices, streaming data
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
#include "matinv2x2sf_inv.h"

#define __OPTIMIZED__ 1

#if HAVE_VFPU
#if !__OPTIMIZED__
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

static complex_float negc(complex_float x)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    z.s.re=-crealf(x);
    z.s.im=-cimagf(x);
    return z.u;
}

static complex_float recipc(complex_float x)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    float32_t d;
    d=crealf(x)*crealf(x) + cimagf(x)*cimagf(x);
    d=1.0f/d;
    z.s.re= crealf(x)*d;
    z.s.im=-cimagf(x)*d;
    return z.u;
}
#endif
/*-------------------------------------------------------------------------
Direct Matrix Inversion For Real Matrices without permutation
Input/output come in stream format

Input/output:
W  - BBE_SIMD_WIDTH/2 for real or BBE_SIMD_WIDTH/4 for complex matrices
X[L*2x2]  matrices 3x3 in the stream order
L         number of matrices

Temporary
pScr            scratch memory
 
Restrictions:
all matrices have to be aligned
L  should be a multiple of W
-------------------------------------------------------------------------*/
void cmatinv2x2f_inv(void * restrict pScr, 
                      complex_float* restrict X, 
                      int L, eLayout layout)
#if !__OPTIMIZED__
{
    int p,l,estride=0,mstride=0,sstride=0;
    NASSERT_ALIGN(X, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    switch(layout)
    {
    case e2x2_stream: estride=L;                sstride=2*L;                mstride=BBE_SIMD_WIDTH/4;    break;
    case e2x2_block:  estride=BBE_SIMD_WIDTH/4; sstride=2*BBE_SIMD_WIDTH/4; mstride=4*BBE_SIMD_WIDTH/4;  break;
    case e3x3_stream: estride=L;                sstride=3*L;                mstride=BBE_SIMD_WIDTH/4;    break;
    case e3x3_block:  estride=BBE_SIMD_WIDTH/4; sstride=3*BBE_SIMD_WIDTH/4; mstride=12*BBE_SIMD_WIDTH/4; break;
    case e4x4_stream: estride=L;                sstride=4*L;                mstride=BBE_SIMD_WIDTH/4;    break;
    case e4x4_block:  estride=BBE_SIMD_WIDTH/4; sstride=4*BBE_SIMD_WIDTH/4; mstride=16*BBE_SIMD_WIDTH/4; break;
    default: NASSERT(0);
    }
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/4),X+=mstride)
    {
        for (p=0; p<(BBE_SIMD_WIDTH/4); p++)
        {
            complex_float a,b,c,d,det,rdet;
            a=X[p+0];    
            b=X[p+estride]; 
            c=X[p+sstride]; 
            d=X[p+estride+sstride];
            det=subc(mulc(a,d),mulc(b,c));
            rdet=recipc(det);
            X[p+0]=mulc( d,rdet);          
            X[p+estride]=negc(mulc(b,rdet)); 
            X[p+sstride]=negc(mulc(c,rdet)); 
            X[p+estride+sstride]= mulc(a,rdet); 
        }
    }
}
#else
{
    xb_vecN_4xcf32 a,b,c,d,det,rdet;
    const xb_vecN_4xcf32 * restrict pArd;
          xb_vecN_4xcf32 * restrict pAwr;
    const xb_vecN_4xcf32 * restrict pRdetrd;
    const xb_vecN_4xcf32 * restrict pRdetrd1;
          xb_vecN_4xcf32 * restrict pRdetwr;
    int l,estride=0,mstride=0,sstride=0;
    NASSERT_ALIGN(X, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    if (L<=0) return;
    switch(layout)
    {
    case e2x2_stream: estride=L;                sstride=2*L;                mstride=BBE_SIMD_WIDTH/4;    break;
    case e2x2_block:  estride=BBE_SIMD_WIDTH/4; sstride=2*BBE_SIMD_WIDTH/4; mstride=4*BBE_SIMD_WIDTH/4;  break;
    case e3x3_stream: estride=L;                sstride=3*L;                mstride=BBE_SIMD_WIDTH/4;    break;
    case e3x3_block:  estride=BBE_SIMD_WIDTH/4; sstride=3*BBE_SIMD_WIDTH/4; mstride=12*BBE_SIMD_WIDTH/4; break;
    case e4x4_stream: estride=L;                sstride=4*L;                mstride=BBE_SIMD_WIDTH/4;    break;
    case e4x4_block:  estride=BBE_SIMD_WIDTH/4; sstride=4*BBE_SIMD_WIDTH/4; mstride=16*BBE_SIMD_WIDTH/4; break;
    default: NASSERT(0);
    }
    estride<<=3;
    sstride<<=3;
    mstride<<=3;

    pRdetwr=(      xb_vecN_4xcf32 *)pScr;
    pArd   =(const xb_vecN_4xcf32 *)X;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        b=BBE_LVN_4XCF32_X (pArd,estride);
        c=BBE_LVN_4XCF32_X (pArd,sstride);
        d=BBE_LVN_4XCF32_X (pArd,sstride+estride);
        BBE_LVN_4XCF32_XP(a,pArd,mstride);
        a=BBE_MULN_4XCF32(a,d);
        b=BBE_MULN_4XCF32(b,c);
        det=BBE_SUBN_4XCF32(a,b);
        BBE_SVN_4XCF32_IP(det,pRdetwr,2*BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
    pRdetrd= (const xb_vecN_4xcf32 *)pScr;
    pRdetrd1=(const xb_vecN_4xcf32 *)pScr;
    pRdetwr=(      xb_vecN_4xcf32 *)pScr;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        xb_vecN_4xcf32 det1;
        xb_vecN_2xf32 c,t8,t1,t2,t5;
        BBE_LVN_4XCF32_IP(det ,pRdetrd,2*BBE_SIMD_WIDTH);
        BBE_LVN_4XCF32_IP(det1,pRdetrd1,2*BBE_SIMD_WIDTH);
        c=BBE_MOVN_2XF32_FROMN_4XCF32(det);
        t8=BBE_MULN_2XF32(c, c);
        BBE_MULMASN_2XF32 (t8, c, c, 0, 3);
        t1=BBE_RECIP0N_2XF32 (t8);
        t2=BBE_CONSTN_2XF32 (1);
        BBE_MULSN_2XF32 (t2, t8, t1);
        BBE_MULANN_2XF32(t1, t1, t2);
        t5=BBE_CONSTN_2XF32 (1);
        BBE_MULSN_2XF32 (t5, t8, t1);
        BBE_MULANN_2XF32(t1, t1, t5);
        c=BBE_MOVN_2XF32_FROMN_4XCF32(det1);
        c=BBE_MULMN_2XF32(t1, c, 2, 12);
        rdet=BBE_MOVN_4XCF32_FROMN_2XF32(c);
        BBE_SVN_4XCF32_IP(rdet,pRdetwr,2*BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
    pArd=(const xb_vecN_4xcf32 *)X;
    pAwr=(      xb_vecN_4xcf32 *)X;
    pRdetrd=(const xb_vecN_4xcf32 *)pScr;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        xb_vecN_4xcf32 t;
        b=BBE_LVN_4XCF32_X (pArd,estride);
        c=BBE_LVN_4XCF32_X (pArd,sstride);
        d=BBE_LVN_4XCF32_X (pArd,sstride+estride);
        BBE_LVN_4XCF32_XP(a,pArd,mstride);
        BBE_LVN_4XCF32_IP(rdet,pRdetrd,2*BBE_SIMD_WIDTH);
        t=BBE_MULMN_4XCF32 (b,rdet, 3, 4); BBE_MULMASN_4XCF32 (t, b,rdet, 2, 11); b=t;
        t=BBE_MULMN_4XCF32 (c,rdet, 3, 4); BBE_MULMASN_4XCF32 (t, c,rdet, 2, 11); c=t;
        BBE_SVN_4XCF32_X (b,pAwr,estride);
        BBE_SVN_4XCF32_X (c,pAwr,sstride);
        BBE_SVN_4XCF32_X (BBE_MULN_4XCF32 (a,rdet),pAwr,sstride+estride);
        BBE_SVN_4XCF32_XP(BBE_MULN_4XCF32 (d,rdet),pAwr,mstride);
    }
}
#endif

/* Return the scratch area size, in bytes. */
size_t cmatinv2x2f_inv_getScratchSize ( int L )
{
  return XT_MAX(L,0)*sizeof(complex_float);
}

#endif
