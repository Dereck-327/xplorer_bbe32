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

/* half-stream DMI without permutations */

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
#include "matinv3x3sf_inv.h"
#include "matinv2x2sf_inv.h"

#define COMBINED_PHASE1 0
#define COMBINED_PHASE2 1
#define COMBINED_PHASE3 1
#define __OPTIMIZED__ 1

#if HAVE_VFPU
#if !__OPTIMIZED__
#include <math.h>
#include <complex.h>
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


// In-place 1x1 submatrix inversion
static void inv1x1x( complex_float* A , int L, int mstride)
#if !__OPTIMIZED__
{
    int l,p;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4,A+=mstride)
    {
        for (p=0; p<BBE_SIMD_WIDTH/4; p++)
        A[p+0]=recipc(A[p+0]);
    }
}
#else
{
    int l;
    const xb_vecN_4xcf32* restrict pArd0;
    const xb_vecN_4xcf32* restrict pArd1;
          xb_vecN_4xcf32* restrict pAwr;
    pArd0=(const xb_vecN_4xcf32 *)A;
    pArd1=(const xb_vecN_4xcf32 *)A;
    pAwr =(      xb_vecN_4xcf32 *)A;
    mstride<<=3;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        xb_vecN_4xcf32 det,det1,rdet;
        xb_vecN_2xf32 c,t8,t1,t2,t5;
        BBE_LVN_4XCF32_XP(det ,pArd0,mstride);
        BBE_LVN_4XCF32_XP(det1,pArd1,mstride);
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
        BBE_SVN_4XCF32_XP(rdet,pAwr,mstride);
    }
}
#endif

#if COMBINED_PHASE1
static void phase1( complex_float * restrict X,
                    complex_float * restrict s, int L, int estride, int mstride, int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT((L%(BBE_SIMD_WIDTH/4))==0);
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(s,2*BBE_SIMD_WIDTH);
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4,X+=mstride,s+=sstride)
    for (p=0; p<BBE_SIMD_WIDTH/4; p++)
    {
        s[p+estride*0]=addc(mulc(X[p+estride*6],X[p+estride*0]),mulc(X[p+estride*7],X[p+estride*3]));
        s[p+estride*1]=addc(mulc(X[p+estride*6],X[p+estride*1]),mulc(X[p+estride*7],X[p+estride*4]));
        X[p+estride*8]=subc(X[p+estride*8],addc(mulc(s[p+estride*0],X[p+estride*2]),mulc(s[p+estride*1],X[p+estride*5])));
    }
}
#else
{
    const xb_vecN_4xcf32 * restrict pX=(const xb_vecN_4xcf32 *)(X);
    const xb_vecN_4xcf32 * restrict pS=(const xb_vecN_4xcf32 *)(s);
          xb_vecN_4xcf32 * restrict pY=(      xb_vecN_4xcf32 *)(X+estride*8);

    int l;
    NASSERT_ALIGN(pX, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pS, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pY, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    estride<<=3;
    mstride<<=3;
    sstride<<=3;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        xb_vecN_4xcf32 x0,x1,x2,x3,x4,x5,x6,x7,x8,s0,s1;
        x1=BBE_LVN_4XCF32_X (pX,1*estride);
        x2=BBE_LVN_4XCF32_X (pX,2*estride);
        x3=BBE_LVN_4XCF32_X (pX,3*estride);
        x4=BBE_LVN_4XCF32_X (pX,4*estride);
        x5=BBE_LVN_4XCF32_X (pX,5*estride);
        x6=BBE_LVN_4XCF32_X (pX,6*estride);
        x7=BBE_LVN_4XCF32_X (pX,7*estride);
        x8=BBE_LVN_4XCF32_X (pX,8*estride);
        BBE_LVN_4XCF32_XP(x0,pX,mstride);
        s0=BBE_MULN_4XCF32(x3,x7);
        s1=BBE_MULN_4XCF32(x4,x7);
        BBE_MULAN_4XCF32(s0,x0,x6);
        BBE_MULAN_4XCF32(s1,x1,x6);
        BBE_MULSN_4XCF32(x8,s0,x2);
        BBE_MULSN_4XCF32(x8,s1,x5);
        BBE_SVN_4XCF32_XP(x8,pY,mstride);
        BBE_SVN_4XCF32_X (s1,pS,1*estride);
        BBE_SVN_4XCF32_XP(s0,pS,sstride);
    }
}
#endif
#else //COMBINED_PHASE1
/* multiply submatrices 1x2 by 2x2 */
static void mul1x2x2x(complex_float * restrict c,
                const complex_float * restrict a,
                const complex_float * restrict b , int L, int estride, int mstride, int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4,c+=sstride,a+=mstride,b+=mstride)
    for (p=0; p<BBE_SIMD_WIDTH/4; p++)
    {
        c[p+estride*0]=addc(mulc(a[p+estride*0],b[p+estride*0]),mulc(a[p+estride*1],b[p+estride*3]));
        c[p+estride*1]=addc(mulc(a[p+estride*0],b[p+estride*1]),mulc(a[p+estride*1],b[p+estride*4]));
    }
}
#else
{
    const xb_vecN_4xcf32 * restrict pA=(const xb_vecN_4xcf32 *)a;
    const xb_vecN_4xcf32 * restrict pB=(const xb_vecN_4xcf32 *)b;
          xb_vecN_4xcf32 * restrict pC=(      xb_vecN_4xcf32 *)c;
    int l;
    NASSERT_ALIGN(pA, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pB, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pC, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    estride<<=3;
    mstride<<=3;
    sstride<<=3;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        xb_vecN_4xcf32 a0,a1,b0,b1,b3,b4,c0,c1;

        a1=BBE_LVN_4XCF32_X(pA,1*estride);
        BBE_LVN_4XCF32_XP(a0,pA,mstride);

        b1=BBE_LVN_4XCF32_X(pB,1*estride);
        b3=BBE_LVN_4XCF32_X(pB,3*estride);
        b4=BBE_LVN_4XCF32_X(pB,4*estride);
        BBE_LVN_4XCF32_XP(b0,pB,mstride);

        c0=BBE_MULN_4XCF32(a0,b0);
        c1=BBE_MULN_4XCF32(a0,b1);
        BBE_MULAN_4XCF32(c0,a1,b3);
        BBE_MULAN_4XCF32(c1,a1,b4);

        BBE_SVN_4XCF32_X (c1,pC,1*estride);
        BBE_SVN_4XCF32_XP(c0,pC,sstride);
    }
}
#endif

// c=c-a*b, submatrices a[1x2], b[2x1]
static void mas1x2x1x(complex_float* restrict c,
                   const complex_float* restrict a, 
                   const complex_float* restrict b, int L, int estride, int mstride, int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4,c+=mstride,a+=sstride,b+=mstride)
    for (p=0; p<BBE_SIMD_WIDTH/4; p++)
    {
        c[p+estride*0]=subc(c[p+estride*0],addc(mulc(a[p+estride*0],b[p+estride*0]),mulc(a[p+estride*1],b[p+estride*3])));
    }
}
#else
{
    const xb_vecN_4xcf32 * restrict pA  =(const xb_vecN_4xcf32 *)a;
    const xb_vecN_4xcf32 * restrict pB  =(const xb_vecN_4xcf32 *)b;
    const xb_vecN_4xcf32 * restrict pCrd=(const xb_vecN_4xcf32 *)c;
          xb_vecN_4xcf32 * restrict pC  =(      xb_vecN_4xcf32 *)c;
    int l;
    NASSERT_ALIGN(pA, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pB, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pC, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    estride<<=3;
    mstride<<=3;
    sstride<<=3;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        xb_vecN_4xcf32 a0,a1,b0,b3,c0;

        BBE_LVN_4XCF32_XP(c0,pCrd,mstride);

        a1=BBE_LVN_4XCF32_X(pA,1*estride);
        BBE_LVN_4XCF32_XP(a0,pA,sstride);

        b3=BBE_LVN_4XCF32_X(pB,3*estride);
        BBE_LVN_4XCF32_XP(b0,pB,mstride);

        BBE_MULSN_4XCF32(c0,a0,b0);
        BBE_MULSN_4XCF32(c0,a1,b3);

        BBE_SVN_4XCF32_XP(c0,pC,mstride);
    }
}
#endif

#endif //COMBINED_PHASE1

#if COMBINED_PHASE2
static void phase2( complex_float * restrict X,
                    complex_float * restrict s, int L, int estride, int mstride, int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT((L%(BBE_SIMD_WIDTH/4))==0);
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(s,2*BBE_SIMD_WIDTH);
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4,X+=mstride,s+=sstride)
    for (p=0; p<BBE_SIMD_WIDTH/4; p++)
    {
        X[p+estride*6]=negc(mulc(X[p+estride*8],s[p+estride*0]));
        X[p+estride*7]=negc(mulc(X[p+estride*8],s[p+estride*1]));
        s[p+estride*0]=addc(mulc(X[p+estride*0],X[p+estride*2]),mulc(X[p+estride*1],X[p+estride*5]));
        s[p+estride*1]=addc(mulc(X[p+estride*3],X[p+estride*2]),mulc(X[p+estride*4],X[p+estride*5]));
    }
}
#else
{
    const xb_vecN_4xcf32 * restrict pX=(const xb_vecN_4xcf32 *)X;
    const xb_vecN_4xcf32 * restrict pS=(const xb_vecN_4xcf32 *)s;
          xb_vecN_4xcf32 * restrict pY=(      xb_vecN_4xcf32 *)(X+6*estride);
          xb_vecN_4xcf32 * restrict pR=(      xb_vecN_4xcf32 *)s;
    int l;
    NASSERT_ALIGN(pX, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pS, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pY, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pR, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    estride<<=3;
    mstride<<=3;
    sstride<<=3;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        xb_vecN_4xcf32 x0,x1,x2,x3,x4,x5,x6,x7,x8,s0,s1;
        x1=BBE_LVN_4XCF32_X (pX,1*estride);
        x2=BBE_LVN_4XCF32_X (pX,2*estride);
        x3=BBE_LVN_4XCF32_X (pX,3*estride);
        x4=BBE_LVN_4XCF32_X (pX,4*estride);
        x5=BBE_LVN_4XCF32_X (pX,5*estride);
        x8=BBE_LVN_4XCF32_X (pX,8*estride);
        BBE_LVN_4XCF32_XP(x0,pX,mstride);
        s1=BBE_LVN_4XCF32_X (pS,1*estride);
        BBE_LVN_4XCF32_XP(s0,pS,sstride);
        x6=BBE_MULMN_4XCF32(x8,s0,3,4); BBE_MULMASN_4XCF32(x6,x8,s0,2,11);
        x7=BBE_MULMN_4XCF32(x8,s1,3,4); BBE_MULMASN_4XCF32(x7,x8,s1,2,11);
        s0=BBE_MULN_4XCF32(x1,x5);
        s1=BBE_MULN_4XCF32(x4,x5);
        BBE_MULAN_4XCF32(s0,x0,x2);
        BBE_MULAN_4XCF32(s1,x3,x2);
        BBE_SVN_4XCF32_X (x7,pY,estride);
        BBE_SVN_4XCF32_XP(x6,pY,mstride);
        BBE_SVN_4XCF32_X (s1,pR,estride);
        BBE_SVN_4XCF32_XP(s0,pR,sstride);
    }
}
#endif
#else //COMBINED_PHASE2

// c=a*b submatrices 1x1 by 1x2 
static void mul1x1x2x( complex_float * restrict c,
                 const complex_float * restrict a,
                 const complex_float * restrict b , int L, int estride, int mstride, int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4,c+=mstride,a+=mstride,b+=sstride)
    for (p=0; p<BBE_SIMD_WIDTH/4; p++)
    {
        c[p+estride*0]=negc(mulc(a[p+estride*0],b[p+estride*0]));
        c[p+estride*1]=negc(mulc(a[p+estride*0],b[p+estride*1]));
    }
}
#else
{
    const xb_vecN_4xcf32 * restrict pA=(const xb_vecN_4xcf32 *)a;
    const xb_vecN_4xcf32 * restrict pB=(const xb_vecN_4xcf32 *)b;
          xb_vecN_4xcf32 * restrict pC=(      xb_vecN_4xcf32 *)c;
    int l;
    NASSERT_ALIGN(pA, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pB, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pC, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    estride<<=3;
    mstride<<=3;
    sstride<<=3;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        xb_vecN_4xcf32 a0,b0,b1,c0,c1;

        BBE_LVN_4XCF32_XP(a0,pA,mstride);

        b1=BBE_LVN_4XCF32_X(pB,1*estride);
        BBE_LVN_4XCF32_XP(b0,pB,sstride);

        c0=BBE_MULMN_4XCF32(a0,b0,3,4); BBE_MULMASN_4XCF32(c0,a0,b0,2,11);
        c1=BBE_MULMN_4XCF32(a0,b1,3,4); BBE_MULMASN_4XCF32(c1,a0,b1,2,11);

        BBE_SVN_4XCF32_X (c1,pC,1*estride);
        BBE_SVN_4XCF32_XP(c0,pC,mstride);
    }
}
#endif

// c=a*b submatrices 2x2 by 2x1
static void mul2x2x1x( complex_float * restrict c,
                 const complex_float * restrict a,
                 const complex_float * restrict b , int L, int estride, int mstride, int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4,c+=sstride,a+=mstride,b+=mstride)
    for (p=0; p<BBE_SIMD_WIDTH/4; p++)
    {
        c[p+estride*0]=addc(mulc(a[p+estride*0],b[p+estride*0]),mulc(a[p+estride*1],b[p+estride*3]));
        c[p+estride*1]=addc(mulc(a[p+estride*3],b[p+estride*0]),mulc(a[p+estride*4],b[p+estride*3]));
    }
}
#else
{
    const xb_vecN_4xcf32 * restrict pA=(const xb_vecN_4xcf32 *)a;
    const xb_vecN_4xcf32 * restrict pB=(const xb_vecN_4xcf32 *)b;
          xb_vecN_4xcf32 * restrict pC=(      xb_vecN_4xcf32 *)c;
    int l;
    NASSERT_ALIGN(pA, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pB, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pC, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    estride<<=3;
    mstride<<=3;
    sstride<<=3;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        xb_vecN_4xcf32 a0,a1,a3,a4,b0,b3,c0,c1;

        a1=BBE_LVN_4XCF32_X(pA,1*estride);
        a3=BBE_LVN_4XCF32_X(pA,3*estride);
        a4=BBE_LVN_4XCF32_X(pA,4*estride);
        BBE_LVN_4XCF32_XP(a0,pA,mstride);

        b3=BBE_LVN_4XCF32_X(pB,3*estride);
        BBE_LVN_4XCF32_XP(b0,pB,mstride);

        c0=BBE_MULN_4XCF32(a0,b0);
        c1=BBE_MULN_4XCF32(a3,b0);
        BBE_MULAN_4XCF32(c0,a1,b3);
        BBE_MULAN_4XCF32(c1,a4,b3);

        BBE_SVN_4XCF32_X (c1,pC,1*estride);
        BBE_SVN_4XCF32_XP(c0,pC,sstride);
    }
}
#endif
#endif//COMBINED_PHASE2

#if COMBINED_PHASE3
static void phase3( complex_float *X, complex_float *s, int L, int estride, int mstride, int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT((L%(BBE_SIMD_WIDTH/4))==0);
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(s,2*BBE_SIMD_WIDTH);
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4,X+=mstride,s+=sstride)
    for (p=0; p<BBE_SIMD_WIDTH/4; p++)
    {
        X[p+estride*0]=subc(X[p+estride*0],mulc(s[p+estride*0],X[p+estride*6]));
        X[p+estride*1]=subc(X[p+estride*1],mulc(s[p+estride*0],X[p+estride*7]));
        X[p+estride*3]=subc(X[p+estride*3],mulc(s[p+estride*1],X[p+estride*6]));
        X[p+estride*4]=subc(X[p+estride*4],mulc(s[p+estride*1],X[p+estride*7]));
        X[p+estride*2]=negc(mulc(s[p+estride*0],X[p+estride*8]));
        X[p+estride*5]=negc(mulc(s[p+estride*1],X[p+estride*8]));
    }
}
#else
{
    const xb_vecN_4xcf32 * restrict pX=(const xb_vecN_4xcf32 *)X;
    const xb_vecN_4xcf32 * restrict pS=(const xb_vecN_4xcf32 *)s;
          xb_vecN_4xcf32 * restrict pY=(      xb_vecN_4xcf32 *)X;
    int l;
    NASSERT_ALIGN(pX, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pS, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pY, (2 * BBE_SIMD_WIDTH));

    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    estride<<=3;
    mstride<<=3;
    sstride<<=3;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        xb_vecN_4xcf32 x0,x1,x2,x3,x4,x5,x6,x7,x8,s0,s1;
        x1=BBE_LVN_4XCF32_X (pX,1*estride);
        x3=BBE_LVN_4XCF32_X (pX,3*estride);
        x4=BBE_LVN_4XCF32_X (pX,4*estride);
        x6=BBE_LVN_4XCF32_X (pX,6*estride);
        x7=BBE_LVN_4XCF32_X (pX,7*estride);
        x8=BBE_LVN_4XCF32_X (pX,8*estride);
        BBE_LVN_4XCF32_XP(x0,pX,mstride);
        s1=BBE_LVN_4XCF32_X (pS,1*estride);
        BBE_LVN_4XCF32_XP(s0,pS,sstride);
        BBE_MULSN_4XCF32(x0,s0,x6);
        BBE_MULSN_4XCF32(x1,s0,x7);
        x2=BBE_MULMN_4XCF32(s0,x8,3,4); BBE_MULMASN_4XCF32(x2,s0,x8,2,11);
        BBE_MULSN_4XCF32(x3,s1,x6);
        BBE_MULSN_4XCF32(x4,s1,x7);
        x5=BBE_MULMN_4XCF32(s1,x8,3,4); BBE_MULMASN_4XCF32(x5,s1,x8,2,11);
        BBE_SVN_4XCF32_X (x1,pY,1*estride);
        BBE_SVN_4XCF32_X (x2,pY,2*estride);
        BBE_SVN_4XCF32_X (x3,pY,3*estride);
        BBE_SVN_4XCF32_X (x4,pY,4*estride);
        BBE_SVN_4XCF32_X (x5,pY,5*estride);
        BBE_SVN_4XCF32_XP(x0,pY,mstride);
    }
}
#endif
#else //COMBINED_PHASE3
// c=a*b submatrices 2x1 by 1x1
static void mul2x1x1x( complex_float * restrict c,
                 const complex_float * restrict a,
                 const complex_float * restrict b , int L, int estride, int mstride, int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4,c+=mstride,a+=sstride,b+=mstride)
    for (p=0; p<BBE_SIMD_WIDTH/4; p++)
    {
        c[p+estride*0]=negc(mulc(a[p+estride*0],b[p+estride*0]));
        c[p+estride*3]=negc(mulc(a[p+estride*1],b[p+estride*0]));
    }
}
#else
{
    const xb_vecN_4xcf32 * restrict pA=(const xb_vecN_4xcf32 *)a;
    const xb_vecN_4xcf32 * restrict pB=(const xb_vecN_4xcf32 *)b;
          xb_vecN_4xcf32 * restrict pC=(      xb_vecN_4xcf32 *)c;
    int l;
    NASSERT_ALIGN(pA, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pB, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pC, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    estride<<=3;
    mstride<<=3;
    sstride<<=3;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        xb_vecN_4xcf32 a0,a1,b0,c0,c3;

        a1=BBE_LVN_4XCF32_X(pA,1*estride);
        BBE_LVN_4XCF32_XP(a0,pA,sstride);

        BBE_LVN_4XCF32_XP(b0,pB,mstride);

        c0=BBE_MULMN_4XCF32(a0,b0,3,4); BBE_MULMASN_4XCF32(c0,a0,b0,2,11);
        c3=BBE_MULMN_4XCF32(a1,b0,3,4); BBE_MULMASN_4XCF32(c3,a1,b0,2,11);

        BBE_SVN_4XCF32_X (c3,pC,3*estride);
        BBE_SVN_4XCF32_XP(c0,pC,mstride);
    }
}
#endif

// c=c-a*b submatrices 2x1 by 1x2
static void mas2x1x2x(complex_float* restrict c,
                  const complex_float* restrict a, 
                  const complex_float* restrict b, int L, int estride, int mstride, int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4,c+=mstride,a+=sstride,b+=mstride)
    for (p=0; p<BBE_SIMD_WIDTH/4; p++)
    {
        c[p+estride*0]=subc(c[p+estride*0],mulc(a[p+estride*0],b[p+estride*0]));
        c[p+estride*1]=subc(c[p+estride*1],mulc(a[p+estride*0],b[p+estride*1]));
        c[p+estride*3]=subc(c[p+estride*3],mulc(a[p+estride*1],b[p+estride*0]));
        c[p+estride*4]=subc(c[p+estride*4],mulc(a[p+estride*1],b[p+estride*1]));
    }
}
#else
{
    const xb_vecN_4xcf32 * restrict pA  =(const xb_vecN_4xcf32 *)a;
    const xb_vecN_4xcf32 * restrict pB  =(const xb_vecN_4xcf32 *)b;
    const xb_vecN_4xcf32 * restrict pCrd=(const xb_vecN_4xcf32 *)c;
          xb_vecN_4xcf32 * restrict pC  =(      xb_vecN_4xcf32 *)c;
    int l;
    NASSERT_ALIGN(pA, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pB, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pC, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    estride<<=3;
    mstride<<=3;
    sstride<<=3;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++)
    {
        xb_vecN_4xcf32 a0,a1,b0,b1,c0,c1,c3,c4;

        c1=BBE_LVN_4XCF32_X (pCrd,1*estride);
        c3=BBE_LVN_4XCF32_X (pCrd,3*estride);
        c4=BBE_LVN_4XCF32_X (pCrd,4*estride);
        BBE_LVN_4XCF32_XP(c0,pCrd,mstride);

        a1=BBE_LVN_4XCF32_X(pA,1*estride);
        BBE_LVN_4XCF32_XP(a0,pA,sstride);

        b1=BBE_LVN_4XCF32_X(pB,1*estride);
        BBE_LVN_4XCF32_XP(b0,pB,mstride);

        BBE_MULSN_4XCF32(c0,a0,b0);
        BBE_MULSN_4XCF32(c1,a0,b1);
        BBE_MULSN_4XCF32(c3,a1,b0);
        BBE_MULSN_4XCF32(c4,a1,b1);

        BBE_SVN_4XCF32_X (c1,pC,1*estride);
        BBE_SVN_4XCF32_X (c3,pC,3*estride);
        BBE_SVN_4XCF32_X (c4,pC,4*estride);
        BBE_SVN_4XCF32_XP(c0,pC,mstride);
    }
}
#endif
#endif//COMBINED_PHASE3

/*-------------------------------------------------------------------------
Direct Matrix Inversion For Real Matrices without permutation
Input/output come in stream format

Matlab formulas:
A=inv(A);
T=C*A;
D=D-T*B;
D=inv(D);
C=-D*T;
T=A*B;
A=A-T*C;
B=-T*D;

Input/output:
W  - BBE_SIMD_WIDTH/2 for real or BBE_SIMD_WIDTH/4 for complex matrices
X[L*3x3]  matrices 3x3 in the stream order
L         number of matrices

Temporary
pScr            scratch memory
 
Restrictions:
all matrices have to be aligned
L  should be a multiple of W
-------------------------------------------------------------------------*/
void cmatinv3x3f_inv(void * restrict pScr, 
                     complex_float* restrict X, 
                     int L, eLayout layout)
{
    int sstride=0,estride=0,mstride=0;
    complex_float *s=(complex_float *)pScr;
    complex_float *temp=s+2*L;
    NASSERT_ALIGN(X, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    if (L<=0) return;

    if (layout==e3x3_stream)
    {
        estride=L;                mstride=BBE_SIMD_WIDTH/4;    sstride=BBE_SIMD_WIDTH/4;
    }
    else
    {
        NASSERT(layout==e3x3_block);  
        estride=BBE_SIMD_WIDTH/4; mstride=12*BBE_SIMD_WIDTH/4; sstride=2*BBE_SIMD_WIDTH/4;
    }

    //inv2x2x  ( temp, X+estride* 0,L,estride,mstride);
    cmatinv2x2f_inv(temp,X,L,layout);
#if COMBINED_PHASE1
    phase1(X,s,L,estride,mstride,sstride);
#else
    mul1x2x2x( s     , X+estride* 6, X+estride* 0 ,L,estride,mstride,sstride);
    mas1x2x1x( X+estride* 8, s     , X+estride* 2 ,L,estride,mstride,sstride);
#endif
    inv1x1x  ( X+estride* 8, L,mstride);
#if COMBINED_PHASE2
    phase2(X,s,L,estride,mstride,sstride);
#else
    mul1x1x2x( X+estride* 6, X+estride* 8, s            ,L,estride,mstride,sstride); 
    mul2x2x1x( s           , X+estride* 0, X+estride* 2 ,L,estride,mstride,sstride);
#endif
#if COMBINED_PHASE3
    phase3( X, s, L,estride,mstride,sstride);
#else
    mas2x1x2x( X+estride* 0, s     , X+estride* 6 ,L,estride,mstride,sstride);
    mul2x1x1x( X+estride* 2, s     , X+estride* 8 ,L,estride,mstride,sstride); 
#endif
}

size_t  cmatinv3x3f_inv_getScratchSize ( int L )
{
    size_t sz=0;
    (void)L;
    sz+=2*L*sizeof(complex_float)+cmatinv2x2f_inv_getScratchSize(L);
    return sz;
}
#endif
