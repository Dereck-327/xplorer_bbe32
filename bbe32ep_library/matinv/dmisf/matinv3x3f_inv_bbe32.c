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
#include <math.h>
#include "matinv3x3sf_inv.h"
#include "matinv2x2sf_inv.h"
#if HAVE_VFPU

#define __OPTIMIZED__ 1

#define COMBINED_PHASE1 0
#define COMBINED_PHASE2 1
#define COMBINED_PHASE3 1


// In-place 1x1 submatrix inversion
static void inv1x1x( float32_t* A , int L,int mstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),A+=mstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        A[p]=1.0f/A[p];
    }
}
#else
{
    xb_vecN_2xf32 a;
    const xb_vecN_2xf32 * restrict pArd;
          xb_vecN_2xf32 * restrict pAwr;
    int l;
    NASSERT_ALIGN(A, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);

    pArd=(const xb_vecN_2xf32 *)A;
    pAwr=(      xb_vecN_2xf32 *)A;
    mstride<<=2;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        BBE_LVN_2XF32_XP(a,pArd,mstride);
        a=BBE_RECIPN_2XF32(a);
        BBE_SVN_2XF32_XP(a,pAwr,mstride);
    }
}
#endif

#if COMBINED_PHASE2
static void phase2( float32_t * restrict X,
                    float32_t * restrict s,
                    int L, int estride,int mstride,int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),X+=mstride,s+=sstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        X[p+estride*6]=-X[p+estride*8]*s[p+estride*0];
        X[p+estride*7]=-X[p+estride*8]*s[p+estride*1];
        s[p+estride*0]=X[p+estride*0]*X[p+estride*2]+X[p+estride*1]*X[p+estride*5];
        s[p+estride*1]=X[p+estride*3]*X[p+estride*2]+X[p+estride*4]*X[p+estride*5];
    }
}
#else
{
    const xb_vecN_2xf32 * restrict pX=(const xb_vecN_2xf32 *)X;
    const xb_vecN_2xf32 * restrict pS=(const xb_vecN_2xf32 *)s;
          xb_vecN_2xf32 * restrict pY=(      xb_vecN_2xf32 *)(X+6*estride);
          xb_vecN_2xf32 * restrict pR=(      xb_vecN_2xf32 *)s;
    int l;
    NASSERT_ALIGN(pX, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pS, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pY, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pR, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    estride<<=2;
    mstride<<=2;
    sstride<<=2;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 x0,x1,x2,x3,x4,x5,x6,x7,x8,s0,s1;
        x1=BBE_LVN_2XF32_X (pX,1*estride);
        x2=BBE_LVN_2XF32_X (pX,2*estride);
        x3=BBE_LVN_2XF32_X (pX,3*estride);
        x4=BBE_LVN_2XF32_X (pX,4*estride);
        x5=BBE_LVN_2XF32_X (pX,5*estride);
        x8=BBE_LVN_2XF32_X (pX,8*estride);
        BBE_LVN_2XF32_XP(x0,pX,mstride);
        s1=BBE_LVN_2XF32_X (pS,estride);
        BBE_LVN_2XF32_XP(s0,pS,sstride);
        x6=BBE_MULMN_2XF32(x8,s0,3,12);
        x7=BBE_MULMN_2XF32(x8,s1,3,12);
        s0=BBE_MULN_2XF32(x1,x5);
        s1=BBE_MULN_2XF32(x4,x5);
        BBE_MULAN_2XF32(s0,x0,x2);
        BBE_MULAN_2XF32(s1,x3,x2);
        BBE_SVN_2XF32_X (x7,pY,1*estride);
        BBE_SVN_2XF32_XP(x6,pY,mstride);
        BBE_SVN_2XF32_X (s1,pR,1*estride);
        BBE_SVN_2XF32_XP(s0,pR,sstride);
    }
}
#endif

#else //COMBINED_PHASE2
// c=-a*b submatrices 1x1 by 1x2 
static void mul1x1x2x( float32_t * restrict c,
               const float32_t * restrict a,
               const float32_t * restrict b , int L, int estride,int mstride,int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),c+=mstride,a+=mstride,b+=sstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        c[p+estride*0]=-a[p+estride*0]*b[p+estride*0];
        c[p+estride*1]=-a[p+estride*0]*b[p+estride*1];
    }
}
#else
{
    const xb_vecN_2xf32 * restrict pA=(const xb_vecN_2xf32 *)a;
    const xb_vecN_2xf32 * restrict pB=(const xb_vecN_2xf32 *)b;
          xb_vecN_2xf32 * restrict pC=(      xb_vecN_2xf32 *)c;
    int l;
    NASSERT_ALIGN(pA, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pB, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pC, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    estride<<=2;
    mstride<<=2;
    sstride<<=2;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 a0,b0,b1,c0,c1;

        BBE_LVN_2XF32_XP(a0,pA,mstride);

        b1=BBE_LVN_2XF32_X(pB,1*estride);
        BBE_LVN_2XF32_XP(b0,pB,sstride);

        c0=BBE_MULMN_2XF32(a0,b0,3,12);
        c1=BBE_MULMN_2XF32(a0,b1,3,12);

        BBE_SVN_2XF32_X (c1,pC,1*estride);
        BBE_SVN_2XF32_XP(c0,pC,mstride);
    }
}
#endif

// c=a*b submatrices 2x2 by 2x1
static void mul2x2x1x( float32_t * restrict c,
               const float32_t * restrict a,
               const float32_t * restrict b , int L, int estride,int mstride,int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),c+=sstride,a+=mstride,b+=mstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        c[p+estride*0]=a[p+estride*0]*b[p+estride*0]+a[p+estride*1]*b[p+estride*3];
        c[p+estride*1]=a[p+estride*3]*b[p+estride*0]+a[p+estride*4]*b[p+estride*3];
    }
}
#else
{
    const xb_vecN_2xf32 * restrict pA=(const xb_vecN_2xf32 *)a;
    const xb_vecN_2xf32 * restrict pB=(const xb_vecN_2xf32 *)b;
          xb_vecN_2xf32 * restrict pC=(      xb_vecN_2xf32 *)c;
    int l;
    NASSERT_ALIGN(pA, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pB, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pC, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    estride<<=2;
    mstride<<=2;
    sstride<<=2;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 a0,a1,a3,a4,b0,b3,c0,c1;

        a1=BBE_LVN_2XF32_X(pA,1*estride);
        a3=BBE_LVN_2XF32_X(pA,3*estride);
        a4=BBE_LVN_2XF32_X(pA,4*estride);
        BBE_LVN_2XF32_XP(a0,pA,mstride);

        b3=BBE_LVN_2XF32_X(pB,3*estride);
        BBE_LVN_2XF32_XP(b0,pB,mstride);

        c0=BBE_MULN_2XF32(a0,b0);
        c1=BBE_MULN_2XF32(a3,b0);
        BBE_MULAN_2XF32(c0,a1,b3);
        BBE_MULAN_2XF32(c1,a4,b3);

        BBE_SVN_2XF32_X (c1,pC,1*estride);
        BBE_SVN_2XF32_XP(c0,pC,sstride);
    }
}
#endif
#endif //COMBINED_PHASE2

#if COMBINED_PHASE1
static void phase1(float32_t *s, float32_t* X,int L, int estride,int mstride,int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),X+=mstride,s+=sstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        s[p+estride*0]=X[p+estride*6]*X[p+estride*0]+X[p+estride*7]*X[p+estride*3];
        s[p+estride*1]=X[p+estride*6]*X[p+estride*1]+X[p+estride*7]*X[p+estride*4];
        X[p+estride*8]-=(s[p+estride*0]*X[p+estride*2]+s[p+estride*1]*X[p+estride*5]);
    }

}
#else
{
    const xb_vecN_2xf32 * restrict pX=(const xb_vecN_2xf32 *)(X);
    const xb_vecN_2xf32 * restrict pS=(const xb_vecN_2xf32 *)(s);
          xb_vecN_2xf32 * restrict pY=(      xb_vecN_2xf32 *)(X+8*estride);

    int l;
    NASSERT_ALIGN(pX, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pS, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pY, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    estride<<=2;
    mstride<<=2;
    sstride<<=2;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 x0,x1,x2,x3,x4,x5,x6,x7,x8,s0,s1;
        x1=BBE_LVN_2XF32_X (pX,1*estride);
        x2=BBE_LVN_2XF32_X (pX,2*estride);
        x3=BBE_LVN_2XF32_X (pX,3*estride);
        x4=BBE_LVN_2XF32_X (pX,4*estride);
        x5=BBE_LVN_2XF32_X (pX,5*estride);
        x6=BBE_LVN_2XF32_X (pX,6*estride);
        x7=BBE_LVN_2XF32_X (pX,7*estride);
        x8=BBE_LVN_2XF32_X (pX,8*estride);
        BBE_LVN_2XF32_XP(x0,pX,mstride);
        s0=BBE_MULN_2XF32(x3,x7);
        s1=BBE_MULN_2XF32(x4,x7);
        BBE_MULAN_2XF32(s0,x0,x6);
        BBE_MULAN_2XF32(s1,x1,x6);
        BBE_MULSN_2XF32(x8,s0,x2);
        BBE_MULSN_2XF32(x8,s1,x5);
        BBE_SVN_2XF32_XP(x8,pY,mstride);
        BBE_SVN_2XF32_X (s1,pS,1*estride);
        BBE_SVN_2XF32_XP(s0,pS,sstride);
    }
}
#endif
#else // COMBINED_PHASE1
/* multiply submatrices 1x2 by 2x2 */
static void mul1x2x2x(float32_t * restrict c,
               const float32_t * restrict a,
               const float32_t * restrict b , int L, int estride,int mstride,int sstride)
#if !__OPTIMIZED
{
    int p,l;
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),c+=sstride,a+=mstride,b+=mstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        c[p+estride*0]=a[p+estride*0]*b[p+estride*0]+a[p+estride*1]*b[p+estride*3];
        c[p+estride*1]=a[p+estride*0]*b[p+estride*1]+a[p+estride*1]*b[p+estride*4];
    }
}
#else
{
    const xb_vecN_2xf32 * restrict pA=(const xb_vecN_2xf32 *)a;
    const xb_vecN_2xf32 * restrict pB=(const xb_vecN_2xf32 *)b;
          xb_vecN_2xf32 * restrict pC=(      xb_vecN_2xf32 *)c;
    int l;
    NASSERT_ALIGN(pA, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pB, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pC, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    estride<<=2;
    mstride<<=2;
    sstride<<=2;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 a0,a1,b0,b1,b3,b4,c0,c1;

        a1=BBE_LVN_2XF32_X(pA,1*estride);
        BBE_LVN_2XF32_XP(a0,pA,mstride);

        b1=BBE_LVN_2XF32_X(pB,1*estride);
        b3=BBE_LVN_2XF32_X(pB,3*estride);
        b4=BBE_LVN_2XF32_X(pB,4*estride);
        BBE_LVN_2XF32_XP(b0,pB,mstride);

        c0=BBE_MULN_2XF32(a0,b0);
        c1=BBE_MULN_2XF32(a0,b1);
        BBE_MULAN_2XF32(c0,a1,b3);
        BBE_MULAN_2XF32(c1,a1,b4);

        BBE_SVN_2XF32_X (c1,pC,1*estride);
        BBE_SVN_2XF32_XP(c0,pC,sstride);
    }
}
#endif
// c=c-a*b, submatrices a[1x2], b[2x1]
static void mas1x2x1x(float32_t* restrict c,
                   const float32_t* restrict a, 
                   const float32_t* restrict b, int L, int estride,int mstride,int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),c+=mstride,a+=sstride,b+=mstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        c[p+estride*0]-=(a[p+estride*0]*b[p+estride*0]+a[p+estride*1]*b[p+estride*3]);
    }
}
#else
{
    const xb_vecN_2xf32 * restrict pA=(const xb_vecN_2xf32 *)a;
    const xb_vecN_2xf32 * restrict pB=(const xb_vecN_2xf32 *)b;
    const xb_vecN_2xf32 * restrict pCrd=(const xb_vecN_2xf32 *)c;
          xb_vecN_2xf32 * restrict pC=(      xb_vecN_2xf32 *)c;
    int l;
    NASSERT_ALIGN(pA, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pB, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pC, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    estride<<=2;
    mstride<<=2;
    sstride<<=2;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 a0,a1,b0,b3,c0;

        BBE_LVN_2XF32_XP(c0,pCrd,mstride);

        a1=BBE_LVN_2XF32_X(pA,1*estride);
        BBE_LVN_2XF32_XP(a0,pA,sstride);

        b3=BBE_LVN_2XF32_X(pB,3*estride);
        BBE_LVN_2XF32_XP(b0,pB,mstride);

        BBE_MULSN_2XF32(c0,a0,b0);
        BBE_MULSN_2XF32(c0,a1,b3);

        BBE_SVN_2XF32_XP(c0,pC,mstride);
    }
}
#endif
#endif // COMBINED_PHASE1


#if COMBINED_PHASE3
static void phase3(float32_t* restrict X,
                   const float32_t* restrict s, 
                   int L, int estride,int mstride,int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),X+=mstride,s+=sstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        X[p+estride*0]-=s[p+estride*0]*X[p+estride*6];
        X[p+estride*1]-=s[p+estride*0]*X[p+estride*7];
        X[p+estride*2]=-s[p+estride*0]*X[p+estride*8];
        X[p+estride*3]-=s[p+estride*1]*X[p+estride*6];
        X[p+estride*4]-=s[p+estride*1]*X[p+estride*7];
        X[p+estride*5]=-s[p+estride*1]*X[p+estride*8];
    }
}
#else
{
    const xb_vecN_2xf32 * restrict pX=(const xb_vecN_2xf32 *)X;
    const xb_vecN_2xf32 * restrict pS=(const xb_vecN_2xf32 *)s;
          xb_vecN_2xf32 * restrict pY=(      xb_vecN_2xf32 *)X;
    int l;
    NASSERT_ALIGN(pX, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pS, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pY, (2 * BBE_SIMD_WIDTH));

    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    estride<<=2;
    mstride<<=2;
    sstride<<=2;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 x0,x1,x2,x3,x4,x5,x6,x7,x8,s0,s1;
        x1=BBE_LVN_2XF32_X (pX,1*estride);
        x3=BBE_LVN_2XF32_X (pX,3*estride);
        x4=BBE_LVN_2XF32_X (pX,4*estride);
        x6=BBE_LVN_2XF32_X (pX,6*estride);
        x7=BBE_LVN_2XF32_X (pX,7*estride);
        x8=BBE_LVN_2XF32_X (pX,8*estride);
        BBE_LVN_2XF32_XP(x0,pX,mstride);
        s1=BBE_LVN_2XF32_X (pS,1*estride);
        BBE_LVN_2XF32_XP(s0,pS,sstride);
        BBE_MULSN_2XF32(x0,s0,x6);
        BBE_MULSN_2XF32(x1,s0,x7);
        x2=BBE_MULMN_2XF32(s0,x8,3,12);
        BBE_MULSN_2XF32(x3,s1,x6);
        BBE_MULSN_2XF32(x4,s1,x7);
        x5=BBE_MULMN_2XF32(s1,x8,3,12);
        BBE_SVN_2XF32_X (x1,pY,1*estride);
        BBE_SVN_2XF32_X (x2,pY,2*estride);
        BBE_SVN_2XF32_X (x3,pY,3*estride);
        BBE_SVN_2XF32_X (x4,pY,4*estride);
        BBE_SVN_2XF32_X (x5,pY,5*estride);
        BBE_SVN_2XF32_XP(x0,pY,mstride);
    }
}
#endif

#else // COMBINED_PHASE3

// c=-a*b submatrices 2x1 by 1x1
static void mul2x1x1x( float32_t * restrict c,
               const float32_t * restrict a,
               const float32_t * restrict b , int L, int estride,int mstride,int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),c+=mstride,a+=sstride,b+=mstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        c[p+estride*0]=-a[p+estride*0]*b[p+estride*0];
        c[p+estride*3]=-a[p+estride*1]*b[p+estride*0];
    }
}
#else
{
    const xb_vecN_2xf32 * restrict pA=(const xb_vecN_2xf32 *)a;
    const xb_vecN_2xf32 * restrict pB=(const xb_vecN_2xf32 *)b;
          xb_vecN_2xf32 * restrict pC=(      xb_vecN_2xf32 *)c;
    int l;
    NASSERT_ALIGN(pA, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pB, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pC, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    estride<<=2;
    mstride<<=2;
    sstride<<=2;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 a0,a1,b0,c0,c3;

        a1=BBE_LVN_2XF32_X(pA,1*estride);
        BBE_LVN_2XF32_XP(a0,pA,sstride);

        BBE_LVN_2XF32_XP(b0,pB,mstride);

        c0=BBE_MULMN_2XF32(a0,b0,3,12);
        c3=BBE_MULMN_2XF32(a1,b0,3,12);

        BBE_SVN_2XF32_X (c3,pC,3*estride);
        BBE_SVN_2XF32_XP(c0,pC,mstride);
    }
}
#endif

// c=c-a*b submatrices 2x1 by 1x2
static void mas2x1x2x(float32_t* restrict c,
                   const float32_t* restrict a, 
                   const float32_t* restrict b, int L, int estride,int mstride,int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),c+=mstride,a+=sstride,b+=mstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        c[p+estride*0]-=a[p+estride*0]*b[p+estride*0];
        c[p+estride*1]-=a[p+estride*0]*b[p+estride*1];
        c[p+estride*3]-=a[p+estride*1]*b[p+estride*0];
        c[p+estride*4]-=a[p+estride*1]*b[p+estride*1];
    }
}
#else
{
    const xb_vecN_2xf32 * restrict pA=(const xb_vecN_2xf32 *)a;
    const xb_vecN_2xf32 * restrict pB=(const xb_vecN_2xf32 *)b;
    const xb_vecN_2xf32 * restrict pCrd=(const xb_vecN_2xf32 *)c;
          xb_vecN_2xf32 * restrict pC=(      xb_vecN_2xf32 *)c;
    int l;
    NASSERT_ALIGN(pA, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pB, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pC, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    estride<<=2;
    mstride<<=2;
    sstride<<=2;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 a0,a1,b0,b1,c0,c1,c3,c4;

        c1=BBE_LVN_2XF32_X (pCrd,1*estride);
        c3=BBE_LVN_2XF32_X (pCrd,3*estride);
        c4=BBE_LVN_2XF32_X (pCrd,4*estride);
        BBE_LVN_2XF32_XP(c0,pCrd,mstride);

        a1=BBE_LVN_2XF32_X(pA,1*estride);
        BBE_LVN_2XF32_XP(a0,pA,sstride);

        b1=BBE_LVN_2XF32_X(pB,1*estride);
        BBE_LVN_2XF32_XP(b0,pB,mstride);

        BBE_MULSN_2XF32(c0,a0,b0);
        BBE_MULSN_2XF32(c1,a0,b1);
        BBE_MULSN_2XF32(c3,a1,b0);
        BBE_MULSN_2XF32(c4,a1,b1);

        BBE_SVN_2XF32_X (c1,pC,1*estride);
        BBE_SVN_2XF32_X (c3,pC,3*estride);
        BBE_SVN_2XF32_X (c4,pC,4*estride);
        BBE_SVN_2XF32_XP(c0,pC,mstride);
    }
}
#endif

#endif //COMBINED_PHASE3
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
void matinv3x3f_inv (void * restrict pScr, 
                     float32_t* restrict X, 
                     int L, eLayout layout)
{
    float32_t *s=(float32_t *)pScr;
    float32_t *temp=s+2*L;
    int estride,mstride,sstride;
    NASSERT_ALIGN(X, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0) return;
    if (layout==e3x3_stream)
    {
        estride=L; mstride=BBE_SIMD_WIDTH/2; sstride=BBE_SIMD_WIDTH/2; 
    }
    else
    {
        NASSERT(layout==e3x3_block); 
        estride=BBE_SIMD_WIDTH/2; mstride=16*BBE_SIMD_WIDTH/2; sstride=2*BBE_SIMD_WIDTH/2;
    }

    matinv2x2f_inv(temp,X+estride* 0,L,layout);
#if COMBINED_PHASE1
    phase1(s,X,L,estride,mstride,sstride);
#else
    mul1x2x2x( s           , X+estride* 6, X+estride* 0 ,L,estride,mstride,sstride);
    mas1x2x1x( X+estride* 8, s           , X+estride* 2 ,L,estride,mstride,sstride);
#endif
    inv1x1x  ( X+estride* 8,L,mstride);       
#if COMBINED_PHASE2
    phase2(X,s,L,estride,mstride,sstride);
#else
    mul1x1x2x( X+estride* 6, X+estride* 8, s            ,L,estride,mstride,sstride); 
    mul2x2x1x( s           , X+estride* 0, X+estride* 2 ,L,estride,mstride,sstride);
#endif
#if COMBINED_PHASE3
    phase3(X,s,L,estride,mstride,sstride);
#else
    mas2x1x2x( X+estride* 0, s, X+estride* 6 ,L,estride,mstride,sstride);
    mul2x1x1x( X+estride* 2, s, X+estride* 8 ,L,estride,mstride,sstride); 
#endif
}

size_t  matinv3x3f_inv_getScratchSize ( int L )
{
    size_t sz=0;
    (void)L;
    sz+=2*L*sizeof(float32_t)+matinv2x2f_inv_getScratchSize(L);
    return sz;
}
#endif
