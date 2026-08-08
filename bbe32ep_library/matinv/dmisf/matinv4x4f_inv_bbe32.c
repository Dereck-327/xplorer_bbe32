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
    Direct inversion of 4x4 floating point matrices, streaming data
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
/* inversion without permutation */

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
#include "matinv2x2sf_inv.h"
#include "matinv4x4sf_inv.h"

#define __OPTIMIZED__ 1

#if HAVE_VFPU

// c=a*b, c - matrix2x2, a,b -submatrices
static void mul2x2mss( float32_t * restrict c,
               const float32_t * restrict a,
               const float32_t * restrict b, int L ,int estride,int mstride,int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT_ALIGN(c,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(a,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(b,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),c+=sstride,a+=mstride,b+=mstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        c[p+estride*0]=a[p+estride*0]*b[p+estride*0]+a[p+estride*1]*b[p+estride*4];
        c[p+estride*1]=a[p+estride*0]*b[p+estride*1]+a[p+estride*1]*b[p+estride*5];
        c[p+estride*2]=a[p+estride*4]*b[p+estride*0]+a[p+estride*5]*b[p+estride*4];
        c[p+estride*3]=a[p+estride*4]*b[p+estride*1]+a[p+estride*5]*b[p+estride*5];
    }
}
#else
{
    const xb_vecN_2xf32* restrict pA;
    const xb_vecN_2xf32* restrict pB;
          xb_vecN_2xf32* restrict pC;
    int l;
    NASSERT_ALIGN(c,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(a,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(b,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    estride<<=2;
    mstride<<=2;
    sstride<<=2;
    pA=(const xb_vecN_2xf32*)a;
    pB=(const xb_vecN_2xf32*)b;
    pC=(      xb_vecN_2xf32*)c;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 a0,a1,a4,a5,b0,b1,b4,b5,c0,c1,c2,c3;

        a1=BBE_LVN_2XF32_X(pA,1*estride);
        a4=BBE_LVN_2XF32_X(pA,4*estride);
        a5=BBE_LVN_2XF32_X(pA,5*estride);
        BBE_LVN_2XF32_XP(a0,pA,mstride);

        b1=BBE_LVN_2XF32_X(pB,1*estride);
        b4=BBE_LVN_2XF32_X(pB,4*estride);
        b5=BBE_LVN_2XF32_X(pB,5*estride);
        BBE_LVN_2XF32_XP(b0,pB,mstride);

        c0=BBE_MULN_2XF32(a1,b4);
        c1=BBE_MULN_2XF32(a1,b5);
        c2=BBE_MULN_2XF32(a5,b4);
        c3=BBE_MULN_2XF32(a5,b5);
        BBE_MULAN_2XF32(c0,a0,b0);
        BBE_MULAN_2XF32(c1,a0,b1);
        BBE_MULAN_2XF32(c2,a4,b0);
        BBE_MULAN_2XF32(c3,a4,b1);

        BBE_SVN_2XF32_X (c1,pC,1*estride);
        BBE_SVN_2XF32_X (c2,pC,2*estride);
        BBE_SVN_2XF32_X (c3,pC,3*estride);
        BBE_SVN_2XF32_XP(c0,pC,sstride);
    }
}
#endif

// c=-a*b b - matrix2x2, c,a -submatrices
static void mul2x2ssm ( float32_t * restrict c,
               const float32_t * restrict a,
               const float32_t * restrict b, int L ,int estride,int mstride,int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT_ALIGN(c,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(a,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(b,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),c+=mstride,a+=mstride,b+=sstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        c[p+estride*0]=-(a[p+estride*0]*b[p+estride*0]+a[p+estride*1]*b[p+estride*2]);
        c[p+estride*1]=-(a[p+estride*0]*b[p+estride*1]+a[p+estride*1]*b[p+estride*3]);
        c[p+estride*4]=-(a[p+estride*4]*b[p+estride*0]+a[p+estride*5]*b[p+estride*2]);
        c[p+estride*5]=-(a[p+estride*4]*b[p+estride*1]+a[p+estride*5]*b[p+estride*3]);
    }
}
#else
{
    const xb_vecN_2xf32* restrict pA;
    const xb_vecN_2xf32* restrict pB;
          xb_vecN_2xf32* restrict pC;
    int l;
    NASSERT_ALIGN(c,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(a,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(b,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    pA=(const xb_vecN_2xf32*)a;
    pB=(const xb_vecN_2xf32*)b;
    pC=(      xb_vecN_2xf32*)c;
    estride<<=2;
    mstride<<=2;
    sstride<<=2;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 a0,a1,a4,a5,b0,b1,b2,b3,c0,c1,c4,c5;

        a1=BBE_LVN_2XF32_X(pA,1*estride);
        a4=BBE_LVN_2XF32_X(pA,4*estride);
        a5=BBE_LVN_2XF32_X(pA,5*estride);
        BBE_LVN_2XF32_XP(a0,pA,mstride);

        b1=BBE_LVN_2XF32_X(pB,1*estride);
        b2=BBE_LVN_2XF32_X(pB,2*estride);
        b3=BBE_LVN_2XF32_X(pB,3*estride);
        BBE_LVN_2XF32_XP(b0,pB,sstride);

        c0=BBE_MULMN_2XF32(a1,b2,3,12);
        c1=BBE_MULMN_2XF32(a1,b3,3,12);
        c4=BBE_MULMN_2XF32(a5,b2,3,12);
        c5=BBE_MULMN_2XF32(a5,b3,3,12);
        BBE_MULSN_2XF32(c0,a0,b0);
        BBE_MULSN_2XF32(c1,a0,b1);
        BBE_MULSN_2XF32(c4,a4,b0);
        BBE_MULSN_2XF32(c5,a4,b1);

        BBE_SVN_2XF32_X (c1,pC,1*estride);
        BBE_SVN_2XF32_X (c4,pC,4*estride);
        BBE_SVN_2XF32_X (c5,pC,5*estride);
        BBE_SVN_2XF32_XP(c0,pC,mstride);
    }
}
#endif

// c=-a*b, a - matrix, c,b - submatrices
static void mul2x2sms( float32_t * restrict c,
               const float32_t * restrict a,
               const float32_t * restrict b, int L ,int estride,int mstride,int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT_ALIGN(c,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(a,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(b,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),c+=mstride,a+=sstride,b+=mstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        c[p+estride*0]=-(a[p+estride*0]*b[p+estride*0]+a[p+estride*1]*b[p+estride*4]);
        c[p+estride*1]=-(a[p+estride*0]*b[p+estride*1]+a[p+estride*1]*b[p+estride*5]);
        c[p+estride*4]=-(a[p+estride*2]*b[p+estride*0]+a[p+estride*3]*b[p+estride*4]);
        c[p+estride*5]=-(a[p+estride*2]*b[p+estride*1]+a[p+estride*3]*b[p+estride*5]);
    }
}
#else
{
    const xb_vecN_2xf32* restrict pA;
    const xb_vecN_2xf32* restrict pB;
          xb_vecN_2xf32* restrict pC;
    int l;
    NASSERT_ALIGN(c,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(a,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(b,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    pA=(const xb_vecN_2xf32*)a;
    pB=(const xb_vecN_2xf32*)b;
    pC=(      xb_vecN_2xf32*)c;
    estride<<=2;
    mstride<<=2;
    sstride<<=2;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 a0,a1,a2,a3,b0,b1,b4,b5,c0,c1,c4,c5;

        a1=BBE_LVN_2XF32_X(pA,1*estride);
        a2=BBE_LVN_2XF32_X(pA,2*estride);
        a3=BBE_LVN_2XF32_X(pA,3*estride);
        BBE_LVN_2XF32_XP(a0,pA,sstride);

        b1=BBE_LVN_2XF32_X(pB,1*estride);
        b4=BBE_LVN_2XF32_X(pB,4*estride);
        b5=BBE_LVN_2XF32_X(pB,5*estride);
        BBE_LVN_2XF32_XP(b0,pB,mstride);

        c0=BBE_MULMN_2XF32(a1,b4,3,12);
        c1=BBE_MULMN_2XF32(a1,b5,3,12);
        c4=BBE_MULMN_2XF32(a3,b4,3,12);
        c5=BBE_MULMN_2XF32(a3,b5,3,12);
        BBE_MULSN_2XF32(c0,a0,b0);
        BBE_MULSN_2XF32(c1,a0,b1);
        BBE_MULSN_2XF32(c4,a2,b0);
        BBE_MULSN_2XF32(c5,a2,b1);

        BBE_SVN_2XF32_X (c1,pC,1*estride);
        BBE_SVN_2XF32_X (c4,pC,4*estride);
        BBE_SVN_2XF32_X (c5,pC,5*estride);
        BBE_SVN_2XF32_XP(c0,pC,mstride);
    }
}
#endif

// c=c-a*b, a - matrix, c,b - submatrices
static void mas2x2sms(float32_t* restrict c,const float32_t* restrict a, const float32_t* restrict b, int L, int estride,int mstride,int sstride)
#if !__OPTIMIZED__
{
    int p,l;
    NASSERT_ALIGN(c,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(a,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(b,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),c+=mstride,a+=sstride,b+=mstride)
    for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
    {
        c[p+estride*0]-=(a[p+estride*0]*b[p+estride*0]+a[p+estride*1]*b[p+estride*4]);
        c[p+estride*1]-=(a[p+estride*0]*b[p+estride*1]+a[p+estride*1]*b[p+estride*5]);
        c[p+estride*4]-=(a[p+estride*2]*b[p+estride*0]+a[p+estride*3]*b[p+estride*4]);
        c[p+estride*5]-=(a[p+estride*2]*b[p+estride*1]+a[p+estride*3]*b[p+estride*5]);
    }
}
#else
{
    const xb_vecN_2xf32* restrict pA;
    const xb_vecN_2xf32* restrict pB;
    const xb_vecN_2xf32* restrict pCrd;
          xb_vecN_2xf32* restrict pCwr;
    int l;
    NASSERT_ALIGN(c,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(a,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(b,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    pA  =(const xb_vecN_2xf32*)a;
    pB  =(const xb_vecN_2xf32*)b;
    pCrd=(const xb_vecN_2xf32*)c;
    pCwr=(      xb_vecN_2xf32*)c;
    estride<<=2;
    mstride<<=2;
    sstride<<=2;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 a0,a1,a2,a3,b0,b1,b4,b5,c0,c1,c4,c5;
        a1=BBE_LVN_2XF32_X(pA,1*estride);
        a2=BBE_LVN_2XF32_X(pA,2*estride);
        a3=BBE_LVN_2XF32_X(pA,3*estride);
        BBE_LVN_2XF32_XP(a0,pA,sstride);

        b1=BBE_LVN_2XF32_X(pB,1*estride);
        b4=BBE_LVN_2XF32_X(pB,4*estride);
        b5=BBE_LVN_2XF32_X(pB,5*estride);
        BBE_LVN_2XF32_XP(b0,pB,mstride);

        c1=BBE_LVN_2XF32_X(pCrd,1*estride);
        c4=BBE_LVN_2XF32_X(pCrd,4*estride);
        c5=BBE_LVN_2XF32_X(pCrd,5*estride);
        BBE_LVN_2XF32_XP(c0,pCrd,mstride);

        BBE_MULSN_2XF32(c0,a1,b4); BBE_MULSN_2XF32(c0,a0,b0); 
        BBE_MULSN_2XF32(c1,a1,b5); BBE_MULSN_2XF32(c1,a0,b1); 
        BBE_MULSN_2XF32(c4,a3,b4); BBE_MULSN_2XF32(c4,a2,b0); 
        BBE_MULSN_2XF32(c5,a3,b5); BBE_MULSN_2XF32(c5,a2,b1); 

        BBE_SVN_2XF32_X (c1,pCwr,1*estride);
        BBE_SVN_2XF32_X (c4,pCwr,4*estride);
        BBE_SVN_2XF32_X (c5,pCwr,5*estride);
        BBE_SVN_2XF32_XP(c0,pCwr,mstride);
    }
}
#endif

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
X[L*4x4]  matrices 4x4 in the stream order
L         number of matrices

Temporary
pScr            scratch memory
 
Restrictions:
all matrices have to be aligned
L  should be a multiple of W
-------------------------------------------------------------------------*/
void  matinv4x4f_inv(void * restrict pScr, 
                     float32_t* restrict X, 
                     int L,
                     eLayout layout)
{
    float32_t *s=(float32_t *)pScr; // [4*L]
    float32_t *tmp=s+4*L;   // [L]
    int estride,mstride,sstride;
    NASSERT_ALIGN(X, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0) return;

    if (layout==e4x4_stream)
    {
        estride=L; mstride=BBE_SIMD_WIDTH/2; sstride=BBE_SIMD_WIDTH/2; 
    }
    else
    {
        NASSERT(layout==e4x4_block);
        estride=BBE_SIMD_WIDTH/2; mstride=16*BBE_SIMD_WIDTH/2; sstride=4*BBE_SIMD_WIDTH/2;
    }

    matinv2x2f_inv(tmp, X+estride* 0 ,L,layout);
    mul2x2mss( s           , X+estride* 8, X+estride* 0,L,estride,mstride,sstride);
    mas2x2sms( X+estride*10, s           , X+estride* 2,L,estride,mstride,sstride);
    matinv2x2f_inv   ( tmp         , X+estride*10              ,L,layout);
    mul2x2ssm( X+estride* 8, X+estride*10, s           ,L,estride,mstride,sstride);
    mul2x2mss( s           , X+estride* 0, X+estride* 2,L,estride,mstride,sstride);
    mas2x2sms( X+estride* 0, s           , X+estride* 8,L,estride,mstride,sstride);
    mul2x2sms( X+estride* 2, s           , X+estride*10,L,estride,mstride,sstride);
}

/* Return the scratch area size, in bytes. */
size_t matinv4x4f_inv_getScratchSize ( int L )
{
    size_t sz=0;
    (void)L;
    sz+=4*L*sizeof(float32_t);
    sz+=matinv2x2f_inv_getScratchSize(L);
    return sz;
}
#endif
