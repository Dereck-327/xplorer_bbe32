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
#include <math.h>
#include "matinv2x2sf_inv.h"
#if HAVE_VFPU

#define __OPTIMIZED__ 1

#define SPLIT_CYCLE 2

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
void matinv2x2f_inv (void * restrict pScr, 
                     float32_t* restrict X, 
                     int L, eLayout layout)
#if !__OPTIMIZED__
{
    int p,l,estride=0,mstride=0,sstride=0;
    float32_t a,b,c,d,det,rdet;
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    switch(layout)
    {
    case e2x2_stream: estride=L;                sstride=2*L;                mstride=BBE_SIMD_WIDTH/2;    break;
    case e2x2_block:  estride=BBE_SIMD_WIDTH/2; sstride=2*BBE_SIMD_WIDTH/2; mstride=4*BBE_SIMD_WIDTH/2;  break;
    case e3x3_stream: estride=L;                sstride=3*L;                mstride=BBE_SIMD_WIDTH/2;    break;
    case e3x3_block:  estride=BBE_SIMD_WIDTH/2; sstride=3*BBE_SIMD_WIDTH/2; mstride=16*BBE_SIMD_WIDTH/2; break;
    case e4x4_stream: estride=L;                sstride=4*L;                mstride=BBE_SIMD_WIDTH/2;    break;
    case e4x4_block:  estride=BBE_SIMD_WIDTH/2; sstride=4*BBE_SIMD_WIDTH/2; mstride=16*BBE_SIMD_WIDTH/2; break;
    default: NASSERT(0);
    }
    for (l=0; l<L; l+=(BBE_SIMD_WIDTH/2),X+=mstride)
    {
        for (p=0; p<(BBE_SIMD_WIDTH/2); p++)
        {
            a=X[p+0              ];    
            b=X[p+estride        ]; 
            c=X[p+sstride        ]; 
            d=X[p+estride+sstride];
            det=a*d-b*c;
            rdet=1.0f/det;
            X[p+0              ]= d*rdet; 
            X[p+estride        ]=-b*rdet; 
            X[p+sstride        ]=-c*rdet; 
            X[p+estride+sstride]= a*rdet; 
        }
    }
}
#else
{
    xb_vecN_2xf32 a,b,c,d,det,rdet;
    const xb_vecN_2xf32 * restrict pArd;
          xb_vecN_2xf32 * restrict pAwr;
    const xb_vecN_2xf32 * restrict pRdetrd;
          xb_vecN_2xf32 * restrict pRdetwr;
    int l,estride=0,mstride=0,sstride=0;
    NASSERT_ALIGN(X,    (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0) return;
    switch(layout)
    {
    case e2x2_stream: estride=L;                sstride=2*L;                mstride=BBE_SIMD_WIDTH/2;    break;
    case e2x2_block:  estride=BBE_SIMD_WIDTH/2; sstride=2*BBE_SIMD_WIDTH/2; mstride=4*BBE_SIMD_WIDTH/2;  break;
    case e3x3_stream: estride=L;                sstride=3*L;                mstride=BBE_SIMD_WIDTH/2;    break;
    case e3x3_block:  estride=BBE_SIMD_WIDTH/2; sstride=3*BBE_SIMD_WIDTH/2; mstride=16*BBE_SIMD_WIDTH/2; break;
    case e4x4_stream: estride=L;                sstride=4*L;                mstride=BBE_SIMD_WIDTH/2;    break;
    case e4x4_block:  estride=BBE_SIMD_WIDTH/2; sstride=4*BBE_SIMD_WIDTH/2; mstride=16*BBE_SIMD_WIDTH/2; break;
    default: NASSERT(0);
    }
    estride<<=2;
    sstride<<=2;
    mstride<<=2;
#if SPLIT_CYCLE==0
    pArd=(const xb_vecN_2xf32 *)X;
    pAwr=(      xb_vecN_2xf32 *)X;
    (void)pRdetrd; (void)pRdetwr;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 a,b,c,d,det,rdet;
        b=BBE_LVN_2XF32_X (pArd,estride);
        c=BBE_LVN_2XF32_X (pArd,sstride);
        d=BBE_LVN_2XF32_X (pArd,sstride+estride);
        BBE_LVN_2XF32_IP(a,pArd,mstride);
        det=BBE_MULN_2XF32(a,d);
        BBE_MULSN_2XF32(det,b,c);
        rdet=BBE_RECIPN_2XF32(det);
        BBE_SVN_2XF32_X (BBE_MULMN_2XF32(b,rdet,3,12),pAwr,estride);
        BBE_SVN_2XF32_X (BBE_MULMN_2XF32(c,rdet,3,12),pAwr,sstride);
        BBE_SVN_2XF32_X (BBE_MULN_2XF32 (a,rdet     ),pAwr,sstride+estride);
        BBE_SVN_2XF32_IP(BBE_MULN_2XF32 (d,rdet     ),pAwr,mstride);
    }
#elif SPLIT_CYCLE==1
    pRdetwr=(      xb_vecN_2xf32 *)pScr;
    pArd=(const xb_vecN_2xf32 *)X;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        b=BBE_LVN_2XF32_X(pArd,estride);
        c=BBE_LVN_2XF32_X(pArd,sstride);
        d=BBE_LVN_2XF32_X(pArd,sstride+estride);
        BBE_LVN_2XF32_IP(a,pArd,mstride);
        det=BBE_MULN_2XF32(a,d);
        BBE_MULSN_2XF32(det,b,c);
        rdet=BBE_RECIPN_2XF32(det);
        BBE_SVN_2XF32_IP(rdet,pRdetwr,2*BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
    pRdetrd=(const xb_vecN_2xf32 *)pScr;
    pArd=(const xb_vecN_2xf32 *)X;
    pAwr=(      xb_vecN_2xf32 *)X;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        b=BBE_LVN_2XF32_X(pArd,estride);
        c=BBE_LVN_2XF32_X(pArd,sstride);
        d=BBE_LVN_2XF32_X(pArd,sstride+estride);
        BBE_LVN_2XF32_IP(a,pArd,mstride);
        BBE_LVN_2XF32_IP(rdet,pRdetrd,2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_X (BBE_MULMN_2XF32(b,rdet,3,12),pAwr,estride);
        BBE_SVN_2XF32_X (BBE_MULMN_2XF32(c,rdet,3,12),pAwr,sstride);
        BBE_SVN_2XF32_X (BBE_MULN_2XF32 (a,rdet     ),pAwr,sstride+estride);
        BBE_SVN_2XF32_IP(BBE_MULN_2XF32 (d,rdet     ),pAwr,mstride);
    }
#elif SPLIT_CYCLE==2

    pRdetwr=(      xb_vecN_2xf32 *)pScr;
    pArd   =(const xb_vecN_2xf32 *)X;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        b=BBE_LVN_2XF32_X(pArd,estride);
        c=BBE_LVN_2XF32_X(pArd,sstride);
        d=BBE_LVN_2XF32_X(pArd,sstride+estride);
        BBE_LVN_2XF32_XP(a,pArd,mstride);
        det=BBE_MULN_2XF32(a,d);
        BBE_MULSN_2XF32(det,b,c);
        BBE_SVN_2XF32_IP(det,pRdetwr,2*BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
    pRdetrd=(const xb_vecN_2xf32 *)pScr;
    pRdetwr=(      xb_vecN_2xf32 *)pScr;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        BBE_LVN_2XF32_IP(det,pRdetrd,2*BBE_SIMD_WIDTH);
        rdet=BBE_RECIPN_2XF32(det);
        BBE_SVN_2XF32_IP(rdet,pRdetwr,2*BBE_SIMD_WIDTH);
    }
    __Pragma("no_reorder")
    pArd=(const xb_vecN_2xf32 *)X;
    pAwr=(      xb_vecN_2xf32 *)X;
    pRdetrd=(const xb_vecN_2xf32 *)pScr;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        b=BBE_LVN_2XF32_X(pArd,estride);
        c=BBE_LVN_2XF32_X(pArd,sstride);
        d=BBE_LVN_2XF32_X(pArd,sstride+estride);
        BBE_LVN_2XF32_XP(a,pArd,mstride);
        BBE_LVN_2XF32_IP(rdet,pRdetrd,2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32_X (BBE_MULMN_2XF32(b,rdet,3,12),pAwr,estride);
        BBE_SVN_2XF32_X (BBE_MULMN_2XF32(c,rdet,3,12),pAwr,sstride);
        BBE_SVN_2XF32_X (BBE_MULN_2XF32 (a,rdet     ),pAwr,sstride+estride);
        BBE_SVN_2XF32_XP(BBE_MULN_2XF32 (d,rdet     ),pAwr,mstride);
    }
#else
#error 
#endif
}

#endif

size_t  matinv2x2f_inv_getScratchSize ( int L )
#if SPLIT_CYCLE
{
    size_t sz=0;
    (void)L;
    sz=L>0? L*sizeof(float32_t):0;
    return sz;
}
#else
{
    (void)L;
    return 0;
}
#endif

#endif // HAVE_VFPU
