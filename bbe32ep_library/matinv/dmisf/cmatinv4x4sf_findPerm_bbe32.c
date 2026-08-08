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
#define TRUNCATED_SEARCH 1
#define NONINTERLEAVED_FINDUL 0 // 0 - use search UL after interleaving, 1 - before
#if 0
#include <math.h>
#include <complex.h>
// determinant of 2x2
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

#if NONINTERLEAVED_FINDUL
/*------------------------------------------
find position of upper left corner
Look for an element of maximum absolute value.
Input:
X[16*L]  -L matrices
Output:
pos[L]   - index of position of element with 
maximum absolute value
------------------------------------------*/
static void findUL(int16_t *pos,const complex_float *X, int L)
#if 0
{
    int l;
    int        n,maxIx;
    float32_t  w,maxW;
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);

    for ( l=0; l<L; l++ )
    {
        for ( maxW=0.f, maxIx=0, n=0; n<16; n++ )
        {
            complex_float t=X[l+L*n];
            w=fabsf(crealf(t)*crealf(t)+cimagf(t)*cimagf(t));
            if ( maxW < w) { maxW = w, maxIx = n; }
        }
        pos[l]=maxIx;
    }
}
#else
{
    vselN_2 sel2_of_8;  // 0,1,8,9,xxxxx
    xb_vecN_2xc16* restrict pPos;
    valign aP;
    vboolN_2 _true=~BBE_LTRN_2I(0);
    int l;
    xb_vecN_2xc16 maxIx,ix;
    xb_vecN_2xf32  a,w,maxW;
    const xb_vecN_2xf32* restrict pX;
    // form special shuffle mask from 2 existing selects
    {
        xb_vecNx16 t=BBE_SEQNX16();
        t=BBE_SELNX16I(t,t,BBE_SELI_16B_EXTRACT_2_OF_4_OFF_0);
        t=BBE_SELNX16I(t,t,BBE_SELI_16B_EXTRACT_2_OF_4_OFF_0);
        sel2_of_8=BBE_MOVVSELN_2NX16(t,0);
    }

    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    aP=BBE_ZALIGN();
    pPos=(xb_vecN_2xc16* )pos;
    pX=(const xb_vecN_2xf32*)X;
    for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-2)); l++ )
    {
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); maxW=a;maxIx = 0;        ix=1;
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); /*BBE_ADDN_2XC16T(ix,ix,1,_true);*/ ix=BBE_ADDN_2XC16(ix,1);
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); /*BBE_ADDN_2XC16T(ix,ix,1,_true);*/ ix=BBE_ADDN_2XC16(ix,1);
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); /*BBE_ADDN_2XC16T(ix,ix,1,_true);*/ ix=BBE_ADDN_2XC16(ix,1);
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); /*BBE_ADDN_2XC16T(ix,ix,1,_true);*/ ix=BBE_ADDN_2XC16(ix,1);
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); /*BBE_ADDN_2XC16T(ix,ix,1,_true);*/ ix=BBE_ADDN_2XC16(ix,1);
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); /*BBE_ADDN_2XC16T(ix,ix,1,_true);*/ ix=BBE_ADDN_2XC16(ix,1);
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); /*BBE_ADDN_2XC16T(ix,ix,1,_true);*/ ix=BBE_ADDN_2XC16(ix,1);
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); /*BBE_ADDN_2XC16T(ix,ix,1,_true);*/ ix=BBE_ADDN_2XC16(ix,1);
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); /*BBE_ADDN_2XC16T(ix,ix,1,_true);*/ ix=BBE_ADDN_2XC16(ix,1);
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); /*BBE_ADDN_2XC16T(ix,ix,1,_true);*/ ix=BBE_ADDN_2XC16(ix,1);
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); BBE_ADDN_2XC16T(ix,ix,1,_true); /*ix=BBE_ADDN_2XC16(ix,1);*/
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); /*BBE_ADDN_2XC16T(ix,ix,1,_true);*/ ix=BBE_ADDN_2XC16(ix,1);
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); /*BBE_ADDN_2XC16T(ix,ix,1,_true);*/ ix=BBE_ADDN_2XC16(ix,1);
        BBE_LVN_2XF32_XP(w,pX,L*sizeof(complex_float)); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); maxW=BBE_MAXN_2XF32(w,maxW); /*BBE_ADDN_2XC16T(ix,ix,1,_true);*/ ix=BBE_ADDN_2XC16(ix,1);
        BBE_LVN_2XF32_XP(w,pX,-15*L*sizeof(complex_float)+2*BBE_SIMD_WIDTH); a=BBE_MULMN_2XF32(w ,w ,0,0); BBE_MULMASN_2XF32(a,w,w,0,15); w=a;  maxIx=BBE_MOVN_2XC16T(maxIx,ix,BBE_OGTN_2XF32(w,maxW)); 
        maxIx=BBE_SHFLN_2XC16(maxIx,sel2_of_8);
        BBE_SAVN_2XC16_XP(maxIx,aP,pPos,BBE_SIMD_WIDTH/2);
        _true= _true & _true;   /* dummy operator to prevent compiler precompute too many ix outside the loop */
    }
    BBE_SAN_2XC16POS_FP(aP,pPos);
}
#endif
#else /* NONINTERLEAVED_FINDUL-*/
#endif /* NONINTERLEAVED_FINDUL-*/

/*-------------------------------------------------------------------------
    find permutation for L matrices written in stream order
    and return permuted matrices
    Input:
    X[L*16]  input matrices
    L        number of matrices
    Output:
    permIx[L] permutation indices
-------------------------------------------------------------------------*/
void cmatinv4x4sf_findAndPerm(int16_t *permIx,complex_float *X, int L)
{
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);

#if NONINTERLEAVED_FINDUL
    findUL(permIx,X, L);
    cmatinv4x4sf_csb4x4_inplace(X,L);
#else
    cmatinv4x4f_csb4x4_inplace(X,L,e4x4_stream);
    cmatinv4x4f_findUL(permIx,X, L,e4x4_stream);
#endif
    cmatinv4x4f_truncatedSearch(permIx,X,L,e4x4_stream);
    cmatinv4x4f_permute(X,X,permIx,matinv4x4_fwd_perm_tbl_bbe32,L,e4x4_stream);
    cmatinv4x4f_cbs4x4_inplace(X,L,e4x4_stream);
}
#endif
