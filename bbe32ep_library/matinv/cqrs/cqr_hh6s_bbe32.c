/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
/*          Copyright (C) 2009-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
  NatureDSP_Baseband library. QR-based matrix decomposition and inversion for streaming order
    helper routines for Housholder algorithm, len_x==6
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "cqr_common.h"

#if !(HAVE_VSAMATH && HAVE_NSAENX40 && 1)
DISCARD_FUN(void, cqrHouseholder6,( void* pScr,const int16_t *   A,int16_t *    v ,int16_t*    Fi_out,const int M,const int N,const int L))

DISCARD_FUN(void,cqrUpdateR6,(  int16_t* restrict R,
const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
int m,int M,int N, int L))
#else
/*
    find Housholder vector by column and transform that column
    len_x==6
*/
void cqrHouseholder6(
    void* pScr,          /* scratch (4*L bytes)        */
    const int16_t *   A, /* input matrix               */
    int16_t *    v,     /* output (housholder vector) */
    int16_t*    Fi_out,  /* output phase rotator       */
    const int M,         /* number of rows             */
    const int N,         /* number of columns          */
    const int L          /* number of matrices         */
    )
{
    int l;
    xb_vecNx16* restrict _pwrs;   // write pointer to the scratch data
    const xb_vecNx16* restrict _prds;   // read pointer to the scratch data
    const xb_vecNx16* restrict px;
    xb_vecNx16* restrict pFi;
    xb_vecNx16 *restrict pv;
    const xb_vecNx16 *px_v16;
    const int V_Lstep = (SIZE_OF_V(M, N) + SIZE_OF_FI(M));
    xb_vecNx16 x, t, mant, k_norm, norm_x, fi;
    vsaN exp;
    xb_vecNx40 a;

    NASSERT_ALIGN(v, BBE_SIMD_WIDTH * 2);
    NASSERT(M > 1);

    //-----------------------------
    // compute diagonal element
    //-----------------------------
    pFi = (xb_vecNx16*)Fi_out;
    px = (xb_vecNx16*)A;
    for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
    {
        BBE_LVNX16_IP(norm_x, px, 2 * BBE_SIMD_WIDTH);
        a = BBE_MULNX16J(norm_x, norm_x);
        a = BBE_ADDNX40(a, a);
        exp = BBE_NSAENX40(a);
        a = BBE_SLLNX40(a, exp);
        BBE_RSQRTLUNX40_0(a, x, t, a);
        BBE_MULUUSNX16(a, t, x);
        a = BBE_SRAINX40(a, 23);
        mant = BBE_PACKLNX40(a);
        exp = BBE_SUBSR1SAVSN(16 + 4, exp);
        x = BBE_MOVVVS(exp);
        mant = BBE_SHFLNX16I(mant, BBE_SHFLI_DUPLICATE_1_EVEN);
        x = BBE_SHFLNX16I(x, BBE_SHFLI_DUPLICATE_1_EVEN);
        exp = BBE_MOVVSV(x, 0);
        a = BBE_MULUSRNX16(mant, norm_x, exp);
        norm_x = BBE_PACKVNX40(a, exp);
        BBE_SVNX16_XP(norm_x, pFi, V_Lstep * 2 * BBE_SIMD_WIDTH);
    }
    //-----------------------------
    // compute rotation vector
    //-----------------------------
    _pwrs = (xb_vecNx16*)pScr;
    px = (xb_vecNx16*)A;
    for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
    {
        px_v16 = (const xb_vecNx16*)px;
        BBE_LVNX16_XP(x, px_v16, 4 * N*L);
        a = BBE_MULNX16J(x, x);
        BBE_LVNX16_XP(x, px_v16, 4 * N*L);
        BBE_MULANX16J(a, x, x);
        BBE_LVNX16_XP(x, px_v16, 4 * N*L);
        BBE_MULANX16J(a, x, x);
        BBE_LVNX16_XP(x, px_v16, 4 * N*L);
        BBE_MULANX16J(a, x, x);
        BBE_LVNX16_XP(x, px_v16, 4 * N*L);
        BBE_MULANX16J(a, x, x);
        BBE_LVNX16_XP(x, px_v16, 4 * N*L);
        BBE_MULANX16J(a, x, x);
        a = BBE_ADDNX40(a, a);
        exp = BBE_NSAENX40(a);
        a = BBE_SLLNX40(a, exp);
        BBE_RSQRTLUNX40_0(a, x, t, a);
        BBE_MULUUSNX16(a, t, x);
        a = BBE_SRAINX40(a, 23);
        mant = BBE_PACKLNX40(a);
        exp = BBE_SUBSR1SAVSN(16 + 4, exp);
        x = BBE_MOVVVS(exp);
        x = BBE_SELNX16I(x, mant, BBE_SELI_INTERLEAVE_1_EVEN);
        BBE_SVNX16_IP(x, _pwrs, 2 * BBE_SIMD_WIDTH);

        px += (BBE_SIMD_WIDTH * 2) / sizeof(*px);
    }
    __Pragma("no_reorder")
        //-------------------------------------
        // affine projection of original column
        //-------------------------------------
        _prds = (const xb_vecNx16*)pScr;
    pFi = (xb_vecNx16*)Fi_out;
    pv = (xb_vecNx16*)(v);
    px = (xb_vecNx16*)A;
    for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
    {
        BBE_LVNX16_XP(fi, pFi, V_Lstep * 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x, _prds, 2 * BBE_SIMD_WIDTH);
        mant = BBE_SHFLNX16I(x, BBE_SHFLI_DUPLICATE_1_EVEN);
        x = BBE_SHFLNX16I(x, BBE_SHFLI_DUPLICATE_1_ODD);
        exp = BBE_MOVVSV(x, 0);

        px_v16 = (const xb_vecNx16*)px;
        BBE_LVNX16_XP(x, px_v16, 4 * N*L);
        a = BBE_MULUSRNX16(mant, x, exp);
        norm_x = BBE_PACKVNX40(a, exp);
        norm_x = BBE_ADDNX16(norm_x, fi);
        {
            vsaN c_vec;
            a = BBE_MULNX16J(norm_x, fi);
            a = BBE_ADDNX40(a, a);
            c_vec = BBE_NSAENX40(a);
            a = BBE_SLLNX40(a, c_vec);
            c_vec = BBE_SUBSR1SAVSN(28, c_vec);
            BBE_RSQRTLUNX40_0(a, t, x, a);
            BBE_MULUUSNX16(a, x, t);
            k_norm = BBE_PACKVNX40(a, c_vec);
            k_norm = BBE_SHFLNX16I(k_norm, BBE_SHFLI_DUPLICATE_1_EVEN);
        }
        norm_x = BBE_MULNX16PACKQ(k_norm, norm_x);
        BBE_SVNX16_IP(norm_x, pv, BBE_SIMD_WIDTH * 2);
        BBE_LVNX16_XP(x, px_v16, 4 * N*L);
        a = BBE_MULUSRNX16(mant, x, exp);
        x = BBE_PACKVNX40(a, exp);
        x = BBE_MULNX16PACKQ(x, k_norm);
        BBE_SVNX16_IP(x, pv, BBE_SIMD_WIDTH * 2);
        BBE_LVNX16_XP(x, px_v16, 4 * N*L);
        a = BBE_MULUSRNX16(mant, x, exp);
        x = BBE_PACKVNX40(a, exp);
        x = BBE_MULNX16PACKQ(x, k_norm);
        BBE_SVNX16_IP(x, pv, BBE_SIMD_WIDTH * 2);
        BBE_LVNX16_XP(x, px_v16, 4 * N*L);
        a = BBE_MULUSRNX16(mant, x, exp);
        x = BBE_PACKVNX40(a, exp);
        x = BBE_MULNX16PACKQ(x, k_norm);
        BBE_SVNX16_IP(x, pv, BBE_SIMD_WIDTH * 2);
        BBE_LVNX16_XP(x, px_v16, 4 * N*L);
        a = BBE_MULUSRNX16(mant, x, exp);
        x = BBE_PACKVNX40(a, exp);
        x = BBE_MULNX16PACKQ(x, k_norm);
        BBE_SVNX16_IP(x, pv, BBE_SIMD_WIDTH * 2);
        BBE_LVNX16_XP(x, px_v16, 4 * N*L);
        a = BBE_MULUSRNX16(mant, x, exp);
        x = BBE_PACKVNX40(a, exp);
        x = BBE_MULNX16PACKQ(x, k_norm);
        BBE_SVNX16_IP(x, pv, BBE_SIMD_WIDTH * 2);
        pv += (V_Lstep - 6)*(BBE_SIMD_WIDTH * 2) / sizeof(*pv);
        px += (BBE_SIMD_WIDTH * 2) / sizeof(*px);
    }
}

/*
    update columns of matrix R by Housholder vectors
    len_x==6
*/
void cqrUpdateR6(int16_t* restrict R,
    const int16_t * restrict V,     /*i  pointer to Householder vector, Q14 */
    int m, int M, int N, int L)
{
    const xb_vecNx16 * restrict pr;
    xb_vecNx16 * restrict pw;
    const xb_vecNx16 * restrict pv;
    vsaN q14 = BBE_MOVVSA32(14);
    xb_vecNx40 acc;
    xb_vecNx16 r0, r1, r2, r3, r4, r5, v0, v1, v2, v3, v4, v5, vr, _0x4000 = 0x4000;
    int l, j, sizeV = (SIZE_OF_V(M, N) + SIZE_OF_FI(M));
    int strideV;
    pr = (const xb_vecNx16*)(R);
    pw = (xb_vecNx16*)(R);
    pv = (const xb_vecNx16*)(V);
    for (j = l = 0; l < (N - m)*((L >> (LOG2_BBE_SIMD_WIDTH - 1))); l++)
    {
        j = BBE_ADDMOD16U(j, ((L >> (LOG2_BBE_SIMD_WIDTH - 1)) << 16) | 1);
        strideV = 2 * BBE_SIMD_WIDTH*(sizeV - 5);
        XT_MOVEQZ(strideV, 2 * BBE_SIMD_WIDTH*(sizeV - 5) - 4 * L*sizeV, j);
        BBE_LVNX16_XP(r0, pr, 4 * N*L);
        BBE_LVNX16_XP(r1, pr, 4 * N*L);
        BBE_LVNX16_XP(r2, pr, 4 * N*L);
        BBE_LVNX16_XP(r3, pr, 4 * N*L);
        BBE_LVNX16_XP(r4, pr, 4 * N*L);
        BBE_LVNX16_XP(r5, pr, -5 * 4 * N*L);
        BBE_LVNX16_IP(v0, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v1, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v2, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v3, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v4, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(v5, pv, strideV);
        acc = BBE_MULRNX16J(r0, v0, q14);
        BBE_MULANX16J(acc, r1, v1);
        BBE_MULANX16J(acc, r2, v2);
        BBE_MULANX16J(acc, r3, v3);
        BBE_MULANX16J(acc, r4, v4);
        BBE_MULANX16J(acc, r5, v5);
        vr = BBE_PACKVNX40(acc, q14);
        BBE_LVNX16_XP(r0, pr, 4 * N*L);
        BBE_LVNX16_XP(r1, pr, 4 * N*L);
        BBE_LVNX16_XP(r2, pr, 4 * N*L);
        BBE_LVNX16_XP(r3, pr, 4 * N*L);
        BBE_LVNX16_XP(r4, pr, 4 * N*L);
        BBE_LVNX16_XP(r5, pr, 2 * BBE_SIMD_WIDTH - 5 * 4 * N*L);
        acc = BBE_MULRNX16(r0, _0x4000, q14); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r1, _0x4000, q14); BBE_MULSNX16C(acc, v1, vr); r1 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r2, _0x4000, q14); BBE_MULSNX16C(acc, v2, vr); r2 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r3, _0x4000, q14); BBE_MULSNX16C(acc, v3, vr); r3 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r4, _0x4000, q14); BBE_MULSNX16C(acc, v4, vr); r4 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r5, _0x4000, q14); BBE_MULSNX16C(acc, v5, vr); r5 = BBE_PACKVNX40(acc, q14);
        BBE_SVNX16_XP(r0, pw, 4 * N*L);
        BBE_SVNX16_XP(r1, pw, 4 * N*L);
        BBE_SVNX16_XP(r2, pw, 4 * N*L);
        BBE_SVNX16_XP(r3, pw, 4 * N*L);
        BBE_SVNX16_XP(r4, pw, 4 * N*L);
        BBE_SVNX16_XP(r5, pw, 2 * BBE_SIMD_WIDTH - 5 * 4 * N*L);
    }
}
#endif
