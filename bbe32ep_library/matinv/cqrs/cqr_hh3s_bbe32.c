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
    helper routines for Housholder algorithm, len_x==3
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
DISCARD_FUN(void,cqrHouseholder3,( 
    const int16_t *   A,
    int16_t *    v ,         
    int16_t*    Fi_out,                                       
    const int M,        
    const int N,        
    const int L         
    ))

    DISCARD_FUN(void,cqrUpdateR3,(  int16_t* restrict R,
    const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
    int m,int M,int N, int L))
#else

/*
    find Housholder vector by column and transform that column
    len_x==3
*/
void cqrHouseholder3(
    const int16_t *   A,/* input matrix               */
    int16_t *    v,     /* output (housholder vector) */
    int16_t*    Fi_out, /* output phase rotator       */
    const int M,        /* number of rows             */
    const int N,        /* number of columns          */
    const int L         /* number of matrices         */
    )
{
    int l;
    const xb_vecNx16* restrict px;
    xb_vecNx16* restrict pFi;
    xb_vecNx16 *restrict pv;
    const int V_Lstep = (SIZE_OF_V(M, N) + SIZE_OF_FI(M));
    xb_vecNx16 x, t, mant, k_norm, norm_x, fi, x0, x1, x2;
    vsaN exp;
    xb_vecNx40 a;

    NASSERT_ALIGN(v, BBE_SIMD_WIDTH * 2);
    NASSERT(M > 1);

    pFi = (xb_vecNx16*)Fi_out;
    px = (xb_vecNx16*)A;
    pv = (xb_vecNx16*)(v);
    for (l = 0; l < ((L >> (LOG2_BBE_SIMD_WIDTH - 1))); l++)
    {
        vsaN c_vec;
        //-----------------------------
        // compute diagonal element
        //-----------------------------
        BBE_LVNX16_IP(x0, px, 0);
        a = BBE_MULNX16J(x0, x0);
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
        a = BBE_MULUSRNX16(mant, x0, exp);
        fi = BBE_PACKVNX40(a, exp);
        BBE_SVNX16_XP(fi, pFi, V_Lstep * 2 * BBE_SIMD_WIDTH);
        //-----------------------------
        // compute rotation vector
        //-----------------------------
        BBE_LVNX16_XP(x0, px, 4 * N*L); a = BBE_MULNX16J(x0, x0);
        BBE_LVNX16_XP(x1, px, 4 * N*L); BBE_MULANX16J(a, x1, x1);
        BBE_LVNX16_XP(x2, px, -2 * 4 * N*L); BBE_MULANX16J(a, x2, x2);

        a = BBE_ADDNX40(a, a);
        exp = BBE_NSAENX40(a);
        a = BBE_SLLNX40(a, exp);
        BBE_RSQRTLUNX40_0(a, x, t, a);
        BBE_MULUUSNX16(a, t, x);
        a = BBE_SRAINX40(a, 23);
        mant = BBE_PACKLNX40(a);
        exp = BBE_SUBSR1SAVSN(16 + 4, exp);
        x = BBE_MOVVVS(exp);
        x = BBE_SHFLNX16I(x, BBE_SHFLI_DUPLICATE_1_EVEN);
        exp = BBE_MOVVSV(x, 0);
        mant = BBE_SHFLNX16I(mant, BBE_SHFLI_DUPLICATE_1_EVEN);
        //-------------------------------------
        // affine projection of original column
        //-------------------------------------
        BBE_LVNX16_XP(x0, px, 4 * N*L);
        a = BBE_MULUSRNX16(mant, x0, exp);
        norm_x = BBE_PACKVNX40(a, exp);
        norm_x = BBE_ADDNX16(norm_x, fi);
        a = BBE_MULNX16J(norm_x, fi);
        a = BBE_ADDNX40(a, a);
        c_vec = BBE_NSAENX40(a);
        a = BBE_SLLNX40(a, c_vec);
        c_vec = BBE_SUBSR1SAVSN(28, c_vec);
        BBE_RSQRTLUNX40_0(a, t, x, a);
        BBE_MULUUSNX16(a, x, t);
        k_norm = BBE_PACKVNX40(a, c_vec);
        k_norm = BBE_SHFLNX16I(k_norm, BBE_SHFLI_DUPLICATE_1_EVEN);

        norm_x = BBE_MULNX16PACKQ(k_norm, norm_x);
        BBE_SVNX16_IP(norm_x, pv, BBE_SIMD_WIDTH * 2);
        BBE_LVNX16_XP(x1, px, 4 * N*L); a = BBE_MULUSRNX16(mant, x1, exp); x = BBE_PACKVNX40(a, exp); x = BBE_MULNX16PACKQ(x, k_norm); BBE_SVNX16_IP(x, pv, BBE_SIMD_WIDTH * 2);
        BBE_LVNX16_XP(x2, px, -2 * 4 * N*L + BBE_SIMD_WIDTH * 2); a = BBE_MULUSRNX16(mant, x2, exp); x = BBE_PACKVNX40(a, exp); x = BBE_MULNX16PACKQ(x, k_norm); BBE_SVNX16_XP(x, pv, (V_Lstep - 3 + 1)*BBE_SIMD_WIDTH * 2);
    }
}

/*
    update columns of matrix R by Housholder vectors
    len_x==3
*/
void cqrUpdateR3(int16_t* restrict R,
    const int16_t * restrict V,         /*i  pointer to Householder vector, Q14 */
    int m, int M, int N, int L)
{
    int strideV;
    const xb_vecNx16 * restrict pr;
    xb_vecNx16 * restrict pw;
    const xb_vecNx16 * restrict pv;
    vsaN q14 = BBE_MOVVSA32(14);
    xb_vecNx40 acc;
    xb_vecNx16 r0, r1, r2, v0, v1, v2, vr, _0x4000 = 0x4000;
    int l, j, sizeV = (SIZE_OF_V(M, N) + SIZE_OF_FI(M));
    pv = (const xb_vecNx16*)(V);
    pr = (const xb_vecNx16*)(R);
    pw = (xb_vecNx16*)(R);
    for (j = l = 0; l < (N - m)*((L >> (LOG2_BBE_SIMD_WIDTH - 1))); l++)
    {
        j = BBE_ADDMOD16U(j, ((L >> (LOG2_BBE_SIMD_WIDTH - 1)) << 16) | 1);
        strideV = 2 * BBE_SIMD_WIDTH*(sizeV - 2);
        XT_MOVEQZ(strideV, 2 * BBE_SIMD_WIDTH*(sizeV - 2) - 4 * L*sizeV, j);
        BBE_LVNX16_XP(r0, pr, 4 * N*L);
        BBE_LVNX16_XP(r1, pr, 4 * N*L);
        BBE_LVNX16_XP(r2, pr, 2 * BBE_SIMD_WIDTH - 2 * 4 * N*L);
        BBE_LVNX16_IP(v0, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v1, pv, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(v2, pv, strideV);
        acc = BBE_MULRNX16J(r0, v0, q14);
        BBE_MULANX16J(acc, r1, v1);
        BBE_MULANX16J(acc, r2, v2);
        vr = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r0, _0x4000, q14); BBE_MULSNX16C(acc, v0, vr); r0 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r1, _0x4000, q14); BBE_MULSNX16C(acc, v1, vr); r1 = BBE_PACKVNX40(acc, q14);
        acc = BBE_MULRNX16(r2, _0x4000, q14); BBE_MULSNX16C(acc, v2, vr); r2 = BBE_PACKVNX40(acc, q14);
        BBE_SVNX16_XP(r0, pw, 4 * N*L);
        BBE_SVNX16_XP(r1, pw, 4 * N*L);
        BBE_SVNX16_XP(r2, pw, 2 * BBE_SIMD_WIDTH - 2 * 4 * N*L);
    }
}
#endif
