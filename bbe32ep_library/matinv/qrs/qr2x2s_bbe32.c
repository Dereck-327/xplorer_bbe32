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
    These functions apply QR decomposition procedure to the sequence of complex
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
#include "qr_common.h"

#if !(HAVE_VSAMATH && HAVE_NSAENX40 && 1)
DISCARD_FUN(void, qr2x2s, (void* pScr,
    int16_t * restrict Q, 
    int16_t * restrict R, 
    int L))
#else

#define CURRENT_M 2

#undef SIZE_OF_FI
#undef SIZE_OF_V

#define SIZE_OF_V(__M) ((__M-1)*(__M+2)/2)
#define SIZE_OF_FI(__M) __M


inline_ void chouseholder_2x16(xb_vecNx16 *pr0,
    xb_vecNx16 *pr1,
    xb_vecNx16 *pv0,
    xb_vecNx16 *pv1,
    xb_vecNx16 *pfi)
{
    xb_vecNx16 tmp0, tmp1;
    vsaN c_vec, v_exp;
    xb_vecNx40 acc0;
    xb_vecNx16 k_norm, mant;
    xb_vecNx16 norm_x, Fi;
    xb_vecNx16  x1, x2;

    x1 = *pr0;

    acc0 = BBE_MULNX16(x1, x1);
    acc0 = BBE_ADDNX40(acc0, acc0);
    c_vec = BBE_NSAENX40(acc0);
    acc0 = BBE_SLLNX40(acc0, c_vec);
    BBE_RSQRTLUNX40_0(acc0, tmp0, tmp1, acc0);
    BBE_RSQRTLUNX40_1(acc0, tmp0, tmp1, acc0);
    BBE_MULUUSNX16(acc0, tmp1, tmp0);
    acc0 = BBE_SRAINX40(acc0, 24);
    tmp0 = BBE_PACKLNX40(acc0);

    c_vec = BBE_SUBSR1SAVSN(18 + 1, c_vec);
    acc0 = BBE_MULUSRNX16(tmp0, x1, c_vec);
    norm_x = BBE_PACKVNX40(acc0, c_vec);

    Fi = norm_x;
    *pfi = norm_x;
    /* Build Householder's vector   v = x/sqrt(x'*x) + e1*x(1)/abs(x(1) */
    x2 = *pr1;

    acc0 = BBE_MULNX16(x1, x1);//acc0 = BBE_MAGINX16C(x1, x1);
    BBE_MULANX16(acc0, x2, x2);//BBE_MAGIANX16C(acc0,x2, x2);
    acc0 = BBE_ADDNX40(acc0, acc0);
    c_vec = BBE_NSAENX40(acc0);
    acc0 = BBE_SLLNX40(acc0, c_vec);
    BBE_RSQRTLUNX40_0(acc0, tmp0, tmp1, acc0);
    BBE_RSQRTLUNX40_1(acc0, tmp0, tmp1, acc0);
    BBE_MULUUSNX16(acc0, tmp1, tmp0);
    acc0 = BBE_SRAINX40(acc0, 24);
    tmp0 = BBE_PACKLNX40(acc0);
    mant = tmp0;//mant= BBE_SHFLNX16I(tmp0, BBE_SHFLI_DUPLICATE_1_EVEN);
    v_exp = BBE_SUBSR1SAVSN(18 + 1, c_vec);
    acc0 = BBE_MULUSRNX16(mant, x1, v_exp);
    norm_x = BBE_PACKVNX40(acc0, v_exp);

    x1 = *pr0;
    norm_x = BBE_ADDNX16(norm_x, Fi);
    acc0 = BBE_MULNX16(norm_x, Fi);
    acc0 = BBE_ADDNX40(acc0, acc0);
    c_vec = BBE_NSAENX40(acc0);
    acc0 = BBE_SLLNX40(acc0, c_vec);
    BBE_RSQRTLUNX40_0(acc0, tmp0, tmp1, acc0);
    BBE_RSQRTLUNX40_1(acc0, tmp0, tmp1, acc0);
    BBE_MULUUSNX16(acc0, tmp1, tmp0);
    c_vec = BBE_SUBSR1SAVSN(28, c_vec);
    k_norm = BBE_PACKVNX40(acc0, c_vec);
    //k_norm = BBE_SHFLNX16I(k_norm, BBE_SHFLI_DUPLICATE_1_EVEN);

    acc0 = BBE_MULNX16(norm_x, k_norm);
    norm_x = BBE_PACKQNX40(acc0);

    *pv0 = norm_x;
    x1 = *pr1;
    acc0 = BBE_MULUSRNX16(mant, x1, v_exp);
    norm_x = BBE_PACKVNX40(acc0, v_exp);
    norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
    *pv1 = norm_x;
}

/*-------------------------------------------------------------------------
These functions apply QR decomposition procedure to the sequence of complex 
matrices written in a streaming order. The transformation is done in-place 
so the result replaces the original input.
Rotation matrix Q is calculated in Q15 fixed point representation. Fixed 
point representation of upper-diagonal matrix R is the same as of input. 

Functions return nonzero if overflow is detected 

NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Input:
R[M*N][L][C]   input matrices
C              1 for real, 2 for complex data
Output:
Q[M*M][L][C]   output rotation matrices (L matrices of size MxM)
R[M*N][L][C]   output upper triangular matrices (L matrices of size MxN)

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Matrix sizes M,N,L must be greater than 1
4. Scratch memory must be aligned on 32-byte boundary. Its size  (in bytes)
   is defined by cqr2x2s_getScratchSize, qr2x2s_getScratchSize
---------------------------------------------------------------------------*/

void  qr2x2s (void* pScr, int16_t * restrict Q, int16_t * restrict R, int L)
{
    int l;
    const int q = 14;
    xb_vecNx16  ar[CURRENT_M*CURRENT_M];
    xb_vecNx16  av[CURRENT_M];

    const xb_vecNx16 * restrict pRRr;
    xb_vecNx16 * restrict pRRw;

#ifdef COMPILER_XTENSA

#pragma ymemory( pRRr )
#pragma ymemory( pRRw )

#endif

    xb_vecNx16 * restrict pQ00 = (xb_vecNx16 *)(Q + 0 * L);
    xb_vecNx16 * restrict pQ01 = (xb_vecNx16 *)(Q + 1 * L);
    xb_vecNx16 * restrict pQ10 = (xb_vecNx16 *)(Q + 2 * L);
    xb_vecNx16 * restrict pQ11 = (xb_vecNx16 *)(Q + 3 * L);
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Q, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L%BBE_SIMD_WIDTH == 0);

    pRRr = (const xb_vecNx16 *)R;
    pRRw = (xb_vecNx16 *)R;

    vsaN  sh14 = BBE_MOVVSA32(q);
    vsaN  sh13 = BBE_MOVVSA32(q - 1);
    xb_vecNx16 one_Q14 = BBE_MOVVA16(1 << 14);
    xb_vecNx16 _0x3fff = BBE_MOVVA16(0x3fff);
    xb_vecNx16 _0x0000 = BBE_MOVVA16(0x0000);

    for (l = 0; l < L; l += BBE_SIMD_WIDTH)
    {
        xb_vecNx40  acc_Nx40, x_vec, a_vec;
        xb_vecNx16 r_Nx16, vr_Nx16_0, vr_Nx16_1;
        xb_vecNx16 manN, k_norm_Nx16, v1_Nx16;
        xb_vecNx16 b_vec, xn_vec, cn_vec;
        vsaN sh_vsaN, c_vec;

        BBE_LVNX16_XP(ar[0], pRRr, L*sizeof(int16_t)); // 0
        BBE_LVNX16_XP(ar[1], pRRr, L*sizeof(int16_t)); // 1L
        BBE_LVNX16_XP(ar[2], pRRr, L*sizeof(int16_t)); // 2L
        BBE_LVNX16_XP(ar[3], pRRr, 2 * BBE_SIMD_WIDTH - 3 * L*sizeof(int16_t)); // 3L

        x_vec = BBE_MULNX16(ar[0], ar[0]);
        BBE_MULANX16(x_vec, ar[2], ar[2]);

        {
            x_vec = BBE_ADDNX40(x_vec, x_vec);
            c_vec = BBE_NSAENX40(x_vec);
            x_vec = BBE_SLLNX40(x_vec, c_vec);
            BBE_RSQRTLUNX40_0(a_vec, b_vec, cn_vec, x_vec);
            BBE_RSQRTLUNX40_1(a_vec, b_vec, cn_vec, x_vec);
            BBE_MULUUSNX16(a_vec, cn_vec, b_vec);
            a_vec = BBE_SRAINX40(a_vec, 23);
            manN = BBE_PACKLNX40(a_vec);
            sh_vsaN = BBE_SUBSR1SAVSN(20, c_vec);
        }

        acc_Nx40 = BBE_MULUSRNX16(manN, ar[0], sh_vsaN);
        k_norm_Nx16 = BBE_PACKVNX40(acc_Nx40, sh_vsaN);

        r_Nx16 = BBE_SRAINX16(ar[0], 15);
        r_Nx16 = BBE_SLLINX16(r_Nx16, 15);
        r_Nx16 = BBE_ADDNX16(r_Nx16, one_Q14);
        v1_Nx16 = BBE_ADDSNX16(r_Nx16, k_norm_Nx16);
        x_vec = BBE_MULNX16(v1_Nx16, r_Nx16);

        {
            x_vec = BBE_ADDNX40(x_vec, x_vec);
            c_vec = BBE_NSAENX40(x_vec);
            x_vec = BBE_SLLNX40(x_vec, c_vec);
            c_vec = BBE_SUBSR1SAVSN(28, c_vec);
            BBE_RSQRTLUNX40_0(a_vec, b_vec, cn_vec, x_vec);
            BBE_RSQRTLUNX40_1(a_vec, b_vec, cn_vec, x_vec);
            BBE_MULUUSNX16(a_vec, cn_vec, b_vec);
            k_norm_Nx16 = BBE_PACKVNX40(a_vec, c_vec);
        }

        av[0] = BBE_MULNX16PACKQ(k_norm_Nx16, v1_Nx16);

        acc_Nx40 = BBE_MULUSNX16(manN, ar[2]);
        r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh_vsaN);
        av[1] = BBE_MULNX16PACKQ(r_Nx16, k_norm_Nx16);
        {
            acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14);
            BBE_MULANX16(acc_Nx40, ar[2], av[1]);
            vr_Nx16_0 = BBE_PACKVNX40(acc_Nx40, sh14);

            acc_Nx40 = BBE_MULRNX16(ar[1], av[0], sh14);
            BBE_MULANX16(acc_Nx40, ar[3], av[1]);
            vr_Nx16_1 = BBE_PACKVNX40(acc_Nx40, sh14);

            acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14);
            BBE_MULSNX16(acc_Nx40, vr_Nx16_0, av[0]);
            r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
            BBE_SVNX16_XP(r_Nx16, pRRw, L*sizeof(int16_t));

            acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14);
            BBE_MULSNX16(acc_Nx40, vr_Nx16_1, av[0]);
            r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
            BBE_SVNX16_XP(r_Nx16, pRRw, L*sizeof(int16_t));

            acc_Nx40 = BBE_MULRNX16(ar[2], one_Q14, sh14);
            BBE_MULSNX16(acc_Nx40, vr_Nx16_0, av[1]);
            r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
            BBE_SVNX16_XP(r_Nx16, pRRw, L*sizeof(int16_t));

            acc_Nx40 = BBE_MULRNX16(ar[3], one_Q14, sh14);
            BBE_MULSNX16(acc_Nx40, vr_Nx16_1, av[1]);
            r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
            BBE_SVNX16_XP(r_Nx16, pRRw, BBE_SIMD_WIDTH * 2 - 3 * L*sizeof(int16_t));
        }
        {
            acc_Nx40 = BBE_MULRNX16(_0x3fff, one_Q14, sh14);
            BBE_MULSNX16(acc_Nx40, av[0], av[0]);
            xn_vec = BBE_PACKVNX40(acc_Nx40, sh13);
            BBE_SVNX16_IP(xn_vec, pQ00, BBE_SIMD_WIDTH * 2); //0L

            acc_Nx40 = BBE_MULRNX16(_0x3fff, one_Q14, sh14);
            BBE_MULSNX16(acc_Nx40, av[1], av[1]);
            xn_vec = BBE_PACKVNX40(acc_Nx40, sh13);
            BBE_SVNX16_IP(xn_vec, pQ11, BBE_SIMD_WIDTH * 2); //3L

            acc_Nx40 = BBE_MULRNX16(_0x0000, one_Q14, sh14);
            BBE_MULSNX16(acc_Nx40, av[1], av[0]);
            xn_vec = BBE_PACKVNX40(acc_Nx40, sh13);
            BBE_SVNX16_IP(xn_vec, pQ01, BBE_SIMD_WIDTH * 2); //1L

            acc_Nx40 = BBE_MULRNX16(_0x0000, one_Q14, sh14);
            BBE_MULSNX16(acc_Nx40, av[0], av[1]);
            xn_vec = BBE_PACKVNX40(acc_Nx40, sh13);
            BBE_SVNX16_IP(xn_vec, pQ10, BBE_SIMD_WIDTH * 2); //2L
        }

    } //for( l=0; l<L; l+=BBE_SIMD_WIDTH )
} /* qr2x2s() */
#endif

size_t  qr2x2s_getScratchSize (int M, int N,int L)
{
    (void)M;
    (void)N;
    (void)L;
    return 0;
} /* qr2x2s_getScratchSize() */
