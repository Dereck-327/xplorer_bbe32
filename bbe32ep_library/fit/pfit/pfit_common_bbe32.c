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
  NatureDSP_Baseband library. Fitting and Interpolation Routines
    Polynomial Fitting and Interpolation for Real Data
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fit.h"
#include "pfit_common.h"
#include <limits.h>

/*-------------------------------------------------------------------------
Polynomial Fitting and Interpolation

Description: the pfit functions fit (in least squares sense) a degree N 
polynomial to input data sampled at a points grid of length M, and use
that polynomial to interpolate data at arbitrary query points. Namely,
the pfitN_grid() functions compute the Vandermonde matrix for the sample
points grid and perform the Cholesky decomposition of that matrix, the
pfitN_process() functions calculate the least squares solution for the
polynomial coefficients. Finally, the pfitN_eval() functions evaluate
the polynomial at query points.

Please refer to the NatureDSP Baseband Library Reference for full details
on these functions.

Representation:
pfit_gridN,      16-bit fixed-point data. Parameter specifications denote
pfit_processN,   fixed-point format for various data items
pfit_evalN       
pfitf_gridN,     IEEE-754 Std single precision floating-point data
pfitf_processN,  
pfitf_evalN     

Note:
Number of fractional bits specidied for various input/output arguments below apply
for the fixed-point variant

Parameters:
Input:
N                     Degree of polynomial, 1..6
M                     Number of sample points
P                     Number of query points
maxIter               Number of least squares solution enhancement iterations. Right 
                      choice depends on required accuracy, the ad-hoc value is (N+1)/2
x[M]                  Sample points grid, Q15 or floating point
y[M]                  Sampled data values, Q15 or floating point
xi[P]                 Query points, Q15 or floating point
M'=(M+7)&(~7), N'=8   for floating point API
M'=(M+15)&(~15),N'=16 for fixed-point API

Intermediate:
V[M'*8]               Vandermonde matrix, Q15 or floating point
R[N'*8]               Upper triangular Cholesky factor of matrix V, Q11 or floating point
Output:
yi[P]                 Data values interpolated at query points, Q15 or floating point
p[N+1]                Polynomial coefficients, Q8.23 or floating point
Temporary:
pScr                  Scratch memory area. To determine the scratch area size required by
                      a function pfitN_<fun>, use the respective helper function 
                      pfit_<fun>_getScratchSize(M,N)

Restrictions:
M>N                   The number of sample points must exceed the degree of polynomial
x,y,xi,yi,V,R,p,pScr  Must not overlap
V,R,pScr              Aligned on 32-byte boundary
---------------------------------------------------------------------------*/

/* Return the scratch area size, in bytes. */
size_t pfit_process_getScratchSize ( int M, int N )
{
    int _M = M + (-M & (BBE_SIMD_WIDTH - 1));

    return (3 *  8 * sizeof(int32_t) + /* yy,xx,pp */
            8 * 16 * sizeof(int16_t) + /* Rt       */
                16 * sizeof(uint16_t)+ /* D        */
                _M * sizeof(int16_t)); /* b        */
} /* pfit_process_getScratchSize() */

/* Return the scratch area size, in bytes. */
size_t pfit_grid_getScratchSize ( int M, int N )
{
    return 0;
} /* pfit_grid_getScratchSize() */


#if HAVE_VSAMATH && HAVE_NSAENX40 && HAVE_RECIP && HAVE_RSQRT

/*-------------------------------------------------------------------------
These functions apply Cholesky decomposition procedure to the sequence of
real matrices written in a streaming order.
R=pfit_chol(A'*A+sigma2*I)
Fixed point representation of upper-diagonal matrix R is the same as of input.

Functions return nonzero if overflow is detected

NOTE:
Data layout for matrices is selected as for other matrices written in a
streaming order.

Input:
A[M*N]  input real matrix, Q15
N       Matrix dimension (number of columns for MxN)
M       Matrix dimension (number of rows for MxN)
sigma2  noise estimate, Q31

Output:
R[N*N]  output real upper-triangle matrices , Q11

Scratch:
t[M+N]

Reference Matlab code:
function R=chol2(A,sigma2)
sz=size(A); M=sz(1); N=sz(2);
R=zeros(N,N);
for m=1:N
Rm=R(:,m);  % take m-th column of original and decomposing matrix
Am=A(:,m);
Amm=Am'*Am+sigma2;
Rmm=Rm'*Rm;
Rmm=sqrt(real(Amm-Rmm));
x(1,1:m)=[zeros(1,m-1) Rmm];
for k=m+1:N
Akm=A(:,k)'*Am;
Rkm=R(:,k)'*Rm;
x(1,k)=(Akm-Rkm)/Rmm;
end
R(m,:)=conj(x);
end
---------------------------------------------------------------------------*/
void pfit_chol(int16_t* tmp,
    int16_t * restrict R,
    const int16_t * restrict A,
    int32_t sigma2,
    int M, int N)
{
    xb_vecNx40 wA, wB;
    xb_int40 tmp40;
    const xb_vecNx16* restrict pAm;
    const xb_vecNx16* restrict pAk;
    xb_vecNx16* pRm;
    xb_vecNx16* pRk;
    xb_vecNx16 t, rk, rm;
    vsaN vq, sh;
    valign vR;
    int m;
    int Dmant, Dexp; // mantissa and exponent of reciprocal of main diagonal
    int n, k, q;
    const int Astride = ((M + (BBE_SIMD_WIDTH - 1)) & ~(BBE_SIMD_WIDTH - 1));
    // compute q 
    wA = BBE_MOVWA32(M << 1);
    vq = BBE_NSAENX40(wA);
    vq = BBE_SUBSAVSN(38, vq);
    t = BBE_MOVVVS(vq); q = BBE_MOVAV16(t);

    NASSERT(q >= 2 && q <= 10); // M in range 2...256
    NASSERT_ALIGN32(A);
    NASSERT_ALIGN32(R);
    NASSERT_ALIGN32(tmp);
    NASSERT(M >= N);
    NASSERT(M>0 && N>0);

    // clean R
    t = BBE_ZERONX16();
    BBE_SVNX16_I(t, (xb_vecNx16*)R, 0 * BBE_SIMD_WIDTH);
    BBE_SVNX16_I(t, (xb_vecNx16*)R, 2 * BBE_SIMD_WIDTH);
    BBE_SVNX16_I(t, (xb_vecNx16*)R, 4 * BBE_SIMD_WIDTH);
    BBE_SVNX16_I(t, (xb_vecNx16*)R, 6 * BBE_SIMD_WIDTH);
    BBE_SVNX16_I(t, (xb_vecNx16*)R, 8 * BBE_SIMD_WIDTH);
    BBE_SVNX16_I(t, (xb_vecNx16*)R, 10 * BBE_SIMD_WIDTH);
    BBE_SVNX16_I(t, (xb_vecNx16*)R, 12 * BBE_SIMD_WIDTH);
    BBE_SVNX16_I(t, (xb_vecNx16*)R, 14 * BBE_SIMD_WIDTH);
    pRm = (xb_vecNx16*)(R);
    pAm = (const xb_vecNx16*)(A);
#ifdef COMPILER_XTENSA
#pragma loop_count min=1
#endif
    for (m = 0; m<N; m++)
    {
        pRk = pRm; vR = BBE_LA_PP(pRm); BBE_LAVNX16_XP(rm, vR, pRm, RSTRIDE << 1);
        /* take colunms of A and R and calculate diagonal elements */
        wA = BBE_ZERONX40();
#ifdef COMPILER_XTENSA
#pragma loop_count min=1
#endif
        for (k = 0; k < ((M + (BBE_SIMD_WIDTH - 1)) >> LOG2_BBE_SIMD_WIDTH); k++)
        {
            BBE_LVNX16_IP(t, pAm, 2 * BBE_SIMD_WIDTH);
            BBE_MULANX16(wA, t, t);
        }
        pAm = (const xb_vecNx16*)XT_ADDX2(-Astride, (uintptr_t)pAm);   // move pAm pointer back

        wB = BBE_MULNX16(rm, rm);
        wB = BBE_SLSINX40(wB, 8);
        wA = BBE_SUBNX40(wA, wB);
        tmp40 = BBE_RADDNX40(wA);
        wA = BBE_MOVNX40_FROM40(tmp40);
        wA = BBE_SRAINX40(wA, 1);
        /* first calculate 1/sqrt() */
        {
            xb_vecNx40 a_vec;
            xb_vecNx16   b_vec, xn_vec, cn_vec;
            vsaN c_vec;
            c_vec = BBE_NSAENX40(wA);
            wA = BBE_SLLNX40(wA, c_vec);
            c_vec = BBE_SUBSR1SAVSN(23 - q, c_vec);
            BBE_RSQRTLUNX40_0(a_vec, b_vec, cn_vec, wA);
            BBE_MULUUSNX16(a_vec, cn_vec, b_vec);
            a_vec = BBE_SRAINX40(a_vec, 23);
            xn_vec = BBE_PACKLNX40(a_vec);
            Dmant = BBE_MOVAV16(xn_vec);
            xn_vec = BBE_MOVVVS(c_vec);
            Dexp = BBE_MOVAV16(xn_vec);
        }
        // compute elements in row from m to N
        vR = BBE_LA_PP(pRk);
        pAk = pAm;
        sh = BBE_MOVVSA32(Dexp);
        for (k = m; k<N; k++)
        {
            xb_vecNx16 ak, am, hi, lo;
            wA = BBE_ZERONX40();
#ifdef COMPILER_XTENSA
#pragma loop_count min=1
#endif
            for (n = 0; n<((M + (BBE_SIMD_WIDTH - 1)) >> LOG2_BBE_SIMD_WIDTH); n++)
            {
                BBE_LVNX16_XP(am, pAm, 2 * BBE_SIMD_WIDTH);
                BBE_LVNX16_XP(ak, pAk, 2 * BBE_SIMD_WIDTH);
                BBE_MULANX16(wA, ak, am);
            }
            pAm = (const xb_vecNx16*)XT_ADDX2(-Astride, (uintptr_t)pAm);

            BBE_LAVNX16_XP(rk, vR, pRk, RSTRIDE << 1);
            wB = BBE_MULNX16(rm, rk);
            wB = BBE_SLSINX40(wB, 8);
            wA = BBE_SUBNX40(wA, wB);
            tmp40 = BBE_RADDNX40(wA);
            wA = BBE_MOVNX40_FROM40(tmp40);
            wA = BBE_SRANX40(wA, vq);

            // normalize to reciprocal of diagonal
            TAKEHILO3(wA, hi, lo);
            t = BBE_MOVVA16(Dmant);
            wA = BBE_MULUUNX16(lo, t);
            wA = BBE_SRAINX40(wA, 16);
            BBE_MULUSANX16(wA, t, hi);
            t = BBE_PACKVNX40(wA, sh);
            BBE_SSNX16_X(t, (int16_t*)pRk, (m - RSTRIDE) << 1);
        }
        // next row
        pAm = (const xb_vecNx16*)XT_ADDX2(Astride, (uintptr_t)pAm);
    }
}


void transposeR(int16_t * Rt, const int16_t * R)
{
    xb_vecNx16 r0, r1, r2, r3, r4, r5, r6, r7;

    r0 = BBE_LVNX16_I((xb_vecNx16*)R, 0 * BBE_SIMD_WIDTH * 2);
    r1 = BBE_LVNX16_I((xb_vecNx16*)R, 1 * BBE_SIMD_WIDTH * 2);
    r2 = BBE_LVNX16_I((xb_vecNx16*)R, 2 * BBE_SIMD_WIDTH * 2);
    r3 = BBE_LVNX16_I((xb_vecNx16*)R, 3 * BBE_SIMD_WIDTH * 2);
    r4 = BBE_LVNX16_I((xb_vecNx16*)R, 4 * BBE_SIMD_WIDTH * 2);
    r5 = BBE_LVNX16_I((xb_vecNx16*)R, 5 * BBE_SIMD_WIDTH * 2);
    r6 = BBE_LVNX16_I((xb_vecNx16*)R, 6 * BBE_SIMD_WIDTH * 2);
    r7 = BBE_LVNX16_I((xb_vecNx16*)R, 7 * BBE_SIMD_WIDTH * 2);

    BBE_DSELNX16I(r1, r0, r1, r0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(r3, r2, r3, r2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(r5, r4, r5, r4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(r7, r6, r7, r6, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(r2, r0, r2, r0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(r6, r4, r6, r4, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(r3, r1, r3, r1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(r7, r5, r7, r5, BBE_DSELI_DEINTERLEAVE_1);

    BBE_DSELNX16I(r4, r0, r4, r0, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(r5, r1, r5, r1, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(r6, r2, r6, r2, BBE_DSELI_DEINTERLEAVE_1);
    BBE_DSELNX16I(r7, r3, r7, r3, BBE_DSELI_DEINTERLEAVE_1);

    r0 = BBE_SELNX16I(BBE_ZERONX16(), r0, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    r1 = BBE_SELNX16I(BBE_ZERONX16(), r1, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    r2 = BBE_SELNX16I(BBE_ZERONX16(), r2, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    r3 = BBE_SELNX16I(BBE_ZERONX16(), r3, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    r4 = BBE_SELNX16I(BBE_ZERONX16(), r4, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    r5 = BBE_SELNX16I(BBE_ZERONX16(), r5, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    r6 = BBE_SELNX16I(BBE_ZERONX16(), r6, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    r7 = BBE_SELNX16I(BBE_ZERONX16(), r7, BBE_SELI_EXTRACT_1_OF_2_OFF_0);

    BBE_SVNX16_I(r0, (xb_vecNx16*)Rt, 0 * BBE_SIMD_WIDTH * 2);
    BBE_SVNX16_I(r1, (xb_vecNx16*)Rt, 1 * BBE_SIMD_WIDTH * 2);
    BBE_SVNX16_I(r2, (xb_vecNx16*)Rt, 2 * BBE_SIMD_WIDTH * 2);
    BBE_SVNX16_I(r3, (xb_vecNx16*)Rt, 3 * BBE_SIMD_WIDTH * 2);
    BBE_SVNX16_I(r4, (xb_vecNx16*)Rt, 4 * BBE_SIMD_WIDTH * 2);
    BBE_SVNX16_I(r5, (xb_vecNx16*)Rt, 5 * BBE_SIMD_WIDTH * 2);
    BBE_SVNX16_I(r6, (xb_vecNx16*)Rt, 6 * BBE_SIMD_WIDTH * 2);
    BBE_SVNX16_I(r7, (xb_vecNx16*)Rt, 7 * BBE_SIMD_WIDTH * 2);
}

/*
calculate main diagonal and its q
Input:
N        matrix size
R[8*8]  R matrix, Q11
Output:
D[8]    reciprocal of main diagonals
Returns D representation
*/
static int calcD(uint16_t * restrict D, const int16_t* R, int N)
{
    xb_vecNx16 x0, x1, b;
    xb_vecNx40 w, a;
    xb_int16 temp16;
    vsaN nsa;
    vboolN goodidx;
    int q;
    NASSERT_ALIGN32(R);
    NASSERT_ALIGN32(D);

    /* Load 8 diagonal elements. */
    BBE_LSNX16_XP(x1, R, (RSTRIDE + 1) * 2);      x0 = BBE_REPNX16(x1, 0);
    BBE_LSNX16_XP(x1, R, (RSTRIDE + 1) * 2);      x0 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_RIGHT_1);
    BBE_LSNX16_XP(x1, R, (RSTRIDE + 1) * 2);      x0 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_RIGHT_1);
    BBE_LSNX16_XP(x1, R, (RSTRIDE + 1) * 2);      x0 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_RIGHT_1);
    BBE_LSNX16_XP(x1, R, (RSTRIDE + 1) * 2);      x0 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_RIGHT_1);
    BBE_LSNX16_XP(x1, R, (RSTRIDE + 1) * 2);      x0 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_RIGHT_1);
    BBE_LSNX16_XP(x1, R, (RSTRIDE + 1) * 2);      x0 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_RIGHT_1);
    BBE_LSNX16_XP(x1, R, -7 * (RSTRIDE + 1) * 2); x0 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_RIGHT_1);
    x0 = BBE_SELNX16I(x0, x0, BBE_SELI_ROTATE_RIGHT_8);
    goodidx = BBE_LTRN(N);
    x1 = BBE_MOVVA16(0x7FFF);
    x0 = BBE_MOVNX16T(x0, x1, goodidx);   // mask unneeded data with 0x7fff

    nsa = BBE_NSANX16(x0);
    x1 = BBE_MOVVVS(nsa);
    temp16 = BBE_RMAXNX16(x1);
    q = xb_int16_rtor_int16(temp16);

    w = BBE_UNPKSNX16(x0);
    nsa = BBE_ADDSAVSN(24, nsa);
    w = BBE_SLSNX40(w, nsa);
    BBE_RECIPLUNX40_0(a, x0, b, w);
    BBE_RECIPLUNX40_1(a, x0, b, w);
    BBE_MULUSANX16(a, b, x0);
    /*nsa = BBE_SUBSAVSN(22 + q, nsa);
    a = BBE_SRANX40(a, nsa);
    x0 = BBE_PACKHNX40(a);*/
    /*nsa = BBE_SUBSAVSN(46 + q, nsa);
    a = BBE_SRANX40(a, nsa);
    x0 = BBE_PACKLNX40(a);*/
    nsa = BBE_SUBSAVSN(47 + q, nsa);
    a = BBE_RNDADJNX40(a, nsa);
    x0 = BBE_PACKVNX40(a, nsa);
    x0 = BBE_SLLINX16(x0, 1);

    BBE_SVNX16_I(x0, (xb_vecNx16*)D, 0);
    return 18 - q;
}

/*-----------------------------
Cholesky forward recursion:
y=R'\(A'*B)
Input:
A[M'*N]  Q15
B[M]     Q15
R[8*8]  Q11
M'=(M+15)&~15
Output:
Y[N]     Q27
-----------------------------*/
static void pfit_cholfwd(
    int32_t* restrict y,
    const int16_t* restrict R,
    const uint16_t* restrict D, int qD,
    const int16_t* restrict A,
    const int16_t* restrict B,
    int M, int N)
{
    const xb_vecNx16* restrict pR = (const xb_vecNx16*)R;
    const xb_vecNx16* restrict pAn;
    const xb_vecNx16* restrict pB;
    const int16_t* restrict pD = (const int16_t*)D;
    xb_int40 tmp40;
    xb_vecNx40 wA, wB, wC;
    xb_vecNx16 Rn, Y, lo, hi, Dn;
    valign vB, vR;
    vsaN vq4, vq, sha;

    int n;
    int m;//,q;
    const int Astride = (M + 15)&~15;
    vboolN bmask;
    xb_vecNx16 zero = BBE_ZERONX16(), mask;    // rotating mask to set y[n] in the Y

    //    q=40-nsa1x40(M);
    //    NASSERT(q>=2);
    NASSERT(qD >= 3 && qD <= 18);
    NASSERT(M > 1);
    NASSERT(N >= 2 && N <= 7);
    NASSERT_ALIGN32(R);
    NASSERT_ALIGN32(A);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(D);

    vR = BBE_LA_PP(pR);
    // set mask: first: 110000..., 2-nd: 00110000..., 3-rd: 00001100... and so on
    mask = BBE_MOVVINT16(1); lo = BBE_ZERONX16();
    bmask = BBE_LTRNI(2);
    mask = BBE_MOVNX16T(mask, lo, bmask);

    Y = BBE_ZERONX16();
    wA = BBE_MOVWA32(M);
    vq = BBE_NSANX40(wA);
    vq4 = BBE_SUBSAVSN(43, vq);
    sha = BBE_SUBSAVSN(43 + 12 - qD, vq);
    vq = BBE_ADDSAVSN(-35, vq);

    pAn = (const xb_vecNx16*)A;
    for (n = 0; n < N; n++)
    {
        // calculate A(:,n)'*B-Rn'*Y, 
        wB = BBE_ZERONX40();
        pB = (const xb_vecNx16*)B;
        vB = BBE_LA_PP(pB);
#ifdef COMPILER_XTENSA
#pragma loop_count min=1
#endif
        for (m = 0; m < (Astride >> 4); m++)
        {
            xb_vecNx16 a, b;
            BBE_LVNX16_IP(a, pAn, 2 * BBE_SIMD_WIDTH);
            BBE_LANX16_IP(b, vB, pB);
            BBE_MULANX16(wB, a, b);
        }
        wB = BBE_SRANX40(wB, vq4);

        BBE_LAVNX16_XP(Rn, vR, pR, 2 * N);
        lo = BBE_SELNX16I(Y, Y, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
        hi = BBE_SELNX16I(Y, Y, BBE_SELI_EXTRACT_1_OF_2_OFF_1);
        wC = BBE_MULUSNX16(lo, Rn);
        wC = BBE_SRAINX40(wC, 16);
        BBE_MULANX16(wC, hi, Rn);
        wC = BBE_SLANX40(wC, vq);
        wB = BBE_SUBNX40(wB, wC);
        tmp40 = BBE_RADDNX40(wB);
        wB = BBE_MOVNX40_FROM40(tmp40);
        BBE_LAVNX16_XP(Rn, vR, pR, 2 * (RSTRIDE - N)); // go to the next row

        // scale by the reciprocal of main diagonal
        TAKEHILO3(wB, hi, lo);
        BBE_LSNX16_IP(Dn, pD, 2);
        wA = BBE_MULUUNX16(lo, Dn);
        wA = BBE_SRAINX40(wA, 16);
        BBE_MULUSANX16(wA, Dn, hi);
        wA = BBE_SLANX40(wA, sha); //->Q27
        // copy one element back to y
        {
            xb_vecNx16 t;
            t = BBE_MOVVWL(wA);
            t = BBE_REPNX16C(t, 0);
            bmask = BBE_NEQNX16(mask, zero);
            Y = BBE_MOVNX16T(t, Y, bmask);
            mask = BBE_SELNX16I(mask, mask, BBE_SELI_ROTATE_LEFT_2);
        }
    }
    // save Y
    BBE_SVNX16_I(Y, (xb_vecNx16*)y, 0);
}

/*-----------------------------
Cholesky backward recursion:
x=R\y
Input:
R[8*8] Q11
y[N]    Q27
Output:
x[N]    Q23
-----------------------------*/
static void  pfit_cholbkw(
    int32_t* restrict x,
    const int16_t* restrict R,
    const uint16_t* restrict D, int qD,
    const int32_t* restrict y,
    int N)
{
    const xb_vecNx16 * pR;
    vsaN _qD22 = BBE_MOVVSA32(22 - qD);
    xb_vecNx40 wB, wA;
    xb_int40 i40;
    xb_vecNx16 t, xlo, xhi, hi, lo;
    vboolN bN;    // true for first N elements
    vboolN bm, bt; // true for m-th element only
    int m;

    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(R);
    NASSERT(N >= 2 && N <= 7);
    NASSERT(qD >= 3 && qD <= 18);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(D);

    pR = (xb_vecNx16*)((uintptr_t)R + (N - 1) * 2 * BBE_SIMD_WIDTH);
    bN = BBE_LTRN(N);
    xlo = xhi = BBE_ZERONX16(); // clean output

    for (m = N - 1; m >= 0; m--)
    {
        // calculate y(m,:)-R(m,:)*X, 1xP
        BBE_LVNX16_IP(t, pR, -2 * BBE_SIMD_WIDTH);
        wB = BBE_MULUSNX16(xlo, t);
        wB = BBE_SRAINX40(wB, 16);
        BBE_MULANX16T(wB, xhi, t, bN);
        i40 = BBE_RADDNX40(wB);
        wB = BBE_MOVNX40_FROM40(i40);
        wA = BBE_MOVWA32(y[m] >> 11);
        wB = BBE_SRAINX40(wB, 2);
        wB = BBE_SUBNX40(wA, wB); //  B_=A_-B_

        wB = BBE_REPNX40(wB, 0);
        TAKEHILO3(wB, hi, lo);
        // scale to reciprocal of diagonal
        t = BBE_MOVVA16(D[m]);
        wB = BBE_MULUUNX16(lo, t);
        wB = BBE_SRLINX40(wB, 16);
        BBE_MULUSANX16(wB, t, hi);
        wB = BBE_SLSNX40(wB, _qD22);
        TAKEHILO3(wB, hi, lo);

        // create selector for m-th element and save m-th element 
        bt = BBE_LTRN(m);
        bm = BBE_LTRN(m + 1);
        bm = BBE_XORB(bt, bm);
        xlo = BBE_MOVNX16T(lo, xlo, bm);
        xhi = BBE_MOVNX16T(hi, xhi, bm);
    }
    t = BBE_SELNX16I(xhi, xlo, BBE_SELI_INTERLEAVE_1_LO);
    BBE_SVNX16_I(t, (xb_vecNx16*)x, 0);
}

void pfit_process
(
void * restrict           pScr,
int32_t * restrict        p,
const int16_t * restrict  V,
const int16_t * restrict  R,
const int16_t * restrict  y,
int                       M,
int                       N,
int                       maxIter,
fndiscr                   discr)
{
    int32_t *yy; // [N+1]
    int32_t *xx; // [N+1]
    int32_t *pp; // [N+1]
    int16_t *Rt; // [8*16]
    uint16_t *D; // [N+1]
    int16_t *b;  // [M]
    int qD, iter;

    NASSERT_ALIGN32(pScr);
    NASSERT(N >= 1 && N <= 6);
    NASSERT(M >= N);

    { /* allocate scratch */
        uintptr_t a = (uintptr_t)pScr;
        /* allocate longer data first */
        yy = (int32_t*)a; a += (8)*sizeof(int32_t);
        xx = (int32_t*)a; a += (8)*sizeof(int32_t);
        pp = (int32_t*)a; a += (8)*sizeof(int32_t);
        Rt = (int16_t*)a; a += (8 * 16)*sizeof(int16_t);
        /* next allocate shorter */
        D = (uint16_t*)a; a += (16)*sizeof(uint16_t);
        b = (int16_t*)a;  a += ((M + 15)&(~15))*sizeof(int16_t);

        NASSERT((uint8_t*)a - (uint8_t*)pScr <= (int)pfit_process_getScratchSize(M, N));
    }
    N = N + 1;

    NASSERT_ALIGN32(xx);
    NASSERT_ALIGN32(yy);
    NASSERT_ALIGN32(D);
    NASSERT_ALIGN32(b);

    qD = calcD(D, R, N); /* calculate 1/diag(R) */
    transposeR(Rt, R);

    pfit_cholfwd(yy, R, D, qD, V, y, M, N);
    pfit_cholbkw(pp, Rt, D, qD, yy, N);

    /* calculate discrepancy : v=y-V*p */
    for (iter = 0; iter<maxIter; iter++)
    {
        discr(b, y, V, pp, M, N);
        pfit_cholfwd(yy, R, D, qD, V, b, M, N);
        pfit_cholbkw(xx, Rt, D, qD, yy, N);
        // add resulted xx to p
        {
            xb_vecNx16 t;
            xb_vecNx40 wp, wx;
            t = BBE_LVNX16_I((const xb_vecNx16*)xx, 0); wx = BBE_MOVSWVL(t);
            t = BBE_LVNX16_I((const xb_vecNx16*)pp, 0); wp = BBE_MOVSWVL(t);
            wp = BBE_ADDNX40(wp, wx);
            t = BBE_MOVSVWL(wp);
            BBE_SVNX16_I(t, (xb_vecNx16*)pp, 0);
        }
    }
    //------------------------------
    // final copying from pp to p
    //------------------------------
    {
        xb_vecNx16* P;
        valign vp;
        xb_vecNx16 t;
        t = BBE_LVNX16_I((xb_vecNx16*)pp, 0);
        vp = BBE_ZALIGN();
        P = (xb_vecNx16*)p;
        BBE_SAVNX16_XP(t, vp, P, 4 * N);
        BBE_SAPOS_FP(vp, P);
    }
}

#endif
