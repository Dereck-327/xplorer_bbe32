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
  NatureDSP_Baseband library. FIR filters and Related Functions
    Blockwise Adaptive LMS Algorithm for Complex Data
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fir.h"

/*-------------------------------------------------------------------------
Blockwise Adaptive LMS Algorithm for Complex Data

Blockwise LMS algorithm performs filtering of complex input samples, 
computation of error over a block of reference samples and makes blockwise
update of IR to minimize the error output.
Algorithm includes FIR filtering, calculation of correlation between the 
error output and reference signal and IR taps update based on that 
correlation.
NOTES: 
1.  For N=1 this algorithm is equivalent to standard LMS, however a bigger 
    block size reduces the computational overhead keeping adaptation 
    properties.
2.  Right selection of N depends on the change rate of impulse response. 
    However, on static or slow varying channels convergence rate depends 
    on selected mu and M, but not on N.
3.  Computation of filter output is done using only higher 16-bit words of 
    IR coefficients (in fact converting them to Q15). However, when it 
    performs coefficient update, it uses full Q30 accuracy. 
4.  Functions use BBE32EP 40-bit accumulators for all computational steps. 
    Thus, the user needs to avoid very long blocks to prevent overflows 
    during IR update step. Typically, N should not exceed 128 to guarantee 
    the 40-bit range.

Parameters:
Temporary:
pScr          Scratch memory area of FIR_BLMS_SCRATCH(M,N) bytes
Input:
h[M]          Complex impulse response, Q30
r[N]          Reference (near end) complex data vector, Q15. First in time 
              value is in r[0].
x[(N+M-1)]    Input (far end) complex data vector, Q15. First in time value 
              is in x[0].
norm          Normalization factor: power of signal multiplied by N, Q31
mu            Adaptation coefficient (LMS step), Q15
N             Length of data block
M             Length of h
Output
e[N]          Estimated error, Q15
h[M]          Updated impulse response, Q30

Restrictions:
pScr,e,h,r,x  Must not overlap
pScr,e,h,r,x  Aligned on 32-byte boundary
N             Multiple of 8 
M             Multiple of 8 (applies to fir_blms() function only)
-------------------------------------------------------------------------*/

#ifdef BBE_SUBSAVSN
#define SUBSAVSN(y,x,e)   y=BBE_SUBSAVSN(x,e) 
#else
#define SUBSAVSN(y,x,e)              \
{                                    \
    xb_vecNx16 res;                  \
    res = BBE_SUBNX16(x, BBE_MOVVVS(e)); \
    y = BBE_MOVVSV(res, 0);             \
}
#endif

void fir_blms8(void * pScr, complex_fract16 * restrict e,
                            complex_fract32 * restrict h,
                      const complex_fract16 * restrict r,
                      const complex_fract16 * restrict x,
                      int32_t norm, int16_t mu, int N)
{
          xb_vecNx16* restrict pE;
    const xb_vecNx16* restrict pHr;
          xb_vecNx16* restrict pHw;
    const xb_vecNx16* restrict pX;
    const xb_vecNx16* restrict pR;
    int n;

    if (N <= 0) return;
    NASSERT(N % 8 == 0);
    NASSERT_ALIGN32(e);
    NASSERT_ALIGN32(h);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(r);

    __Pragma("no_reorder");
        // convert IR from Q30 to Q15 and save in reversed order
    {
        xb_vecNx40 A;
        xb_vecNx16 H, H_;

        pHw = (xb_vecNx16*)(pScr);
        pHr = (const xb_vecNx16*)(h + (8 * 8 - 2 * BBE_SIMD_WIDTH) / sizeof(complex_fract32));

        BBE_LVNX16_IP(H, pHr, -2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(H_,pHr, -2 * BBE_SIMD_WIDTH);
        A = BBE_MOVSWV(H, H_);
        H = BBE_PACKQNX40(A);
        H = BBE_SHFLNX16I(H, BBE_SHFLI_REVERSE_2);
        H = BBE_CONJSNX16C(H);
        BBE_SVNX16_IP(H, pHw, 2 * BBE_SIMD_WIDTH);
    }

    __Pragma("no_reorder");
        // compute error
    {
        xb_vecNx40 A;
        xb_vecNx16 X0, X1, H0, H, X, R0;

        pE = (xb_vecNx16*)e;
        pHr = (const xb_vecNx16*)pScr;
        pX = (const xb_vecNx16*)x;
        pR = (const xb_vecNx16*)r;

        BBE_LVNX16_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
        {
            BBE_LVNX16_IP(R0, pR, 2 * BBE_SIMD_WIDTH);
            A = BBE_UNPKQNX16(R0); //unpack with fractional format conversion from Q15 to Q9.30           

            BBE_LVNX16_IP(X1, pX, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(H0, pHr, 0);
            X = X0;                                            H = BBE_REPNX16C(H0, 0); BBE_MULSNX16C(A, X, H);
            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_2); H = BBE_REPNX16C(H0, 1); BBE_MULSNX16C(A, X, H);
            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_4); H = BBE_REPNX16C(H0, 2); BBE_MULSNX16C(A, X, H);
            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_6); H = BBE_REPNX16C(H0, 3); BBE_MULSNX16C(A, X, H);

            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_8);  H = BBE_REPNX16C(H0, 4); BBE_MULSNX16C(A, X, H);
            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_10); H = BBE_REPNX16C(H0, 5); BBE_MULSNX16C(A, X, H);
            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_12); H = BBE_REPNX16C(H0, 6); BBE_MULSNX16C(A, X, H);
            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_14); H = BBE_REPNX16C(H0, 7); BBE_MULSNX16C(A, X, H);
            X0 = X1;

            X = BBE_PACKQNX40(A);
            BBE_SVNX16_IP(X, pE, 2 * BBE_SIMD_WIDTH);
        }
    }

    __Pragma("no_reorder");
        // compute correlation 
    {
        xb_vecNx16 E0, E;
        xb_vecNx16 X0, X1, X;
        xb_vecNx40 A;

        pHw = (xb_vecNx16*)pScr;
        pE = (xb_vecNx16*)e;
        pX = (const xb_vecNx16*)x;

        BBE_LVNX16_IP(X0, pX, 2 * BBE_SIMD_WIDTH);

        A = 0;

        for (n = 0; n < N / (BBE_SIMD_WIDTH / 2); n++)
        {
            BBE_LVNX16_IP(X1, pX, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(E0, pE, 2 * BBE_SIMD_WIDTH);

            X = X0;                                             E = BBE_REPNX16C(E0, 0); BBE_MULANX16J(A, X, E);
            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_2);  E = BBE_REPNX16C(E0, 1); BBE_MULANX16J(A, X, E);
            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_4);  E = BBE_REPNX16C(E0, 2); BBE_MULANX16J(A, X, E);
            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_6);  E = BBE_REPNX16C(E0, 3); BBE_MULANX16J(A, X, E);

            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_8 ); E = BBE_REPNX16C(E0, 4); BBE_MULANX16J(A, X, E);
            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_10); E = BBE_REPNX16C(E0, 5); BBE_MULANX16J(A, X, E);
            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_12); E = BBE_REPNX16C(E0, 6); BBE_MULANX16J(A, X, E);
            X = BBE_SELNX16I(X1, X0, BBE_SELI_ROTATE_RIGHT_14); E = BBE_REPNX16C(E0, 7); BBE_MULANX16J(A, X, E);

            X0 = X1;
        }
        X0 = BBE_MOVSVWL(A);
        BBE_SVNX16_IP(X0, pHw, 2 * BBE_SIMD_WIDTH);

        X0 = BBE_MOVSVWH(A);
        BBE_SVNX16_IP(X0, pHw, 2 * BBE_SIMD_WIDTH);
    }

    __Pragma("no_reorder");
        // update IR
    {
        xb_vecNx16 X0, X1, X, bman;
        xb_vecNx40 A, B;
        vsaN a_e, sh;
        int b_exp;

        // compute normalization factor
        {
            int norm_exp, mu_exp;
            int16_t norm_n;
            int32_t rl, mu_n;

            mu_exp = XT_NSA(mu) - 16;
            mu_n = ((int32_t)mu) << mu_exp;
            mu_n <<= 14;
            // Q(15+mu_exp) <- Q15 + mu_exp
            norm_exp = XT_NSA(norm);
            // Q(15+norm_exp) <- Q31 + norm_exp - 16
            norm_n = (int16_t)((norm << norm_exp) >> 16);
            // Q(14 + mu_exp-norm_exp) <- Q(15+mu_exp+14)/Q(15+norm_exp)
            rl = mu_n / norm_n;
            bman = BBE_MOVVA16(rl);
            b_exp = mu_exp - norm_exp;
        }

        pHr = (const xb_vecNx16*)(h + (8 * 8 - 2 * BBE_SIMD_WIDTH) / sizeof(complex_fract32));
        pHw = (xb_vecNx16*)(pHr);
        pX = (const xb_vecNx16*)(pScr);

        BBE_LVNX16_IP(X0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(X1, pX, 2 * BBE_SIMD_WIDTH);

        A = BBE_MOVSWV(X1, X0);
        a_e = BBE_NSANX40(A);
        SUBSAVSN(sh, 24, a_e);
        X = BBE_PACKVNX40(A, sh);

        A = BBE_MULNX16(X, bman);
        SUBSAVSN(sh, 10 - b_exp, a_e);
        A = BBE_SLSNX40(A, sh);

        // reverse result
        X0 = BBE_MOVSVWL(A);
        X1 = BBE_MOVSVWH(A);
        X0 = BBE_SHFLNX16I(X0, BBE_SHFLI_REVERSE_4);
        X1 = BBE_SHFLNX16I(X1, BBE_SHFLI_REVERSE_4);
        A = BBE_MOVSWV(X0, X1);

        BBE_LVNX16_IP(X1, pHr, -2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(X0, pHr, -2 * BBE_SIMD_WIDTH);
        B = BBE_MOVSWV(X1, X0);
        A = BBE_ADDNX40(A, B);
        X0 = BBE_MOVSVWL(A);
        X1 = BBE_MOVSVWH(A);
        BBE_SVNX16_IP(X1, pHw, -2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(X0, pHw, -2 * BBE_SIMD_WIDTH);
    }
} /* fir_blms8() */
