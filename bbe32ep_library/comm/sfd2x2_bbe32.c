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
  NatureDSP_Baseband library. Communications
    Alamouti SFD
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_comm.h"

/*-------------------------------------------------------------------------
Alamouti SFD

Function calculates maximum likelihood decision by Alamouti formulas
for given Tx/Rx diversity. See 
A Simple Transmit Diversity Technique for Wireless Communications. Siavash 
M. Alamouti, IEEE JOURNAL ON SELECT AREAS IN COMMUNICATIONS, VOL. 16, NO. 
8, OCTOBER 1998
LTE - The UMTS Long Term Evolution: From Theory to Practice. Stefania Sesia,
Issam Toufik, Matthew Baker

Input:
r0,r1,r2,r3[2*N]  Input signals, complex vectors of N elements, 
                  fixed point presentation Qr
h[2*2][N][2]      Impulse response, sequence of N 2x2 complex matrices
                  stored in the streaming order, fixed point 
                  presentation Qh
N                 Size of vectors
rsh               Additional right shift amount for output data

Output:
s0,s1[2*N]        Output signals, complex vectors of N elements, 
                  fixed point presentation Qr+Qh-rsh

Restrictions:
Arrays must not overlap
Arrays must be aligned on 32-byte boundary
N must be a multiple of 8
-------------------------------------------------------------------------*/

void sfd2x2 ( int16_t * restrict s0,
              int16_t * restrict s1,
        const int16_t * restrict r0,
        const int16_t * restrict r1,
        const int16_t * restrict r2,
        const int16_t * restrict r3,
        const int16_t * restrict h,
        int N, int rsh )
{
    int n;
    xb_vecNx16 *  restrict S0 = (xb_vecNx16  *)s0;
    xb_vecNx16 *  restrict S1 = (xb_vecNx16  *)s1;
    xb_vecNx16 *  restrict R0 = (xb_vecNx16  *)r0;
    xb_vecNx16 *  restrict R1 = (xb_vecNx16  *)r1;
    xb_vecNx16 *  restrict R2 = (xb_vecNx16  *)r2;
    xb_vecNx16 *  restrict R3 = (xb_vecNx16  *)r3;
    xb_vecNx16 *  restrict H0 = (xb_vecNx16  *)h;
    xb_vecNx16    y0, y1, x0, x1, x2, x3, h0, h1, h2, h3;
    vsaN        shv;
    xb_vecNx40  A10, A32;
    /*
    Alamouti formula:
    s0=conj(h0).*r0+h1.*conj(r1)+conj(h2).*r2+h3.*conj(r3);
    s1=conj(h1).*r0-h0.*conj(r1)+conj(h3).*r2-h2.*conj(r3);
    */
    NASSERT_ALIGN(s0, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(s1, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(r0, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(r1, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(r2, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(r3, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(h, 2 * BBE_SIMD_WIDTH);
    NASSERT((N % (BBE_SIMD_WIDTH / 2)) == 0);
    NASSERT(rsh >= 0);
    shv = BBE_MOVVSA32(rsh + 16);

    for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
    {
        BBE_LVNX16_IP(x0, R0, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1, R1, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x2, R2, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x3, R3, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(h0, H0, 4 * N);
        BBE_LVNX16_XP(h1, H0, 4 * N);
        BBE_LVNX16_XP(h2, H0, 4 * N);
        BBE_LVNX16_XP(h3, H0, 2 * (-6 * N + BBE_SIMD_WIDTH));

        A10 = BBE_MULRNX16J(x0, h0, shv);
        BBE_MULANX16J(A10, h1, x1);
        BBE_MULANX16J(A10, x2, h2);
        BBE_MULANX16J(A10, h3, x3);
        A32 = BBE_MULRNX16J(x0, h1, shv);
        BBE_MULSNX16J(A32, h0, x1);
        BBE_MULANX16J(A32, x2, h3);
        BBE_MULSNX16J(A32, h2, x3);

        y0 = BBE_PACKVNX40(A10, shv);
        y1 = BBE_PACKVNX40(A32, shv);
        BBE_SVNX16_IP(y0, S0, 2 * BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(y1, S1, 2 * BBE_SIMD_WIDTH);
    }
} /* sfd2x2() */
