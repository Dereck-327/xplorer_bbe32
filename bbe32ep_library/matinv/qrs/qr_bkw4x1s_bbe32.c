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
    cqr_bkwNxPs/qr_bkwNxPs
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"

#if !(HAVE_VSAMATH && HAVE_DIV && 1)
DISCARD_FUN(void, qr_bkw4x1s, (void* pScr, int16_t* restrict B, const int16_t* restrict R, int L))
#else
#define CURRENT_M 4
#define CURRENT_P 1

/*-------------------------------------------------------------------------
cqr_bkwNxPs/qr_bkwNxPs

Last stage of solving a set of L complex-valued linear problems A*X=B
through the QR decomposition by Householder reflections: back substitution
process for L systems of complex-valued linear equations R*X=QB, where R is
an MxM upper triangular matrix, X is an MxP matrix of unknowns, QB is an MxP
matrix resulting from Householder reflections being applied to the right
hand matrix B of the original linear problem: QB=Q'*B.

Fixed-point representation for output data is a function of fixed-point
format of input data: FPP(X) = FPP(QB)-FPP(R)+10, where FPP(x) stands for
the Fixed-Point Position of data item x.

Data transform is performed in-place.

NOTE:
1. Data layout for matrices is selected as for other matrices written 
   in a stream order. So, shorter dimension of output matrix B (NxP 
   instead of MxP as on input) does not require special management - 
   remaining (M-N)*P*L elements are kept unchanged

Input
B[M*P][L]  Matrices QB=Q'*B (L matrices of size MxP)
R[M*N][L]  upper triangular matrices R (L matrices of size MxN)
Output:
B[N*P][L]  Matrices X (L matrices of size NxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(N,P,L)
4. Matrix sizes N,L must be greater than 1
---------------------------------------------------------------------------*/

void  qr_bkw4x1s (void* pScr, int16_t* restrict B, const int16_t* restrict R, int L)
{
    const xb_vecNx40 z0 = BBE_MOVWA32(0x20000000);
    const int q = 10;
    int i;

    xb_vecNx16 * restrict  pr11;
    xb_vecNx16 * restrict  pr22;
    xb_vecNx16 * restrict  pr33;
    xb_vecNx16 * restrict  pr44;

    xb_vecNx16 * restrict  pb11;
    xb_vecNx16 * restrict  pb22;
    xb_vecNx16 * restrict  pb33;
    xb_vecNx16 * restrict  pb44;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT((L%BBE_SIMD_WIDTH) == 0 && L > 0);
    (void)pScr;

    pr11 = (xb_vecNx16 *)(R + 0 * (CURRENT_M + 1)*L);// 0 *L
    pr22 = (xb_vecNx16 *)(R + 1 * (CURRENT_M + 1)*L);// 5 *L
    pr33 = (xb_vecNx16 *)(R + 2 * (CURRENT_M + 1)*L);// 10*L
    pr44 = (xb_vecNx16 *)(R + 3 * (CURRENT_M + 1)*L);// 15*L

    pb11 = (xb_vecNx16 *)(B + 0 * L);
    pb22 = (xb_vecNx16 *)(B + 1 * L);
    pb33 = (xb_vecNx16 *)(B + 2 * L);
    pb44 = (xb_vecNx16 *)(B + 3 * L);

    __Pragma("loop_count min=1");
    for (i = 0; i < (L>>LOG2_BBE_SIMD_WIDTH); i++)
    {
        vsaN e;
        vsaN sh11 = BBE_MOVVSA32(q + 1);
        vsaN ex_vsaN;

        xb_vecNx16 dp;
        xb_vecNx16 x0, f0, r0, b, x, x3, x2, x1;
        xb_vecNx40 acc0;

        //m=3;
        {
            BBE_LVNX16_XP(x0, pr44, BBE_SIMD_WIDTH*sizeof(int16_t)); // 15L
            e = BBE_NSANX16(x0);
            x0 = BBE_SLANX16(x0, e);
            f0 = BBE_QUONX32(z0, x0);
            ex_vsaN = BBE_SUBSAVSN(15 + 4 - 1, e);

            //dp = BBE_ZERONX16();
            b = BBE_LVNX16_X(pb44, 0);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            //BBE_MULSNX16(acc0, f0, dp);
            acc0 = BBE_RNDADJNX40(acc0, ex_vsaN);
            x3 = BBE_PACKVNX40(acc0, ex_vsaN);
            BBE_SVNX16_IP(x3, pb44, BBE_SIMD_WIDTH*sizeof(int16_t));
        }
        //m=2;
        {
            r0 = BBE_LVNX16_X(pr11, (3 * L + 2 * L*CURRENT_M)*sizeof(int16_t));
            acc0 = BBE_MULRNX16(x3, r0, sh11);

            BBE_LVNX16_XP(x0, pr33, BBE_SIMD_WIDTH*sizeof(int16_t));
            e = BBE_NSANX16(x0);
            x0 = BBE_SLANX16(x0, e);
            f0 = BBE_QUONX32(z0, x0);
            ex_vsaN = BBE_SUBSAVSN(15 + 4 - 1, e);

            dp = BBE_PACKVNX40(acc0, sh11);
            b = BBE_LVNX16_X(pb33, 0);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            acc0 = BBE_RNDADJNX40(acc0, ex_vsaN);
            x2 = BBE_PACKVNX40(acc0, ex_vsaN);
            BBE_SVNX16_IP(x2, pb33, BBE_SIMD_WIDTH*sizeof(int16_t));
        }
        //m=1;
        {
            r0 = BBE_LVNX16_X(pr11, (3 * L + 1 * L*CURRENT_M)*sizeof(int16_t));
            acc0 = BBE_MULRNX16(x3, r0, sh11);

            r0 = BBE_LVNX16_X(pr11, (2 * L + 1 * L*CURRENT_M)*sizeof(int16_t));
            BBE_MULANX16(acc0, x2, r0);

            BBE_LVNX16_XP(x0, pr22, BBE_SIMD_WIDTH*sizeof(int16_t));
            e = BBE_NSANX16(x0);
            x0 = BBE_SLANX16(x0, e);
            f0 = BBE_QUONX32(z0, x0);
            ex_vsaN = BBE_SUBSAVSN(15 + 4 - 1, e);

            dp = BBE_PACKVNX40(acc0, sh11);
            b = BBE_LVNX16_X(pb22, 0);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            acc0 = BBE_RNDADJNX40(acc0, ex_vsaN);
            x1 = BBE_PACKVNX40(acc0, ex_vsaN);
            BBE_SVNX16_IP(x1, pb22, BBE_SIMD_WIDTH*sizeof(int16_t));
        }
        //m=0;
        {
            r0 = BBE_LVNX16_X(pr11, (3 * L + 0 * L*CURRENT_M)*sizeof(int16_t));
            acc0 = BBE_MULRNX16(x3, r0, sh11);

            r0 = BBE_LVNX16_X(pr11, (2 * L + 0 * L*CURRENT_M)*sizeof(int16_t));
            BBE_MULANX16(acc0, x2, r0);

            r0 = BBE_LVNX16_X(pr11, (1 * L + 0 * L*CURRENT_M)*sizeof(int16_t));
            BBE_MULANX16(acc0, x1, r0);

            BBE_LVNX16_XP(x0, pr11, BBE_SIMD_WIDTH*sizeof(int16_t));
            e = BBE_NSANX16(x0);
            x0 = BBE_SLANX16(x0, e);
            f0 = BBE_QUONX32(z0, x0);
            ex_vsaN = BBE_SUBSAVSN(15 + 4 - 1, e);

            dp = BBE_PACKVNX40(acc0, sh11);
            b = BBE_LVNX16_X(pb11, 0);
            acc0 = BBE_MULNX16(f0, b);
            acc0 = BBE_SRAINX40(acc0, 1);
            BBE_MULSNX16(acc0, f0, dp);
            acc0 = BBE_RNDADJNX40(acc0, ex_vsaN);
            x = BBE_PACKVNX40(acc0, ex_vsaN);
            BBE_SVNX16_IP(x, pb11, BBE_SIMD_WIDTH*sizeof(int16_t));
        }
    }//for(i=0; i<L; i+=BBE_SIMD_WIDTH) 
} /* qr_bkw4x1s() */
#endif

size_t  qr_bkw4x1s_getScratchSize (int N, int P, int L)
{
    (void)N;(void)P;(void)L;
    return 0;
} /* qr_bkw4x1s_getScratchSize() */
