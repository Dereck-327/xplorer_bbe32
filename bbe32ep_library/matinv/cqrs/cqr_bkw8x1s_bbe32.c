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

#if HAVE_VSAMATH && HAVE_DIV && 1
#define CURRENT_M 8
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

void cqr_bkw8x1s (void* pScr, complex_fract16* restrict _B, const complex_fract16* restrict _R, int L)
{
          int16_t* restrict B=(      int16_t*)_B;
    const int16_t* restrict R=(const int16_t*)_R;
    const int q = 10;
    int i;

    xb_vecNx16 * restrict pb0  = (xb_vecNx16 *)(B); 
    xb_vecNx16 * restrict pb4  = (xb_vecNx16 *)(B + 4*(2*L)); 
    const xb_vecNx16 * restrict pR27 = (const xb_vecNx16 *)(R+2*L*(2*CURRENT_M+7)); 
    const xb_vecNx16 * restrict pR67 = (const xb_vecNx16 *)(R+2*L*(6*CURRENT_M+7)); 
    const xb_vecNx16 * restrict pr77 = (const xb_vecNx16 *)(R+7*(CURRENT_M+1)*L*2);

#ifdef COMPILER_XTENSA
#pragma ymemory( pr77 )
#pragma ymemory( pR27 ) 
#pragma ymemory( pR67 )
#endif

#define INVDIAG1(f0,ex,f0sav,ex0sav,pr,princ1,princ2)            \
{                                                                \
    const xb_vecNx16 c18= BBE_MOVVA16(30-4-q+3-1);              \
    vsaN nsa_vsa;                                               \
    xb_vecNx40 z0;                                              \
    xb_vecNx16 z4;                                              \
    z0 = BBE_MOVWA32(0x2000000);                                \
    z4 = BBE_MOVVINT16(4);                                      \
    BBE_LVNX16_XP(f0, pr, princ1);                              \
    BBE_LVNX16_XP(ex, pr, princ2);                              \
    f0 = BBE_SELNX16I ( ex, f0,  BBE_SELI_EXTRACT_1_OF_2_OFF_0);\
    nsa_vsa = BBE_NSANX16(f0);                                  \
    ex = BBE_MOVVVS(nsa_vsa);                                   \
    ex0sav = BBE_SUBNX16(c18, ex);                              \
    ex = BBE_SUBNX16(ex,z4);                                    \
    nsa_vsa = BBE_MOVVSV(ex, 0);                                \
    f0 =  BBE_SLANX16(f0,nsa_vsa);                              \
    f0sav = BBE_QUONX32(z0,f0);                                 \
    f0 = BBE_SHFLNX16I(f0sav, BBE_SHFLI_DOUBLE_1_LO);           \
    ex = BBE_SHFLNX16I(ex0sav, BBE_SHFLI_DOUBLE_1_LO);          \
}

#define INVDIAG2(f0,ex,f0sav,ex0sav)                            \
{                                                               \
    f0 = BBE_SHFLNX16I(f0sav, BBE_SHFLI_DOUBLE_1_HI);           \
    ex = BBE_SHFLNX16I(ex0sav, BBE_SHFLI_DOUBLE_1_HI);          \
}

    (void)pScr;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);;
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);;
    NASSERT(L%(BBE_SIMD_WIDTH/2) == 0 && L > 0);

    {
        xb_vecNx16
            r01, r02, r03, r04, r05, r06, r07,
            r12, r13, r14, r15, r16, r17,
            r23, r24, r25, r26, r27, 
            r34, r35, r36, r37,
            r45, r46, r47, 
            r56, r57,
            r67; 
        const int o1 = 1*2*L*sizeof(int16_t); 
        const int o2 = 2*2*L*sizeof(int16_t); 
        const int o3 = 3*2*L*sizeof(int16_t); 
        const int stepR = -(2*L*(CURRENT_M+1))*sizeof(int16_t);
        const int last_step_pR67 = sizeof(*pR67)-stepR*6;
        const int last_step_pR27 = sizeof(*pR27)-stepR*2;
        for(i=0; i<L; i+=(BBE_SIMD_WIDTH/2))
        {
            vsaN _11=BBE_MOVVSA32(11);
            vsaN ex_vsa;
            xb_vecNx16 dp;
            xb_vecNx16 ex, f0, b,f0sav,exsav; 
            xb_vecNx40 acc0;
            xb_vecNx16 b7,b6,b5,b4,b3,b2,b1,b0;

            // m=7
            INVDIAG1(f0,ex,f0sav,exsav,pr77,-(CURRENT_M+1)*L*4,-(CURRENT_M+1)*L*4);
            ex_vsa = BBE_MOVVSV(ex, 0); 
            ex_vsa=BBE_ADDSAVSN(1,ex_vsa);
            b = BBE_LVNX16_X(pb4, o3); 
            acc0 = BBE_MULRNX16(  f0, b,ex_vsa); 
            b7 = BBE_PACKVNX40(acc0, ex_vsa); 
            BBE_SVNX16_X(b7, pb4, o3); 
            // m=6
            INVDIAG2(f0,ex,f0sav,exsav);
            ex_vsa = BBE_MOVVSV(ex, 0); 
            BBE_LVNX16_XP(r67, pR67, stepR);
            acc0 = BBE_MULRNX16C(b7, r67,_11); 
            dp  = BBE_PACKVNX40(acc0,_11);
            b = BBE_LVNX16_X(pb4, o2); 
            b=BBE_SRAINX16(b,1);
            b = BBE_SUBNX16(b, dp);
            acc0 = BBE_MULRNX16(  f0, b,ex_vsa); 
            b6 = BBE_PACKVNX40(acc0, ex_vsa); 
            BBE_SVNX16_X(b6, pb4, o2); 
            // m=5
            INVDIAG1(f0,ex,f0sav,exsav,pr77,-(CURRENT_M+1)*L*4,-(CURRENT_M+1)*L*4);
            ex_vsa = BBE_MOVVSV(ex, 0); 
            r57 = BBE_LVNX16_X( pR67, o1);
            BBE_LVNX16_XP(r56, pR67, stepR);
            acc0 = BBE_MULRNX16C( b7, r57,_11); 
            BBE_MULANX16C(acc0, b6, r56); 
            dp = BBE_PACKVNX40(acc0,_11); 
            b = BBE_LVNX16_X(pb4, o1); 
            b=BBE_SRAINX16(b,1);
            b = BBE_SUBNX16(b, dp);
            acc0 = BBE_MULRNX16(  f0, b,ex_vsa); 
            b5 = BBE_PACKVNX40(acc0, ex_vsa); 
            BBE_SVNX16_X(b5, pb4, o1); 
            // m=4
            INVDIAG2(f0,ex,f0sav,exsav);
            ex_vsa = BBE_MOVVSV(ex, 0); 
            r47 = BBE_LVNX16_X(pR67, o2);
            r46 = BBE_LVNX16_X(pR67, o1);
            BBE_LVNX16_XP(r45, pR67, stepR);
            acc0 = BBE_MULRNX16C( b7, r47,_11); 
            BBE_MULANX16C(acc0, b6, r46); 
            BBE_MULANX16C(acc0, b5, r45); 
            dp =  BBE_PACKVNX40(acc0,_11); 
            b = BBE_LVNX16_I(pb4, 0); 
            b=BBE_SRAINX16(b,1);
            b = BBE_SUBNX16(b, dp);
            acc0 = BBE_MULRNX16(  f0, b,ex_vsa); 
            b4 = BBE_PACKVNX40(acc0, ex_vsa); 
            BBE_SVNX16_IP(b4, pb4, sizeof(*pb4)); 
            // m=3
            INVDIAG1(f0,ex,f0sav,exsav,pr77,-(CURRENT_M+1)*L*4,-(CURRENT_M+1)*L*4);
            ex_vsa = BBE_MOVVSV(ex, 0); 
            r37 = BBE_LVNX16_X(pR67, o3);
            r36 = BBE_LVNX16_X(pR67, o2);
            r35 = BBE_LVNX16_X(pR67, o1);
            BBE_LVNX16_XP(r34, pR67, stepR);
            acc0 = BBE_MULRNX16C( b7, r37,_11); 
            BBE_MULANX16C(acc0, b6, r36); 
            BBE_MULANX16C(acc0, b5, r35); 
            BBE_MULANX16C(acc0, b4, r34); 
            dp = BBE_PACKVNX40(acc0,_11);
            b = BBE_LVNX16_X(pb0, o3);
            b=BBE_SRAINX16(b,1);
            b = BBE_SUBNX16(b, dp);
            acc0 = BBE_MULRNX16(  f0, b,ex_vsa); 
            b3 = BBE_PACKVNX40(acc0, ex_vsa); 
            BBE_SVNX16_X(b3, pb0, o3);
            // m=2
            INVDIAG2(f0,ex,f0sav,exsav);
            ex_vsa = BBE_MOVVSV(ex, 0); 
            BBE_LVNX16_XP(r27, pR27, stepR);
            r26 = BBE_LVNX16_X(pR67, o3);
            r25 = BBE_LVNX16_X(pR67, o2);
            r24 = BBE_LVNX16_X(pR67, o1);
            BBE_LVNX16_XP(r23, pR67, stepR);
            acc0 = BBE_MULRNX16C(  b7, r27,_11); 
            BBE_MULANX16C(acc0, b6, r26); 
            BBE_MULANX16C(acc0, b5, r25); 
            BBE_MULANX16C(acc0, b4, r24); 
            BBE_MULANX16C(acc0, b3, r23); 
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
            dp = BBE_PACKVNX40(acc0,_11);
            b = BBE_LVNX16_X(pb0, o2);
            b=BBE_SRAINX16(b,1);
            b = BBE_SUBNX16(b, dp);
            acc0 = BBE_MULRNX16(  f0, b,ex_vsa); 
            b2 = BBE_PACKVNX40(acc0, ex_vsa); 
            BBE_SVNX16_X(b2, pb0, o2);
            // m=1
            INVDIAG1(f0,ex,f0sav,exsav,pr77,-(CURRENT_M+1)*L*4,2*BBE_SIMD_WIDTH+(CURRENT_M-1)*(CURRENT_M+1)*L*4);
            ex_vsa = BBE_MOVVSV(ex, 0); 
            r17 = BBE_LVNX16_X( pR27, o1);
            BBE_LVNX16_XP(r16, pR27, stepR);
            r15 = BBE_LVNX16_X(pR67, o3);
            r14 = BBE_LVNX16_X(pR67, o2);
            r13 = BBE_LVNX16_X(pR67, o1);
            BBE_LVNX16_XP(r12, pR67, stepR);
            acc0 = BBE_MULRNX16C( b7, r17,_11); 
            BBE_MULANX16C(acc0, b6, r16); 
            BBE_MULANX16C(acc0, b5, r15); 
            BBE_MULANX16C(acc0, b4, r14); 
            BBE_MULANX16C(acc0, b3, r13); 
            BBE_MULANX16C(acc0, b2, r12); 
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
            dp = BBE_PACKVNX40(acc0,_11);
            b = BBE_LVNX16_X(pb0, o1);
            b=BBE_SRAINX16(b,1);
            b = BBE_SUBNX16(b, dp);
            acc0 = BBE_MULRNX16(  f0, b,ex_vsa); 
            b1 = BBE_PACKVNX40(acc0, ex_vsa); 
            BBE_SVNX16_X(b1, pb0, o1);
            // m=0
            INVDIAG2(f0,ex,f0sav,exsav);
            ex_vsa = BBE_MOVVSV(ex, 0); 
            r07 = BBE_LVNX16_X( pR27, o2);
            r06 = BBE_LVNX16_X( pR27, o1);
            BBE_LVNX16_XP(r05, pR27, last_step_pR27);
            r04 = BBE_LVNX16_X(pR67, o3);
            r03 = BBE_LVNX16_X(pR67, o2);
            r02 = BBE_LVNX16_X(pR67, o1);
            BBE_LVNX16_XP(r01, pR67, last_step_pR67);
            acc0=BBE_MULRNX16C(  b7, r07,_11); 
            BBE_MULANX16C(acc0, b6, r06); 
            BBE_MULANX16C(acc0, b5, r05); 
            BBE_MULANX16C(acc0, b4, r04); 
            BBE_MULANX16C(acc0, b3, r03); 
            BBE_MULANX16C(acc0, b2, r02); 
            BBE_MULANX16C(acc0, b1, r01); 
            dp = BBE_PACKVNX40(acc0,_11);
            b = BBE_LVNX16_I(pb0, 0);
            b=BBE_SRAINX16(b,1);
            b = BBE_SUBNX16(b, dp);
            acc0 = BBE_MULRNX16(  f0, b,ex_vsa); 
            b0 = BBE_PACKVNX40(acc0, ex_vsa); 
            BBE_SVNX16_IP(b0, pb0, sizeof(*pb0));
        }
    }
} /* cqr_bkw8x1s() */
#else
DISCARD_FUN(void,cqr_bkw8x1s,(void* pScr, complex_fract16* restrict B,const complex_fract16* restrict R, int L))
#endif

size_t cqr_bkw8x1s_getScratchSize (int N, int P, int L)
{
    (void)P; (void)L;
    return 0;
} /* cqr_bkw8x1s_getScratchSize() */
