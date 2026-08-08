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
/* Common utility declarations. */
#include "cqr_common.h"

#if HAVE_VSAMATH && HAVE_DIV && 1

#define CURRENT_M 4
#define CURRENT_P 4
#if 0
#define UPDATEB(b,f0,shift)          \
{                                    \
    xb_vecNx16 t;                    \
    xb_vecNx40 acc0;                 \
    t = BBE_SUBSNX16(b, dp);         \
    acc0 = BBE_MULNX16( f0, t);      \
    b = BBE_PACKVNX40(acc0, shift);  \
}
#else
#define UPDATEB(b,f0,shift)          \
{                                    \
    xb_vecNx40 acc0;                 \
    acc0 = BBE_MULRNX16(f0, b,shift);\
    BBE_MULSNX16(acc0,f0,dp);        \
    b = BBE_PACKVNX40(acc0, shift);  \
}
#endif

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

void cqr_bkw4x4s (void* pScr, complex_fract16* restrict _B, const complex_fract16* restrict _R, int L)
{
          int16_t* restrict B=(      int16_t*)_B;
    const int16_t* restrict R=(const int16_t*)_R;
    const int q = 10;
    int i;
    int16_t *QB = B;

    const int16_t *start_R = R;
    xb_vecNx16 * restrict  pB0;
    xb_vecNx16 * restrict  pB1;
    xb_vecNx16 * restrict  pB2;
    xb_vecNx16 * restrict  pB3;

    const xb_vecNx16 * restrict  pr12;
    const xb_vecNx16 * restrict  pr23;
    const xb_vecNx16 * restrict  pr34;
    const xb_vecNx16 * restrict  pr44;

    xb_vecNx40 z0  = BBE_MOVWA32(0x2000000);
                xb_vecNx16 z4  = BBE_MOVVINT16(4);
                xb_vecNx16 c19 = BBE_MOVVA16(30-4-q+3);
                xb_vecNx16 nsa;
                vsaN nsa_vsa;
#define SHRB(b) b=BBE_SRAINX16(b,1);

#ifdef COMPILER_XTENSA
#pragma ymemory( pr12 )
#pragma ymemory( pr23 )
#pragma ymemory( pr34 )
#endif

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/2)==0);

    pr44 = (xb_vecNx16 *)(R+3*(CURRENT_M+1)*L*2);

#define INVDIAG1(f0,ex,f0sav,exsav,pr,princ1,princ2)              \
    {                                                             \
    xb_vecNx40 z0  = BBE_MOVWA32(0x2000000);                      \
    xb_vecNx16 z4  = BBE_MOVVINT16(4);                            \
    xb_vecNx16 c19 = BBE_MOVVA16(30-4-q+3);                       \
    xb_vecNx16 nsa;                                               \
    vsaN nsa_vsa;                                                 \
    BBE_LVNX16_XP(f0, pr, princ1);                                \
    BBE_LVNX16_XP(ex, pr, princ2);                                \
    f0 = BBE_SELNX16I ( ex, f0,  BBE_SELI_EXTRACT_1_OF_2_OFF_0);  \
    nsa_vsa = BBE_NSANX16(f0);                                    \
    nsa = BBE_MOVVVS(nsa_vsa);                                    \
    exsav = BBE_SUBNX16(c19, nsa);                                \
    nsa = BBE_SUBNX16(nsa,z4);                                    \
    nsa_vsa = BBE_MOVVSV(nsa, 0);                                 \
    f0 =  BBE_SLANX16(f0,nsa_vsa);                                \
    f0sav = BBE_QUONX32(z0,f0);                                   \
    f0 = BBE_SHFLNX16I(f0sav, BBE_SHFLI_DOUBLE_1_LO);             \
    ex = BBE_SHFLNX16I(exsav, BBE_SHFLI_DOUBLE_1_LO);             \
    }

#define INVDIAG2(f0,ex,f0sav,exsav)                              \
    {                                                            \
    f0 = BBE_SHFLNX16I(f0sav, BBE_SHFLI_DOUBLE_1_HI);            \
    ex = BBE_SHFLNX16I(exsav, BBE_SHFLI_DOUBLE_1_HI);            \
    }

    R = start_R;
    {

        const int start_m = CURRENT_M-1;

        pr34 = (xb_vecNx16 *) (R + (start_m-1)*CURRENT_M*L*2/*row = start_m-1*/ +start_m*L*2 /*column = start_m*/);
        pr23 = pr34 - (1*CURRENT_M+1)*L*2*sizeof(int16_t)/sizeof(*pr34);
        pr12 = pr34 - (2*CURRENT_M+2)*L*2*sizeof(int16_t)/sizeof(*pr34);

        pB0= (xb_vecNx16 *)(QB +start_m*CURRENT_P*L*2 +0*L*2);
        pB1= (xb_vecNx16 *)(QB +start_m*CURRENT_P*L*2 +1*L*2);
        pB2= (xb_vecNx16 *)(QB +start_m*CURRENT_P*L*2 +2*L*2);
        pB3= (xb_vecNx16 *)(QB +start_m*CURRENT_P*L*2 +3*L*2);
//#ifdef COMPILER_XTENSA
//#pragma no_reorder
//#endif
        for(i=0; i<L; i+=(BBE_SIMD_WIDTH/2))
        {
            vsaN _11=BBE_MOVVSA32(11);
            xb_vecNx16 ex,f0sav,exsav;
            xb_vecNx16 dp,b,b3;
            xb_vecNx40 acc0  ;
            xb_vecNx40 acc00 ;
            xb_vecNx16 inv_r, inv_r00;
            vsaN       inv_re, inv_r11e;
            xb_vecNx16 r34,r14,r13,r12,r24,r23;

            /**************cols 0...3 row 3 *********************/
            //INVDIAG1(inv_r,ex,f0sav,exsav,pr44,-(CURRENT_M+1)*L*4,-(CURRENT_M+1)*L*4);
            BBE_LVNX16_XP(inv_r, pr44, -(CURRENT_M+1)*L*4);
            BBE_LVNX16_XP(ex, pr44, -(CURRENT_M+1)*L*4);
            inv_r = BBE_SELNX16I ( ex, inv_r,  BBE_SELI_EXTRACT_1_OF_2_OFF_0);
            nsa_vsa = BBE_NSANX16(inv_r);
            nsa = BBE_MOVVVS(nsa_vsa);
            exsav = BBE_SUBNX16(c19, nsa);
            nsa = BBE_SUBNX16(nsa,z4);
            nsa_vsa = BBE_MOVVSV(nsa, 0);
            inv_r =  BBE_SLANX16(inv_r,nsa_vsa);
            f0sav = BBE_QUONX32(z0,inv_r);
            inv_r = BBE_SHFLNX16I(f0sav, BBE_SHFLI_DOUBLE_1_LO);
            ex = BBE_SHFLNX16I(exsav, BBE_SHFLI_DOUBLE_1_LO);

            inv_re = BBE_MOVVSV(ex, 0);
            b = BBE_LVNX16_I(pB0, 0);
            acc0 = BBE_MULRNX16( inv_r, b, inv_re);
            b3 = BBE_PACKVNX40(acc0, inv_re);
            BBE_SVNX16_XP(b3, pB0, -CURRENT_P*L*2*sizeof(int16_t));
            b = BBE_LVNX16_I(pB1, 0);
            acc0 = BBE_MULRNX16(inv_r, b, inv_re);
            b3 = BBE_PACKVNX40(acc0, inv_re);
            BBE_SVNX16_XP(b3, pB1, -CURRENT_P*L*2*sizeof(int16_t));
            b = BBE_LVNX16_I(pB2, 0);
            acc0 = BBE_MULRNX16( inv_r, b, inv_re);
            b3 = BBE_PACKVNX40(acc0, inv_re);
            BBE_SVNX16_XP(b3, pB2, -CURRENT_P*L*2*sizeof(int16_t));
            b = BBE_LVNX16_I(pB3, 0);
            acc0 = BBE_MULRNX16(inv_r, b, inv_re);
            b3 = BBE_PACKVNX40(acc0, inv_re);
            BBE_SVNX16_XP(b3, pB3, -CURRENT_P*L*2*sizeof(int16_t));
            /************** cols 0...3 row 2 *********************/
            //INVDIAG2(inv_r,ex,f0sav,exsav);
            inv_r = BBE_SHFLNX16I(f0sav, BBE_SHFLI_DOUBLE_1_HI);
            ex = BBE_SHFLNX16I(exsav, BBE_SHFLI_DOUBLE_1_HI);
            inv_re = BBE_MOVVSV(ex, 0);
            inv_re=BBE_ADDSAVSN(-1,inv_re);
            BBE_LVNX16_XP( r34, pr34,   -1*CURRENT_M*L*4);
            BBE_LVNX16_XP( r24, pr34,   -1*CURRENT_M*L*4);
            BBE_LVNX16_XP( r14, pr34, (2*BBE_SIMD_WIDTH)+2*CURRENT_M*L*4);
            b3 = BBE_LVNX16_X (pB0, 1*CURRENT_P*L*2*sizeof(int16_t));
            acc0 = BBE_MULRNX16C( b3, r34,_11);
            dp = BBE_PACKVNX40 (acc0,_11);
            b = BBE_LVNX16_I(pB0, 0);
            SHRB(b);
            UPDATEB(b,inv_r,inv_re);
            BBE_SVNX16_IP(b, pB0, 0);
            b3 = BBE_LVNX16_X (pB1, 1*CURRENT_P*L*2*sizeof(int16_t));
            acc0 = BBE_MULRNX16C(b3, r34,_11);
            dp = BBE_PACKVNX40 (acc0,_11);
            b = BBE_LVNX16_I(pB1, 0);
            SHRB(b);
            UPDATEB(b,inv_r,inv_re);
            BBE_SVNX16_IP(b, pB1, 0);
            b3 = BBE_LVNX16_X (pB2, 1*CURRENT_P*L*2*sizeof(int16_t));
            acc0 = BBE_MULRNX16C(b3, r34,_11);
            dp = BBE_PACKVNX40 (acc0,_11);
            b = BBE_LVNX16_I(pB2, 0);
            SHRB(b);
            UPDATEB(b,inv_r,inv_re);
            BBE_SVNX16_IP(b, pB2, 0);
            b3 = BBE_LVNX16_X (pB3, 1*CURRENT_P*L*2*sizeof(int16_t));
            acc0 = BBE_MULRNX16C( b3, r34,_11);
            dp = BBE_PACKVNX40 ( acc0,_11);
            b = BBE_LVNX16_I(pB3, 0);
            SHRB(b);
            UPDATEB(b,inv_r,inv_re);
            BBE_SVNX16_IP(b, pB3, 0);

            /************** col 0 row 1 *********************/
            //INVDIAG1(inv_r,ex,f0sav,exsav,pr44,-(CURRENT_M+1)*L*4,2*BBE_SIMD_WIDTH+(CURRENT_M-1)*(CURRENT_M+1)*L*4);
            BBE_LVNX16_XP(inv_r, pr44, -(CURRENT_M+1)*L*4);
            BBE_LVNX16_XP(ex, pr44, 2*BBE_SIMD_WIDTH+(CURRENT_M-1)*(CURRENT_M+1)*L*4);
            inv_r = BBE_SELNX16I ( ex, inv_r,  BBE_SELI_EXTRACT_1_OF_2_OFF_0);
            nsa_vsa = BBE_NSANX16(inv_r);
            nsa = BBE_MOVVVS(nsa_vsa);
            exsav = BBE_SUBNX16(c19, nsa);
            nsa = BBE_SUBNX16(nsa,z4);
            nsa_vsa = BBE_MOVVSV(nsa, 0);
            inv_r =  BBE_SLANX16(inv_r,nsa_vsa);
            f0sav = BBE_QUONX32(z0,inv_r);
            inv_r = BBE_SHFLNX16I(f0sav, BBE_SHFLI_DOUBLE_1_LO);
            ex = BBE_SHFLNX16I(exsav, BBE_SHFLI_DOUBLE_1_LO);

            inv_r11e = BBE_MOVVSV(ex, 0);
            inv_r11e=BBE_ADDSAVSN(-1,inv_r11e);
            //INVDIAG2(inv_r00,ex,f0sav,exsav);
            inv_r00 = BBE_SHFLNX16I(f0sav, BBE_SHFLI_DOUBLE_1_HI);
            ex = BBE_SHFLNX16I(exsav, BBE_SHFLI_DOUBLE_1_HI);

            inv_re = BBE_MOVVSV(ex, 0);
            inv_re=BBE_ADDSAVSN(-1,inv_re);
            BBE_LVNX16_XP( r23, pr23,   -1*CURRENT_M*L*4);
            BBE_LVNX16_XP( r13, pr23, (2*BBE_SIMD_WIDTH)+1*CURRENT_M*L*4);
            BBE_LVNX16_IP( r12, pr12, sizeof(*pr34));

            b3 = BBE_LVNX16_X (pB0, 1*CURRENT_P*L*2*sizeof(int16_t));
            BBE_LVNX16_XP (b, pB0,  -CURRENT_P*L*2*sizeof(int16_t));
            acc0 = BBE_MULRNX16C( b3, r24,_11);
            BBE_MULANX16C(acc0, b, r23);
            dp = BBE_PACKVNX40 (acc0,_11);
            acc00 = BBE_MULRNX16C( b3, r14,_11);
            BBE_MULANX16C(acc00, b, r13);
            b = BBE_LVNX16_I(pB0, 0);
            SHRB(b);
            UPDATEB(b,inv_r,inv_r11e);
            BBE_SVNX16_XP(b, pB0, -CURRENT_P*L*2*sizeof(int16_t));
            /************** col 0 row 0 *********************/
            BBE_MULANX16C(acc00, b, r12);
            dp = BBE_PACKVNX40 (acc00,_11);
            b = BBE_LVNX16_I(pB0, 0);
            SHRB(b);
            UPDATEB(b,inv_r00,inv_re);
            BBE_SVNX16_XP(b, pB0, ((CURRENT_M-1)*CURRENT_P*L*2+BBE_SIMD_WIDTH)*sizeof(int16_t));
            /************** col 1 row 1 *********************/
            b3 = BBE_LVNX16_X (pB1, 1*CURRENT_P*L*2*sizeof(int16_t));
            BBE_LVNX16_XP (b, pB1,  -CURRENT_P*L*2*sizeof(int16_t));
            acc0 = BBE_MULRNX16C(b3, r24,_11);
            BBE_MULANX16C(acc0, b, r23);
            dp = BBE_PACKVNX40 (acc0,_11);
            acc00 = BBE_MULRNX16C(b3, r14,_11);
            BBE_MULANX16C(acc00, b, r13);
            b = BBE_LVNX16_I(pB1, 0);
            SHRB(b);
            UPDATEB(b,inv_r,inv_r11e);
            BBE_SVNX16_XP(b, pB1, -CURRENT_P*L*2*sizeof(int16_t));
            /************** col 1 row 0 *********************/
            BBE_MULANX16C(acc00, b, r12);
            dp = BBE_PACKVNX40 (acc00,_11);
            b = BBE_LVNX16_I(pB1, 0);
            SHRB(b);
            UPDATEB(b,inv_r00,inv_re);
            BBE_SVNX16_XP(b, pB1, ((CURRENT_M-1)*CURRENT_P*L*2+BBE_SIMD_WIDTH) * sizeof(int16_t));
            /************** col 2 row 1 *********************/
            b3 = BBE_LVNX16_X (pB2, 1*CURRENT_P*L*2*sizeof(int16_t));
            BBE_LVNX16_XP (b, pB2,  -CURRENT_P*L*2*sizeof(int16_t));
            acc0 = BBE_MULRNX16C(b3, r24,_11);
            BBE_MULANX16C(acc0, b, r23);
            dp = BBE_PACKVNX40 (acc0,_11);
            acc00 = BBE_MULRNX16C(b3, r14,_11);
            BBE_MULANX16C(acc00, b, r13);
            b = BBE_LVNX16_I(pB2, 0);
            SHRB(b);
            UPDATEB(b,inv_r,inv_r11e);
            BBE_SVNX16_XP(b, pB2, -CURRENT_P*L*2*sizeof(int16_t));
            /************** col 2 row 0 *********************/
            BBE_MULANX16C(acc00, b, r12);
            dp = BBE_PACKVNX40 (acc00,_11);
            b = BBE_LVNX16_I(pB2, 0);
            SHRB(b);
            UPDATEB(b,inv_r00,inv_re);
            BBE_SVNX16_XP(b, pB2, ((CURRENT_M-1)*CURRENT_P*L*2+BBE_SIMD_WIDTH)*sizeof(int16_t));
            /************** col 3 row 1 *********************/
            b3 = BBE_LVNX16_X (pB3, 1*CURRENT_P*L*2*sizeof(int16_t));
            BBE_LVNX16_XP (b, pB3,  -CURRENT_P*L*2*sizeof(int16_t));
            acc0 = BBE_MULRNX16C( b3, r24,_11);
            BBE_MULANX16C(acc0, b, r23);
            dp = BBE_PACKVNX40 (acc0,_11);
            acc00 = BBE_MULRNX16C( b3, r14,_11);
            BBE_MULANX16C(acc00, b, r13);
            b = BBE_LVNX16_I(pB3, 0);
            SHRB(b);
            UPDATEB(b,inv_r,inv_r11e);
            BBE_SVNX16_XP(b, pB3, -CURRENT_P*L*2*sizeof(int16_t));
            /************** col 3 row 0 *********************/
            BBE_MULANX16C(acc00, b, r12);
            dp = BBE_PACKVNX40 ( acc00,_11);
            b = BBE_LVNX16_I(pB3, 0);
            SHRB(b);
            UPDATEB(b,inv_r00,inv_re);
            BBE_SVNX16_XP(b, pB3, ((CURRENT_M-1)*CURRENT_P*L*2+BBE_SIMD_WIDTH)*sizeof(int16_t));
        }
    }
} /* cqr_bkw4x4s() */
#else
DISCARD_FUN(void,cqr_bkw4x4s,(void* pScr, complex_fract16* restrict B,const complex_fract16* restrict R, int L))
#endif

size_t cqr_bkw4x4s_getScratchSize (int N, int P, int L)
{
    (void)N; (void)P; (void)L;
    return 0;
} /* cqr_bkw4x4s_getScratchSize() */
