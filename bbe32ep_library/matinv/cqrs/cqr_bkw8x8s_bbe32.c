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

#if (HAVE_DIV && 1)

#define CURRENT_M 8
#define CURRENT_P 8

#undef __Pragma
#define __Pragma(x)

#define SOLVERXP_MCYCLE_BODY(_iii_)                       \
{                                                         \
                                                          \
    vsaN _11=BBE_MOVVSA32(11);                            \
    vsaN ex_vsa;                                          \
    BBE_LVNX16_IP(f0, pf, (2*BBE_SIMD_WIDTH));            \
    BBE_LVNX16_IP(ex, pf, (2*BBE_SIMD_WIDTH));            \
                                                          \
    ex_vsa = BBE_MOVVSV(ex, 0);                           \
    /*ex_vsa=BBE_ADDSAVSN(-1,ex_vsa);*/                 \
    if(_iii_==7)                                          \
    {                                                     \
        acc0 = 0;                            \
    }                                                     \
    else if(_iii_==6)                                     \
    {                                                     \
            BBE_LVNX16_XP(r0, pR0, R_last_step);          \
            acc0 = BBE_MULRNX16C( Xreg[7], r0,_11);       \
            R_last_step+=R_last_step_inc;                 \
    }                                                     \
    else if(_iii_==5)                                     \
    {                                                     \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            acc0 = BBE_MULRNX16C( Xreg[7], r0,_11);       \
            BBE_LVNX16_XP(r0, pR0, R_last_step);          \
            BBE_MULANX16C(acc0, Xreg[6], r0);             \
            R_last_step+=R_last_step_inc;                 \
    }                                                     \
    else if(_iii_==4)                                     \
    {                                                     \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            acc0 = BBE_MULRNX16C( Xreg[7], r0,_11);       \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[6], r0);             \
            BBE_LVNX16_XP(r0, pR0, R_last_step);          \
            BBE_MULANX16C(acc0, Xreg[5], r0);             \
            R_last_step+=R_last_step_inc;                 \
    }                                                     \
    else if(_iii_==3)                                     \
    {                                                     \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            acc0 = BBE_MULRNX16C( Xreg[7], r0,_11);       \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[6], r0);             \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[5], r0);             \
            BBE_LVNX16_XP(r0, pR0, R_last_step);          \
            BBE_MULANX16C(acc0, Xreg[4], r0);             \
            R_last_step+=R_last_step_inc;                 \
    }                                                     \
    else if(_iii_==2)                                     \
    {                                                     \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            acc0 = BBE_MULRNX16C( Xreg[7], r0,_11);       \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[6], r0);             \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[5], r0);             \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[4], r0);             \
            BBE_LVNX16_XP(r0, pR0, R_last_step);          \
            BBE_MULANX16C(acc0, Xreg[3], r0);             \
            R_last_step+=R_last_step_inc;                 \
    }                                                     \
    else if(_iii_==1)                                     \
    {                                                     \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            acc0 = BBE_MULRNX16C( Xreg[7], r0,_11);       \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[6], r0);             \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[5], r0);             \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[4], r0);             \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[3], r0);             \
            BBE_LVNX16_XP(r0, pR0, R_last_step);          \
            BBE_MULANX16C(acc0, Xreg[2], r0);             \
            R_last_step+=R_last_step_inc;                 \
    }                                                     \
    else if(_iii_==0)                                     \
    {                                                     \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            acc0 = BBE_MULRNX16C( Xreg[7], r0,_11);       \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[6], r0);             \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[5], r0);             \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[4], r0);             \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[3], r0);             \
            BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t)); \
            BBE_MULANX16C(acc0, Xreg[2], r0);             \
            BBE_LVNX16_XP(r0, pR0, R_last_step);          \
            BBE_MULANX16C(acc0, Xreg[1], r0);             \
            R_last_step+=R_last_step_inc;                 \
    }                                                     \
    dp = BBE_PACKVNX40(acc0,_11);                         \
    BBE_LVNX16_XP(b, pB, -CURRENT_P*L*2*sizeof(int16_t)); \
    b=BBE_SRAINX16(b,1);                                  \
    b = BBE_SUBNX16(b, dp);                               \
    acc0 = BBE_MULRNX16(  f0, b,ex_vsa);                  \
    x = BBE_PACKVNX40(acc0, ex_vsa);                      \
    Xreg[_iii_] = x;                                      \
    BBE_SVNX16_XP(x, pX, -CURRENT_P*L*2*sizeof(int16_t)); \
}

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
void cqr_bkw8x8s (void* pScr, complex_fract16* restrict _B, const complex_fract16* restrict _R, int L)
{
          int16_t* restrict B=(      int16_t*)_B;
    const int16_t* restrict R=(const int16_t*)_R;
    const int q = 10;
    int i; 
    int16_t *invDiagR = (int16_t *)pScr; 
    int16_t *QB = B; 
    int16_t *X = QB; 
    const xb_vecNx16 c0 = BBE_MOVVA16(30-4-q-1); 
    xb_vecNx40 z0;
    xb_vecNx16 z4, z3;

    xb_vecNx16 * restrict pf ; 
    xb_vecNx16 * restrict pr00 = (xb_vecNx16 *)R; 
    xb_vecNx16 * restrict pr11 = (xb_vecNx16 *)(R+1*(CURRENT_M+1)*L*2); 
    xb_vecNx16 * restrict pr22 = (xb_vecNx16 *)(R+2*(CURRENT_M+1)*L*2);
    xb_vecNx16 * restrict pr33 = (xb_vecNx16 *)(R+3*(CURRENT_M+1)*L*2);
    xb_vecNx16 * restrict pr44 = (xb_vecNx16 *)(R+4*(CURRENT_M+1)*L*2); 
    xb_vecNx16 * restrict pr55 = (xb_vecNx16 *)(R+5*(CURRENT_M+1)*L*2);
    xb_vecNx16 * restrict pr66 = (xb_vecNx16 *)(R+6*(CURRENT_M+1)*L*2); 
    xb_vecNx16 * restrict pr77 = (xb_vecNx16 *)(R+7*(CURRENT_M+1)*L*2);
    xb_vecNx16 * restrict pR0,  *pR00; 
    xb_vecNx16 * restrict pX; 
    xb_vecNx16 * restrict pB; 

    z0 = BBE_MOVWA32(0x2000000);
    z4 = BBE_MOVVINT16(4);
    z3 = BBE_MOVVINT16(3);
#ifdef COMPILER_XTENSA
#pragma ymemory( pr00 ) 
#pragma ymemory( pr11 )
#pragma ymemory( pr22 ) 
#pragma ymemory( pr33 )
#pragma ymemory( pr44 ) 
#pragma ymemory( pr55 )
#pragma ymemory( pr66 ) 
#pragma ymemory( pr77 )
#pragma ymemory( pr66 ) 
#pragma ymemory( pr77 )
#pragma ymemory( pR0 ) 
#pragma ymemory( pR00 )
#endif
    pf = (xb_vecNx16 *) invDiagR ;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);;
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);;
    NASSERT(L%(BBE_SIMD_WIDTH/2) == 0 && L > 0);

    for(i=0; i<L; i+=(BBE_SIMD_WIDTH/2))
    {
        xb_vecNx16 f0,nsa,ex;
        vsaN nsa_vsa; 
        xb_vecNx16 x0, x1, x2, x3; 

        BBE_LVNX16_IP(x0, pr77, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_IP(x1, pr66, (2*BBE_SIMD_WIDTH));

        x0 = BBE_SELNX16I ( x1, x0,  BBE_SELI_EXTRACT_1_OF_2_OFF_0); 

        {
            nsa_vsa = BBE_NSANX16(x0);
            nsa = BBE_MOVVVS(nsa_vsa); 
            ex = BBE_SUBNX16(nsa,z3);
            nsa = BBE_SUBNX16(nsa,z4);
            nsa_vsa = BBE_MOVVSV(nsa, 0);
            x0 =  BBE_SLANX16(x0,nsa_vsa);
            f0 = BBE_QUONX32(z0,x0);
            ex = BBE_SUBNX16(c0, ex);
            x0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_LO); 
            x1 = BBE_SHFLNX16I(ex, BBE_SHFLI_DOUBLE_1_LO); 
            x2 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_HI); 
            x3 = BBE_SHFLNX16I(ex, BBE_SHFLI_DOUBLE_1_HI); 

            BBE_SVNX16_IP(x0, pf,(2*BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(x1, pf,(2*BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(x2, pf,(2*BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(x3, pf,(2*BBE_SIMD_WIDTH));
        }
        BBE_LVNX16_IP(x2, pr55, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_IP(x3, pr44, (2*BBE_SIMD_WIDTH));
        x2 = BBE_SELNX16I ( x3, x2,  BBE_SELI_EXTRACT_1_OF_2_OFF_0); 

        {
            nsa_vsa = BBE_NSANX16(x2);
            nsa = BBE_MOVVVS(nsa_vsa); 

            ex = BBE_SUBNX16(nsa,z3);
            nsa = BBE_SUBNX16(nsa,z4);
            nsa_vsa = BBE_MOVVSV(nsa, 0);
            x2 =  BBE_SLANX16(x2,nsa_vsa);
            f0 = BBE_QUONX32(z0,x2);
            ex = BBE_SUBNX16(c0, ex);
            x0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_LO);
            x1 = BBE_SHFLNX16I(ex, BBE_SHFLI_DOUBLE_1_LO);
            x2 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_HI);
            x3 = BBE_SHFLNX16I(ex, BBE_SHFLI_DOUBLE_1_HI);

            BBE_SVNX16_IP(x0, pf,(2*BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(x1, pf,(2*BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(x2, pf,(2*BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(x3, pf,(2*BBE_SIMD_WIDTH));
        }

        BBE_LVNX16_IP(x0, pr33, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_IP(x1, pr22, (2*BBE_SIMD_WIDTH));

        x0 = BBE_SELNX16I ( x1, x0,  BBE_SELI_EXTRACT_1_OF_2_OFF_0); 

        {
            nsa_vsa = BBE_NSANX16(x0);
            nsa = BBE_MOVVVS(nsa_vsa); 
            ex = BBE_SUBNX16(nsa,z3);
            nsa = BBE_SUBNX16(nsa,z4);
            nsa_vsa = BBE_MOVVSV(nsa, 0);
            x0 =  BBE_SLANX16(x0,nsa_vsa);
            f0 = BBE_QUONX32(z0,x0);
            ex = BBE_SUBNX16(c0, ex);
            x0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_LO); 
            x1 = BBE_SHFLNX16I(ex, BBE_SHFLI_DOUBLE_1_LO); 
            x2 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_HI); 
            x3 = BBE_SHFLNX16I(ex, BBE_SHFLI_DOUBLE_1_HI); 
            BBE_SVNX16_IP(x0, pf,(2*BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(x1, pf,(2*BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(x2, pf,(2*BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(x3, pf,(2*BBE_SIMD_WIDTH));
        }
        BBE_LVNX16_IP(x2, pr11, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_IP(x3, pr00, (2*BBE_SIMD_WIDTH));
        x2 = BBE_SELNX16I ( x3, x2,  BBE_SELI_EXTRACT_1_OF_2_OFF_0); 

        {
            nsa_vsa = BBE_NSANX16(x2);
            nsa = BBE_MOVVVS(nsa_vsa); 
            ex = BBE_SUBNX16(nsa,z3);
            nsa = BBE_SUBNX16(nsa,z4);
            nsa_vsa = BBE_MOVVSV(nsa, 0);
            x2 =  BBE_SLANX16(x2,nsa_vsa);
            f0 = BBE_QUONX32(z0,x2);
            ex = BBE_SUBNX16(c0, ex);
            x0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_LO);
            x1 = BBE_SHFLNX16I(ex, BBE_SHFLI_DOUBLE_1_LO);
            x2 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_HI);
            x3 = BBE_SHFLNX16I(ex, BBE_SHFLI_DOUBLE_1_HI);
            BBE_SVNX16_IP(x0, pf,(2*BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(x1, pf,(2*BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(x2, pf,(2*BBE_SIMD_WIDTH));
            BBE_SVNX16_IP(x3, pf,(2*BBE_SIMD_WIDTH));
        }
    }

    for(i=0; i<L; i+=(BBE_SIMD_WIDTH/2))
    {
        const int end_m = CURRENT_M-1-7;
        int k;

        const int R_last_step0 = -CURRENT_M*L*2*sizeof(int16_t); 
        const int R_last_step_inc =L*2*sizeof(int16_t) ; 
        const int X_last_step  = (((CURRENT_M-1)-end_m+1) * CURRENT_P-1)*L*2 *sizeof(int16_t); 
        int R_last_step; 
        pR00 = (xb_vecNx16 *) (R + ((CURRENT_M-1)-1)*CURRENT_M*L*2/*row = (CURRENT_M-1)-1*/ +(CURRENT_M-1)*L*2 /*column = (CURRENT_M-1)*/); 
        pX = (xb_vecNx16 *)(X  +(CURRENT_M-1)*CURRENT_P*L*2+(CURRENT_P-1)*L*2); 
        pB = (xb_vecNx16 *)(QB +(CURRENT_M-1)*CURRENT_P*L*2+(CURRENT_P-1)*L*2); 

        for(k=CURRENT_P-1; k>=0; k--)
        {
            xb_vecNx16 dp;
            xb_vecNx16 r0, ex, f0, b, x; 
            xb_vecNx40 acc0;
            xb_vecNx16 Xreg[8]; 
            pf = (xb_vecNx16*)invDiagR; 
            pR0 = pR00; 
            R_last_step = R_last_step0;

            SOLVERXP_MCYCLE_BODY(7);
            SOLVERXP_MCYCLE_BODY(6);
            SOLVERXP_MCYCLE_BODY(5);
            SOLVERXP_MCYCLE_BODY(4);
            SOLVERXP_MCYCLE_BODY(3);
            SOLVERXP_MCYCLE_BODY(2);
            SOLVERXP_MCYCLE_BODY(1);
            SOLVERXP_MCYCLE_BODY(0);
            pX += X_last_step/sizeof(*pX); 
            pB += X_last_step/sizeof(*pB);
        } 

        invDiagR+=2*BBE_SIMD_WIDTH*CURRENT_M;
        X+=BBE_SIMD_WIDTH;
        R+=BBE_SIMD_WIDTH;
        QB+=BBE_SIMD_WIDTH;
    }

} /* cqr_bkw8x8s() */

size_t cqr_bkw8x8s_getScratchSize (int N, int P, int L)
{
    (void)P;
    return (L)*4*(N)*sizeof(int16_t);
} /* cqr_bkw8x8s_getScratchSize() */

#else
DISCARD_FUN(void, cqr_bkw8x8s,  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L))
size_t cqr_bkw8x8s_getScratchSize (int N, int P, int L) { (void)N;(void)P;(void)L; return 0; }
#endif
