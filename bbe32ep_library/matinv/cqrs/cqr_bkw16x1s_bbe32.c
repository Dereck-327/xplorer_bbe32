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

#define CURRENT_M 16
#define CURRENT_P 1

#define SOLVERXP_MCYCLE_BODY_V2(_iii_, shiftb)                   \
{                                                        \
    vsaN _11=BBE_MOVVSA32(11);                           \
    vsaN ex_vsa;                                         \
    BBE_LVNX16_IP(ex, pe, -(2*BBE_SIMD_WIDTH));          \
    BBE_LVNX16_IP(f0, pf, -(2*BBE_SIMD_WIDTH));          \
    ex_vsa = BBE_MOVVSV(ex, 0);                          \
    ex_vsa=BBE_ADDSAVSN(-1,ex_vsa);                    \
   if(_iii_==7)                                          \
   {                                                     \
    acc0 = 0;                               \
   }                                                     \
   else if(_iii_==6)                                     \
   {                                                     \
     BBE_LVNX16_XP(r0, pR0, R_last_step);                \
     acc0=BBE_MULRNX16C(Xreg[7], r0,_11);                \
     R_last_step+=R_last_step_inc;                       \
   }                                                     \
   else if(_iii_==5)                                     \
   {                                                     \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     acc0=BBE_MULRNX16C(Xreg[7], r0,_11);                \
     BBE_LVNX16_XP(r0, pR0, R_last_step);                \
     BBE_MULANX16C(acc0, Xreg[6], r0);                   \
     R_last_step+=R_last_step_inc;                       \
   }                                                     \
   else if(_iii_==4)                                     \
   {                                                     \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     acc0=BBE_MULRNX16C(Xreg[7], r0,_11);                \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[6], r0);                   \
     BBE_LVNX16_XP(r0, pR0, R_last_step);                \
     BBE_MULANX16C(acc0, Xreg[5], r0);                   \
     R_last_step+=R_last_step_inc;                       \
   }                                                     \
   else if(_iii_==3)                                     \
   {                                                     \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     acc0=BBE_MULRNX16C(Xreg[7], r0,_11);                \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[6], r0);                   \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[5], r0);                   \
     BBE_LVNX16_XP(r0, pR0, R_last_step);                \
     BBE_MULANX16C(acc0, Xreg[4], r0);                   \
     R_last_step+=R_last_step_inc;                       \
   }                                                     \
   else if(_iii_==2)                                     \
   {                                                     \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     acc0=BBE_MULRNX16C(Xreg[7], r0,_11);                \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[6], r0);                   \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[5], r0);                   \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[4], r0);                   \
     BBE_LVNX16_XP(r0, pR0, R_last_step);                \
     BBE_MULANX16C(acc0, Xreg[3], r0);                   \
     R_last_step+=R_last_step_inc;                       \
   }                                                     \
   else if(_iii_==1)                                     \
   {                                                     \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     acc0=BBE_MULRNX16C(Xreg[7], r0,_11);                \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[6], r0);                   \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[5], r0);                   \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[4], r0);                   \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[3], r0);                   \
     BBE_LVNX16_XP(r0, pR0, R_last_step);                \
     BBE_MULANX16C(acc0, Xreg[2], r0);                   \
     R_last_step+=R_last_step_inc;                       \
   }                                                     \
   else if(_iii_==0)                                     \
   {                                                     \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     acc0=BBE_MULRNX16C(Xreg[7], r0,_11);                \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[6], r0);                   \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[5], r0);                   \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[4], r0);                   \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[3], r0);                   \
     BBE_LVNX16_XP(r0, pR0, -L*2*sizeof(int16_t));       \
     BBE_MULANX16C(acc0, Xreg[2], r0);                   \
     BBE_LVNX16_XP(r0, pR0, R_last_step);                \
     BBE_MULANX16C(acc0, Xreg[1], r0);                   \
     R_last_step+=R_last_step_inc;                       \
   }                                                     \
   dp = BBE_PACKVNX40(acc0, _11);                        \
   BBE_LVNX16_XP(b, pB, -CURRENT_P*L*2*sizeof(int16_t)); \
    if(shiftb) b=BBE_SRAINX16(b, shiftb);                              \
   b = BBE_SUBNX16(b, dp);                               \
   acc0 = BBE_MULRNX16(  f0, b, ex_vsa);                 \
   x = BBE_PACKVNX40(acc0, ex_vsa);                      \
   Xreg[_iii_] = x;                                      \
   BBE_SVNX16_XP(x, pX, -CURRENT_P*L*2*sizeof(int16_t)); \
}

inline_ void SolveRxP_M16_P16x8(
                               int16_t *X, 
                               const int16_t *R, 
                               int16_t *QB,
                               int L,
                               int start_m, 
                               int16_t *pe_end, 
                               int16_t *pf_end
                               )
                               ATTRIBUTE_ALWAYS_INLINE; 

inline_ void SolveRxP_M16_P16x8_last8(
                               int16_t *X,
                               const int16_t *R,
                               int16_t *QB,
                               int L,
                               int start_m, 
                               int16_t *pe_end,
                               int16_t *pf_end

                               )
                                ATTRIBUTE_ALWAYS_INLINE; 

inline_ void SolveRxP_M16_P16x8(
                               int16_t *X, 
                               const int16_t *R, 
                               int16_t *QB,
                               int L,
                               int start_m, 
                               int16_t *pe_end, 
                               int16_t *pf_end
                               )
{
    const int end_m = start_m-7;

    xb_vecNx16 * pR0,  *pR00; 
    xb_vecNx16 * pX; 
    xb_vecNx16 * pB; 

    int k;
    const int R_last_step0 = -CURRENT_M*L*2*sizeof(int16_t); 
    const int R_last_step_inc =L*2*sizeof(int16_t) ; 
    const int X_last_step  = ((start_m-end_m+1) * CURRENT_P-1)*L*2 *sizeof(int16_t); 
    int R_last_step; 
    pR00 = (xb_vecNx16 *) (R + (start_m-1)*CURRENT_M*L*2/*row = start_m-1*/ +start_m*L*2 /*column = start_m*/); 
    pX = (xb_vecNx16 *)(X  +start_m*CURRENT_P*L*2+(CURRENT_P-1)*L*2); 
    pB = (xb_vecNx16 *)(QB +start_m*CURRENT_P*L*2+(CURRENT_P-1)*L*2); 

    for(k=CURRENT_P-1; k>=0; k--)
    {
        xb_vecNx16  *pe, *pf; 
        xb_vecNx16 dp;
        xb_vecNx16 r0, ex, f0, b, x; 
        xb_vecNx40 acc0;//, acc1;
        xb_vecNx16 Xreg[8]; 

        pR0 = pR00; 
        pe = (xb_vecNx16  *)pe_end; 
        pf = (xb_vecNx16  *)pf_end;
        R_last_step = R_last_step0;

        SOLVERXP_MCYCLE_BODY_V2(7,1);
        SOLVERXP_MCYCLE_BODY_V2(6,1);
        SOLVERXP_MCYCLE_BODY_V2(5,1);
        SOLVERXP_MCYCLE_BODY_V2(4,1);
        SOLVERXP_MCYCLE_BODY_V2(3,1);
        SOLVERXP_MCYCLE_BODY_V2(2,1);
        SOLVERXP_MCYCLE_BODY_V2(1,1);
        SOLVERXP_MCYCLE_BODY_V2(0,1);
        pX += X_last_step/sizeof(*pX);
        pB += X_last_step/sizeof(*pB);
    } //for(k=CURRENT_P-1; k>=0; k--)
}


inline_ void SolveRxP_M16_P16x8_last8(
                               int16_t *X,
                               const int16_t *R,
                               int16_t *QB,
                               int L,
                               int start_m, 
                               int16_t *pe_end,
                               int16_t *pf_end

                               )
{
    const int end_m = start_m-7;

    xb_vecNx16 * pR0,  *pR00;
    xb_vecNx16 * pX;
    xb_vecNx16 * pB;

    int k;
    const int R_last_step0 = -CURRENT_M*L*2*sizeof(int16_t);
    const int R_last_step_inc =L*2*sizeof(int16_t) ;
    const int X_last_step  = ((start_m-end_m+1) * CURRENT_P-1)*L*2 *sizeof(int16_t);
    int R_last_step;
    pR00 = (xb_vecNx16 *) (R + (start_m-1)*CURRENT_M*L*2/*row = start_m-1*/ +start_m*L*2 /*column = start_m*/);
    pX = (xb_vecNx16 *)(X  +start_m*CURRENT_P*L*2+(CURRENT_P-1)*L*2);
    pB = (xb_vecNx16 *)(QB +start_m*CURRENT_P*L*2+(CURRENT_P-1)*L*2);

   
    for(k=CURRENT_P-1; k>=0; k--)
    {
        xb_vecNx16  *pe, *pf;
        xb_vecNx16 dp;
        xb_vecNx16 r0, ex, f0, b, x;
        xb_vecNx40 acc0;
        xb_vecNx16 Xreg[8];

        pR0 = pR00;
        pe = (xb_vecNx16  *)pe_end;
        pf = (xb_vecNx16  *)pf_end;
        R_last_step = R_last_step0;

        SOLVERXP_MCYCLE_BODY_V2(7,0);
        SOLVERXP_MCYCLE_BODY_V2(6,0);
        SOLVERXP_MCYCLE_BODY_V2(5,0);
        SOLVERXP_MCYCLE_BODY_V2(4,0);
        SOLVERXP_MCYCLE_BODY_V2(3,0);
        SOLVERXP_MCYCLE_BODY_V2(2,0);
        SOLVERXP_MCYCLE_BODY_V2(1,0);
        SOLVERXP_MCYCLE_BODY_V2(0,0);

        pX += X_last_step/sizeof(*pX);
        pB += X_last_step/sizeof(*pB);
    } //for(k=CURRENT_P-1; k>=0; k--)
}

 void UpdateB_M16_P16x8(int16_t *X, 
                              const int16_t *R,
                              int16_t *QB,
                              int _M, 
                              int _P,
                              int L,
                              int q,
                              int _start_m, 
                              int _end_m)
{
    vsaN _q=BBE_MOVVSA32(q+1);
    const int start_m = 7; 
    const int end_m = 0; 
    int k;
    int m;
    xb_vecNx40 acc0;
    xb_vecNx16 b; 

    xb_vecNx16 * restrict pX;
    xb_vecNx16 * restrict pR;
    xb_vecNx16 * restrict pB;

    pX = (xb_vecNx16 *) (X +(start_m+1)*CURRENT_P*L*2 + (CURRENT_P-1)*L*2 +L*2*CURRENT_P*(CURRENT_M-start_m-2)); 
    pR = (xb_vecNx16 *) (R + start_m*CURRENT_M*L*2 +(start_m+1)*L*2 +L*2*(CURRENT_M-1-start_m-1));
    pB = (xb_vecNx16 *)(QB+start_m*CURRENT_P*L*2+(CURRENT_P-1)*L*2);

    for(k=CURRENT_P-1; k>=0; k--)
    {
        xb_vecNx16 dp;
        xb_vecNx16 r0, x0;
        xb_vecNx16 x7, x6, x5, x4, x3, x2, x1 ;

        pR = (xb_vecNx16 *) (R + start_m*CURRENT_M*L*2 +(start_m+1)*L*2 +L*2*(CURRENT_M-1-start_m-1));
//        pe = pe_end; 
  //      pf = pf_end;

        x7=BBE_LVNX16_X( pX, -L*2*CURRENT_P*0*sizeof(int16_t)); 
        x6=BBE_LVNX16_X( pX, -L*2*CURRENT_P*1*sizeof(int16_t)); 
        x5=BBE_LVNX16_X( pX, -L*2*CURRENT_P*2*sizeof(int16_t)); 
        x4=BBE_LVNX16_X( pX, -L*2*CURRENT_P*3*sizeof(int16_t)); 
        x3=BBE_LVNX16_X( pX, -L*2*CURRENT_P*4*sizeof(int16_t)); 
        x2=BBE_LVNX16_X( pX, -L*2*CURRENT_P*5*sizeof(int16_t)); 
        x1=BBE_LVNX16_X( pX, -L*2*CURRENT_P*6*sizeof(int16_t)); 
        x0=BBE_LVNX16_X( pX, -L*2*CURRENT_P*7*sizeof(int16_t)); 

        for(m=start_m; m>=end_m; m--)
        {
            BBE_LVNX16_XP(r0, pR, -L*2*sizeof(int16_t));
                acc0 = BBE_MULRNX16C(x7, r0, _q);
            BBE_LVNX16_XP(r0, pR, -L*2*sizeof(int16_t));
            BBE_MULANX16C(acc0, x6, r0); 
            BBE_LVNX16_XP(r0, pR, -L*2*sizeof(int16_t));
            BBE_MULANX16C(acc0, x5, r0); 
            BBE_LVNX16_XP(r0, pR, -L*2*sizeof(int16_t));
            BBE_MULANX16C(acc0, x4, r0); 
            BBE_LVNX16_XP(r0, pR, -L*2*sizeof(int16_t));
            BBE_MULANX16C(acc0, x3, r0); 
            BBE_LVNX16_XP(r0, pR, -L*2*sizeof(int16_t));
            BBE_MULANX16C(acc0, x2, r0); 
            BBE_LVNX16_XP(r0, pR, -L*2*sizeof(int16_t));
            BBE_MULANX16C(acc0, x1, r0); 
            BBE_LVNX16_XP(r0, pR, -L*2*sizeof(int16_t));
            BBE_MULANX16C(acc0, x0, r0); 
            pR +=  -(-(1+start_m) +CURRENT_M)*L*2*sizeof(int16_t)/sizeof(*pR);
            dp = BBE_PACKVNX40(acc0, _q);
            b = BBE_LVNX16_I( pB, 0); 
            b = BBE_SRAINX16(b, 1); 
            b = BBE_SUBNX16(b, dp);
            BBE_SVNX16_XP(b, pB, -CURRENT_P*L*2*sizeof(int16_t)); 

        }
        pB += (CURRENT_P*L*2*(start_m-end_m+1)-2*L)*sizeof(int16_t)/sizeof(*pB);
        pX -= L*2*sizeof(int16_t)/sizeof(*pX);

    } //for(k=CURRENT_P-1; k>=0; k--)
    __Pragma("no_reorder");
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

void cqr_bkw16x1s (void* pScr, complex_fract16* restrict _B, const complex_fract16* restrict _R, int L)
{
          int16_t* restrict B=(      int16_t*)_B;
    const int16_t* restrict R=(const int16_t*)_R;
    const int q = 10;
    int i; 
    int16_t *invDiagR = (int16_t *)pScr; 
    int16_t *exp = CURRENT_M*BBE_SIMD_WIDTH + (int16_t*)pScr; 
    int16_t *QB = B; 
    int16_t *X = QB; 

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/2)==0);

    for(i=0; i<L; i+=(BBE_SIMD_WIDTH/2))
    {
        int k;
        const xb_vecNx16 * px = (xb_vecNx16 *) R ;
        vsaN shft;
        xb_vecNx16 * pf = (xb_vecNx16 *) invDiagR ;
        xb_vecNx16 * pe = (xb_vecNx16 *) exp ;
        xb_vecNx16 x0, c0 = BBE_MOVVA16(30-4-q);
        xb_vecNx40 z0;
        xb_vecNx16 f0,nsa,z4,ex,z3;
        NASSERT_ALIGN(invDiagR,2*BBE_SIMD_WIDTH);
        NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
        NASSERT_ALIGN(exp,2*BBE_SIMD_WIDTH);

        z0 = BBE_MOVWA32(0x2000000);
        z4 = BBE_MOVVINT16(4);
        z3 = BBE_MOVVINT16(3);

        for (k=0; k<CURRENT_M; k++)
        {
            BBE_LVNX16_XP(x0, px, (CURRENT_M+1)*L*2*sizeof(int16_t));
            shft = BBE_NSANX16(x0);
            nsa = BBE_MOVVVS(shft);
            ex = BBE_SUBNX16(nsa,z3);
            nsa = BBE_SUBNX16(nsa,z4);
            shft = BBE_MOVVSV(nsa, 0);
            x0 = BBE_SLANX16(x0,shft);
            f0 = BBE_QUONX32(z0,x0);

            f0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DUPLICATE_1_EVEN);
            ex = BBE_SHFLNX16I(ex, BBE_SHFLI_DUPLICATE_1_EVEN);
            ex = BBE_SUBNX16(c0, ex); 
            BBE_SVNX16_IP(f0,pf,2*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(ex,pe,2*BBE_SIMD_WIDTH);
        }

        SolveRxP_M16_P16x8(X, R, QB,L, CURRENT_M-1, exp +(CURRENT_M-1)*BBE_SIMD_WIDTH, invDiagR +(CURRENT_M-1)*BBE_SIMD_WIDTH);
        UpdateB_M16_P16x8(X, R, QB, CURRENT_M, CURRENT_P, L, 10, CURRENT_M/2-1, 0); 
        SolveRxP_M16_P16x8_last8(X, R, QB,L, CURRENT_M/2-1, exp +(CURRENT_M/2-1)*BBE_SIMD_WIDTH, invDiagR +(CURRENT_M/2-1)*BBE_SIMD_WIDTH);

        X+=BBE_SIMD_WIDTH;
        R+=BBE_SIMD_WIDTH;
        QB+=BBE_SIMD_WIDTH;
    }
} /* cqr_bkw16x1s() */

size_t cqr_bkw16x1s_getScratchSize (int N, int P, int L)
{
    (void)P;
    return (L)*4*(N)*sizeof(int16_t);
} /* cqr_bkw16x1s_getScratchSize() */
#else
DISCARD_FUN(void,cqr_bkw16x1s ,(void* pScr, complex_fract16* restrict B,const complex_fract16* restrict R, int L))
size_t cqr_bkw16x1s_getScratchSize (int N, int P, int L) { (void)N;(void)P;(void)L; return 0; }
#endif
