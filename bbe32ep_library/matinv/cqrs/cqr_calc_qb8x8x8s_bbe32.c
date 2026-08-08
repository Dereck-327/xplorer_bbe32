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
    cqr_calc_qbMxNxPs/qr_calc_qbMxNxPs
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"

#define CURRENT_M 8
#define CURRENT_P 8

#define UPDATE(_M,_t,stride)           \
if(_M>_t)                              \
{                                      \
    BBE_LVNX16_IC(v0, pV0);            \
    BBE_LVNX16_XP(r0, pX0,stride);     \
    w0 = BBE_MULRNX16(r0,_0x4000,_14); \
    BBE_MULSNX16C(w0, v0, vR0);        \
    r0 = BBE_PACKVNX40(w0,_14);        \
    BBE_SVNX16_XP(r0, pY0, stride);    \
    BBE_LVNX16_XP(r1, pX1,stride);     \
    w0 = BBE_MULRNX16(r1,_0x4000,_14); \
    BBE_MULSNX16C(w0, v0, vR1);        \
    r1 = BBE_PACKVNX40(w0,_14);        \
    BBE_SVNX16_XP(r1, pY1, stride);    \
}


#define __UpdateMtxNonInplace8x8(X,Y,vR,V,L,P,_M,row_stride) \
{                                                            \
    int i;                                                   \
    vsaN _14=BBE_MOVVSA32(14);                               \
    xb_vecNx40 w0;                                           \
    xb_vecNx16 _0x4000=BBE_MOVPINT16(16);                  \
    xb_vecNx16 r0, v0, vR0, vR1, r1;                         \
    pX0 = (xb_vecNx16 *)(X);                                 \
    pY0 = (xb_vecNx16 *)(Y);                                 \
    pV0 = (xb_vecNx16 *)(V);                                 \
    pY1 = (xb_vecNx16 *)(Y+P/2*2*L);                         \
    pX1 = (xb_vecNx16 *)(X+P/2*2*L);                         \
    pVR  = (xb_vecNx16 *)vR;                                 \
    pVR1 = (xb_vecNx16 *)(vR+P/2*2*L);                       \
    NASSERT_ALIGN(vR,2*BBE_SIMD_WIDTH);                       \
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);                                      \
    NASSERT(L%(BBE_SIMD_WIDTH/2) == 0 && L > 0);                             \
    NASSERT(_M > 1 && P>=1);                                 \
    WUR_CBEGIN( (unsigned) V);                               \
    WUR_CEND  ( (unsigned) (V + _M*2*L));                    \
    __Pragma("loop_count factor=4, min=4")                   \
    for(i=0; i < ((P/2)*L)/(BBE_SIMD_WIDTH/2); i++)          \
        {                                                        \
        BBE_LVNX16_IP(vR0, pVR,  2*BBE_SIMD_WIDTH);          \
        BBE_LVNX16_IP(vR1, pVR1, 2*BBE_SIMD_WIDTH);          \
        UPDATE(_M,1,row_stride)                              \
        UPDATE(_M,2,row_stride)                              \
        UPDATE(_M,3,row_stride)                              \
        UPDATE(_M,4,row_stride)                              \
        UPDATE(_M,5,row_stride)                              \
        UPDATE(_M,6,row_stride)                              \
        UPDATE(_M,7,row_stride)                              \
        UPDATE(_M,0,-row_stride*(_M-1) + 2*BBE_SIMD_WIDTH)   \
        }                                                        \
    __Pragma("no_reorder");                                  \
}


/*
    Calculate product  vR = v' * R ,
    v - Householder's vectors, LxMx1,
    R - matrix MxPxL,
    vR -vector Px1xL,
    vR and R in a streaming format
    v - in a blocking format
    */
#if 0
#define __CalcVR8x8(vR,R,V,L,P,M,row_stride)     \
{                                                \
    vsaN _14=BBE_MOVVSA32(14);                   \
    int i;                                       \
    xb_vecNx16 r0, v0, vR0;                      \
    xb_vecNx40 acc0;                             \
    pR0 = (xb_vecNx16 *)R;                       \
    pV0 = (xb_vecNx16 *)V;                       \
    pVR = (xb_vecNx16 *)vR;                      \
    NASSERT_ALIGN(vR,2*BBE_SIMD_WIDTH);;                         \
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);;                          \
    NASSERT_ALIGN(vR,2*BBE_SIMD_WIDTH);;                         \
    NASSERT(L%(BBE_SIMD_WIDTH/2) == 0 && L > 0);                 \
    NASSERT(M > 1 && P>=1);                      \
    WUR_CBEGIN( (unsigned) V);                   \
    WUR_CEND  ( (unsigned) (V + M*2*L));         \
    __Pragma("loop_count factor=8, min=8")       \
    for(i=0; i < (P*L)/(BBE_SIMD_WIDTH/2); i++)  \
        {                                            \
                  BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);acc0 =BBE_MULRNX16J(r0, v0,_14); \
        if(M>2) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>3) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>4) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>5) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>6) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>7) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        BBE_LVNX16_XP(r0, pR0, -row_stride*(M-1) + 2*BBE_SIMD_WIDTH );                                        \
        BBE_LVNX16_IC(v0, pV0);                                                                               \
        BBE_MULANX16J(acc0, r0, v0);                                                                          \
        vR0 = BBE_PACKVNX40(acc0, _14);   /* Q13*/                                                              \
        BBE_SVNX16_IP(vR0, pVR, 2*BBE_SIMD_WIDTH);                                                            \
        }                                                                                                         \
}
#else
#define __CalcVR8x8(vR,R,V,L,P,M,row_stride)     \
{                                                \
    vsaN _14=BBE_MOVVSA32(14);                   \
    int i;                                       \
    xb_vecNx16 r0, v0, vR0;                      \
    xb_vecNx40 acc0;                             \
    pR0 = (xb_vecNx16 *)R;                       \
    pV0 = (xb_vecNx16 *)V;                       \
    pVR = (xb_vecNx16 *)vR;                      \
    NASSERT_ALIGN(vR,2*BBE_SIMD_WIDTH);;                         \
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);;                          \
    NASSERT_ALIGN(vR,2*BBE_SIMD_WIDTH);;                         \
    NASSERT(L%(BBE_SIMD_WIDTH/2) == 0 && L > 0);                 \
    NASSERT(M > 1 && P>=1);                      \
    WUR_CBEGIN( (unsigned) V);                   \
    WUR_CEND  ( (unsigned) (V + M*2*L));         \
    __Pragma("loop_count factor=4, min=4")       \
    for(i=0; i < (P*L)/(BBE_SIMD_WIDTH/2)/2; i++)  \
        {                                            \
                  BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);acc0 =BBE_MULRNX16J(r0, v0,_14); \
        if(M>2) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>3) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>4) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>5) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>6) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>7) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        BBE_LVNX16_XP(r0, pR0, -row_stride*(M-1) + 2*BBE_SIMD_WIDTH );                                        \
        BBE_LVNX16_IC(v0, pV0);                                                                               \
        BBE_MULANX16J(acc0, r0, v0);                                                                          \
        vR0 = BBE_PACKVNX40(acc0, _14);   /* Q13*/                                                              \
        BBE_SVNX16_IP(vR0, pVR, 2*BBE_SIMD_WIDTH);                                                            \
                                                                                                              \
                  BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);acc0 =BBE_MULRNX16J(r0, v0,_14); \
        if(M>2) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>3) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>4) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>5) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>6) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>7) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        BBE_LVNX16_XP(r0, pR0, -row_stride*(M-1) + 2*BBE_SIMD_WIDTH );                                        \
        BBE_LVNX16_IC(v0, pV0);                                                                               \
        BBE_MULANX16J(acc0, r0, v0);                                                                          \
        vR0 = BBE_PACKVNX40(acc0, _14);   /* Q13*/                                                              \
        BBE_SVNX16_IP(vR0, pVR, 2*BBE_SIMD_WIDTH);                                                            \
        }                                                                                                         \
    __Pragma("no_reorder");                                                                                     \
}
#endif
#if defined COMPILER_XTENSA
inline_ void RotateRows( int16_t * R, /*(io)*/
    const int16_t * Fi,/* (i)*/                          
    int L, 
    const int M, 
    const int P) __attribute__((always_inline));
#endif



inline_ void RotateRows(int16_t * R, /*(io)*/
    const int16_t * Fi,/* (i)*/
    int L,
    const int M,
    const int P)
{
    int i, j, l;
    xb_vecNx16* pFi0 = (xb_vecNx16 *)Fi;
    xb_vecNx16* pFi1 = (xb_vecNx16 *)(Fi + 1 * 2 * L*M / 4);
    xb_vecNx16* pFi2 = (xb_vecNx16 *)(Fi + 2 * 2 * L*M / 4);
    xb_vecNx16* pFi3 = (xb_vecNx16 *)(Fi + 3 * 2 * L*M / 4);

    xb_vecNx16* pR0 = (xb_vecNx16 *)R, *pR1, *pR2, *pR3;
    xb_vecNx16 r0, f0, f1, f2, f3, r1, r2, r3;

    NASSERT_ALIGN(Fi, 2 * BBE_SIMD_WIDTH);;
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);;
#if 0 
    for(i=0; i<M; i++)
    {
        for(l=0; l<L; l+=BBE_SIMD_WIDTH/2)
        {

            pR0  =  (xb_vecNx16 *)(R + i*2*L*P + 2*l);
            pFi =  (xb_vecNx16 *)(Fi + i*2*L + 2*l ); 

            f0 = BBE_LVNX16_I(pFi, 0);
            f0 = BBE_SLSINX16(f0, 1); 

            for(j=0; j<P; j++)
            {
                BBE_LVNX16_XP(r0, pR0, 0); 
                r0 = BBE_MULNX16JPACKQ(r0, f0); 
                BBE_SVNX16_XP(r0, pR0, 2*L*sizeof(int16_t)); 
            }
        }
    }
    return ;
#endif
    ASSERT(M % 4 == 0);

    pR0 = (xb_vecNx16 *)(R + 0 * 2 * L*P);
    pR1 = (xb_vecNx16 *)(R + (0 + 1 * M / 4) * 2 * L*P);
    pR2 = (xb_vecNx16 *)(R + (0 + 2 * M / 4) * 2 * L*P);
    pR3 = (xb_vecNx16 *)(R + (0 + 3 * M / 4) * 2 * L*P);

#ifdef COMPILER_XTENSA
#pragma ymemory(pR0)
#pragma ymemory(pR1)
#pragma ymemory(pR2)
#pragma ymemory(pR3)
#endif

    for (i = 0; i < M / 4; i++)
    {
        //  pR0  =  (xb_vecNx16 *)(R +         i*2*L*P);
        //  pR1  =  (xb_vecNx16 *)(R + (i+1*M/4)*2*L*P);
        //  pR2  =  (xb_vecNx16 *)(R + (i+2*M/4)*2*L*P);
        //  pR3  =  (xb_vecNx16 *)(R + (i+3*M/4)*2*L*P);

        for (l = 0; l < L; l += BBE_SIMD_WIDTH / 2)
        {

            BBE_LVNX16_IP(f0, pFi0, sizeof(*pFi0));
            BBE_LVNX16_IP(f1, pFi1, sizeof(*pFi0));
            BBE_LVNX16_IP(f2, pFi2, sizeof(*pFi0));
            BBE_LVNX16_IP(f3, pFi3, sizeof(*pFi0));

            f0 = BBE_SLSINX16(f0, 1);
            f1 = BBE_SLSINX16(f1, 1);
            f2 = BBE_SLSINX16(f2, 1);
            f3 = BBE_SLSINX16(f3, 1);

            for (j = 0; j < P - 1; j++)
            {
                BBE_LVNX16_XP(r0, pR0, 0);
                BBE_LVNX16_XP(r1, pR1, 0);
                BBE_LVNX16_XP(r2, pR2, 0);
                BBE_LVNX16_XP(r3, pR3, 0);

                r0 = BBE_MULNX16JPACKQ(r0, f0);
                r1 = BBE_MULNX16JPACKQ(r1, f1);
                r2 = BBE_MULNX16JPACKQ(r2, f2);
                r3 = BBE_MULNX16JPACKQ(r3, f3);

                BBE_SVNX16_XP(r0, pR0, 2 * L*sizeof(int16_t));
                BBE_SVNX16_XP(r1, pR1, 2 * L*sizeof(int16_t));
                BBE_SVNX16_XP(r2, pR2, 2 * L*sizeof(int16_t));
                BBE_SVNX16_XP(r3, pR3, 2 * L*sizeof(int16_t));
            }

            BBE_LVNX16_XP(r0, pR0, 0);
            BBE_LVNX16_XP(r1, pR1, 0);
            BBE_LVNX16_XP(r2, pR2, 0);
            BBE_LVNX16_XP(r3, pR3, 0);

            r0 = BBE_MULNX16JPACKQ(r0, f0);
            r1 = BBE_MULNX16JPACKQ(r1, f1);
            r2 = BBE_MULNX16JPACKQ(r2, f2);
            r3 = BBE_MULNX16JPACKQ(r3, f3);

            BBE_SVNX16_XP(r0, pR0, sizeof(*pR0) - (P - 1) * 2 * L*sizeof(int16_t));
            BBE_SVNX16_XP(r1, pR1, sizeof(*pR0) - (P - 1) * 2 * L*sizeof(int16_t));
            BBE_SVNX16_XP(r2, pR2, sizeof(*pR0) - (P - 1) * 2 * L*sizeof(int16_t));
            BBE_SVNX16_XP(r3, pR3, sizeof(*pR0) - (P - 1) * 2 * L*sizeof(int16_t));
        }
        // Adjust pointers to begin of a next rows
        pR0 = (xb_vecNx16*)((P - 1) * 2 * L*sizeof(int16_t) + (uintptr_t)pR0);
        pR1 = (xb_vecNx16*)((P - 1) * 2 * L*sizeof(int16_t) + (uintptr_t)pR1);
        pR2 = (xb_vecNx16*)((P - 1) * 2 * L*sizeof(int16_t) + (uintptr_t)pR2);
        pR3 = (xb_vecNx16*)((P - 1) * 2 * L*sizeof(int16_t) + (uintptr_t)pR3);
    }
}

/*-------------------------------------------------------------------------
cqr_calc_qbMxNxPs/qr_calc_qbMxNxPs

These functions apply Householder reflections to L MxP matrices B in the
course of solving a set of complex-valued linear problems A*X=B through
the QR decomposition of matrices A: A*X=B, A=Q*R => Q*R*X=B => R*X=Q'*B.
Instead of direct multiplication of each matrix B by conjugate transpose
of the corresponding matrix Q, we use a set of vectors V to perform a
sequence of Householder  reflections (see QR decomposition routines
cqr_build_rMxN/qr_build_rMxN).

Fixed point representation of output matrices is the same as for input
matrices.

Data transform is performed in-place.

NOTE:
Data layout for matrices is selected as for other matrices written in a stream 
order. 

Input:
B[M*P]L]                   Matrices B (L matrices of size MxP)
V[((M*N+((N-1)*N)/2+M)*L]  L sets of Householder vectors
Output:
B[M*P][L]                  Matrices Q'*B (L matrices of size MxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,P,L)
4. M must be greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
---------------------------------------------------------------------------*/
void cqr_calc_qb8x8x8s (void *pScr,complex_fract16 *_B, const complex_fract16 *_V, int L)
{
          int16_t *B=(      int16_t *)_B;
    const int16_t *V=(const int16_t *)_V;
    xb_vecNx16 * restrict pR0;
    xb_vecNx16 * restrict pX0;
    xb_vecNx16 * restrict pY0;
    xb_vecNx16 * restrict pV0;
    xb_vecNx16 * restrict pY1;
    xb_vecNx16 * restrict pX1;
    xb_vecNx16 * restrict pVR;
    xb_vecNx16 * restrict pVR1;
    const int16_t *v = V;
    int16_t *startB = B;
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);;
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT(L % (BBE_SIMD_WIDTH / 2) == 0 && L > 0);

    __CalcVR8x8(((int16_t*)pScr), B, v, L, CURRENT_P, 8, (CURRENT_P * 4 * L));
    __UpdateMtxNonInplace8x8(B, B, ((int16_t*)pScr), v, L, CURRENT_P, 8, (CURRENT_P * 4 * L));
    //__Pragma("no_reorder");
    v += 2 * L * 8;
    B += CURRENT_P * 2 * L;

    __CalcVR8x8(((int16_t*)pScr), B, v, L, CURRENT_P, 7, (CURRENT_P * 4 * L));
    __UpdateMtxNonInplace8x8(B, B, ((int16_t*)pScr), v, L, CURRENT_P, 7, (CURRENT_P * 4 * L));
    //__Pragma("no_reorder");
    v += 2 * L * 7;
    B += CURRENT_P * 2 * L;

    __CalcVR8x8(((int16_t*)pScr), B, v, L, CURRENT_P, 6, (CURRENT_P * 4 * L));
    __UpdateMtxNonInplace8x8(B, B, ((int16_t*)pScr), v, L, CURRENT_P, 6, (CURRENT_P * 4 * L));
    //__Pragma("no_reorder");
    v += 2 * L * 6;
    B += CURRENT_P * 2 * L;

    __CalcVR8x8(((int16_t*)pScr), B, v, L, CURRENT_P, 5, (CURRENT_P * 4 * L));
    __UpdateMtxNonInplace8x8(B, B, ((int16_t*)pScr), v, L, CURRENT_P, 5, (CURRENT_P * 4 * L));
    //__Pragma("no_reorder");
    v += 2 * L * 5;
    B += CURRENT_P * 2 * L;

    __CalcVR8x8(((int16_t*)pScr), B, v, L, CURRENT_P, 4, (CURRENT_P * 4 * L));
    __UpdateMtxNonInplace8x8(B, B, ((int16_t*)pScr), v, L, CURRENT_P, 4, (CURRENT_P * 4 * L));
    //__Pragma("no_reorder");
    v += 2 * L * 4;
    B += CURRENT_P * 2 * L;

    __CalcVR8x8(((int16_t*)pScr), B, v, L, CURRENT_P, 3, (CURRENT_P * 4 * L));
    __UpdateMtxNonInplace8x8(B, B, ((int16_t*)pScr), v, L, CURRENT_P, 3, (CURRENT_P * 4 * L));
    //__Pragma("no_reorder");
    v += 2 * L * 3;
    B += CURRENT_P * 2 * L;

    __CalcVR8x8(((int16_t*)pScr), B, v, L, CURRENT_P, 2, (CURRENT_P * 4 * L));
    __UpdateMtxNonInplace8x8(B, B, ((int16_t*)pScr), v, L, CURRENT_P, 2, (CURRENT_P * 4 * L));
    // __Pragma("no_reorder");
    v += 2 * L * 2;
    B += CURRENT_P * 2 * L;

    B = startB;
    RotateRows(B, v, L, CURRENT_M, CURRENT_P);
} /* cqr_calc_qb8x8x8s() */

size_t cqr_calc_qb8x8x8s_getScratchSize (int M, int P, int L)
{
    (void)M; (void)P; (void)L;
    return M * 2 * L*sizeof(int16_t);
} /* cqr_calc_qb8x8x8s_getScratchSize() */
