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
/* Common utility declarations. */
#include "cqr_common.h"

#define CURRENT_M 16
#define CURRENT_P 16

static int LeftRotateBy_conj_FIx8_16x16(
                                      int16_t * restrict R, /*(io)*/
                                      int16_t * restrict Fi,/* (i)*/
                                      int L)
{
    int i; 
    xb_vecNx16* restrict pFi = (xb_vecNx16 *)Fi; 
    xb_vecNx16* restrict pRrow; 
    xb_vecNx16* restrict pRrow_wr; 
    vsaN _14=BBE_MOVVSA32(14);

    NASSERT_ALIGN(Fi,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);

    for(i=0; i<CURRENT_M; i++)
    {
        xb_vecNx16  fi, r;
        xb_vecNx40 acc0;
        pRrow  =  (xb_vecNx16 *)(R + i*2*L*CURRENT_M);
        pRrow_wr  =  (xb_vecNx16 *)(R + i*2*L*CURRENT_M);

        BBE_LVNX16_IP( fi, pFi, sizeof(*pFi) );

        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
        BBE_LVNX16_XP( r, pRrow, 4*L );acc0 = BBE_MULRNX16J(r, fi,_14); r = BBE_PACKVNX40(acc0 ,_14); BBE_SVNX16_XP( r,  pRrow_wr, 4*L );
    }
    return 0; 
}

#define UPDATE_QB_16x16(R__,v__,m__,reload)                          \
{                                                                    \
    int  j;                                                          \
    const int v_back_step = -(CURRENT_M-m__-1)*2*BBE_SIMD_WIDTH;     \
    const int R_last_step = -CURRENT_P*L*4*(CURRENT_M-m__-1) +  L*4; \
    vsaN _14=BBE_MOVVSA32(14);                                       \
    xb_vecNx16 _0x4000=BBE_MOVPINT16(16);                            \
    xb_vecNx40 acc0;                                                 \
    xb_vecNx16  v0, vR0;                                             \
    xb_vecNx16  ar[16];                                              \
    NASSERT_ALIGN(R__,2*BBE_SIMD_WIDTH);                             \
    NASSERT_ALIGN(v__,2*BBE_SIMD_WIDTH);                             \
    NASSERT( (CURRENT_M-m__)>=2);                                    \
    pRrd__ = (xb_vecNx16 *)(R__ + CURRENT_P*L*2*m__ );               \
    pRwr__ = (xb_vecNx16 *)(R__ + CURRENT_P*L*2*m__ );               \
    pv__ = (xb_vecNx16 *)v__;                                        \
    for(j=0; j<CURRENT_P; j++)                                       \
    {                                                                \
        pRrd0  = pRrd__;                                             \
                      BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[ 0], pRrd__, CURRENT_P*L*4); acc0 =BBE_MULRNX16J(ar[ 0], v0,_14);\
        if( m__<14) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[ 1], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[ 1], v0); }\
        if( m__<13) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[ 2], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[ 2], v0); }\
        if( m__<12) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[ 3], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[ 3], v0); }\
        if( m__<11) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[ 4], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[ 4], v0); }\
        if( m__<10) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[ 5], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[ 5], v0); }\
        if( m__< 9) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[ 6], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[ 6], v0); }\
        if( m__< 8) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[ 7], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[ 7], v0); }\
        if( m__< 7) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[ 8], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[ 8], v0); }\
        if( m__< 6) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[ 9], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[ 9], v0); }\
        if( m__< 5) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[10], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[10], v0); }\
        if( m__< 4) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[11], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[11], v0); }\
        if( m__< 3) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[12], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[12], v0); }\
        if( m__< 2) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[13], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[13], v0); }\
        if( m__< 1) { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); BBE_LVNX16_XP(ar[14], pRrd__, CURRENT_P*L*4); BBE_MULANX16J(acc0, ar[14], v0); }\
                      BBE_LVNX16_XP(v0, pv__, v_back_step);      BBE_LVNX16_XP(ar[15], pRrd__, R_last_step);   BBE_MULANX16J(acc0, ar[15], v0); \
        vR0 = BBE_PACKVNX40(acc0, _14);\
        if(reload>0) {            BBE_LVNX16_XP(ar[ 0], pRrd0, CURRENT_P*L*4); } \
                  { BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 0],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[ 0]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[ 0], pRwr__, CURRENT_P*L*4);}\
        if(reload>1) { if(m__<14) BBE_LVNX16_XP(ar[ 1], pRrd0, CURRENT_P*L*4); } \
        if(m__<14){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 1],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[ 1]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[ 1], pRwr__, CURRENT_P*L*4);}\
        if(reload>2) { if(m__<13) BBE_LVNX16_XP(ar[ 2], pRrd0, CURRENT_P*L*4); } \
        if(m__<13){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 2],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[ 2]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[ 2], pRwr__, CURRENT_P*L*4);}\
        if(reload>3) { if(m__<12) BBE_LVNX16_XP(ar[ 3], pRrd0, CURRENT_P*L*4); } \
        if(m__<12){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 3],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[ 3]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[ 3], pRwr__, CURRENT_P*L*4);}\
        if(reload>4) { if(m__<11) BBE_LVNX16_XP(ar[ 4], pRrd0, CURRENT_P*L*4); } \
        if(m__<11){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 4],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[ 4]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[ 4], pRwr__, CURRENT_P*L*4);}\
        if(reload>5) { if(m__<10) BBE_LVNX16_XP(ar[ 5], pRrd0, CURRENT_P*L*4); } \
        if(m__<10){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 5],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[ 5]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[ 5], pRwr__, CURRENT_P*L*4);}\
        if(reload>6) { if(m__< 9) BBE_LVNX16_XP(ar[ 6], pRrd0, CURRENT_P*L*4); } \
        if(m__< 9){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 6],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[ 6]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[ 6], pRwr__, CURRENT_P*L*4);}\
        if(reload>7) { if(m__< 8) BBE_LVNX16_XP(ar[ 7], pRrd0, CURRENT_P*L*4); } \
        if(m__< 8){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 7],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[ 7]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[ 7], pRwr__, CURRENT_P*L*4);}\
        if(m__< 7){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 8],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[ 8]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[ 8], pRwr__, CURRENT_P*L*4);}\
        if(m__< 6){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 9],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[ 9]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[ 9], pRwr__, CURRENT_P*L*4);}\
        if(m__< 5){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[10],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[10]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[10], pRwr__, CURRENT_P*L*4);}\
        if(m__< 4){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[11],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[11]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[11], pRwr__, CURRENT_P*L*4);}\
        if(m__< 3){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[12],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[12]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[12], pRwr__, CURRENT_P*L*4);}\
        if(m__< 2){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[13],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[13]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[13], pRwr__, CURRENT_P*L*4);}\
        if(m__< 1){ BBE_LVNX16_IP(v0, pv__, 2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[14],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[14]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[14], pRwr__, CURRENT_P*L*4);}\
                    BBE_LVNX16_XP(v0, pv__, v_back_step);      acc0=BBE_MULRNX16(ar[15],_0x4000,_14); BBE_MULSNX16C(acc0,vR0, v0); ar[15]=BBE_PACKVNX40(acc0,_14);BBE_SVNX16_XP(ar[15],  pRwr__,  R_last_step );\
    }\
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
void cqr_calc_qb16x16x16s (void *pScr,complex_fract16 *_B, const complex_fract16 *_V, int L)
{
          int16_t *B=(      int16_t *)_B;
    const int16_t *V=(const int16_t *)_V;
    xb_vecNx16 * restrict pv__;
    xb_vecNx16 * restrict pRrd__;
    xb_vecNx16 * restrict pRrd0;
    xb_vecNx16 * restrict pRwr__;
    int l;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V,2*BBE_SIMD_WIDTH);
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/2)==0);

    for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2))
    {
        const int16_t *v = V; 
        UPDATE_QB_16x16(B, v,  0,7);    v += BBE_SIMD_WIDTH*(CURRENT_M-  0); 
        UPDATE_QB_16x16(B, v,  1,7);    v += BBE_SIMD_WIDTH*(CURRENT_M-  1); 
        UPDATE_QB_16x16(B, v,  2,7);    v += BBE_SIMD_WIDTH*(CURRENT_M-  2); 
        UPDATE_QB_16x16(B, v,  3,6);    v += BBE_SIMD_WIDTH*(CURRENT_M-  3); 
        UPDATE_QB_16x16(B, v,  4,2);    v += BBE_SIMD_WIDTH*(CURRENT_M-  4); 
        UPDATE_QB_16x16(B, v,  5,1);    v += BBE_SIMD_WIDTH*(CURRENT_M-  5); 
        UPDATE_QB_16x16(B, v,  6,1);    v += BBE_SIMD_WIDTH*(CURRENT_M-  6); 
        UPDATE_QB_16x16(B, v,  7,0);    v += BBE_SIMD_WIDTH*(CURRENT_M-  7); 
        UPDATE_QB_16x16(B, v,  8,0);    v += BBE_SIMD_WIDTH*(CURRENT_M-  8); 
        UPDATE_QB_16x16(B, v,  9,0);    v += BBE_SIMD_WIDTH*(CURRENT_M-  9); 
        UPDATE_QB_16x16(B, v, 10,0);    v += BBE_SIMD_WIDTH*(CURRENT_M- 10); 
        UPDATE_QB_16x16(B, v, 11,0);    v += BBE_SIMD_WIDTH*(CURRENT_M- 11); 
        UPDATE_QB_16x16(B, v, 12,0);    v += BBE_SIMD_WIDTH*(CURRENT_M- 12); 
        UPDATE_QB_16x16(B, v, 13,0);    v += BBE_SIMD_WIDTH*(CURRENT_M- 13); 
        UPDATE_QB_16x16(B, v, 14,0);    v += BBE_SIMD_WIDTH*(CURRENT_M- 14); 
        LeftRotateBy_conj_FIx8_16x16( B, (int16_t*)v,L); 
        B += BBE_SIMD_WIDTH;
        V += BBE_SIMD_WIDTH*(SIZE_OF_V(CURRENT_M, CURRENT_M)+SIZE_OF_FI(CURRENT_M)); 
    }

} /* cqr_calc_qb16x16x16s() */

size_t cqr_calc_qb16x16x16s_getScratchSize (int M, int P, int L)
{
    (void)M;(void)P;(void)L;
    return 0;
} /* cqr_calc_qb16x16x16s_getScratchSize() */
