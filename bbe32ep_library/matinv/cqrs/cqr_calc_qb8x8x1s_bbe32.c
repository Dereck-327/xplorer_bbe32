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

#define CURRENT_M 8
#define CURRENT_P 1

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
void cqr_calc_qb8x8x1s (void *pScr,complex_fract16 *_B, const complex_fract16 *_V, int L)
{
          int16_t *B=(      int16_t *)_B;
    const int16_t *V=(const int16_t *)_V;
    vsaN _14 = BBE_MOVVSA32(14);
    const int M = CURRENT_M;
    const int v_inc = 2 * BBE_SIMD_WIDTH - 4 * L;

    int16_t *start_B = B;
    const int16_t *start_V = V;
    int l;
    xb_vecNx16 *  restrict  pb0;
    xb_vecNx16 *  restrict  pb1;
    xb_vecNx16 *  restrict  pb2;
    xb_vecNx16 *  restrict  pb3;
    xb_vecNx16 *  restrict  pb4;
    xb_vecNx16 *  restrict  pb5;
    xb_vecNx16 *  restrict  pb6;
    xb_vecNx16 *  restrict  pb7;

    xb_vecNx16 *  restrict   _pv = (xb_vecNx16 *)V;
    xb_vecNx16 *  restrict  _pb12 = (xb_vecNx16 *)(B + (CURRENT_M - 4) * 2 * L);
    xb_vecNx16 *  restrict  _pb13 = (xb_vecNx16 *)(B + (CURRENT_M - 3) * 2 * L);
    xb_vecNx16 *  restrict  _pb14 = (xb_vecNx16 *)(B + (CURRENT_M - 2) * 2 * L);
    xb_vecNx16 *  restrict  _pb15 = (xb_vecNx16 *)(B + (CURRENT_M - 1) * 2 * L);

    const int     _o_4 = -4 * 4 * L;   // offset 4 rows
    int     _v_step = 2 * BBE_SIMD_WIDTH;


#ifdef COMPILER_XTENSA
#pragma ymemory( _pb12 ) 
#pragma ymemory( _pb13 )
#pragma ymemory( _pb14 ) 
#pragma ymemory( _pb15 )
#endif
    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);;
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && L % (BBE_SIMD_WIDTH / 2) == 0);


    _v_step = -7 * 2 * BBE_SIMD_WIDTH + 4 * L*M;
    for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2), V += CURRENT_M*BBE_SIMD_WIDTH)
    {
        xb_vecNx16 v0, vR0, ar[16], _0x4000 = BBE_MOVVA16(0x4000);
        xb_vecNx40 acc0;
        _pv = (xb_vecNx16 *)(V);

        BBE_LVNX16_IP(v0, _pv, 0 * 2 * BBE_SIMD_WIDTH); ar[8] = BBE_LVNX16_X(_pb12, _o_4);  acc0 = BBE_MULRNX16J(ar[8], v0, _14);
        v0 = BBE_LVNX16_I(_pv, 1 * 2 * BBE_SIMD_WIDTH); ar[9] = BBE_LVNX16_X(_pb13, _o_4);  BBE_MULANX16J(acc0, ar[9], v0);
        v0 = BBE_LVNX16_I(_pv, 2 * 2 * BBE_SIMD_WIDTH); ar[10] = BBE_LVNX16_X(_pb14, _o_4);  BBE_MULANX16J(acc0, ar[10], v0);
        v0 = BBE_LVNX16_I(_pv, 3 * 2 * BBE_SIMD_WIDTH); ar[11] = BBE_LVNX16_X(_pb15, _o_4);  BBE_MULANX16J(acc0, ar[11], v0);
        v0 = BBE_LVNX16_I(_pv, 4 * 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[12], _pb12, 0);  BBE_MULANX16J(acc0, ar[12], v0);
        v0 = BBE_LVNX16_I(_pv, 5 * 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[13], _pb13, 0);  BBE_MULANX16J(acc0, ar[13], v0);
        v0 = BBE_LVNX16_I(_pv, 6 * 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[14], _pb14, 0);  BBE_MULANX16J(acc0, ar[14], v0);
        v0 = BBE_LVNX16_I(_pv, 7 * 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[15], _pb15, 0);  BBE_MULANX16J(acc0, ar[15], v0);
        vR0 = BBE_PACKVNX40(acc0, _14);

        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH);                                acc0 = BBE_MULRNX16(ar[8], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[8] = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_X(ar[8], _pb12, _o_4);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH);                                acc0 = BBE_MULRNX16(ar[9], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[9] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH);                                acc0 = BBE_MULRNX16(ar[10], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[10] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH);                                acc0 = BBE_MULRNX16(ar[11], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[11] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[12], _pb12, 0); acc0 = BBE_MULRNX16(ar[12], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[12] = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_IP(ar[12], _pb12, 0);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[13], _pb13, 0); acc0 = BBE_MULRNX16(ar[13], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[13] = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_IP(ar[13], _pb13, 0);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[14], _pb14, 0); acc0 = BBE_MULRNX16(ar[14], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[14] = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_IP(ar[14], _pb14, 0);
        BBE_LVNX16_XP(v0, _pv, _v_step); BBE_LVNX16_IP(ar[15], _pb15, 0); acc0 = BBE_MULRNX16(ar[15], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[15] = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_IP(ar[15], _pb15, 0);

        BBE_LVNX16_IP(v0, _pv, 0 * 2 * BBE_SIMD_WIDTH);                                acc0 = BBE_MULRNX16J(ar[9], v0, _14);
        v0 = BBE_LVNX16_I(_pv, 1 * 2 * BBE_SIMD_WIDTH);                                BBE_MULANX16J(acc0, ar[10], v0);
        v0 = BBE_LVNX16_I(_pv, 2 * 2 * BBE_SIMD_WIDTH);                                BBE_MULANX16J(acc0, ar[11], v0);
        v0 = BBE_LVNX16_I(_pv, 3 * 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[12], _pb12, 0); BBE_MULANX16J(acc0, ar[12], v0);
        v0 = BBE_LVNX16_I(_pv, 4 * 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[13], _pb13, 0); BBE_MULANX16J(acc0, ar[13], v0);
        v0 = BBE_LVNX16_I(_pv, 5 * 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[14], _pb14, 0); BBE_MULANX16J(acc0, ar[14], v0);
        v0 = BBE_LVNX16_I(_pv, 6 * 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[15], _pb15, 0); BBE_MULANX16J(acc0, ar[15], v0);
        vR0 = BBE_PACKVNX40(acc0, _14);
        _v_step += v_inc;

        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[9], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[9] = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_X(ar[9], _pb13, _o_4);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[10], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[10] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[11], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[11] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[12], _pb12, 0); acc0 = BBE_MULRNX16(ar[12], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[12] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[13], _pb13, 0); acc0 = BBE_MULRNX16(ar[13], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[13] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); BBE_LVNX16_IP(ar[14], _pb14, 0); acc0 = BBE_MULRNX16(ar[14], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[14] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_XP(v0, _pv, _v_step); BBE_LVNX16_IP(ar[15], _pb15, 0); acc0 = BBE_MULRNX16(ar[15], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[15] = BBE_PACKVNX40(acc0, _14);

        BBE_LVNX16_IP(v0, _pv, 0 * 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16J(ar[10], v0, _14);
        v0 = BBE_LVNX16_I(_pv, 1 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[11], v0);
        v0 = BBE_LVNX16_I(_pv, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[12], v0);
        v0 = BBE_LVNX16_I(_pv, 3 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[13], v0);
        v0 = BBE_LVNX16_I(_pv, 4 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[14], v0);
        v0 = BBE_LVNX16_I(_pv, 5 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[15], v0);
        vR0 = BBE_PACKVNX40(acc0, _14);
        _v_step += v_inc;

        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[10], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[10] = BBE_PACKVNX40(acc0, _14);  BBE_SVNX16_X(ar[10], _pb14, _o_4);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[11], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[11] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[12], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[12] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[13], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[13] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[14], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[14] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_XP(v0, _pv, _v_step); acc0 = BBE_MULRNX16(ar[15], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[15] = BBE_PACKVNX40(acc0, _14);

        BBE_LVNX16_IP(v0, _pv, 0 * 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16J(ar[11], v0, _14);
        v0 = BBE_LVNX16_I(_pv, 1 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[12], v0);
        v0 = BBE_LVNX16_I(_pv, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[13], v0);
        v0 = BBE_LVNX16_I(_pv, 3 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[14], v0);
        v0 = BBE_LVNX16_I(_pv, 4 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[15], v0);
        vR0 = BBE_PACKVNX40(acc0, _14);
        _v_step += v_inc;

        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[11], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[11] = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_X(ar[11], _pb15, _o_4);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[12], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[12] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[13], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[13] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[14], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[14] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_XP(v0, _pv, _v_step); acc0 = BBE_MULRNX16(ar[15], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[15] = BBE_PACKVNX40(acc0, _14);

        BBE_LVNX16_IP(v0, _pv, 0 * 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16J(ar[12], v0, _14);
        v0 = BBE_LVNX16_I(_pv, 1 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[13], v0);
        v0 = BBE_LVNX16_I(_pv, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[14], v0);
        v0 = BBE_LVNX16_I(_pv, 3 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[15], v0);
        vR0 = BBE_PACKVNX40(acc0, _14);
        _v_step += v_inc;

        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[12], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[12] = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_IP(ar[12], _pb12, (2 * BBE_SIMD_WIDTH));
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[13], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[13] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[14], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[14] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_XP(v0, _pv, _v_step); acc0 = BBE_MULRNX16(ar[15], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[15] = BBE_PACKVNX40(acc0, _14);

        BBE_LVNX16_IP(v0, _pv, 0 * 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16J(ar[13], v0, _14);
        v0 = BBE_LVNX16_I(_pv, 1 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[14], v0);
        v0 = BBE_LVNX16_I(_pv, 2 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[15], v0);
        vR0 = BBE_PACKVNX40(acc0, _14);
        _v_step += v_inc;

        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[13], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[13] = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_IP(ar[13], _pb13, (2 * BBE_SIMD_WIDTH));
        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[14], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[14] = BBE_PACKVNX40(acc0, _14);
        BBE_LVNX16_XP(v0, _pv, _v_step); acc0 = BBE_MULRNX16(ar[15], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[15] = BBE_PACKVNX40(acc0, _14);

        v0 = BBE_LVNX16_I(_pv, 0 * 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16J(ar[14], v0, _14);
        v0 = BBE_LVNX16_I(_pv, 1 * 2 * BBE_SIMD_WIDTH); BBE_MULANX16J(acc0, ar[15], v0);
        vR0 = BBE_PACKVNX40(acc0, _14);
        _v_step += v_inc;

        BBE_LVNX16_IP(v0, _pv, 2 * BBE_SIMD_WIDTH); acc0 = BBE_MULRNX16(ar[14], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[14] = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_IP(ar[14], _pb14, (2 * BBE_SIMD_WIDTH));
        BBE_LVNX16_XP(v0, _pv, _v_step); acc0 = BBE_MULRNX16(ar[15], _0x4000, _14); BBE_MULSNX16C(acc0, vR0, v0); ar[15] = BBE_PACKVNX40(acc0, _14); BBE_SVNX16_IP(ar[15], _pb15, (2 * BBE_SIMD_WIDTH));
        _v_step -= 6 * v_inc + 2 * BBE_SIMD_WIDTH;
    }
    V = start_V;
    // final rotation
    {
        pb0 = (xb_vecNx16 *)(start_B + 2 * L * 0);
        pb1 = (xb_vecNx16 *)(start_B + 2 * L * 1);
        pb2 = (xb_vecNx16 *)(start_B + 2 * L * 2);
        pb3 = (xb_vecNx16 *)(start_B + 2 * L * 3);
        pb4 = (xb_vecNx16 *)(start_B + 2 * L * 4);
        pb5 = (xb_vecNx16 *)(start_B + 2 * L * 5);
        pb6 = (xb_vecNx16 *)(start_B + 2 * L * 6);
        pb7 = (xb_vecNx16 *)(start_B + 2 * L * 7);
        _pv = (xb_vecNx16 *)(V + BBE_SIMD_WIDTH*(SIZE_OF_V(CURRENT_M, CURRENT_M)));


        for (l = 0; l < L; l += (BBE_SIMD_WIDTH / 2))
        {
            xb_vecNx40 acc0;
            xb_vecNx16 b0, b1, b2, b3, b4, b5, b6, b7;
            xb_vecNx16 v0, v1, v2, v3, v4, v5, v6, v7;
            _pv = (xb_vecNx16 *)(V + 2 * l + 2 * L*(M + M - 1 + M - 2 + M - 3 + M - 4 + M - 5 + M - 6));

            BBE_LVNX16_IP(b0, pb0, 0);
            BBE_LVNX16_IP(b1, pb1, 0);
            BBE_LVNX16_IP(b2, pb2, 0);
            BBE_LVNX16_IP(b3, pb3, 0);
            BBE_LVNX16_IP(b4, pb4, 0);
            BBE_LVNX16_IP(b5, pb5, 0);
            BBE_LVNX16_IP(b6, pb6, 0);
            BBE_LVNX16_IP(b7, pb7, 0);

            BBE_LVNX16_XP(v0, _pv, 2 * L*sizeof(int16_t));
            BBE_LVNX16_XP(v1, _pv, 2 * L*sizeof(int16_t));
            BBE_LVNX16_XP(v2, _pv, 2 * L*sizeof(int16_t));
            BBE_LVNX16_XP(v3, _pv, 2 * L*sizeof(int16_t));
            BBE_LVNX16_XP(v4, _pv, 2 * L*sizeof(int16_t));
            BBE_LVNX16_XP(v5, _pv, 2 * L*sizeof(int16_t));
            BBE_LVNX16_XP(v6, _pv, 2 * L*sizeof(int16_t));
            BBE_LVNX16_XP(v7, _pv, 2 * L*sizeof(int16_t));

            acc0 = BBE_MULRNX16J(b0, v0, _14);    v0 = BBE_PACKVNX40(acc0, _14);  BBE_SVNX16_IP(v0, pb0, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16J(b1, v1, _14);    v1 = BBE_PACKVNX40(acc0, _14);  BBE_SVNX16_IP(v1, pb1, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16J(b2, v2, _14);    v2 = BBE_PACKVNX40(acc0, _14);  BBE_SVNX16_IP(v2, pb2, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16J(b3, v3, _14);    v3 = BBE_PACKVNX40(acc0, _14);  BBE_SVNX16_IP(v3, pb3, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16J(b4, v4, _14);    v4 = BBE_PACKVNX40(acc0, _14);  BBE_SVNX16_IP(v4, pb4, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16J(b5, v5, _14);    v5 = BBE_PACKVNX40(acc0, _14);  BBE_SVNX16_IP(v5, pb5, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16J(b6, v6, _14);    v6 = BBE_PACKVNX40(acc0, _14);  BBE_SVNX16_IP(v6, pb6, 2 * BBE_SIMD_WIDTH);
            acc0 = BBE_MULRNX16J(b7, v7, _14);    v7 = BBE_PACKVNX40(acc0, _14);  BBE_SVNX16_IP(v7, pb7, 2 * BBE_SIMD_WIDTH);
        }
    }
} /* cqr_calc_qb8x8x1s() */

size_t cqr_calc_qb8x8x1s_getScratchSize (int M, int P, int L)
{
    (void)M; (void)P; (void)L;
    return 0;
} /* cqr_calc_qb8x8x1s_getScratchSize() */
