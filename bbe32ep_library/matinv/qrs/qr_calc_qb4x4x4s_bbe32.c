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

#define CURRENT_M 4
#define CURRENT_P 4

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
void  qr_calc_qb4x4x4s (void *pScr, int16_t *B, const int16_t *V , int L)
{
#define MANUAL_UNROLL
#undef MANUAL_UNROLL
    int l;

    xb_vecNx16 *pR;
    xb_vecNx16 *pV;

    xb_vecNx40  acc_Nx40;
    xb_vecNx16 r_Nx16, vr_Nx16;
    vsaN  sh14 = BBE_MOVVSA32(14);
    xb_vecNx16 one_Q14 = BBE_MOVVA16(1 << 14);

    xb_vecNx16 ar[CURRENT_M];
    xb_vecNx16 av[CURRENT_M];

    pR = (xb_vecNx16*)(B);
    pV = (xb_vecNx16*)(V);

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT(L%BBE_SIMD_WIDTH == 0);
    (void)pScr;

    for (l = 0; l < L; l += BBE_SIMD_WIDTH)
    {
        {
#ifdef MANUAL_UNROLL
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 0*L*CURRENT_P + 0*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH*2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14); 
                BBE_LVNX16_XP(ar[1], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                BBE_LVNX16_XP(ar[2], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[2], pV, BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[2], av[2]);
                BBE_LVNX16_XP(ar[3], pR, -3*L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_XP(av[3], pV, -3*BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[3], av[3]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[2], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[2]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[3], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[3]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, 1*L*sizeof(int16_t) - 3*L*CURRENT_P*sizeof(int16_t)); 
            }
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 0*L*CURRENT_P + 1*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH*2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14); 
                BBE_LVNX16_XP(ar[1], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                BBE_LVNX16_XP(ar[2], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[2], pV, BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[2], av[2]);
                BBE_LVNX16_XP(ar[3], pR, -3*L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_XP(av[3], pV, -3*BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[3], av[3]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[2], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[2]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[3], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[3]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, 1*L*sizeof(int16_t) - 3*L*CURRENT_P*sizeof(int16_t));
            }
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 0*L*CURRENT_P + 2*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH*2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14); 
                BBE_LVNX16_XP(ar[1], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                BBE_LVNX16_XP(ar[2], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[2], pV, BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[2], av[2]);
                BBE_LVNX16_XP(ar[3], pR, -3*L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_XP(av[3], pV, -3*BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[3], av[3]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[2], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[2]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[3], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[3]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, 1*L*sizeof(int16_t) - 3*L*CURRENT_P*sizeof(int16_t));
            }
#else
            for (int n = 0; n < 3; n++)
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 0*L*CURRENT_P + 0*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH * 2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14);
                BBE_LVNX16_XP(ar[1], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, BBE_SIMD_WIDTH * 2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                BBE_LVNX16_XP(ar[2], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[2], pV, BBE_SIMD_WIDTH * 2); BBE_MULANX16(acc_Nx40, ar[2], av[2]);
                BBE_LVNX16_XP(ar[3], pR, -3 * L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_XP(av[3], pV, -3 * BBE_SIMD_WIDTH * 2); BBE_MULANX16(acc_Nx40, ar[3], av[3]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t));
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t));
                acc_Nx40 = BBE_MULRNX16(ar[2], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[2]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t));
                acc_Nx40 = BBE_MULRNX16(ar[3], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[3]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, 1 * L*sizeof(int16_t) - 3 * L*CURRENT_P*sizeof(int16_t));
            }
#endif // MANUAL_UNROLL
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 0*L*CURRENT_P + 3*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH * 2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14);
                BBE_LVNX16_XP(ar[1], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, BBE_SIMD_WIDTH * 2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                BBE_LVNX16_XP(ar[2], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[2], pV, BBE_SIMD_WIDTH * 2); BBE_MULANX16(acc_Nx40, ar[2], av[2]);
                BBE_LVNX16_XP(ar[3], pR, -3 * L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[3], pV, BBE_SIMD_WIDTH*(CURRENT_M - 0)*sizeof(int16_t) - 3 * BBE_SIMD_WIDTH*sizeof(int16_t)); BBE_MULANX16(acc_Nx40, ar[3], av[3]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t));
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t));
                acc_Nx40 = BBE_MULRNX16(ar[2], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[2]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t));
                acc_Nx40 = BBE_MULRNX16(ar[3], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[3]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, 1 * L*sizeof(int16_t) - 3 * L*CURRENT_P*sizeof(int16_t));
            }
            //v += BBE_SIMD_WIDTH*(CURRENT_M-0); 
        }
        {
#ifdef MANUAL_UNROLL
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 1*L*CURRENT_P + 0*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH*2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14); 
                BBE_LVNX16_XP(ar[1], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                BBE_LVNX16_XP(ar[2], pR, -2*L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[2], pV, -2*BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[2], av[2]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[2], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[2]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, 1*L*sizeof(int16_t) - 2*L*CURRENT_P*sizeof(int16_t)); 
            }
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 1*L*CURRENT_P + 1*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH*2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14); 
                BBE_LVNX16_XP(ar[1], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                BBE_LVNX16_XP(ar[2], pR, -2*L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[2], pV, -2*BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[2], av[2]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[2], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[2]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, 1*L*sizeof(int16_t) - 2*L*CURRENT_P*sizeof(int16_t)); 
            }
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 1*L*CURRENT_P + 2*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH*2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14); 
                BBE_LVNX16_XP(ar[1], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                BBE_LVNX16_XP(ar[2], pR, -2*L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[2], pV, -2*BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[2], av[2]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[2], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[2]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, 1*L*sizeof(int16_t) - 2*L*CURRENT_P*sizeof(int16_t)); 
            }
#else // MANUAL_UNROLL
            for (int n = 0; n < 3; n++)
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 1*L*CURRENT_P + 0*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH * 2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14);
                BBE_LVNX16_XP(ar[1], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, BBE_SIMD_WIDTH * 2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                BBE_LVNX16_XP(ar[2], pR, -2 * L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[2], pV, -2 * BBE_SIMD_WIDTH * 2); BBE_MULANX16(acc_Nx40, ar[2], av[2]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t));
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t));
                acc_Nx40 = BBE_MULRNX16(ar[2], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[2]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, 1 * L*sizeof(int16_t) - 2 * L*CURRENT_P*sizeof(int16_t));
            }
#endif
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 1*L*CURRENT_P + 3*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH * 2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14);
                BBE_LVNX16_XP(ar[1], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, BBE_SIMD_WIDTH * 2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                BBE_LVNX16_XP(ar[2], pR, -2 * L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[2], pV, BBE_SIMD_WIDTH*(CURRENT_M - 1)*sizeof(int16_t) - 2 * BBE_SIMD_WIDTH*sizeof(int16_t)); BBE_MULANX16(acc_Nx40, ar[2], av[2]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t));
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t));
                acc_Nx40 = BBE_MULRNX16(ar[2], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[2]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, 1 * L*sizeof(int16_t) - 2 * L*CURRENT_P*sizeof(int16_t));
            }
            //v += BBE_SIMD_WIDTH*(CURRENT_M-1); 
        }
        {
#ifdef MANUAL_UNROLL
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 2*L*CURRENT_P + 0*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH*2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14); 
                BBE_LVNX16_XP(ar[1], pR, -1*L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, -1*BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, 1*L*sizeof(int16_t) - 1*L*CURRENT_P*sizeof(int16_t)); 
            }
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 2*L*CURRENT_P + 1*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH*2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14); 
                BBE_LVNX16_XP(ar[1], pR, -1*L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, -1*BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, 1*L*sizeof(int16_t) - 1*L*CURRENT_P*sizeof(int16_t));
            }
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 2*L*CURRENT_P + 2*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH*2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14); 
                BBE_LVNX16_XP(ar[1], pR, -1*L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, -1*BBE_SIMD_WIDTH*2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t)); 
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14); 
                BBE_SVNX16_XP(r_Nx16, pR, 1*L*sizeof(int16_t) - 1*L*CURRENT_P*sizeof(int16_t));
            }
#else 
            for (int n = 0; n < 3; n++)
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 2*L*CURRENT_P + 1*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH * 2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14);
                BBE_LVNX16_XP(ar[1], pR, -1 * L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, -1 * BBE_SIMD_WIDTH * 2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t));
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, 1 * L*sizeof(int16_t) - 1 * L*CURRENT_P*sizeof(int16_t));
            }
#endif// MANUAL_UNROLL
            {
                //xb_vecNx16 *pR = (xb_vecNx16*)(B + 2*L*CURRENT_P + 3*L); 
                //xb_vecNx16 *pV = (xb_vecNx16*)(v);
                BBE_LVNX16_XP(ar[0], pR, L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[0], pV, BBE_SIMD_WIDTH * 2); acc_Nx40 = BBE_MULRNX16(ar[0], av[0], sh14);
                BBE_LVNX16_XP(ar[1], pR, -1 * L*CURRENT_P*sizeof(int16_t)); BBE_LVNX16_IP(av[1], pV, BBE_SIMD_WIDTH*(CURRENT_M - 2)*sizeof(int16_t) - 1 * BBE_SIMD_WIDTH * 2); BBE_MULANX16(acc_Nx40, ar[1], av[1]);
                vr_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);

                acc_Nx40 = BBE_MULRNX16(ar[0], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[0]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, L*CURRENT_P*sizeof(int16_t));
                acc_Nx40 = BBE_MULRNX16(ar[1], one_Q14, sh14); BBE_MULSNX16(acc_Nx40, vr_Nx16, av[1]); r_Nx16 = BBE_PACKVNX40(acc_Nx40, sh14);
                BBE_SVNX16_XP(r_Nx16, pR, BBE_SIMD_WIDTH*sizeof(int16_t) - 3 * L*sizeof(int16_t) - 2 * L*CURRENT_P*sizeof(int16_t) - L*CURRENT_P*sizeof(int16_t));
            }
            //v += BBE_SIMD_WIDTH*(CURRENT_M-2); 
        }
    }//for(l=0; l<L; l+=BBE_SIMD_WIDTH)

} /* qr_calc_qb4x4x4s() */

size_t  qr_calc_qb4x4x4s_getScratchSize (int M, int P, int L)
{
    (void)M;
    (void)P;
    (void)L;
    return 0;
} /* qr_calc_qb4x4x4s_getScratchSize() */
