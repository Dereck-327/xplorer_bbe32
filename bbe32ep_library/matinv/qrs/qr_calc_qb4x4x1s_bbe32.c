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
#define CURRENT_P 1

#define USE_RWPTR 0

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
void  qr_calc_qb4x4x1s (void *pScr, int16_t *B, const int16_t *V , int L)
{
    vsaN _14=BBE_MOVVSA32(14);
    int l;

    xb_vecNx16 *  restrict   _pv = (xb_vecNx16 *)V; 
    const xb_vecNx16 *  restrict  _pbr=(const xb_vecNx16 *)B; 
          xb_vecNx16 *  restrict  _pbw=(xb_vecNx16 *)B; 
    const int     _v_last_step = BBE_SIMD_WIDTH*sizeof(int16_t); 

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V,2*BBE_SIMD_WIDTH);
    NASSERT(L>0 && L%BBE_SIMD_WIDTH==0);
    (void)pScr;

    for( l=0; l<L; l+=BBE_SIMD_WIDTH )
    {
        xb_vecNx16 _0x4000=BBE_MOVPINT16(16);
        xb_vecNx16 v0, v1, v2, v3, vR0,vR1,vR2, ar0,ar1,ar2,ar3; 
        xb_vecNx40 acc0; 

        BBE_LVNX16_IP(v0, _pv,  sizeof(*_pv)); 
        BBE_LVNX16_IP(v1, _pv,  sizeof(*_pv)); 
        BBE_LVNX16_IP(v2, _pv,  sizeof(*_pv)); 
        BBE_LVNX16_IP(v3, _pv,  sizeof(*_pv)); 

        BBE_LVNX16_XP(ar0, _pbr, 2*L); 
        BBE_LVNX16_XP(ar1, _pbr, 2*L); 
        BBE_LVNX16_XP(ar2, _pbr, 2*L); 
        BBE_LVNX16_XP(ar3, _pbr, 2*BBE_SIMD_WIDTH-3*2*L); 

        acc0 = BBE_MULRNX16(ar0, v0,_14);
        BBE_MULANX16(acc0,  ar1, v1);
        BBE_MULANX16(acc0,  ar2, v2); 
        BBE_MULANX16(acc0,  ar3, v3); 
        vR0 = BBE_PACKVNX40(acc0,_14);
        acc0=BBE_MULRNX16(ar0,_0x4000,_14); BBE_MULSNX16(acc0,vR0,v0); ar0=BBE_PACKVNX40(acc0,_14);
        acc0=BBE_MULRNX16(ar1,_0x4000,_14); BBE_MULSNX16(acc0,vR0,v1); ar1=BBE_PACKVNX40(acc0,_14);
        acc0=BBE_MULRNX16(ar2,_0x4000,_14); BBE_MULSNX16(acc0,vR0,v2); ar2=BBE_PACKVNX40(acc0,_14);
        acc0=BBE_MULRNX16(ar3,_0x4000,_14); BBE_MULSNX16(acc0,vR0,v3); ar3=BBE_PACKVNX40(acc0,_14);

        BBE_LVNX16_IP(v1, _pv,  sizeof(*_pv)); 
        BBE_LVNX16_IP(v2, _pv,  sizeof(*_pv)); 
        BBE_LVNX16_XP(v3, _pv,  sizeof(*_pv)); 
        acc0 = BBE_MULNX16( ar1, v1);
        BBE_MULANX16(acc0,  ar2, v2); 
        BBE_MULANX16(acc0,  ar3, v3); 
        vR1 = BBE_PACKVNX40(acc0,_14);
        acc0=BBE_MULRNX16(ar1,_0x4000,_14); BBE_MULSNX16(acc0,vR1,v1); ar1=BBE_PACKVNX40(acc0,_14);
        acc0=BBE_MULRNX16(ar2,_0x4000,_14); BBE_MULSNX16(acc0,vR1,v2); ar2=BBE_PACKVNX40(acc0,_14);
        acc0=BBE_MULRNX16(ar3,_0x4000,_14); BBE_MULSNX16(acc0,vR1,v3); ar3=BBE_PACKVNX40(acc0,_14);

        BBE_LVNX16_IP(v2, _pv,  sizeof(*_pv)); 
        BBE_LVNX16_XP(v3, _pv,  _v_last_step); 
        acc0 = BBE_MULNX16(  ar2, v2); 
        BBE_MULANX16(acc0,  ar3, v3); 
        vR2 = BBE_PACKVNX40(acc0,_14);
        acc0=BBE_MULRNX16(ar2,_0x4000,_14); BBE_MULSNX16(acc0,vR2,v2); ar2=BBE_PACKVNX40(acc0,_14);
        acc0=BBE_MULRNX16(ar3,_0x4000,_14); BBE_MULSNX16(acc0,vR2,v3); ar3=BBE_PACKVNX40(acc0,_14);

        BBE_SVNX16_XP(ar0, _pbw, 2*L);
        BBE_SVNX16_XP(ar1, _pbw, 2*L);
        BBE_SVNX16_XP(ar2, _pbw, 2*L);
        BBE_SVNX16_XP(ar3, _pbw, 2*BBE_SIMD_WIDTH-3*2*L);        
    }

} /* qr_calc_qb4x4x1s() */

size_t  qr_calc_qb4x4x1s_getScratchSize (int M, int P, int L)
{
    (void)M;
    (void)P;
    (void)L;
    return 0;
} /* qr_calc_qb4x4x1s_getScratchSize() */
