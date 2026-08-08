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
  NatureDSP_Baseband library. Apply the QR decomposition to the matrix of normal equations system
    Update right side of equations for QR process for block ordered matrices.
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
#include "cqrn_common.h"

#if HAVE_CQRN

#define sz_i16 sizeof(int16_t)

/*-------------------------------------------------------------------------
Update right side of equations for QR process for block ordered matrices.
Matrix sizes SB,SV are selected as usual for complex block ordered matrix 
sequencies, i.e. total size is rounded up to the closest bigger multiple of 
BBE_SIMD_WIDTH/2==8 elements or, if it is lesser, to the closest bigger 
multiple of degree of 2. 
SB=size(M*P)
SV=size(((2*M-N+1)*N/2+N)*L)
Scratch size in bytes is defined by cqr_calc_qbmxnn_getScratchSize(M,N,P,L)
functions

Input:
 M, N, P      dimensional parameters
 L            Number of matrices
Input/output:
 B[L][SB]     On input it is the sequence of L complex matrices B. 
              At the end of the process, matrices Z replace input
              matrices A. In a case of non-square matrices (N!=M), 
              only N*P elements of each output matrix will be valid.
Input:
 V[SV]        Sequence of L Housholder rotation vectors 

Restrictions:
1. B, V, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. N <= M
---------------------------------------------------------------------------*/
void cqr_calc_qb8x8x1n (void *pScr,
                          complex_fract16* _B,
                    const complex_fract16* _V,
                    int L)
{
        int16_t* B=(      int16_t*)_B;
  const int16_t* V=(const int16_t*)_V;
/*
Reference code:
% compute Q'B matrix
% input:
% V  - sequence of Housholder vectors  [(2*M-N+1)*N/2,1]
% Fi - common rotation diagonal matrix [Nx1]
% R  - upper triangle decomposition
function [B] = cqr_calcQB(B,V,Fi)
[M, P] = size(B); 
[N, t] = size(Fi);
Z=zeros(M,P);
im=1;
for m=1:N
    v=V(im:im+M-m);
    im=im+(M-m+1);
    Bm=B(m:end,:);
    Bm=(Bm-2*v*v'*Bm);
    B(m:end,:)=Bm;
end
B=diag([Fi;ones(M-N,1)])'*B;
*/
    const xb_vecNx16 * restrict B_rd;
          xb_vecNx16 * restrict B_wr;
    const xb_vecNx16 * restrict V_rd;
    xb_vecNx40 A;
    xb_vecNx16 _0;
    vsaN _14;
    const int SB=2*8;
    int m, l;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B   ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V   ,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);

    /* scale down B by 1 bit right */
    cqrnScaleB(B,SB*L);

    /* apply Housholder vectors */
    V_rd = (const xb_vecNx16*)V;
    _0=0;
    _14=BBE_MOVVSA32(14);
    for (m=0; m<8; m++)
    {
        vsaN _11=BBE_MOVVSA32(11);
        xb_vecNx16 _0x4000=0x4000;
        xb_vecNx16 V0,B0,Z0;
        xb_c40 Im;
        valign aV;

        vboolN mask;
        xb_vecNx16 shuffle;

        Z0=BBE_SEQNX16(); B0=BBE_MOVVA16((m)*2-16); Z0=BBE_SUBNX16(Z0,B0);
        shuffle=BBE_MOVVSELNX16(Z0,0);
        mask = BBE_LTRN(m*2-16);

        B_rd=(const xb_vecNx16 * )B;
        B_wr=(      xb_vecNx16 * )B;
        aV=BBE_LA_PP(V_rd);
        __Pragma("loop_count min=1"); 
        for (l=0; l<L; l++)
        {
            BBE_LAVNX16_XP(V0,aV,V_rd ,4*(8-m));
            V0 = BBE_SHFLNX16(V0,shuffle);       
            BBE_LVNX16F_IP(B0, B_rd, 2*BBE_SIMD_WIDTH, mask); 

            A=BBE_MULRNX16J(B0,V0,_11); Im=BBE_RADDNX40C(A); A=BBE_MOVNX40_FROMC40(Im);
            Z0=BBE_PACKVNX40(A,_14);
            Z0=BBE_REPNX16C(Z0,0);
            A=BBE_MULRNX16(B0,_0x4000,_14);
            BBE_MULSNX16C(A,Z0,V0);
            Z0=BBE_PACKVNX40(A,_14);
            BBE_SVNX16F_IP(Z0, B_wr, 2*BBE_SIMD_WIDTH, mask);
        }
    }

    /* Rotate matrix B by diagonal matrix Fi' */
    cqrnRotateB8(pScr,B,V+2*36*L,L);
} /* cqr_calc_qb8x8x1n() */

size_t cqr_calc_qb8x8x1n_getScratchSize (int M, int N,int P,int L)
{
    size_t Bsize;
    NASSERT(M==8 && N==8 && P==1 && L>0);
    (void)M;(void)P;(void)L;
    Bsize = N*2*sz_i16*BBE_SIMD_WIDTH/2;
    return Bsize;
} /* cqr_calc_qb8x8x1n_getScratchSize() */

#else
DISCARD_FUN(void,cqr_calc_qb8x8x1n,(void *pScr,
                          complex_fract16* B,
                    const complex_fract16* V,
                    int L))
size_t cqr_calc_qb8x8x1n_getScratchSize(int M, int N,int P,int L)
{
    (void)M;(void)N;(void)P;(void)L;
    return 0;
}
#endif
