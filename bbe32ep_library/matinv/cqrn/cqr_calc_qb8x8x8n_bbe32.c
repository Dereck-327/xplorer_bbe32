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

void cqr_calc_qb8x8x8n (void *pScr,
                          complex_fract16* _B,
                    const complex_fract16* _V,
                    int L)
{
        int16_t* B=(      int16_t*)_B;
  const int16_t* V=(const int16_t*)_V;
    int16_t* Z=(int16_t*)pScr;
    int m;
    int SB=2*8*8;
    const int16_t* pV;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B   ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V   ,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    /* scale down B by 1 bit right */
    cqrnScaleB(B,SB*L);

    /* apply Housholder vectors */
    for (pV=V,m=0; m<8; m++)
    {
        cqrnUpdateR8(Z,B+2*m*8,pV,SB,(8-m),8,8,L); 
        pV+=2*(8-m)*L;
    }

    /* Rotate matrix B by diagonal matrix Fi' */
    cqrnRotateR8(B,pV,L,1);
} /* cqr_calc_qb8x8x8n() */

size_t cqr_calc_qb8x8x8n_getScratchSize (int M, int N,int P,int L)
{
    NASSERT(M==8 && N==8 && P==8 && L>0);
    (void)N;(void)M;
    return 2*P*L*sizeof(int16_t);
} /* cqr_calc_qb8x8x8n_getScratchSize() */
#else
DISCARD_FUN(void,cqr_calc_qb8x8x8n,(void *pScr,
                          complex_fract16* B,
                    const complex_fract16* V,
                    int L))
size_t cqr_calc_qb8x8x8n_getScratchSize(int M, int N,int P,int L)
{
    (void)M;(void)N;(void)P;(void)L;
    return 0;
}
#endif
