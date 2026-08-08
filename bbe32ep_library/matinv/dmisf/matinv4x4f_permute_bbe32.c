/* ------------------------------------------------------------------------ */
/* Copyright (c) 2017 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
    NatureDSP_Baseband library. Direct Matrix Inversion
    Direct inversion of 4x4 floating point matrices 
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "matinv4x4Tbl.h"
#include "matinv4x4_common.h"
#if HAVE_VFPU
/*-------------------------------------------------------------------------
    inplace forward/backward permutation
    Input/output:
    X[16*L]  - matrices in the intermediate format (after conversion by
               cmatinv4x4sf_csb4x4_inplace)
    Input:
    permIx[L] - permutation indexes
    permTbl[36*32] permutation table (matinv4x4_bkw_perm_tbl_bbe32/
                   matinv4x4_fwd_perm_tbl_bbe32)
-------------------------------------------------------------------------*/
void matinv4x4f_permute(
                               float32_t * restrict Y,
                         const float32_t * restrict X,
                         const int16_t *permIx,
                         const int16_t *permTbl,
                         int L,
                         eLayout layout)
{
    const xb_vecN_2xf32* restrict pXrd;
          xb_vecN_2xf32* restrict pXwr;
    const xb_vecNx16   * restrict pTbl;
    int l,inc0,inc1,inc2;
    pXrd=(const xb_vecN_2xf32*)X;
    pXwr=(      xb_vecN_2xf32*)Y;
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    if(layout==e4x4_block)
    {
        inc0 =inc1=4*BBE_SIMD_WIDTH; inc2=2*BBE_SIMD_WIDTH; 
    }
    else
    {
        NASSERT(layout==e4x4_stream);
        inc0 =L*2*sizeof(float32_t); inc1=(-14*L+8)*sizeof(float32_t); inc2=L*sizeof(float32_t); 
    }
    for ( l=0; l<L; l++ )
    {
        xb_vecN_2xf32 X0,X1,A,B;
        xb_vecNx16 permA,permB;
        int addx,maxIx;
        maxIx=permIx[0];
        addx = inc0;
        XT_MOVEQZ(addx, inc1, (l&7)^7);
        pTbl=(const xb_vecNx16*)(permTbl+(maxIx<<5));
        BBE_LVNX16_IP(permA,pTbl,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(permB,pTbl,2*BBE_SIMD_WIDTH);
        X1=BBE_LVN_2XF32_X(pXrd,inc2);
        BBE_LVN_2XF32_XP(X0,pXrd,addx);
        A=BBE_SELN_2XF32(X1,X0,BBE_MOVVSELN_2NX16(permA,0));
        B=BBE_SELN_2XF32(X1,X0,BBE_MOVVSELN_2NX16(permB,0));
        BBE_SVN_2XF32_X(B,pXwr,inc2);
        BBE_SVN_2XF32_XP(A,pXwr,addx);
        permIx++;
    }
}

#endif
