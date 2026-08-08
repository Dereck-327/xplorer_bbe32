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
    Direct inversion of 3x3 floating point matrices 
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "matinv3x3Tbl.h"
#include "matinv3x3_common.h"
#if HAVE_VFPU
#define __OPTIMIZED__ 1
/*-------------------------------------------------------------------------
    inplace forward/backward permutation
    Input/output:
    X[9*L]  - matrices in the intermediate format (after conversion by
              matinv3x3sf_rsb3x3_inplace)
    Input:
    permIx[L] - permutation indexes
    permTbl[8*32] permutation table (matinv3x3_bkw_perm_tbl_bbe32/
                  matinv3x3_fwd_perm_tbl_bbe32)
-------------------------------------------------------------------------*/
void matinv3x3f_permute(
                               float32_t * restrict Y,
                         const float32_t * restrict X,
                         const int16_t *permIx,
                         const int16_t *permTbl,
                         int L, eLayout layout)
#if !__OPTIMIZED__
{
    int n,k,j,l;
    NASSERT_ALIGN(X,BBE_SIMD_WIDTH*2);
    NASSERT(L>0);
    for ( l=0; l<L; l++ )
    {
        float32_t x[9],t[9];
        const int16_t* tbl;
       int ix;
       ix=permIx[l];
       tbl=&permTbl[32*ix];
       // load from intermediate format
       if (layout==e3x3_stream)
       {
           for (n=0; n<8 ; n++)
           {
               x[n ]= X[(l&~7)+n+L*(l&7)];
           }
           x[8]=X[8*L+l];
       }
       else
       {
           NASSERT(layout==e3x3_block);
           for (n=0; n<9 ; n++) x[n]=X[16*l+n];
       }

        for (k=0;k<9;k++) { j=tbl[(k*2)]>>1; t[k]=x[j]; }
        for (k=0;k<9;k++) { x[k]=t[k];  }

       // save back
       if (layout==e3x3_stream)
       {
           for (n=0; n<8 ; n++)
           {
               Y[(l&~7)+n+L*(l&7)]=x[n ];
           }
           Y[8*L+l]=x[8];
       }
       else
       {
           NASSERT(layout==e3x3_block);
           for (n=0; n<9 ; n++) Y[16*l+n]=x[n];
       }
    }
}
#else
{
    const xb_vecN_2xf32* restrict pXrd;
          xb_vecN_2xf32* restrict pXwr;
    const xtfloat* restrict pXrd8=NULL;
          xtfloat* restrict pXwr8;

    const xb_vecNx16   * restrict pTbl;
    int l,estride,mstride,estride8;

    if (layout==e3x3_block)
    {
        estride = 4*BBE_SIMD_WIDTH;
        mstride =16*2*BBE_SIMD_WIDTH;
        pXrd8=( xtfloat*)(X+8);
        pXwr8=( xtfloat*)(Y+8);
        estride8=4*BBE_SIMD_WIDTH;
    }
    else
    {
        NASSERT(layout==e3x3_stream);
        estride = L*sizeof(float32_t);
        mstride = 2*BBE_SIMD_WIDTH;
        pXrd8=(     xtfloat*)(X+8*L);
        pXwr8=(     xtfloat*)(Y+8*L);
        estride8=sizeof(float32_t);
    }

    pXrd=(const xb_vecN_2xf32*)X;
    pXwr=(      xb_vecN_2xf32*)Y;
    for ( l=0; l<L; l++ )
    {
        xb_vecN_2xf32 X0,X1,A,B;
        xb_vecNx16 permA,permB;
        int addx,maxIx;
        maxIx=permIx[0];
        addx = estride;
        XT_MOVEQZ(addx, -7*estride+mstride, (l&7)^7);
        pTbl=(const xb_vecNx16*)(permTbl+(maxIx<<5));
        BBE_LVNX16_IP(permA,pTbl,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(permB,pTbl,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_XP(X0,pXrd,addx);
        BBE_LSN_2XF32_XP(X1,pXrd8,estride8);
        A=BBE_SELN_2XF32(X1,X0,BBE_MOVVSELN_2NX16(permA,0));
        B=BBE_SELN_2XF32(X1,X0,BBE_MOVVSELN_2NX16(permB,0));
        BBE_SVN_2XF32_XP(A,pXwr,addx);
        BBE_SSN_2XF32_XP(B,pXwr8,estride8);
        permIx++;
    }
}
#endif

#endif
