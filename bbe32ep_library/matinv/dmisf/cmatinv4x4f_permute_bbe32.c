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
#include "matinv4x4_common.h"

#define __OPTIMIZED__ 1

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
void cmatinv4x4f_permute(
                               complex_float * restrict Y,
                         const complex_float * restrict X,
                         const int16_t *permIx,
                         const int16_t *permTbl,
                         int L,eLayout layout)
#if !__OPTIMIZED__
{
    int l;
    NASSERT_ALIGN(X,BBE_SIMD_WIDTH*2);
    NASSERT(L>0);
    for ( l=0; l<L; l++ )
    {
       int ix,j,k,n;
       const int16_t * tbl;
       ix=permIx[l];
       complex_float t[16],x[16];
       if (layout==e4x4_stream)
       {
           for (n=0; n<4 ; n++)
           {
               x[n+ 0]= X[(l&~3)+n+L*((l&3)*4+0)];
               x[n+ 4]= X[(l&~3)+n+L*((l&3)*4+1)];
               x[n+ 8]= X[(l&~3)+n+L*((l&3)*4+2)];
               x[n+12]= X[(l&~3)+n+L*((l&3)*4+3)];
           }
       }
       else
       {
           NASSERT(layout==e4x4_block);
           for (n=0; n<16; n++) x[n]=X[l*16+n];
       }
       tbl=&permTbl[32*ix];
       for (k=0;k<16;k++) { j=tbl[k*2]>>1; t[k]=x[j]; }
       for (k=0;k<16;k++) { x[k]=t[k];  }
       if (layout==e4x4_stream)
       {
           for (n=0; n<4 ; n++)
           {
               Y[(l&~3)+n+L*((l&3)*4+0)]= x[n+ 0];
               Y[(l&~3)+n+L*((l&3)*4+1)]= x[n+ 4];
               Y[(l&~3)+n+L*((l&3)*4+2)]= x[n+ 8];
               Y[(l&~3)+n+L*((l&3)*4+3)]= x[n+12];
           }
       }
       else
       {
           NASSERT(layout==e4x4_block);
           for (n=0; n<16; n++) Y[l*16+n]=x[n];
       }
    }
}
#elif 0
{
    const xb_vecN_2xf32* restrict pXrd;
          xb_vecN_2xf32* restrict pXwr;
    const xb_vecNx16   * restrict pTbl;
    int l,addx,estride,mstride;
    NASSERT_ALIGN(X,BBE_SIMD_WIDTH*2);
    NASSERT_ALIGN(permTbl,BBE_SIMD_WIDTH*2);
    NASSERT(L>0);

    switch(layout)
    {
    case e4x4_block:  estride=2*BBE_SIMD_WIDTH;  mstride=16*2*BBE_SIMD_WIDTH; break;
    case e4x4_stream: estride=L<<3;  mstride=2*BBE_SIMD_WIDTH; break;
    }
    pXrd=(const xb_vecN_2xf32*)X;
    pXwr=(      xb_vecN_2xf32*)Y;
    for ( l=0; l<L; l++ )
    {
        xb_vecN_2xf32 X0,X1,X2,X3,reX0,reX1,imX0,imX1,A,B,C,D;
        xb_vecNx16 permA,permB;
        int maxIx;

        addx = 4*estride;
        XT_MOVEQZ(addx, -12*estride+mstride, (l&3)^3);
        maxIx=permIx[0];permIx++;
        pTbl=(const xb_vecNx16*)(permTbl+(maxIx<<5));
        BBE_LVNX16_IP(permA,pTbl,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(permB,pTbl,2*BBE_SIMD_WIDTH);
        X1=BBE_LVN_2XF32_X (pXrd,1*estride);
        X2=BBE_LVN_2XF32_X (pXrd,2*estride);
        X3=BBE_LVN_2XF32_X (pXrd,3*estride);
        BBE_LVN_2XF32_XP(X0,pXrd,addx);

        BBE_DSELN_2XF32I(imX0,reX0,X1,X0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(imX1,reX1,X3,X2,BBE_DSELI_DEINTERLEAVE_2);

        A=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permA,0));
        B=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permB,0));
        C=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permA,0));
        D=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permB,0));

        BBE_DSELN_2XF32I(C,A,C,A,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(D,B,D,B,BBE_DSELI_INTERLEAVE_2);

        BBE_SVN_2XF32_X (C,pXwr,1*estride);
        BBE_SVN_2XF32_X (B,pXwr,2*estride);
        BBE_SVN_2XF32_X (D,pXwr,3*estride);
        BBE_SVN_2XF32_XP(A,pXwr,addx);
    }
}
#else
{
    const xb_vecN_2xf32* restrict pXrd;
          xb_vecN_2xf32* restrict pXwr;
    const xb_vecNx16   * restrict pTbl;
    int l,estride,mstride;
    NASSERT_ALIGN(X,BBE_SIMD_WIDTH*2);
    NASSERT_ALIGN(permTbl,BBE_SIMD_WIDTH*2);
    NASSERT(L>0);

    if(layout==e4x4_block)
    {
        estride=2*BBE_SIMD_WIDTH;  mstride=16*2*BBE_SIMD_WIDTH; 
    }
    else
    {
        NASSERT(layout==e4x4_stream);
        estride=L<<3;  mstride=2*BBE_SIMD_WIDTH; 
    }
    pXrd=(const xb_vecN_2xf32*)X;
    pXwr=(      xb_vecN_2xf32*)Y;
    for ( l=0; l<(L>>2); l++ )
    {
        xb_vecN_2xf32 X0,X1,X2,X3,reX0,reX1,imX0,imX1,A,B,C,D;
        xb_vecNx16 permA,permB;
        int maxIx;

        maxIx=permIx[0];permIx++;
        pTbl=(const xb_vecNx16*)(permTbl+(maxIx<<5));
        BBE_LVNX16_IP(permA,pTbl,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(permB,pTbl,2*BBE_SIMD_WIDTH);
        X1=BBE_LVN_2XF32_X (pXrd,1*estride);
        X2=BBE_LVN_2XF32_X (pXrd,2*estride);
        X3=BBE_LVN_2XF32_X (pXrd,3*estride);
        BBE_LVN_2XF32_XP(X0,pXrd,4*estride);

        BBE_DSELN_2XF32I(imX0,reX0,X1,X0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(imX1,reX1,X3,X2,BBE_DSELI_DEINTERLEAVE_2);

        A=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permA,0));
        B=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permB,0));
        C=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permA,0));
        D=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permB,0));

        BBE_DSELN_2XF32I(C,A,C,A,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(D,B,D,B,BBE_DSELI_INTERLEAVE_2);

        BBE_SVN_2XF32_X (C,pXwr,1*estride);
        BBE_SVN_2XF32_X (B,pXwr,2*estride);
        BBE_SVN_2XF32_X (D,pXwr,3*estride);
        BBE_SVN_2XF32_XP(A,pXwr,4*estride);

        maxIx=permIx[0];permIx++;
        pTbl=(const xb_vecNx16*)(permTbl+(maxIx<<5));
        BBE_LVNX16_IP(permA,pTbl,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(permB,pTbl,2*BBE_SIMD_WIDTH);
        X1=BBE_LVN_2XF32_X (pXrd,1*estride);
        X2=BBE_LVN_2XF32_X (pXrd,2*estride);
        X3=BBE_LVN_2XF32_X (pXrd,3*estride);
        BBE_LVN_2XF32_XP(X0,pXrd,4*estride);

        BBE_DSELN_2XF32I(imX0,reX0,X1,X0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(imX1,reX1,X3,X2,BBE_DSELI_DEINTERLEAVE_2);

        A=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permA,0));
        B=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permB,0));
        C=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permA,0));
        D=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permB,0));

        BBE_DSELN_2XF32I(C,A,C,A,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(D,B,D,B,BBE_DSELI_INTERLEAVE_2);

        BBE_SVN_2XF32_X (C,pXwr,1*estride);
        BBE_SVN_2XF32_X (B,pXwr,2*estride);
        BBE_SVN_2XF32_X (D,pXwr,3*estride);
        BBE_SVN_2XF32_XP(A,pXwr,4*estride);

        maxIx=permIx[0];permIx++;
        pTbl=(const xb_vecNx16*)(permTbl+(maxIx<<5));
        BBE_LVNX16_IP(permA,pTbl,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(permB,pTbl,2*BBE_SIMD_WIDTH);
        X1=BBE_LVN_2XF32_X (pXrd,1*estride);
        X2=BBE_LVN_2XF32_X (pXrd,2*estride);
        X3=BBE_LVN_2XF32_X (pXrd,3*estride);
        BBE_LVN_2XF32_XP(X0,pXrd,4*estride);

        BBE_DSELN_2XF32I(imX0,reX0,X1,X0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(imX1,reX1,X3,X2,BBE_DSELI_DEINTERLEAVE_2);

        A=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permA,0));
        B=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permB,0));
        C=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permA,0));
        D=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permB,0));

        BBE_DSELN_2XF32I(C,A,C,A,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(D,B,D,B,BBE_DSELI_INTERLEAVE_2);

        BBE_SVN_2XF32_X (C,pXwr,1*estride);
        BBE_SVN_2XF32_X (B,pXwr,2*estride);
        BBE_SVN_2XF32_X (D,pXwr,3*estride);
        BBE_SVN_2XF32_XP(A,pXwr,4*estride);

        maxIx=permIx[0];permIx++;
        pTbl=(const xb_vecNx16*)(permTbl+(maxIx<<5));
        BBE_LVNX16_IP(permA,pTbl,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(permB,pTbl,2*BBE_SIMD_WIDTH);
        X1=BBE_LVN_2XF32_X (pXrd,1*estride);
        X2=BBE_LVN_2XF32_X (pXrd,2*estride);
        X3=BBE_LVN_2XF32_X (pXrd,3*estride);
        BBE_LVN_2XF32_XP(X0,pXrd,-12*estride+mstride);

        BBE_DSELN_2XF32I(imX0,reX0,X1,X0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(imX1,reX1,X3,X2,BBE_DSELI_DEINTERLEAVE_2);

        A=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permA,0));
        B=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permB,0));
        C=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permA,0));
        D=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permB,0));

        BBE_DSELN_2XF32I(C,A,C,A,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(D,B,D,B,BBE_DSELI_INTERLEAVE_2);

        BBE_SVN_2XF32_X (C,pXwr,1*estride);
        BBE_SVN_2XF32_X (B,pXwr,2*estride);
        BBE_SVN_2XF32_X (D,pXwr,3*estride);
        BBE_SVN_2XF32_XP(A,pXwr,-12*estride+mstride);
    }

    for ( l=0; l<(L&3); l++ )
    {
        xb_vecN_2xf32 X0,X1,X2,X3,reX0,reX1,imX0,imX1,A,B,C,D;
        xb_vecNx16 permA,permB;
        int maxIx;

        maxIx=permIx[0];permIx++;
        pTbl=(const xb_vecNx16*)(permTbl+(maxIx<<5));
        BBE_LVNX16_IP(permA,pTbl,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(permB,pTbl,2*BBE_SIMD_WIDTH);
        X1=BBE_LVN_2XF32_X (pXrd,1*estride);
        X2=BBE_LVN_2XF32_X (pXrd,2*estride);
        X3=BBE_LVN_2XF32_X (pXrd,3*estride);
        BBE_LVN_2XF32_XP(X0,pXrd,4*estride);

        BBE_DSELN_2XF32I(imX0,reX0,X1,X0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELN_2XF32I(imX1,reX1,X3,X2,BBE_DSELI_DEINTERLEAVE_2);

        A=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permA,0));
        B=BBE_SELN_2XF32(reX1,reX0,BBE_MOVVSELN_2NX16(permB,0));
        C=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permA,0));
        D=BBE_SELN_2XF32(imX1,imX0,BBE_MOVVSELN_2NX16(permB,0));

        BBE_DSELN_2XF32I(C,A,C,A,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(D,B,D,B,BBE_DSELI_INTERLEAVE_2);

        BBE_SVN_2XF32_X (C,pXwr,1*estride);
        BBE_SVN_2XF32_X (B,pXwr,2*estride);
        BBE_SVN_2XF32_X (D,pXwr,3*estride);
        BBE_SVN_2XF32_XP(A,pXwr,4*estride);
    }
}
#endif

#endif
