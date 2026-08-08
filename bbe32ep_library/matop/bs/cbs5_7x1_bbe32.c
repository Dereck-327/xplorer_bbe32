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
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
  NatureDSP_Baseband library. Matrix Operations
    Packed to Streaming Conversion for Real and Complex Matrices
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matop.h"
#include "bs_common.h"

/*   N*M=5,6,7, Sx=8, Sy=MN */
void cbs5_7x1(complex_fract16 * restrict y, const complex_fract16 * restrict x, int MN, int L)
{
    const xb_vecNx16* restrict w=(const xb_vecNx16*)x;
    xb_vecNx16* restrict z=(xb_vecNx16* )y;
    int l,off5,off6;
    xb_vecNx16* restrict z1=z+(L>>(LOG2_BBE_SIMD_WIDTH-1));
    xb_vecNx16* restrict z2=z+(L>>(LOG2_BBE_SIMD_WIDTH-2));
    xb_vecNx16* restrict z3=z1+(L>>(LOG2_BBE_SIMD_WIDTH-2));
    xb_vecNx16* restrict z4=z2+(L>>(LOG2_BBE_SIMD_WIDTH-2));
    xb_vecNx16 X0,X1,X2,X3,X4,X5,X6,X7;
    xb_vecNx16 Y0,Y1,Y2,Y3,Y4,Y5,Y6,Y7;
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT((L&(BBE_SIMD_WIDTH/2-1))==0);
    NASSERT(MN>=5 && MN<=7);
    off5=(MN>5)?(L*(4*5)):0;
    off6=(MN>6)?(L*(4*5)):0;
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        BBE_LVNX16_IP(Y0,w, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(Y1,w, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(Y2,w, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(Y3,w, 2*BBE_SIMD_WIDTH);
        INTLV(X1,X0,Y1,Y0);
        INTLV(X3,X2,Y3,Y2);
        INTLV(Y2,Y0,X2,X0);
        INTLV(Y3,Y1,X3,X1);

        BBE_LVNX16_IP(Y4,w, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(Y5,w, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(Y6,w, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(Y7,w, 2*BBE_SIMD_WIDTH);
        INTLV(X5,X4,Y5,Y4);
        INTLV(X7,X6,Y7,Y6);
        INTLV(Y6,Y4,X6,X4);
        INTLV(Y7,Y5,X7,X5);

        INTLV(X4,X0,Y4,Y0);
        INTLV(X5,X1,Y5,Y1);
        INTLV(X6,X2,Y6,Y2);
        INTLV(X7,X3,Y7,Y3);

        BBE_SVNX16_X(X5,z , off5);
        BBE_SVNX16_X(X6,z1, off6);
        BBE_SVNX16_IP(X0,z , 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(X1,z1, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(X2,z2, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(X3,z3, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(X4,z4, 2*BBE_SIMD_WIDTH);
    }
}
