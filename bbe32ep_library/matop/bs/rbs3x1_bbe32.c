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
#include "NatureDSP_Math.h"
#include "bs_common.h"

/*   N=3,M=1, Sx=4, Sy=3 */
void rbs3x1(int16_t* restrict y, const int16_t* restrict x, int MN,int L)
{
    const xb_vecNx16* restrict X=(const xb_vecNx16*)x;
          xb_vecNx16* restrict Y=(xb_vecNx16* )y;
          xb_vecNx16* restrict y1=Y+(L>>4);
          xb_vecNx16* restrict y2=Y+(L>>3);
       
    xb_vecNx16 X0,X1,X2,X3,Y0,Y1,Y2,Y3;
    int l;
    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT((L&(BBE_SIMD_WIDTH - 1)) == 0);
    NASSERT(MN==3);
    (void)MN;
    if (L <= 0) return;
    for (l = 0; l<(L >> LOG2_BBE_SIMD_WIDTH); l++)
    {
      BBE_LVNX16_IP(X0,X, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X1,X, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X2,X, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(X3,X, 2*BBE_SIMD_WIDTH);

      BBE_DSELNX16I(Y1,Y0,X1,X0,BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(Y3,Y2,X3,X2,BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(X1,X0,Y2,Y0,BBE_DSELI_DEINTERLEAVE_1);
      BBE_DSELNX16I(X3,X2,Y3,Y1,BBE_DSELI_DEINTERLEAVE_1);

      BBE_SVNX16_IP(X0,Y , 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(X1,y1, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(X2,y2, 2*BBE_SIMD_WIDTH);
    }
}
