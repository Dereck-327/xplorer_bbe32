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
    Streaming to Packed Conversion for Real and Complex Matrices
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matop.h"
#include "bs_common.h"

/*   N*M=3, Sx=3, Sy=4 */
void csb3x1(complex_fract16* restrict y, const complex_fract16* restrict x, int MN, int L)
{
  const xb_vecNx16* restrict w=(const xb_vecNx16*)x;
  xb_vecNx16* restrict z=(xb_vecNx16* )y;
  xb_vecNx16 X0,X1,X2,X3,Y0,Y1,Y2,Y3;
  int l;
  const xb_vecNx16* restrict w1=w+(L>>(LOG2_BBE_SIMD_WIDTH-1));
  const xb_vecNx16* restrict w2=w+(L>>(LOG2_BBE_SIMD_WIDTH-2));
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT((L&(BBE_SIMD_WIDTH / 2 - 1)) == 0);
  NASSERT(MN==3);
  (void)MN;
  X3=0;
  for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
  {
    BBE_LVNX16_IP(X0,w , 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1,w1, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X2,w2, 2*BBE_SIMD_WIDTH);
    DEINTLV1(Y3,Y1,X3,X1);
    DEINTLV1(Y2,Y0,X2,X0);
    DEINTLV1(X3,X2,Y3,Y2);
    DEINTLV1(X1,X0,Y1,Y0);
    BBE_SVNX16_IP(X0,z, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X1,z, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X2,z, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(X3,z, 2*BBE_SIMD_WIDTH);
  }
}
