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
  NatureDSP_Baseband library. Vector Operations
    Dot Product of Complex Vectors
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Vector Operations. */
#include "NatureDSP_Baseband_vector.h"

/*-------------------------------------------------------------------------
Dot Product of Complex Vectors

Description: These routines take two complex vectors and calculate their
dot product.

Representation:
cvdot   Signed fixed-point format
        Input vectors x and y are comprised of 32-bit complex data with
        16-bit signed real and imaginary components of arbitrary formats
        Qx and Qy. Complex dot product is computed in a pair of 40-bit 
        accumulators, which are then rounded, shifted to the right by rsh
        bit positions and saturated to form a 64-bit complex result with
        Qx+Qy-rsh fractional bits in 32-bit real and imaginary components.
cvdotf  IEEE-754 Std. single precision floating-point format for real/imaginary
        components of 64-bit input data and dot product result

Parameters:
Input:
x[N]    Input vector
y[N]    Input vector
rsh     Right shift amount (cvdot)
N       Length of vectors
Returned Value:
Dot product result

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 8 (cvdot) or 4 (cvdotf)
rsh>=0  Right shift amount must be non-negative
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(complex_float,cvdotf,( const complex_float   * restrict x,
                         const complex_float   * restrict y,
                         int N ))
#else
complex_float   cvdotf ( const complex_float   * restrict x,
                         const complex_float   * restrict y,
                         int N )
{
  complex_float ALIGN(8)   z;
  int n;
  xb_vecN_2xf32 x0, y0;
  xb_vecN_2xf32 Acc0, Acc1, Acc2, Acc3;
  xb_vecN_2xf32 a0, a1, a2, a3;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  const xb_vecN_2xf32  * restrict pY = (const xb_vecN_2xf32  *)y;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 4) == 0);
   
  MOV_COMPLEXFLOAT_FROM_N_4XCF32(z,BBE_ZERON_4XCF32());
  if (N <= 0) return z;
  Acc0 = BBE_ZERON_2XF32();
  Acc1 = BBE_ZERON_2XF32(); 
  Acc2 = BBE_ZERON_2XF32();
  Acc3 = BBE_ZERON_2XF32();
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    a0 = BBE_MULMN_2XF32(x0, y0, 0, 4);  /* 00 01 00  x.re*y.re     x.re*y.im */
    a1 = BBE_MULMN_2XF32(x0, y0, 1, 11); /* 01 10 11 -x.im*y.im     x.im*y.re */
    Acc0 = BBE_ADDN_2XF32(Acc0, a0); 
    Acc1 = BBE_ADDN_2XF32(Acc1, a1);  

    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    a2 = BBE_MULMN_2XF32(x0, y0, 0, 4);   /* 00 01 00  x.re*y.re     x.re*y.im */
    a3 = BBE_MULMN_2XF32(x0, y0, 1, 11);  /* 01 10 11 -x.im*y.im     x.im*y.re */
    Acc2 = BBE_ADDN_2XF32(Acc2, a2);
    Acc3 = BBE_ADDN_2XF32(Acc3, a3);
  }
  if (N&(BBE_SIMD_WIDTH/4))
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, pY, 2 * BBE_SIMD_WIDTH);
    a0 = BBE_MULMN_2XF32(x0, y0, 0, 4);  /* 00 01 00  x.re*y.re     x.re*y.im */
    a1 = BBE_MULMN_2XF32(x0, y0, 1, 11); /* 01 10 11 -x.im*y.im     x.im*y.re */
    Acc0 = BBE_ADDN_2XF32(Acc0, a0);
    Acc1 = BBE_ADDN_2XF32(Acc1, a1);
  }
  Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);
  Acc2 = BBE_ADDN_2XF32(Acc2, Acc3);
  Acc0 = BBE_ADDN_2XF32(Acc0, Acc2);
  Acc1 = BBE_SELN_2XF32I(Acc0, Acc0, BBE_SELI_ROTATE_RIGHT_8);
  Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);
  Acc1 = BBE_SELN_2XF32I(Acc0, Acc0, BBE_SELI_ROTATE_RIGHT_4);
  Acc0 = BBE_ADDN_2XF32(Acc0, Acc1);
  {
    xb_vecN_4xcf32 t;
    t = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(Acc0));
    MOV_COMPLEXFLOAT_FROM_N_4XCF32(z,t);
  }
  return z;

} /* cvdotf() */
#endif
