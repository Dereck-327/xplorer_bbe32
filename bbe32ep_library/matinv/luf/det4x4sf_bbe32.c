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
 * NatureDSP_Baseband Library API
 * Matrix Decomposition and Inversion Functions

    Compute Determinant from LU decomposition for real matrices (stream ordered)
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"

#if HAVE_VFPU

/*-------------------------------------------------------------------------
Determinant For Stream Ordered Matrices

Description: compute determinant of a real/complex matrix from its LU 
decomposition (see lu<size>sf() and clu<size>sf() functions) by multiplying
together diagonal elements of the upper triangular factor U. This operation
is accomplished for a sequence of LU matrices stored in stream order.

Data format: IEEE-754 Std single precision floating-point

Input:
  N           Matrix size
  L           Number of matrices
Input/Output:
  LU[N*N][L]  Packed L and U factors computed by [c]lu<size>sf()
Output:
  D[L]        Determinant values
Restrictions:
  D,LU        Must not overlap and must be aligned on 32-byte boundary 
  N           Must be greater than 1
  L           Must be a multiple of 8 for real-valued functions, or a mutiple
              of 4 for complex-valued functions.
---------------------------------------------------------------------------*/
void det4x4sf   ( float32_t * restrict D, const float32_t * restrict LU, int L )
#if 0
{
    int l,k;
    NASSERT_ALIGN(D ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(LU,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0 ) return;

    for (l=0; l<L; l++)
    {
        float32_t d=LU[l];
        for (k=1; k<4; k++) d*=LU[l+L*(k+k*4)];
        D[l]=d;
    }
}
#else
{
    const xb_vecN_2xf32 * restrict pLU00=(const xb_vecN_2xf32 *)LU;
    const xb_vecN_2xf32 * restrict pLU11=(const xb_vecN_2xf32 *)(LU+(1+1*4)*L);
    const xb_vecN_2xf32 * restrict pLU22=(const xb_vecN_2xf32 *)(LU+(2+2*4)*L);
    const xb_vecN_2xf32 * restrict pLU33=(const xb_vecN_2xf32 *)(LU+(3+3*4)*L);
          xb_vecN_2xf32 * restrict pD   =(      xb_vecN_2xf32 *)D;
    int l;
    NASSERT_ALIGN(D ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(LU,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    if (L<=0 ) return;

    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecN_2xf32 d,a00,a11,a22,a33;
        BBE_LVN_2XF32_IP(a00,pLU00,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(a11,pLU11,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(a22,pLU22,2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(a33,pLU33,2*BBE_SIMD_WIDTH);
        d=BBE_MULN_2XF32(BBE_MULN_2XF32(a00,a11),BBE_MULN_2XF32(a22,a33));
        BBE_SVN_2XF32_IP(d,pD,2*BBE_SIMD_WIDTH);
    }
}
#endif
#else
DISCARD_FUN(void, det4x4sf, ( float32_t * restrict D, const float32_t * restrict LU, int L ))
#endif
