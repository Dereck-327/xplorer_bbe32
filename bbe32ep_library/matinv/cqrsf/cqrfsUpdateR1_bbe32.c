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
  NatureDSP_Baseband Library API
  Matrix Decomposition and Inversion Functions
  QR decomposition, floating point, complex data, stream format
  C code optimized for BBE32EP with VFPU
  IntegrIT, 2006-2017
*/

#include "NatureDSP_types.h"
#include "common.h"
#include "cqrsf_common.h"

#if HAVE_VFPU

/*-----------------------------------------------------------
    partial update of R matrix - update begins from 
    v         pointer to m-th Householder vector 
    K         current number of column in the matrix R
    M,N       matrix size
    L         number of matrices
    Input/output:
    R         pointer to the beginning of column  
-----------------------------------------------------------*/
// M==1
void cqrfsUpdateR1(
                              float32_t* restrict R,
                        const float32_t* restrict v,
                        int K,int M,int N, int L)
{
#if 0
    int l, n;
    int _2L = 2 * L;
    float32_t v_re, v_im, b_re, b_im, z_re, z_im;
    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(v, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && (L % (BBE_SIMD_WIDTH / 4) == 0));
    NASSERT(M == 1);

    for (n = 0; n < (N - K); n++)
    {
        for (l = 0; l < L; l++)
        {
            // compute v'B first
            v_re = v[l * 2 + 0];
            v_im = v[l * 2 + 1];
            b_re = R[l * 2 + 0];
            b_im = R[l * 2 + 1];
            z_re = (v_re*b_re) + (v_im*b_im);
            z_im = (v_re*b_im) - (v_im*b_re);
            z_re = 2.f*z_re;
            z_im = 2.f*z_im;

            b_re -= (z_re*v_re) - (z_im*v_im);
            b_im -= (z_re*v_im) + (z_im*v_re);
            R[l * 2 + 0] = b_re;
            R[l * 2 + 1] = b_im;
        }
        R += _2L; // next column
    }
#endif // 0

    int l;

    const xb_vecN_2xf32 * restrict pRr;
          xb_vecN_2xf32 * restrict pR;
    const xb_vecN_2xf32 * restrict pV;

    xb_vecN_2xf32 Acc, V0, R0, Z0;

    NASSERT_ALIGN(R, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(v, 2 * BBE_SIMD_WIDTH);
    NASSERT(L > 0 && (L % (BBE_SIMD_WIDTH / 4) == 0));
    NASSERT(M == 1);

    WUR_CBEGIN((uintptr_t)v);
    WUR_CEND((uintptr_t)v + 2 * L * sizeof(float32_t));

    pRr = (const xb_vecN_2xf32 *)R;
    pR = (xb_vecN_2xf32 *)R;
    pV = (const xb_vecN_2xf32 *)v;
    for (l = 0; l < (N - K)*L; l += (BBE_SIMD_WIDTH / 4))
    {
        BBE_LVN_2XF32_IC(V0, pV);
        BBE_LVN_2XF32_IP(R0, pRr, 2 * BBE_SIMD_WIDTH);
        Acc = BBE_MULMN_2XF32(V0, R0, 0, 4);
        BBE_MULMASN_2XF32(Acc, V0, R0, 2, 11);
        Z0 = BBE_ADDN_2XF32(Acc, Acc);
        //Z0 = BBE_MULN_2XF32(Acc, 2.f);
        BBE_MULMASN_2XF32(R0, Z0, V0, 3, 4);
        BBE_MULMASN_2XF32(R0, Z0, V0, 2, 11);
        BBE_SVN_2XF32_IP(R0, pR, 2 * BBE_SIMD_WIDTH);
    }
}
#endif
