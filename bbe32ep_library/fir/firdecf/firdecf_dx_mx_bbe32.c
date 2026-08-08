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
    NatureDSP_Baseband library. FIR filters and Related Functions
    Decimating Block Complex FIR Filter
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
#include "firdecf_common.h"

/*-------------------------------------------------------------------------
Decimating Block Complex FIR Filter

Computes a complex FIR filter (direct-form) with decimation using real IR 
stored in vector h. The complex data input is stored in vector x. The filter
output result is stored in vector y. The filter calculates N output samples
using M coefficients and requires last D*N+M-1 samples in the delay line.

NOTE:
To avoid aliasing, the IR should be synthesized in such a way that filter pass
band is limited by input sample frequency divided by 2*D.

Representation:
firdec   16-bit signed fixed-point format
         Filter coefficients are Q15
         Number of fractional bits for input/output samples is user-difined
firdecf  IEEE-754 Std. single precision floating-point format for filter 
         coefficients and input/output samples

Parameters:
Input:
D        Decimation factor
N        Length of output sample block
M        Length of filter
h[M]     Filter coefficients; h[0] is to be multiplied by the newest 
         sample
x[N*D]   Input complex samples
Output:
y[N]     Output complex samples

Restrictions:
x,y      Must not overlap
x,y      Aligned on 32-byte boundary
N        Multiple of 8 (firdec) or 4 (firdecf)
M        2,4,8 or a positive multiple of 16 for D=2,3,4; or 
         a positive multiple of 16 for D>4
D>1

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for
filter lengths M=2,4,8,16 and 32 and decimation factors D=2,3 and 4, in
any combination.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, firdec[f]_init returns NULL handle.
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN (void,firdecf_proc_DX_MX,( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D ))
#else
void firdecf_proc_DX_MX( complex_float * restrict y,
                   const complex_float * restrict x,
                   const float32_t     * restrict coef,
                         float32_t     * restrict delayLine,
                   int M, int N, int D )
{
    int m, n, NN;
    const int NSamples = 64;
          float32_t     * restrict pY;
    const xb_vecN_2xf32 *          pX;
          xb_vecN_2xf32 * restrict pD;
          xb_vecN_2xf32 * restrict pDr;
    const xb_vecN_2xf32 * restrict pH;
    const xb_vecN_2xf32 * restrict S0;
    const xb_vecN_2xf32 * restrict S1;
    const xb_vecN_2xf32 * restrict S2;
    const xb_vecN_2xf32 * restrict S3;

    xb_vecN_2xf32 c0, c1;
    xb_vecN_2xf32 x0, x1;
    xb_vecN_2xf32 y0, y1, y2, y3;
    xb_vecN_2xf32 t0, t1, t2, t3, t4, t5, t6, t7;
    valign        S0_va, S1_va, S2_va, S3_va;
    uintptr_t     px0, px1, px2, px3;

    NASSERT(N>0 && (N % 4) == 0);
    NASSERT(M>0 && (M % 16) == 0);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(coef);
    NASSERT_ALIGN32(delayLine);

    pY = (      float32_t     *)y;
    pX = (const xb_vecN_2xf32 *)x;
    pD = (      xb_vecN_2xf32 *)delayLine + M / (BBE_SIMD_WIDTH / 4);

    for (; N > 0; N -= NSamples)
    {
        NN = N > NSamples ? NSamples : N;
        for (n = 0; n < (D*NN) / (BBE_SIMD_WIDTH / 4); n++)
        {
            BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(x0, pD, 2 * BBE_SIMD_WIDTH);
        }

        for (n = 0; n < NN / 4; n++)
        {
            y0 = y1 = y2 = y3 = BBE_ZERON_2XF32();
            pH = (xb_vecN_2xf32*)coef;
            px0 = (uintptr_t)(delayLine + 8 * D*n + 2);
            px1 = (uintptr_t)(delayLine + 8 * D*n + 2 * D + 2);
            px2 = XT_ADDX8(2 * D, px0);
            px3 = XT_ADDX8(2 * D, px1);
            S0 = (xb_vecN_2xf32*)(px0);
            S1 = (xb_vecN_2xf32*)(px1);
            S2 = (xb_vecN_2xf32*)(px2);
            S3 = (xb_vecN_2xf32*)(px3);
            S0_va = BBE_LAN_2XF32_PP(S0);
            S1_va = BBE_LAN_2XF32_PP(S1);
            S2_va = BBE_LAN_2XF32_PP(S2);
            S3_va = BBE_LAN_2XF32_PP(S3);

            for (m = 0; m < M / BBE_SIMD_WIDTH; m++)
            {
                BBE_LVN_2XF32_IP(c0, pH, 2 * BBE_SIMD_WIDTH);
                BBE_DSELN_2XF32I(c1, c0, c0, c0, BBE_DSELI_INTERLEAVE_2);

                BBE_LAN_2XF32_IP(x0, S0_va, S0);
                BBE_LAN_2XF32_IP(x1, S0_va, S0);
                BBE_MULAN_2XF32(y0, x0, c0);
                BBE_MULAN_2XF32(y0, x1, c1);

                BBE_LAN_2XF32_IP(x0, S1_va, S1);
                BBE_LAN_2XF32_IP(x1, S1_va, S1);
                BBE_MULAN_2XF32(y1, x0, c0);
                BBE_MULAN_2XF32(y1, x1, c1);

                BBE_LAN_2XF32_IP(x0, S2_va, S2);
                BBE_LAN_2XF32_IP(x1, S2_va, S2);
                BBE_MULAN_2XF32(y2, x0, c0);
                BBE_MULAN_2XF32(y2, x1, c1);

                BBE_LAN_2XF32_IP(x0, S3_va, S3);
                BBE_LAN_2XF32_IP(x1, S3_va, S3);
                BBE_MULAN_2XF32(y3, x0, c0);
                BBE_MULAN_2XF32(y3, x1, c1);

                BBE_LVN_2XF32_IP(c0, pH, 2 * BBE_SIMD_WIDTH);
                BBE_DSELN_2XF32I(c1, c0, c0, c0, BBE_DSELI_INTERLEAVE_2);

                BBE_LAN_2XF32_IP(x0, S0_va, S0);
                BBE_LAN_2XF32_IP(x1, S0_va, S0);
                BBE_MULAN_2XF32(y0, x0, c0);
                BBE_MULAN_2XF32(y0, x1, c1);

                BBE_LAN_2XF32_IP(x0, S1_va, S1);
                BBE_LAN_2XF32_IP(x1, S1_va, S1);
                BBE_MULAN_2XF32(y1, x0, c0);
                BBE_MULAN_2XF32(y1, x1, c1);

                BBE_LAN_2XF32_IP(x0, S2_va, S2);
                BBE_LAN_2XF32_IP(x1, S2_va, S2);
                BBE_MULAN_2XF32(y2, x0, c0);
                BBE_MULAN_2XF32(y2, x1, c1);

                BBE_LAN_2XF32_IP(x0, S3_va, S3);
                BBE_LAN_2XF32_IP(x1, S3_va, S3);
                BBE_MULAN_2XF32(y3, x0, c0);
                BBE_MULAN_2XF32(y3, x1, c1);
            }
                
            t0 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), y0, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
            t1 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), y1, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
            t2 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), y2, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
            t3 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), y3, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
            t4 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), y0, BBE_SELI_EXTRACT_2_OF_4_OFF_2);
            t5 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), y1, BBE_SELI_EXTRACT_2_OF_4_OFF_2);
            t6 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), y2, BBE_SELI_EXTRACT_2_OF_4_OFF_2);
            t7 = BBE_SELN_2XF32I(BBE_ZERON_2XF32(), y3, BBE_SELI_EXTRACT_2_OF_4_OFF_2);
                
            t0 = BBE_MOVN_2XF32_FROMF32(BBE_RADDN_2XF32(t0));
            t1 = BBE_MOVN_2XF32_FROMF32(BBE_RADDN_2XF32(t1));
            t2 = BBE_MOVN_2XF32_FROMF32(BBE_RADDN_2XF32(t2));
            t3 = BBE_MOVN_2XF32_FROMF32(BBE_RADDN_2XF32(t3));
            t4 = BBE_MOVN_2XF32_FROMF32(BBE_RADDN_2XF32(t4));
            t5 = BBE_MOVN_2XF32_FROMF32(BBE_RADDN_2XF32(t5));
            t6 = BBE_MOVN_2XF32_FROMF32(BBE_RADDN_2XF32(t6));
            t7 = BBE_MOVN_2XF32_FROMF32(BBE_RADDN_2XF32(t7));

            BBE_SSN_2XF32_IP(t0, pY, 2 * 2);
            BBE_SSN_2XF32_IP(t4, pY, 2 * 2);
            BBE_SSN_2XF32_IP(t1, pY, 2 * 2);
            BBE_SSN_2XF32_IP(t5, pY, 2 * 2);
            BBE_SSN_2XF32_IP(t2, pY, 2 * 2);
            BBE_SSN_2XF32_IP(t6, pY, 2 * 2);
            BBE_SSN_2XF32_IP(t3, pY, 2 * 2);
            BBE_SSN_2XF32_IP(t7, pY, 2 * 2);
        }

        pDr = (xb_vecN_2xf32*)(delayLine + 2 * D*NN);
        pD = (xb_vecN_2xf32*)(delayLine);
        for (m = 0; m < M / (BBE_SIMD_WIDTH / 4); m++)
        {
            BBE_LVN_2XF32_IP(x0, pDr, 2 * BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(x0, pD, 2 * BBE_SIMD_WIDTH);
        }
    }
}// firdecf_proc_DX_MX()
#endif
const tFirFxdxns firdecf_xd_xm = { &firdecf_alloc_dx_mx, firdecf_proc_DX_MX };
