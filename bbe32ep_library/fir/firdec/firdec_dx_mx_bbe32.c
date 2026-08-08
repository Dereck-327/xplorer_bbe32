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
#include "firdec_common.h"

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

/* processing function M%16==0 */
static void firdec_proc_Dx_Mx(int16_t * restrict y,
                        const int16_t * restrict x,
                        const int16_t * restrict coef,
                              int16_t * restrict delayLine,
                              int M, int NN, int D )
{

    int16_t   * restrict Y;
    const xb_vecNx16 *          X;
    xb_vecNx16 * restrict pD;
    xb_vecNx16 * restrict pDr;
    xb_vecNx16 * restrict pH;
    const xb_vecNx16 * restrict S0;
    const xb_vecNx16 * restrict S1;
    const xb_vecNx16 * restrict S2;
    const xb_vecNx16 * restrict S3;

    int m, n, N;
    const int NSamples = 64;

    xb_vecNx16 x0, x1;
    xb_vecNx16 y0, y1, y2, y3;
    xb_vecNx40 w0, w1, w2, w3;
    xb_c40     b0, b1, b2, b3;
    valign     S0_va, S1_va, S2_va, S3_va;
    xb_vecNx16 cf0, cf1;
    uintptr_t  px0, px1, px2, px3;

    NASSERT(NN>0 && (NN % 8) == 0);
    NASSERT(M>0 && (M % 16) == 0);
    NASSERT_ALIGN(y,         2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x,         2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(coef,      2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(delayLine, 2 * BBE_SIMD_WIDTH);

    Y = (int16_t   *)y;
    X = (const xb_vecNx16*)x;
    pD = (xb_vecNx16*)delayLine + M / (BBE_SIMD_WIDTH / 2);
     //__Pragma("loop_count min=1");
    for (; NN > 0; NN -= NSamples)
    {
        N = NN > NSamples ? NSamples : NN;
        if ((D*N) & (BBE_SIMD_WIDTH / 2))
        {
          BBE_LVNX16_IP(x0, X, 2 * BBE_SIMD_WIDTH);
          BBE_SVNX16_IP(x0, pD, 2 * BBE_SIMD_WIDTH);
        }
        __Pragma("ymemory(X)");
        __Pragma("ymemory(pD)");
        __Pragma("loop_count min=2");
        for (n = 0; n < (D*N) / (BBE_SIMD_WIDTH ); n++)
        {
          BBE_LVNX16_IP(x0, X, 2 * BBE_SIMD_WIDTH);
          BBE_SVNX16_IP(x0, pD, 2 * BBE_SIMD_WIDTH);
          BBE_LVNX16_IP(x0, X, 2 * BBE_SIMD_WIDTH);
          BBE_SVNX16_IP(x0, pD, 2 * BBE_SIMD_WIDTH);
        }

        if (M == BBE_SIMD_WIDTH)
        {
          pH = (xb_vecNx16*)coef;
          BBE_LVNX16_IP(cf0, pH, 0 * 2 * BBE_SIMD_WIDTH);
          BBE_DSELNX16I(cf1, cf0, cf0, cf0, BBE_DSELI_INTERLEAVE_1);

          __Pragma("loop_count min=2");
          __Pragma("loop_count max=16");
          for (n = 0; n < N / 4; n++)
          {

            px0 = (uintptr_t)(delayLine + 8 * D*n + 2);
            px1 = (uintptr_t)(delayLine + 8 * D*n + 2 * D + 2);

            px2 = XT_ADDX8(D, px0);
            px3 = XT_ADDX8(D, px1);
          //
            S0 = (xb_vecNx16*)(px0);
            S1 = (xb_vecNx16*)(px1);
            S2 = (xb_vecNx16*)(px2);
            S3 = (xb_vecNx16*)(px3);


            S0_va = BBE_LANX16_PP(S0);
            S1_va = BBE_LANX16_PP(S1);
            S2_va = BBE_LANX16_PP(S2);
            S3_va = BBE_LANX16_PP(S3);

            //BBE_LVNX16_IP(cf0, pH, 0*2 * BBE_SIMD_WIDTH);
            //BBE_DSELNX16I(cf1, cf0, cf0, cf0, BBE_DSELI_INTERLEAVE_1);

            BBE_LAVNX16_XP(x0, S0_va, S0, 2 * BBE_SIMD_WIDTH);
            BBE_LAVNX16_XP(x1, S0_va, S0, 2 * BBE_SIMD_WIDTH);
            w0 = BBE_MULNX16( x0, cf0);
            BBE_MULANX16(w0, x1, cf1);

            //--------------------------------------
            BBE_LAVNX16_XP(x0, S1_va, S1, 2 * BBE_SIMD_WIDTH);
            BBE_LAVNX16_XP(x1, S1_va, S1, 2 * BBE_SIMD_WIDTH);
            w1 = BBE_MULNX16(x0, cf0);
            BBE_MULANX16(w1, x1, cf1);

            //--------------------------------------
            BBE_LAVNX16_XP(x0, S2_va, S2, 2 * BBE_SIMD_WIDTH);
            BBE_LAVNX16_XP(x1, S2_va, S2, 2 * BBE_SIMD_WIDTH);
            w2 = BBE_MULNX16(x0, cf0);
            BBE_MULANX16(w2, x1, cf1);

            //--------------------------------------
            BBE_LAVNX16_XP(x0, S3_va, S3, 2 * BBE_SIMD_WIDTH);
            BBE_LAVNX16_XP(x1, S3_va, S3, 2 * BBE_SIMD_WIDTH);
            w3 = BBE_MULNX16(x0, cf0);
            BBE_MULANX16(w3, x1, cf1);

            b0 = BBE_RADDNX40C(w0);
            b1 = BBE_RADDNX40C(w1);
            b2 = BBE_RADDNX40C(w2);
            b3 = BBE_RADDNX40C(w3);

            w0 = BBE_MOVNX40_FROMC40(b0);
            w1 = BBE_MOVNX40_FROMC40(b1);
            w2 = BBE_MOVNX40_FROMC40(b2);
            w3 = BBE_MOVNX40_FROMC40(b3);

            y0 = BBE_PACKQNX40(w0);
            y1 = BBE_PACKQNX40(w1);
            y2 = BBE_PACKQNX40(w2);
            y3 = BBE_PACKQNX40(w3);
            BBE_SPNX16_IP(y0, Y, 2 * 2);
            BBE_SPNX16_IP(y1, Y, 2 * 2);
            BBE_SPNX16_IP(y2, Y, 2 * 2);
            BBE_SPNX16_IP(y3, Y, 2 * 2);

          }
          pDr = (xb_vecNx16*)(delayLine + 2 * D*N);
          pD = (xb_vecNx16*)(delayLine);
          BBE_LVNX16_IP(x0, pDr, 2 * BBE_SIMD_WIDTH);
          BBE_SVNX16_IP(x0, pD, 2 * BBE_SIMD_WIDTH);
          BBE_LVNX16_IP(x0, pDr, 2 * BBE_SIMD_WIDTH);
          BBE_SVNX16_IP(x0, pD, 2 * BBE_SIMD_WIDTH);
        }
        else
        {

          //__Pragma("ymemory(pH)");
         __Pragma("loop_count min=2");
         // __Pragma("concurrent");
          for (n = 0; n < N / 4; n++)
          {
            pH = (xb_vecNx16*)coef;
            px0 = (uintptr_t)(delayLine + 8 * D*n + 2);
            px1 = (uintptr_t)(delayLine + 8 * D*n + 2 * D + 2);

            px2 = XT_ADDX8(D, px0);
            px3 = XT_ADDX8(D, px1);

            S0 = (xb_vecNx16*)(px0);
            S1 = (xb_vecNx16*)(px1);
            S2 = (xb_vecNx16*)(px2);
            S3 = (xb_vecNx16*)(px3);

            S0_va = BBE_LANX16_PP(S0);
            S1_va = BBE_LANX16_PP(S1);
            S2_va = BBE_LANX16_PP(S2);
            S3_va = BBE_LANX16_PP(S3);

            w0 = w1 = w2 = w3 = 0;
            __Pragma("loop_count min=1");
            for (m = 0; m < M / (2*BBE_SIMD_WIDTH); m++)
            {
              BBE_LVNX16_IP(cf0, pH, 2 * BBE_SIMD_WIDTH);
              BBE_DSELNX16I(cf1, cf0, cf0, cf0, BBE_DSELI_INTERLEAVE_1);

              BBE_LAVNX16_XP(x0, S0_va, S0, 2 * BBE_SIMD_WIDTH);
              BBE_LAVNX16_XP(x1, S0_va, S0, 2 * BBE_SIMD_WIDTH);
              BBE_MULANX16(w0, x0, cf0);
              BBE_MULANX16(w0, x1, cf1);

              //--------------------------------------
              BBE_LAVNX16_XP(x0, S1_va, S1, 2 * BBE_SIMD_WIDTH);
              BBE_LAVNX16_XP(x1, S1_va, S1, 2 * BBE_SIMD_WIDTH);
              BBE_MULANX16(w1, x0, cf0);
              BBE_MULANX16(w1, x1, cf1);

              //--------------------------------------
              BBE_LAVNX16_XP(x0, S2_va, S2, 2 * BBE_SIMD_WIDTH);
              BBE_LAVNX16_XP(x1, S2_va, S2, 2 * BBE_SIMD_WIDTH);
              BBE_MULANX16(w2, x0, cf0);
              BBE_MULANX16(w2, x1, cf1);

              //--------------------------------------
              BBE_LAVNX16_XP(x0, S3_va, S3, 2 * BBE_SIMD_WIDTH);
              BBE_LAVNX16_XP(x1, S3_va, S3, 2 * BBE_SIMD_WIDTH);
              BBE_MULANX16(w3, x0, cf0);
              BBE_MULANX16(w3, x1, cf1);

              BBE_LVNX16_IP(cf0, pH, 2 * BBE_SIMD_WIDTH);
              BBE_DSELNX16I(cf1, cf0, cf0, cf0, BBE_DSELI_INTERLEAVE_1);

              BBE_LAVNX16_XP(x0, S0_va, S0, 2 * BBE_SIMD_WIDTH);
              BBE_LAVNX16_XP(x1, S0_va, S0, 2 * BBE_SIMD_WIDTH);
              BBE_MULANX16(w0, x0, cf0);
              BBE_MULANX16(w0, x1, cf1);

              //--------------------------------------
              BBE_LAVNX16_XP(x0, S1_va, S1, 2 * BBE_SIMD_WIDTH);
              BBE_LAVNX16_XP(x1, S1_va, S1, 2 * BBE_SIMD_WIDTH);
              BBE_MULANX16(w1, x0, cf0);
              BBE_MULANX16(w1, x1, cf1);

              //--------------------------------------
              BBE_LAVNX16_XP(x0, S2_va, S2, 2 * BBE_SIMD_WIDTH);
              BBE_LAVNX16_XP(x1, S2_va, S2, 2 * BBE_SIMD_WIDTH);
              BBE_MULANX16(w2, x0, cf0);
              BBE_MULANX16(w2, x1, cf1);

              //--------------------------------------
              BBE_LAVNX16_XP(x0, S3_va, S3, 2 * BBE_SIMD_WIDTH);
              BBE_LAVNX16_XP(x1, S3_va, S3, 2 * BBE_SIMD_WIDTH);
              BBE_MULANX16(w3, x0, cf0);
              BBE_MULANX16(w3, x1, cf1);
            }
            if (M&(2 * BBE_SIMD_WIDTH - 1))
            {
              BBE_LVNX16_IP(cf0, pH, 2 * BBE_SIMD_WIDTH);
              BBE_DSELNX16I(cf1, cf0, cf0, cf0, BBE_DSELI_INTERLEAVE_1);

              BBE_LAVNX16_XP(x0, S0_va, S0, 4 * 8);
              BBE_LAVNX16_XP(x1, S0_va, S0, 4 * 8);
              BBE_MULANX16(w0, x0, cf0);
              BBE_MULANX16(w0, x1, cf1);

              //---------------------------------------
              BBE_LAVNX16_XP(x0, S1_va, S1, 4 * 8);
              BBE_LAVNX16_XP(x1, S1_va, S1, 4 * 8);
              BBE_MULANX16(w1, x0, cf0);
              BBE_MULANX16(w1, x1, cf1);

              //---------------------------------------
              BBE_LAVNX16_XP(x0, S2_va, S2, 4 * 8);
              BBE_LAVNX16_XP(x1, S2_va, S2, 4 * 8);
              BBE_MULANX16(w2, x0, cf0);
              BBE_MULANX16(w2, x1, cf1);

              //---------------------------------------
              BBE_LAVNX16_XP(x0, S3_va, S3, 4 * 8);
              BBE_LAVNX16_XP(x1, S3_va, S3, 4 * 8);
              BBE_MULANX16(w3, x0, cf0);
              BBE_MULANX16(w3, x1, cf1);
            }
            b0 = BBE_RADDNX40C(w0);
            b1 = BBE_RADDNX40C(w1);
            b2 = BBE_RADDNX40C(w2);
            b3 = BBE_RADDNX40C(w3);

            w0 = BBE_MOVNX40_FROMC40(b0);
            w1 = BBE_MOVNX40_FROMC40(b1);
            w2 = BBE_MOVNX40_FROMC40(b2);
            w3 = BBE_MOVNX40_FROMC40(b3);

            y0 = BBE_PACKQNX40(w0);
            y1 = BBE_PACKQNX40(w1);
            y2 = BBE_PACKQNX40(w2);
            y3 = BBE_PACKQNX40(w3);

            BBE_SPNX16_IP(y0, Y, 2 * 2);
            BBE_SPNX16_IP(y1, Y, 2 * 2);
            BBE_SPNX16_IP(y2, Y, 2 * 2);
            BBE_SPNX16_IP(y3, Y, 2 * 2);
          }
          pDr = (xb_vecNx16*)(delayLine + 2 * D*N);
          pD = (xb_vecNx16*)(delayLine);
          __Pragma("ymemory(pD)");
          __Pragma("ymemory(pDr)");
          __Pragma("loop_count min=1");
          for (m = 0; m < M / (BBE_SIMD_WIDTH / 2); m++)
          {
            BBE_LVNX16_IP(x0, pDr, 2 * BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(x0, pD, 2 * BBE_SIMD_WIDTH);
          }
        }
        
    }
}
const tFirFxdxns firdec_dx_8n    ={&firdec_alloc_dx ,firdec_proc_Dx_Mx};
