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
    Interpolating Block Complex FIR Filter
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

#include "firinterp_common.h"

/*-------------------------------------------------------------------------
Interpolating Block Complex FIR Filter

Computes a complex FIR filter (direct-form) with interpolation using real
IR stored in vector h. The complex data input is stored in vector x. The
filter output result is stored in vector y. The filter calculates N*D complex
output samples using M*D coefficients and requires last N+M-1 samples in the
delay line.

Representation:
firinterp   16-bit signed fixed-point format
            Filter coefficients are Q15
            Number of fractional bits for input/output samples is user-difined
firinterpf  IEEE-754 Std. single precision floating-point format for filter 
            coefficients and input/output samples

Parameters:
Input:
D           Interpolation ratio 
N           Length of input sample block
M           Length of subfilter. Total length of filter is M*D
h[M*D]      Filter coefficients; h[0] is to be multiplied by the newest 
            sample,Q15
x[N]        Input complex samples
Output:
y[N*D]      Output complex samples

Restrictions:
x,y         Must not overlap
x,y         Aligned on 32-byte boundary
N           Multiple of 8 (firinterp) or 4 (firinterpf)
M           2,4,8 or a positive multiple of 16 for D=2,3,4,6,12; or 
            a positive multiple of 8 for other D
D>1

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for
subfilter lengths M=2,4,8,16 and 32 and interpolation factors D=2,3 and 4,
in any combination.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, firinterp[f]_init returns NULL handle.
-------------------------------------------------------------------------*/

/* Interleaving function for D=12 */
void firinterp_12d_intlv (int16_t * y, int N)
{
    int n;
    xb_vecNx16 *  restrict pY;
    xb_vecNx16 *  restrict pY0;
    valign y_align;

    pY0 = (xb_vecNx16 *)y;
    pY = (xb_vecNx16 *)(y);
    y_align = BBE_ZALIGN();
    for (n = 0; n<N / (BBE_SIMD_WIDTH / 2); n++)
    {
        xb_vecNx16 y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11;
        xb_vecNx16 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11;
        xb_vecNx16 z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11;
        //y0[0],...,y0[7]
        BBE_LVNX16_XP(y0, pY0, 2 * BBE_SIMD_WIDTH);
        //y1[0],...,y1[7]
        BBE_LVNX16_XP(y1, pY0, 2 * BBE_SIMD_WIDTH);
        //y2[0],...,y2[7]
        BBE_LVNX16_XP(y2, pY0, 2 * BBE_SIMD_WIDTH);
        //y3[0],...,y3[7]
        BBE_LVNX16_XP(y3, pY0, 2 * BBE_SIMD_WIDTH);
        //y4[0],...,y4[7]
        BBE_LVNX16_XP(y4, pY0, 2 * BBE_SIMD_WIDTH);
        //y5[0],...,y5[7]
        BBE_LVNX16_XP(y5, pY0, 2 * BBE_SIMD_WIDTH);
        //y6[0],...,y6[7]
        BBE_LVNX16_XP(y6, pY0, 2 * BBE_SIMD_WIDTH);
        //y7[0],...,y7[7]
        BBE_LVNX16_XP(y7, pY0, 2 * BBE_SIMD_WIDTH);
        //y8[0],...,y8[7]
        BBE_LVNX16_XP(y8, pY0, 2 * BBE_SIMD_WIDTH);
        //y9[0],...,y9[7]
        BBE_LVNX16_XP(y9, pY0, 2 * BBE_SIMD_WIDTH);
        //y10[0],...,y10[7]
        BBE_LVNX16_XP(y10, pY0, 2 * BBE_SIMD_WIDTH);
        //y11[0],...,y11[7]
        BBE_LVNX16_XP(y11, pY0, 2 * BBE_SIMD_WIDTH);

        //y0[0],y1[0],y0[1],y1[1],y0[2],y1[2],...,y0[7],y1[7]
        BBE_DSELNX16I(x1, x0, y1, y0, BBE_DSELI_INTERLEAVE_2);
        //y2[0],y3[0],y2[1],y3[1],y2[2],y3[2],...,y2[7],y3[7]
        BBE_DSELNX16I(x3, x2, y3, y2, BBE_DSELI_INTERLEAVE_2);

        BBE_DSELNX16I(z1, z0, x2, x0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(z3, z2, x3, x1, BBE_DSELI_INTERLEAVE_4);

        //y4[0],y5[0],y4[1],y5[1],y4[2],y5[2],...,y4[7],y5[7]
        BBE_DSELNX16I(x5, x4, y5, y4, BBE_DSELI_INTERLEAVE_2);
        //y6[0],y7[0],y6[1],y7[1],y6[2],y7[2],...,y6[7],y7[7]
        BBE_DSELNX16I(x7, x6, y7, y6, BBE_DSELI_INTERLEAVE_2);

        BBE_DSELNX16I(z5, z4, x6, x4, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(z7, z6, x7, x5, BBE_DSELI_INTERLEAVE_4);

        //y8[0],y9[0],y8[1],y9[1],y8[2],y9[2],...,y8[7],y9[7]
        BBE_DSELNX16I(x9, x8, y9, y8, BBE_DSELI_INTERLEAVE_2);
        //y10[0],y11[0],y10[1],y11[1],y10[2],y11[2],...,y10[7],y11[7]
        BBE_DSELNX16I(x11, x10, y11, y10, BBE_DSELI_INTERLEAVE_2);

        BBE_DSELNX16I(z9, z8, x10, x8, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(z11, z10, x11, x9, BBE_DSELI_INTERLEAVE_4);

        //y0[0],y1[0],...,y7[0]
        y0 = BBE_SELNX16I(z4, z0, BBE_SELI_EXTRACT_LO_HALVES);
        //y0[1],y1[1],..,y7[1]
        y1 = BBE_SELNX16I(z4, z0, BBE_SELI_EXTRACT_HI_HALVES);
        //y0[2],y1[2],..,y7[2],
        y2 = BBE_SELNX16I(z5, z1, BBE_SELI_EXTRACT_LO_HALVES);
        //y0[3],y1[3],...,y7[3],
        y3 = BBE_SELNX16I(z5, z1, BBE_SELI_EXTRACT_HI_HALVES);
        //y0[4],y1[4],..,y7[4]
        y4 = BBE_SELNX16I(z6, z2, BBE_SELI_EXTRACT_LO_HALVES);
        //y0[5],y1[5],..,y7[5]
        y5 = BBE_SELNX16I(z6, z2, BBE_SELI_EXTRACT_HI_HALVES);
        //y0[6],y1[6],..,y7[6]
        y6 = BBE_SELNX16I(z7, z3, BBE_SELI_EXTRACT_LO_HALVES);
        //y0[7],y1[7],..,y7[7]
        y7 = BBE_SELNX16I(z7, z3, BBE_SELI_EXTRACT_HI_HALVES);

        BBE_SAVNX16_XP(y0, y_align, pY, 2 * 16);
        BBE_SAVNX16_XP(z8, y_align, pY, 2 * 8);

        z8 = BBE_SELNX16I(z8, z8, BBE_SELI_ROTATE_RIGHT_8);
        BBE_SAVNX16_XP(y1, y_align, pY, 2 * 16);
        BBE_SAVNX16_XP(z8, y_align, pY, 2 * 8);

        BBE_SAVNX16_XP(y2, y_align, pY, 2 * 16);
        BBE_SAVNX16_XP(z9, y_align, pY, 2 * 8);

        z9 = BBE_SELNX16I(z9, z9, BBE_SELI_ROTATE_RIGHT_8);
        BBE_SAVNX16_XP(y3, y_align, pY, 2 * 16);
        BBE_SAVNX16_XP(z9, y_align, pY, 2 * 8);

        BBE_SAVNX16_XP(y4, y_align, pY, 2 * 16);
        BBE_SAVNX16_XP(z10, y_align, pY, 2 * 8);

        z10 = BBE_SELNX16I(z10, z10, BBE_SELI_ROTATE_RIGHT_8);
        BBE_SAVNX16_XP(y5, y_align, pY, 2 * 16);
        BBE_SAVNX16_XP(z10, y_align, pY, 2 * 8);

        BBE_SAVNX16_XP(y6, y_align, pY, 2 * 16);
        BBE_SAVNX16_XP(z11, y_align, pY, 2 * 8);

        z11 = BBE_SELNX16I(z11, z11, BBE_SELI_ROTATE_RIGHT_8);
        BBE_SAVNX16_XP(y7, y_align, pY, 2 * 16);
        BBE_SAVNX16_XP(z11, y_align, pY, 2 * 8);
    }
    BBE_SAVNX16POS_FP(y_align, pY);
}

