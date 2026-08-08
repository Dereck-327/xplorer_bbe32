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

/* Interleaving function for D=6 */
void firinterp_6d_intlv (int16_t * y, int N)
{

#if 1

  xb_vecNx16 *  restrict pY0;
  xb_vecNx16 *  restrict pY;

  xb_vecNx16 y0, y1, y2, y3, y4, y5;
  xb_vecNx16 x0, x1, x2, x3;
  xb_vecNx16 z0, z1, z2, z3, z4, z5;
  int n;

  pY0 = (      xb_vecNx16 *) y;
  pY  = (      xb_vecNx16 *) y;

  for (n=0;n<N/(BBE_SIMD_WIDTH/2);n++)
  {
   //y0[0],...,y0[7]
    BBE_LVNX16_IP(y0,pY0,2*BBE_SIMD_WIDTH);
    //y1[0],...,y1[7]
    BBE_LVNX16_IP(y1,pY0,2*BBE_SIMD_WIDTH);
    //y2[0],...,y2[7]
    BBE_LVNX16_IP(y2,pY0,2*BBE_SIMD_WIDTH);
    //y3[0],...,y3[7]
    BBE_LVNX16_IP(y3,pY0,2*BBE_SIMD_WIDTH);
    //y4[0],...,y4[7]
    BBE_LVNX16_IP(y4,pY0,2*BBE_SIMD_WIDTH);
    //y5[0],...,y5[7]
    BBE_LVNX16_IP(y5,pY0,2*BBE_SIMD_WIDTH);

    //y0[0],y1[0],y0[1],y1[1],y0[2],y1[2],...,y0[7],y1[7]
    BBE_DSELNX16I(x1,x0,y1,y0,BBE_DSELI_INTERLEAVE_2);
    //y3[0],y4[0],y3[1],y4[1],y3[2],y4[2],...,y3[7],y4[7]
    BBE_DSELNX16I(x3,x2,y4,y3,BBE_DSELI_INTERLEAVE_2);

    BBE_DSELNX16I(y1,y0,x2,x0,BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(y4,y3,x3,x1,BBE_DSELI_INTERLEAVE_4);

    // y2[0],y5[0],...,y2[7],y5[7]
    BBE_DSELNX16I(x1,x0,y5,y2,BBE_DSELI_INTERLEAVE_2);

    BBE_DSELNX16I(z1,z0,x0,y0,BBE_DSELI_INTERLEAVE_C3_STEP_0); 
    BBE_DSELNX16I_H(z1,z2,x0,y1,BBE_DSELI_INTERLEAVE_C3_STEP_1);
    BBE_DSELNX16I(z4,z3,x1,y3,BBE_DSELI_INTERLEAVE_C3_STEP_0); 
    BBE_DSELNX16I_H(z4,z5,x1,y4,BBE_DSELI_INTERLEAVE_C3_STEP_1);

    BBE_SVNX16_IP(z0,pY,2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(z1,pY,2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(z2,pY,2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(z3,pY,2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(z4,pY,2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(z5,pY,2*BBE_SIMD_WIDTH);
  }

#endif

#if 0

  xb_vecNx16 *  restrict pY0;
  xb_vecNx16 *  restrict pY;

  xb_vecNx16 y0, y1, y2, y3, y4, y5;
  xb_vecNx16 x0, x1, x2, x3;
  xb_vecNx16 z0, z1, z2, z3, z4, z5;
  int n;

  pY0 = (      xb_vecNx16 *) y;
  pY  = (      xb_vecNx16 *) y;

  for (n=0;n<N/(BBE_SIMD_WIDTH/2);n++)
  {
   //y0[0],...,y0[7]
    BBE_LVNX16_IP(y0,pY0,2*BBE_SIMD_WIDTH);
    //y1[0],...,y1[7]
    BBE_LVNX16_IP(y1,pY0,2*BBE_SIMD_WIDTH);
    //y2[0],...,y2[7]
    BBE_LVNX16_IP(y2,pY0,2*BBE_SIMD_WIDTH);
    //y3[0],...,y3[7]
    BBE_LVNX16_IP(y3,pY0,2*BBE_SIMD_WIDTH);
    //y4[0],...,y4[7]
    BBE_LVNX16_IP(y4,pY0,2*BBE_SIMD_WIDTH);
    //y5[0],...,y5[7]
    BBE_LVNX16_IP(y5,pY0,2*BBE_SIMD_WIDTH);

    //y0[0],y1[0],y0[1],y1[1],y0[2],y1[2],...,y0[7],y1[7]
    BBE_DSELNX16I(x1,x0,y1,y0,BBE_DSELI_INTERLEAVE_2);
    //y3[0],y4[0],y3[1],y4[1],y3[2],y4[2],...,y3[7],y4[7]
    BBE_DSELNX16I(x3,x2,y4,y3,BBE_DSELI_INTERLEAVE_2);

    BBE_DSELNX16I(y1,y0,x2,x0,BBE_DSELI_INTERLEAVE_4);
    BBE_DSELNX16I(y4,y3,x3,x1,BBE_DSELI_INTERLEAVE_4);

    // y2[0],y5[0],...,y2[7],y5[7]
    BBE_DSELNX16I(x1,x0,y5,y2,BBE_DSELI_INTERLEAVE_2);

    BBE_DSELNX16I(z1,z0,x0,y0,BBE_DSELI_INTERLEAVE_C3_STEP_0); 
    BBE_DSELNX16I_H(z1,z2,x0,y1,BBE_DSELI_INTERLEAVE_C3_STEP_1);
    BBE_DSELNX16I(z4,z3,x1,y3,BBE_DSELI_INTERLEAVE_C3_STEP_0); 
    BBE_DSELNX16I_H(z4,z5,x1,y4,BBE_DSELI_INTERLEAVE_C3_STEP_1);

    BBE_SVNX16_IP(z0,pY,2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(z1,pY,2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(z2,pY,2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(z3,pY,2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(z4,pY,2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(z5,pY,2*BBE_SIMD_WIDTH);

  }

#endif


}

