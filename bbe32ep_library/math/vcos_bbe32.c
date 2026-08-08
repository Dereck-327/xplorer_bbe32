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
  NatureDSP_Baseband library. Vector Mathematics
    Sine/Cosine
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_math.h"
#include "sineTbl.h"

/*-------------------------------------------------------------------------
Sine/Cosine 

Description: These functions compute sine or cosine of input data

Representation:
vsine,vcos,    16-bit signed fixed-point format Q15 for input/output data
ssine,scos     It is assumed that input angular values are normalized by pi.
               That is, Fixed-point functions actually compute sin(pi*x) or
               cos(pi*x)
vsinef,vcosf,  IEEE-754 Std. single precision floating-point format for
ssinef,scosf   input/output data. Input data are treated as angular values
               specified in radians. Floating-point functions limit the
               rangw of allowable input values, see note 2.

Accuracy:
2 LSB - vsine(),vcos(), ssine(), scos(),
2 ULP - vsinef(), vcosf(), ssinef(), scosf(),
3 ULP - vfastsinef(),vfastcosf()

Notes for non-fast versions:
1. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
2. Floating-point functions require that input value belongs to the 
   closed range [-102940.0,102940.0], otherwise the respective result
   is NaN.

Input domain for 'fast' versions vfastsinef(),vfastcosf()
|x|<804.2477
The output value is not defined outside of this range or accuracy is degraded

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vsine,vcos) or 8 (vsinef,vcosf,vfastsinef,vfastcosf)
-------------------------------------------------------------------------*/

void vcos ( int16_t * restrict y,
       const int16_t * restrict x,
       int N )
{
    vselN sel;
    vboolN neg;
    xb_vecNx16 ofs, t;
    xb_vecNx16 x0, y0, /*SrcX,*/ ArgShift;
    xb_vecNx16 t1, t3, t4, tmp_t;
    xb_vecNx40 w0, z0;
    int k;
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    //const xb_vecNx16 * restrict pX1 = (const xb_vecNx16 *)x;
    const xb_vecNx16 * restrict T = (const xb_vecNx16 *)sineTbl;
    xb_vecNx16 * restrict pY = (xb_vecNx16 *)y;

    if (N <= 0) return;
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(N>0 && N % BBE_SIMD_WIDTH == 0);

    ArgShift = BBE_MOVVA16(16384); // shift to converse sin(x) to cos(x)

    t1 = BBE_LVNX16_I(T, 2 * BBE_SIMD_WIDTH);
    t3 = BBE_LVNX16_I(T, 6 * BBE_SIMD_WIDTH);
    t4 = BBE_LVNX16_I(T, 8 * BBE_SIMD_WIDTH);
    z0 = BBE_MOVWA32(0);

    __Pragma("loop_count min=1");
    for (k = 0; k<(N >> LOG2_BBE_SIMD_WIDTH); k++)
    {
        BBE_LVNX16_XP(x0, pX, 2 * BBE_SIMD_WIDTH);

        x0 = BBE_ADDNX16(x0, ArgShift); // converse cosine argument to sine argument
        neg = BBE_LTNX16(x0, 0);

        //x0 = BBE_ABSNX16(x0);
        sel = BBE_MOVVSELNX16(x0, 11);     // interval number (index in sine table)
        //x0 = BBE_ADDNX16(x0, x0); // 
        x0 = BBE_SLLINX16(x0, 1);

        // POLI_SU
        ofs = BBE_POLYNX16_OFF(x0, 12, 1); // offset in interval

        tmp_t = BBE_SHFLNX16(t1, sel);
        w0 = BBE_MULNX16(tmp_t, ofs);
        w0 = BBE_ADDNX40(w0, z0);
        t = BBE_PACKQNX40(w0);

        tmp_t = BBE_SHFLNX16(t3, sel);
        y0 = BBE_ADDNX16(tmp_t, t);

        t = BBE_MULNX16PACKQ(y0, ofs);

        tmp_t = BBE_SHFLNX16(t4, sel);
        y0 = BBE_ADDSNX16(tmp_t, t);

        //BBE_LVNX16_IP(SrcX, pX1, 2 * BBE_SIMD_WIDTH);
        //SrcX = BBE_ADDNX16(SrcX, ArgShift); // converse cosine argument to sine argument
        //neg = BBE_LTNX16(SrcX, 0);
        BBE_NEGNX16T(y0, y0, neg);
        //y0 = BBE_MULSGNNX16(SrcX, y0); // correct sign

        BBE_SVNX16_XP(y0, pY, 2 * BBE_SIMD_WIDTH);
    }
} /* vcos() */
