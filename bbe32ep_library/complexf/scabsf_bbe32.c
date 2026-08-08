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
  NatureDSP_Baseband library. Complex Math functions
    Magnitude Of Complex Number
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Complex Math Functions. */
#include "NatureDSP_Baseband_complex.h"

#if !(HAVE_VFPU)
DISCARD_FUN(float32_t, scabsf, ( complex_float x ))
#else

/*-------------------------------------------------------------------------
Magnitude Of Complex Number

Description: These functions multiply complex input number by its conjugate
and take the square root of the result.

Representation:
vcabs          16-bit signed fixed-point format
               Number of fractional bits for input data Qx is user-defined.
               Fixed-point format for output data is Qz = Qx+sh/2, where sh
               is the shift control argument. If a resulting value is too
               large to be represented in Qz format, then it is saturated
               by 32767.
vcabsf,scabsf  IEEE-754 Std. single precision floating-point format


Input domain for 'fast' version vfastcabsf():
1.1755e-038 < |real(x)*real(x)+ imag(x)*imag(x)| < Inf 
The output value is not defined outside of this range or accuracy is degraded.

Accuracy:
1 LSB - vcabs
2 ULP - vcabsf,scabsf
3 ULP - vfastcabsf

Parameters:
Input:
x[N]  Complex numbers
sh    (vcabs only) bi-directional shift control, an even integer from
      the range [-2,30]. Positive value corresponds to a left shift
N     Length of vectors
Output:
z[N]  Magnitudes

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 16 (vcabs) or 8 (vcabsf,vfastcabsf)
-------------------------------------------------------------------------*/

float32_t scabsf ( complex_float x )
{
    float32_t X,X_re,X_im,Sp,Sm;
    int Exp_re,Exp_im,Exp_abs;

    union
    {
        struct 
        {
            float32_t re, im;
        } s;
        complex_float z;
    } temp;/* used to extract real and imaginary parts */
    temp.z = x;

    X_re=temp.s.re;
    X_im=temp.s.im;

    /* compute normalization factors */
    Exp_re = XT_RFR(X_re)>>23;
    Exp_im = XT_RFR(X_im)>>23;
    Exp_re = Exp_re & 0xFF;
    Exp_im = Exp_im & 0xFF;
    Exp_abs= XT_MAX(Exp_re,Exp_im);
    Exp_abs= XT_MIN(Exp_abs,256-3);
    Exp_abs= XT_MAX(Exp_abs,    1);
    Sp=XT_WFR((256-2-Exp_abs)<<23);
    Sm=XT_WFR((      Exp_abs)<<23);
    /* compute normalized result */
    X_re=XT_MUL_S(X_re,Sp);
    X_im=XT_MUL_S(X_im,Sp);
    X=XT_MUL_S(X_re,X_re);
    XT_MADD_S(X,X_im,X_im);
    /* square root and denormalization */
    {
        float32_t z, x_red, x_adj, r, r_err, t0, t1;
        float32_t zero, half;

        zero = XT_CONST_S(0);
        half = XT_CONST_S(3);

        x_red = XT_NEXP01_S(X);/* x with reduced exponent range */

        /* Initial rsqrt approximation with exponent range reduction */
        r = XT_SQRT0_S(X);
        /* compute approximation error */
        t0 = XT_MUL_S(r, r);
        t1 = x_red; XT_ADDEXP_S(t1, half);/* -0.5*x */
        r_err = half; XT_MADDN_S(r_err, t0, t1);/* approximation error is (0.5-0.5*x*r*r) */
        XT_MADDN_S(r, r, r_err);/* Second recip sqrt approximation */

        /* Compute reduced range sqrt approximation */
        z = zero;
        XT_MSUBN_S(z, x_red, r);/* z = x*rsqrt(x) */

        /* Make final adjustment and restore range */
        x_adj = XT_MKSADJ_S(X);
        XT_MADDN_S(x_red, z, z);
        XT_ADDEXPM_S(z, x_adj);
        t0 = zero;
        XT_MSUBN_S(t0, half, r);
        XT_ADDEXP_S(t0, x_adj);
        XT_DIVN_S(z, x_red, t0);
        X = z;
    }
    X=XT_MUL_S(X,Sm);

    return X;
} /* scabsf() */

#endif/* !HAVE_VFPU */
