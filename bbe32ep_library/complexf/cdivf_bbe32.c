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
    Division of Complex Numbers
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
DISCARD_FUN(complex_float, cdivf, ( complex_float x, complex_float y ))
#else

/*-------------------------------------------------------------------------
Division of Complex Numbers

Representation:
vdividecx         16-bit signed fixed-point format (Q15)
vdividecxf,cdivf  IEEE-754 Std. single precision floating-point format

Fixed-point function returns the fractional and exponential portion of real
and imaginary components of the division result. Fixed-point format for the
fractional part is 16-bit Q(15-exp), where exp denotes the respective 
exponential value. Full division result can be restored in 48-bit Q31 format
by sign extending the fractional part to 64 bits and shifting it to the left
by 16+exp bit positions

Special cases:
      x   |    y    |  Result |  Extra Conditions    
  --------|---------|---------|---------------------
    +/-0  |  +/-0   |   NaN   |
     x    | +/-inf  |    0    | x is a finite number          (floating-point functions)
   +/-inf |   y     | +/-inf  | y is a finite number (see *)
   +/-inf | +/-inf  |   NaN   |
  --------|---------|---------|---------------------
     0    |   0     |   not   |                               (fixed-point functions)
          |         | defined |

* - quotent of infinite number to finite number may have NaN in its real or 
    imaginary part depending on signs of input arguments

Input domain for 'fast' version vfastdividecxf is limited by usage of direct 
formula without normalization of inputs, so input domain is:
|real(x)*real(y)|<Inf
|imag(x)*imag(y)|<Inf
|real(y)*real(y)|<Inf
|imag(y)*imag(y)|<Inf
1.1755e-038 < |real(x)*real(y) + imag(x)*imag(y)| < Inf
1.1755e-038 < |real(y)*real(y) + imag(y)*imag(y)| < Inf
The output value is not defined outside of this range or accuracy is degraded.

Accuracy:
vdividecx         1.5 LSB of the fractional part
vdividecxf,cdivf  2 ULP
vfastdividecxf    2 ULP

Input:
x[N]      Input vector of dividends
y[N]      Input vector of divisors
N         Length of vectors
Output:
vdividecx:
fract[N]  Fractional part of quotients, Q(15-exp); if non-zero, then
            8192<=|fract|<32768
exp[2*N]  Exponential part of quotients, -14..16. Exponential values for
          real and imaginary components are interleaved, with exponential
          parts of real components going first (at even indices)
vdividecxf,cdivf
z[N]      Quotients
Restrictions:
x,y,z,fract,exp   Aligned on 32-byte boundary
x,y,z,fract,exp   Must not overlap
N                 Multiple of 8 (vdividecx) or 4 (vdividecxf, vfastdividecxf)
---------------------------------------------------------------------------*/

complex_float cdivf ( complex_float x, complex_float y )
{
    float32_t A,B,C,D,XMAX,YMAX,S0,S1,CC,DD,V,W,DEN;
    float32_t absA,absB,absC,absD;
    float32_t zero,one,plusInf,maxDenorm,pow2_23;
    vbool1 FORCEZERO;
    int  EX,EY,SY,SY0,SY1;
    
    union
    {
        struct 
        {
            float32_t re, im;
        } s;
        complex_float z;
    } tempX, tempY;/* used to extract real and imaginary parts */

    zero=XT_CONST_S(0);
    one =XT_CONST_S(1);
    plusInf  =XT_WFR(0x7F800000);
    maxDenorm=XT_WFR(0x007FFFFF);
    pow2_23  =XT_WFR((127+23)<<23);
    /* load, select real/imaginary parts and compute rescale factor */
    tempX.z = x;
    tempY.z = y;
    A=tempX.s.re;
    B=tempX.s.im;
    C=tempY.s.re;
    D=tempY.s.im;

    absA=XT_ABS_S(A);
    absB=XT_ABS_S(B);
    absC=XT_ABS_S(C);
    absD=XT_ABS_S(D);

    XMAX=XT_MAX_S(absA,absB);
    XT_MOVT_S(XMAX,absA,XT_UN_S(absA,absA));
    YMAX=XT_MAX_S(absC,absD);
    XT_MOVT_S(YMAX,absC,XT_UN_S(absC,absC));

    EX = XT_RFR(XMAX)>>23;
    EY = XT_RFR(YMAX)>>23;
    SY = EY+XT_MAX(127-1,XT_MAX(EX,EY));
    SY = ((127-1)*3+(127)*2-SY);

    FORCEZERO= BBE_ANDNOTB1(XT_OEQ_S(YMAX,plusInf), XT_UEQ_S(XMAX,plusInf));

    /* scale c and d by 2^sy */
    SY0=SY>>1;
    SY1=SY-SY0;
    S0=XT_WFR(SY0<<23);
    S1=XT_WFR(SY1<<23);
    CC = XT_MUL_S(C,S0);
    DD = XT_MUL_S(D,S0);
    CC = XT_MUL_S(CC,S1);
    DD = XT_MUL_S(DD,S1);
    
    /* prepare nominator and denominator for division */
    V=XT_ADD_S(XT_MUL_S(A,CC),XT_MUL_S(B,DD));
    W=XT_SUB_S(XT_MUL_S(B,CC),XT_MUL_S(A,DD));
    DEN=XT_ADD_S(XT_MUL_S(C,CC),XT_MUL_S(D,DD));
    XT_MOVT_S(V  ,zero,FORCEZERO);
    XT_MOVT_S(W  ,zero,FORCEZERO);
    XT_MOVT_S(DEN,one ,FORCEZERO);

    /* divide */
    {
    vbool1 DENDENORM,DENZERO;
    float32_t Y;
    DENDENORM = XT_OLE_S(DEN,maxDenorm);
    Y=XT_MUL_S(V  ,pow2_23); XT_MOVT_S(V  ,Y,DENDENORM);
    Y=XT_MUL_S(W  ,pow2_23); XT_MOVT_S(W  ,Y,DENDENORM);
    Y=XT_MUL_S(DEN,pow2_23); XT_MOVT_S(DEN,Y,DENDENORM);
    DENZERO=XT_OEQ_S(DEN,zero);
    /* divide v,w by den */
    Y=XT_RECIP_S(DEN);
    XT_MOVT_S(Y,plusInf,DENZERO);
    V=XT_MUL_S(Y,V);
    W=XT_MUL_S(Y,W);
    }

    tempX.s.re = V;
    tempX.s.im = W;
    x = tempX.z;
    return x;
} /* cdivf() */

#endif/* !HAVE_VFPU */
