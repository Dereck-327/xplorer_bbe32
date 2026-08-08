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
  NatureDSP_Baseband library. Fitting and Interpolation Routines
    Polynomial Fitting and Interpolation for Real Data
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fit.h"
#include "pfit_common.h"

/*-------------------------------------------------------------------------
Polynomial Fitting and Interpolation

Description: the pfit functions fit (in least squares sense) a degree N 
polynomial to input data sampled at a points grid of length M, and use
that polynomial to interpolate data at arbitrary query points. Namely,
the pfitN_grid() functions compute the Vandermonde matrix for the sample
points grid and perform the Cholesky decomposition of that matrix, the
pfitN_process() functions calculate the least squares solution for the
polynomial coefficients. Finally, the pfitN_eval() functions evaluate
the polynomial at query points.

Please refer to the NatureDSP Baseband Library Reference for full details
on these functions.

Representation:
pfit_gridN,      16-bit fixed-point data. Parameter specifications denote
pfit_processN,   fixed-point format for various data items
pfit_evalN       
pfitf_gridN,     IEEE-754 Std single precision floating-point data
pfitf_processN,  
pfitf_evalN     

Note:
Number of fractional bits specidied for various input/output arguments below apply
for the fixed-point variant

Parameters:
Input:
N                     Degree of polynomial, 1..6
M                     Number of sample points
P                     Number of query points
maxIter               Number of least squares solution enhancement iterations. Right 
                      choice depends on required accuracy, the ad-hoc value is (N+1)/2
x[M]                  Sample points grid, Q15 or floating point
y[M]                  Sampled data values, Q15 or floating point
xi[P]                 Query points, Q15 or floating point
M'=(M+7)&(~7), N'=8   for floating point API
M'=(M+15)&(~15),N'=16 for fixed-point API

Intermediate:
V[M'*8]               Vandermonde matrix, Q15 or floating point
R[N'*8]               Upper triangular Cholesky factor of matrix V, Q11 or floating point
Output:
yi[P]                 Data values interpolated at query points, Q15 or floating point
p[N+1]                Polynomial coefficients, Q8.23 or floating point
Temporary:
pScr                  Scratch memory area. To determine the scratch area size required by
                      a function pfitN_<fun>, use the respective helper function 
                      pfit_<fun>_getScratchSize(M,N)

Restrictions:
M>N                   The number of sample points must exceed the degree of polynomial
x,y,xi,yi,V,R,p,pScr  Must not overlap
V,R,pScr              Aligned on 32-byte boundary
---------------------------------------------------------------------------*/
void pfit_eval4(int16_t * restrict yi, const int16_t * restrict xi, const int32_t * restrict p, int P)
{
#if (HAVE_ADVPMUL)
    int n;
    xb_vecNx16 yp, xp, lo, hi;
    xb_mvecNx32  mvec;
    valign apoly;
    xb_vecNx16 vc0, t, zero;
    xb_vecNx40 w;
    xb_vecNx40 Am, Ap;
    valign vx, vy;
    vsaN shr;
    int bytecount;
    int rbytecount;

    const xb_vecNx16 * ppoly = (const xb_vecNx16*)p;
    const xb_vecNx16 * restrict pxi = (const xb_vecNx16 *)xi;
    xb_vecNx16 * restrict pyi = (xb_vecNx16 *)yi;
    vx = BBE_LA_PP(pxi);
    shr = BBE_MOVVSA32(4);

    rbytecount = bytecount = P << 1;
    rbytecount += (P >> 31);  // add zero !!!! to disable compiler's remapping of bytecount to bytecount
    vy = BBE_ZALIGN();
    apoly = BBE_LA_PP(ppoly);
    BBE_LAVNX16_XP(vc0, apoly, ppoly, 5 * 4); // 5 int32_t coefficients
    t = BBE_SEQNX16(); shr = BBE_MOVVSV(t, 0); // shift them by 0,1...4 bits right
    w = BBE_MOVSWV(vc0, vc0);
    w = BBE_SRANX40(w, shr);
    vc0 = BBE_MOVSVWL(w);
    zero = BBE_ZERONX16();

    shr = BBE_MOVVSA32(4);


    BBE_LAVNX16_XP(xp, vx, pxi, rbytecount);
    lo = BBE_REPNX16(vc0, 0); // p[0], low part
    hi = BBE_REPNX16(vc0, 1); // p[0], high part
    Am = BBE_MULUSNX16(lo, xp);
    mvec = BBE_MULMNX16(hi, xp);
    BBE_SRAIWADDMNX40(Am, mvec, 16);
    t = BBE_REPNX16C(vc0, 1); Ap = BBE_MOVSWV(t, t); // p[1]>>1
    Am = BBE_ADDNX40(Am, Ap);

    lo = BBE_MOVSVWXL(Am);
    hi = BBE_MOVSVWXH(Am);
    Am = BBE_MULUSNX16(lo, xp);
    mvec = BBE_MULMNX16(hi, xp);
    BBE_SRAIWADDMNX40(Am, mvec, 16);
    t = BBE_REPNX16C(vc0, 2); Ap = BBE_MOVSWV(t, t); // p[2]>>2
    Am = BBE_ADDNX40(Am, Ap);

    lo = BBE_MOVSVWXL(Am);
    hi = BBE_MOVSVWXH(Am);
    Am = BBE_MULUSNX16(lo, xp);
    mvec = BBE_MULMNX16(hi, xp);
    BBE_SRAIWADDMNX40(Am, mvec, 16);
    t = BBE_REPNX16C(vc0, 3); Ap = BBE_MOVSWV(t, t); // p[3]>>3
    Am = BBE_ADDNX40(Am, Ap);

    lo = BBE_MOVSVWXL(Am);
    hi = BBE_MOVSVWXH(Am);
    Am = BBE_MULUSNX16(lo, xp);
    mvec = BBE_MULMNX16(hi, xp);
    BBE_SRAIWADDMNX40(Am, mvec, 16);
    t = BBE_REPNX16C(vc0, 4); Ap = BBE_MOVSWV(t, t); // p[4]>>4
    Am = BBE_ADDNX40(Am, Ap);
    Am = BBE_RNDSADJNX40(Am, shr);
    yp = BBE_PACKVNX40(Am, shr);

    __Pragma("loop_count min=1");
    for (n = 0; n<((P + BBE_SIMD_WIDTH - 1) >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        BBE_SAVNX16_XP(yp, vy, pyi, bytecount);
        bytecount = XT_ADDX8(-(BBE_SIMD_WIDTH / 4), bytecount);
        rbytecount = XT_ADDX8(-(BBE_SIMD_WIDTH / 4), rbytecount);
        
        BBE_LAVNX16_XP(xp, vx, pxi, rbytecount);
        lo = BBE_REPNX16(vc0, 0); // p[0], low part
        hi = BBE_REPNX16(vc0, 1); // p[0], high part
        Am = BBE_MULUSNX16(lo, xp);
        mvec = BBE_MULMNX16(hi, xp);
        BBE_SRAIWADDMNX40(Am, mvec, 16);
        vc0 = BBE_ADDNX16(vc0, zero); // improves scheduling !!!
        t = BBE_REPNX16C(vc0, 1); Ap = BBE_MOVSWV(t, t); // p[1]>>1
        Am = BBE_ADDNX40(Am, Ap);

        lo = BBE_MOVSVWXL(Am);
        hi = BBE_MOVSVWXH(Am);
        Am = BBE_MULUSNX16(lo, xp);
        mvec = BBE_MULMNX16(hi, xp);
        BBE_SRAIWADDMNX40(Am, mvec, 16);
        vc0 = BBE_ADDNX16(vc0, zero); // improves scheduling !!!
        t = BBE_REPNX16C(vc0, 2); Ap = BBE_MOVSWV(t, t); // p[2]>>2
        Am = BBE_ADDNX40(Am, Ap);

        lo = BBE_MOVSVWXL(Am);
        hi = BBE_MOVSVWXH(Am);
        Am = BBE_MULUSNX16(lo, xp);
        mvec = BBE_MULMNX16(hi, xp);
        BBE_SRAIWADDMNX40(Am, mvec, 16);
        vc0 = BBE_ADDNX16(vc0, zero); // improves scheduling !!!
        t = BBE_REPNX16C(vc0, 3); Ap = BBE_MOVSWV(t, t); // p[3]>>3
        Am = BBE_ADDNX40(Am, Ap);

        lo = BBE_MOVSVWXL(Am);
        hi = BBE_MOVSVWXH(Am);
        Am = BBE_MULUSNX16(lo, xp);
        mvec = BBE_MULMNX16(hi, xp);
        BBE_SRAIWADDMNX40(Am, mvec, 16);
        vc0 = BBE_ADDNX16(vc0, zero); // improves scheduling !!!
        t = BBE_REPNX16C(vc0, 4); Ap = BBE_MOVSWV(t, t); // p[4]>>4
        Am = BBE_ADDNX40(Am, Ap);
        Am = BBE_RNDSADJNX40(Am, shr);
        yp = BBE_PACKVNX40(Am, shr);
    }
    BBE_SAPOS_FP(vy, pyi);
#else
  int n;
  xb_vecNx16 yp, xp, lo, hi;
  xb_vecNx40 Am, Ap;
  valign vx, vy;
  vsaN _8;
  int bytecount;
  int rbytecount;
  const xb_vecNx16 * restrict pxi = (const xb_vecNx16 *)xi;
  xb_vecNx16 * restrict pyi = (xb_vecNx16 *)yi;

  vx = BBE_LA_PP(pxi);
  _8 = BBE_MOVVSA32(8);

  rbytecount = bytecount = P << 1;
  rbytecount += (P >> 31);  // add zero !!!! to disable compiler's remapping of bytecount to bytecount
  vy = BBE_ZALIGN();
  #define ITERATION(n)        \
  {                           \
      Ap=BBE_MOVWA32(p[n]);\
      TAKEHILO3(Am,hi,lo);    \
      Am=BBE_MULUSNX16(lo,xp);\
      Am=BBE_SRAINX40(Am,16); \
      BBE_MULANX16(Am,hi,xp); \
      Am=BBE_SLLINX40(Am,1);  \
      Am=BBE_ADDNX40(Am,Ap);  \
  }

  BBE_LAVNX16_XP(xp, vx, pxi, rbytecount);
  for (n = 0; n<((P + 15) >> 5); n++)
  {
    lo = BBE_MOVVA16C(p[0]);
    hi = BBE_SHFLNX16I(lo, BBE_SHFLI_DUPLICATE_1_ODD);
    lo = BBE_SHFLNX16I(lo, BBE_SHFLI_DUPLICATE_1_EVEN);

    Ap = BBE_MOVWA32(p[1]);
    Am = BBE_MULUSNX16(lo, xp);
    Am = BBE_SRAINX40(Am, 16);
    BBE_MULANX16(Am, hi, xp);
    Am = BBE_SLLINX40(Am, 1);
    Am = BBE_ADDNX40(Am, Ap);

    ITERATION(2)
    ITERATION(3)
    ITERATION(4)

    Am = BBE_RNDSADJNX40(Am, _8);
    yp = BBE_PACKVNX40(Am, _8);
    rbytecount = XT_ADDX8(-4, rbytecount);
    BBE_LAVNX16_XP(xp, vx, pxi, rbytecount);
    BBE_SAVNX16_XP(yp, vy, pyi, bytecount);
    bytecount = XT_ADDX8(-4, bytecount);

    lo = BBE_MOVVA16C(p[0]);
    hi = BBE_SHFLNX16I(lo, BBE_SHFLI_DUPLICATE_1_ODD);
    lo = BBE_SHFLNX16I(lo, BBE_SHFLI_DUPLICATE_1_EVEN);

    Ap = BBE_MOVWA32(p[1]);
    Am = BBE_MULUSNX16(lo, xp);
    Am = BBE_SRAINX40(Am, 16);
    BBE_MULANX16(Am, hi, xp);
    Am = BBE_SLLINX40(Am, 1);
    Am = BBE_ADDNX40(Am, Ap);

    ITERATION(2)
    ITERATION(3)
    ITERATION(4)

    Am = BBE_RNDSADJNX40(Am, _8);
    yp = BBE_PACKVNX40(Am, _8);
    rbytecount = XT_ADDX8(-4, rbytecount);
    BBE_LAVNX16_XP(xp, vx, pxi, rbytecount);
    BBE_SAVNX16_XP(yp, vy, pyi, bytecount);
    bytecount = XT_ADDX8(-4, bytecount);
  }

  if ((P + 15) >> 4)
  {
    lo = BBE_MOVVA16C(p[0]);
    hi = BBE_SHFLNX16I(lo, BBE_SHFLI_DUPLICATE_1_ODD);
    lo = BBE_SHFLNX16I(lo, BBE_SHFLI_DUPLICATE_1_EVEN);

    Ap = BBE_MOVWA32(p[1]);
    Am = BBE_MULUSNX16(lo, xp);
    Am = BBE_SRAINX40(Am, 16);
    BBE_MULANX16(Am, hi, xp);
    Am = BBE_SLLINX40(Am, 1);
    Am = BBE_ADDNX40(Am, Ap);

    ITERATION(2)
    ITERATION(3)
    ITERATION(4)

    Am = BBE_RNDSADJNX40(Am, _8);
    yp = BBE_PACKVNX40(Am, _8);
    rbytecount = XT_ADDX8(-4, rbytecount);
    BBE_LAVNX16_XP(xp, vx, pxi, rbytecount);
    BBE_SAVNX16_XP(yp, vy, pyi, bytecount);
    bytecount = XT_ADDX8(-4, bytecount);
  }

  BBE_SAPOS_FP(vy, pyi);

  #undef ITERATION
#endif
} /* pfit_eval4() */
