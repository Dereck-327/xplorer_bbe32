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

#ifndef PFIT_COMMON_H__
#define PFIT_COMMON_H__

#define RSTRIDE 16

#define TAKEHILO3(A,hi,lo)  \
{   vsaN _16;               \
    _16=BBE_MOVVSA32(16);   \
    lo=BBE_PACKLNX40(A);    \
    hi=BBE_PACKVNX40(A,_16);\
}

#define BBE_MULNX16PACKQ_SAT(Z,X,Y) \
{                                   \
    xb_vecNx40 T;                   \
    T=BBE_MULNX16(X,Y);             \
    Z=BBE_PACKQNX40(T);             \
}

// forward declarations
void pfit_chol
(
  int16_t * restrict        t,
  int16_t * restrict        R, 
  const int16_t * restrict  A, 
  int32_t                   sigma2,
  int                       M,
  int                       N
);

typedef void (*fndiscr)
(
  int16_t * restrict        f,
  const int16_t * restrict  b,
  const int16_t * restrict  A,
  const int32_t * restrict  x,
  int                       M,
  int                       N
);

void pfit_process
(
  void * restrict           pScr,
  int32_t * restrict        p,
  const int16_t * restrict  V,
  const int16_t * restrict  R,
  const int16_t * restrict  y,
  int                       M,
  int                       N,
  int                       maxIter,
  fndiscr                   discr
);

#endif
