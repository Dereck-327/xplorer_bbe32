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
  NatureDSP_Baseband library. Math functions
    Arctangent
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

#include <errno.h>
/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Inline functions for floating-point exceptions and environment control. */
#include "__fenv.h"
/* Math Functions. */
#include "NatureDSP_Baseband_math.h"
/* pi, pi/2 values */
#include "pif_tbl.h"
/* atan polynomial coeff table */
#include "atanf_tbl.h"
#include "inff_tbl.h"
/*-------------------------------------------------------------------------
Arctangent 

Description: These functions compute the principal value of arctangent.

Representation:
vatan16,satan16  16-bit signed fixed-point format
                 Input data are Q15. Functions compute atan(x)/(pi/4)
                 and output results in Q15 format.
vatanf,satanf    IEEE-754 Std. single precision floating-point format
vfastatanf       Functions compute atan(x) and output results in radians

Special cases:
    Input | Result 
   -------+--------
    +inf  |  pi/2  (floating-point functions)
    -inf  | -pi/2  

Accuracy:
1 LSB for fixed point fixed point functions
1 ULP for vatanf(), satanf()
2 ULP for vfastatanf()

Notes:
1. These functions are much faster than full-quadrant arctangent atan2,
   so they are preferable when the full phase is not required.
2. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.

Input domain for 'fast' version vfastatanf():
|x|<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
z[N]  Results

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 16 (vatan16) or 8 (vatanf,vfastatanf)
-------------------------------------------------------------------------*/
#if HAVE_VFPU
float32_t satanf ( float32_t x )
{
  float32_t zero, one, half;
  float32_t y, x2, z;
  vbool1 big, verybig;
  int lesshalf;
  const union ufloat32uint32* p;
  xb_int32v ux;

  xb_int32v SCF; /* Floating-point Status and Control Register values. */
  
  /* check for NaN on input */  
  if (BBE_MOVAB1(XT_UN_S(x, x))) 
  {
    __Pragma("frequency_hint never");
    errno = EDOM; 
    x = XT_ADD_S(x, x);
    return x;
  }

  SCF = BBE_MOVVSCF(); /* Sample floating-point exception flags. */

  /* range reduction */
  zero =(float32_t)XT_CONST_S(0);
  one  =(float32_t)XT_CONST_S(1);
  half =(float32_t)XT_CONST_S(3);

  ux = BBE_MOV32_FROMF32(x); 
  ux = BBE_OPERATOR_AND32(ux, 0x80000000);/* take sign */
  x=XT_ABS_S(x);
  big=XT_OLT_S(one,x);

  if (BBE_MOVAB1(big))
  {
    verybig=XT_OLE_S(realmaxf.f,x);
    XT_MOVT_S(x,realmaxf.f,verybig);
    x = BBE_MOVAB1(verybig) ? zero : (float32_t)XT_RECIP_S(x);
  }

  lesshalf = BBE_MOVAB1(XT_OLT_S(x, half));
  /* approximate atan(x)/x-1 */

  p=(lesshalf) ? atanftbl1a:atanftbl2a;
  z=x;
  x = (lesshalf) ? x: x - 1.f ;
  /* Apply a combination of Estrin's rule and Horner's method */
  {
    float32_t cf0,cf1,cf2,cf3;

    cf0 = p[1].f; XT_MADD_S( cf0, x, p[0].f );
    cf1 = p[3].f; XT_MADD_S( cf1, x, p[2].f );
    cf2 = p[5].f; XT_MADD_S( cf2, x, p[4].f );
    cf3 = p[7].f; XT_MADD_S( cf3, x, p[6].f );

    x2 = XT_MUL_S(x, x);    

    y = cf0;
    XT_MADD_S( cf1, y, x2 ); y = cf1;
    XT_MADD_S( cf2, y, x2 ); y = cf2;
    XT_MADD_S( cf3, y, x2 ); y = cf3;
  }
  x=z;
  /* convert result to true atan(x) */
  XT_MADD_S(x,x,y);
  x2 = XT_SUB_S(pi2m1f.f, x);
  x2 = XT_ADD_S(1.f, x2);
  XT_MOVT_S(x, x2, big);
  y=BBE_MOVF32_FROM32(BBE_OPERATOR_XOR32(BBE_MOV32_FROMF32(x),ux)); /* apply sign */

  BBE_MOVSCFV(SCF);
  return y;

} /* satanf() */
#else
DISCARD_FUN(float32_t,satanf,( float32_t x ))
#endif
