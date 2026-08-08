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
    NatureDSP_Baseband library. Direct Matrix Inversion
    Direct inversion of 2x2 complex floating point matrices, streaming data
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"

#include "matinv2x2sf_inv.h"

#if HAVE_VFPU

/*-------------------------------------------------------------------------
Direct Matrix Inversion For Complex Matrices

Description: perform in-place inversion of 2x2 and 4x4 complex matrices. 
2x2 matrices are inverted by Cramer's rule. For 3x3, 4x4 matrices we employ 
the blockwise inversion algorithm encompassed with a suboptimal row/column
permutation that gains better conditioning of the block structure.

Data format and order options:
  Suffix   Data Order                 Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Notes:
1. In general, accuracy of a matrix inversion algorithm implementation is a
   function of input matrix condition number. Thus it is user's responsibility
   to qualify the reliability of numeric results. Refer to NatureDSP Baseband 
   Library Reference for details. 
2. For blockwise inversion of 4x4 fixed-point matrices, it is reasonable to
   limit the dynamic range of input data by 11..13 significant bits. This
   measure reduces the possibility of an overflow at internal computations.

Parameters:
Temporary:
pScr      Scratch memory area. To determine the scratch area size required by
          a function <fun>, use the respective helper function 
          <fun>_getScratchSize()
Input:
L         Number of matrices
qA        Number of fractional bits for fixed-point input/output data
Input/Output:
A[L][SA]  Sequence of L NxN complex input/result matrices. SA is the number
          of data elements occupied by a single NxN matrix in a block 
          (stream) ordered sequence, see function specifications.

Restrictions:
pScr,A    Aligned on 32-byte boundary

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/
void cmatinv2x2sf(void * restrict pScr, complex_float * restrict A, int L)
#if 0
{
    int l;
    NASSERT_ALIGN(A, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    for (l=0; l<L; l++,A++)
    {
        complex_float a,b,c,d,det,rdet;
        a=A[0*L]; b=A[1*L]; c=A[2*L]; d=A[3*L];
        det=subc(mulc(a,d),mulc(b,c));
        rdet=recipc(det);
        A[0*L]=mulc( d,rdet);      A[1*L]=negc(mulc(b,rdet)); 
        A[2*L]=negc(mulc(c,rdet)); A[3*L]= mulc(a,rdet); 
    }
}
#else
{
    cmatinv2x2f_inv(pScr,A,L,e2x2_stream);
}
#endif

/* Return the scratch area size, in bytes. */
size_t cmatinv2x2sf_getScratchSize ( int L )
{
  return cmatinv2x2f_inv_getScratchSize(L);
}
#else
DISCARD_FUN(void, cmatinv2x2sf,(void * pScr, complex_float * A, int L))
size_t cmatinv2x2sf_getScratchSize ( int L )
{
  return (0);
}
#endif
