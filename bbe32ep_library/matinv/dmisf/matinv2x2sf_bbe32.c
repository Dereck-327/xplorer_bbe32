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
    Direct inversion of 2x2 floating point matrices, streaming data
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
Direct Matrix Inversion For Real Matrices

Description: perform in-place inversion of 2x2 and 4x4 real matrices. 
2x2 matrices are inverted by Cramer's rule. For 3x3, 4x4 matrices we employ 
the blockwise inversion algorithm encompassed with a suboptimal row/column
permutation that gains better conditioning of the block structure.

Data format and order options:
  Suffix   Data Order                Data Format   
    nf       Block     IEEE-754 Std single precision floating-point
    sf       Stream    IEEE-754 Std single precision floating-point

Note:
In general, accuracy of a matrix inversion algorithm implementation is a
function of input matrix condition number. Thus it is user's responsibility
to qualify the reliability of numeric results. Refer to NatureDSP Baseband 
Library Reference for details. 

Parameters:
Temporary:
pScr      Scratch memory area. To determine the scratch area size required by
          a function <fun>, use the respective helper function 
          <fun>_getScratchSize()
Input:
L         Number of matrices
Input/Output:
A[L][SA]  Sequence of L NxN input/result matrices. SA is the number
          of data elements occupied by a single NxN matrix in a block 
          (stream) ordered sequence, see function specifications.

Restrictions:
pScr,A    Aligned on 32-byte boundary

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/
/* Stream Order, Floating-Point, 2x2, SA=4
   Restrictions: 
     L must be a multiple of 8
*/
void matinv2x2sf(void * restrict pScr, float32_t* restrict A, int L)
{
    matinv2x2f_inv (pScr, A, L, e2x2_stream);
}

/* Return the scratch area size, in bytes. */
size_t matinv2x2sf_getScratchSize ( int L )
{
    return matinv2x2f_inv_getScratchSize(L);
}
#else
DISCARD_FUN(void, matinv2x2sf,(void * restrict pScr, float32_t* restrict A, int L))
size_t matinv2x2sf_getScratchSize ( int L )
{
  return (0);
}
#endif
