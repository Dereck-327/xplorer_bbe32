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
    Direct inversion of 3x3 floating point matrices, streaming data
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Fixed point arithmetics. */
#include "NatureDSP_Math.h"
/* Common utility declarations. */
#include "common.h"
#include <math.h>
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
#include "matinv3x3Tbl.h"
#include "matinv3x3_common.h"
#include "matinv3x3sf_inv.h"
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
void matinv3x3sf(void * restrict pScr, float32_t* restrict A, int L)
{
    int16_t* permIx;
    NASSERT_ALIGN(A, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0) return;

    // allocate data in the scratch
    {
        permIx=(int16_t*)(pScr);
        pScr=permIx+L;
        pScr=(void*)((((uintptr_t)pScr)+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1));
    }

    matinv3x3f_rsb3x3_inplace(A,L);
    matinv3x3f_search(permIx,A,L,e3x3_stream);
    matinv3x3f_permute(A,A,permIx,matinv3x3_fwd_perm_tbl_bbe32,L,e3x3_stream);
    matinv3x3f_rsb3x3_inplace(A,L);

    matinv3x3f_inv(pScr, A, L,e3x3_stream);
    matinv3x3f_rsb3x3_inplace(A,L);
    matinv3x3f_permute(A,A,permIx,matinv3x3_bkw_perm_tbl_bbe32,L,e3x3_stream);
    matinv3x3f_rsb3x3_inplace(A,L);

}

/* Return the scratch area size, in bytes. */
size_t matinv3x3sf_getScratchSize ( int L )
{
    size_t sz=0;
    sz+=L*sizeof(int16_t);  /* permIx */
    sz=(sz+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    sz+=matinv3x3f_inv_getScratchSize(L);
    return sz;
}
#else
DISCARD_FUN(void, matinv3x3sf,(void * restrict pScr, float32_t* restrict A, int L))
size_t matinv3x3sf_getScratchSize ( int L )
{
  return (0);
}
#endif
