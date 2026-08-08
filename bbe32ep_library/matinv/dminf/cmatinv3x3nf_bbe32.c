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
    Direct inversion of 3x3 complex floating point matrices 
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
#include "matinv3x3_common.h"
#include "matinv3x3sf_inv.h"
#include "matinv3x3Tbl.h"

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
void cmatinv3x3nf(void * restrict pScr, complex_float * restrict _A, int L)
{
    complex_float *X;
    int16_t *permIx;
    int L0;
    NASSERT_ALIGN(_A, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr, (2 * BBE_SIMD_WIDTH));
    if (L<=0) return;
    L0=(L)&~(BBE_SIMD_WIDTH/4-1);

    /* allocate data in the scratch */
    {
        X=(complex_float*)pScr;
        permIx=(int16_t*)(X+12*BBE_SIMD_WIDTH/4);
        pScr=permIx+L;
        pScr=(void*)((((uintptr_t)pScr)+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1));
        {
            /* clean X to prevent ferret warnings */
            xb_vecN_2xf32 *pX=(xb_vecN_2xf32 *)X;
            BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pX,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pX,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pX,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pX,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pX,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pX,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pX,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pX,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pX,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pX,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pX,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pX,2*BBE_SIMD_WIDTH);
        }
    }

    cmatinv3x3f_search(permIx,_A,L,e3x3_block);
    if (L0)
    {
        cmatinv3x3f_permute(_A,_A,permIx,matinv3x3_fwd_perm_tbl_bbe32,L0,e3x3_block);
        cmatinv3x3f_cbs3x3_inplace(_A,L0,e3x3_block);
        cmatinv3x3f_inv(pScr,_A,L0,e3x3_block);
        cmatinv3x3f_csb3x3_inplace(_A,L0,e3x3_block);
        cmatinv3x3f_permute(_A,_A,permIx,matinv3x3_bkw_perm_tbl_bbe32,L0,e3x3_block);

    }
    if (L0!=L)
    {
        _A+=12*L0;
        permIx+=L0;
        L-=L0;
        cmatinv3x3f_permute(X,_A,permIx,matinv3x3_fwd_perm_tbl_bbe32,L,e3x3_block);
        cmatinv3x3f_cbs3x3_inplace(X,BBE_SIMD_WIDTH/4,e3x3_block);
        cmatinv3x3f_inv(pScr,X,BBE_SIMD_WIDTH/4,e3x3_block);
        cmatinv3x3f_csb3x3_inplace(X,BBE_SIMD_WIDTH/4,e3x3_block);
        cmatinv3x3f_permute(_A,X,permIx,matinv3x3_bkw_perm_tbl_bbe32,L,e3x3_block);
    }
}

/* Return the scratch area size, in bytes. */
size_t cmatinv3x3nf_getScratchSize ( int L )
{
    size_t sz=0;
    if (L<=0) return 0;
    L=(L+BBE_SIMD_WIDTH/4-1)&~(BBE_SIMD_WIDTH/4-1);
    sz+=12*(BBE_SIMD_WIDTH/4)*sizeof(complex_float); /* X */
    sz+=L*sizeof(int16_t);  /* permIx */
    sz=(sz+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    sz+=cmatinv3x3f_inv_getScratchSize(XT_MAX(BBE_SIMD_WIDTH/4,L)); /* for matinv3x3hfs */
    return sz;
    return (0);
}
#else
DISCARD_FUN(void, cmatinv3x3nf,(void * pScr, complex_float * A, int L))
size_t cmatinv3x3nf_getScratchSize ( int L )
{
    (void)L;
    return (0);
}
#endif
