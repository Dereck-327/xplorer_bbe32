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
    Direct inversion of 4x4 floating point matrices 
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
#include "matinv4x4_common.h"
#include "matinv4x4Tbl.h"
#include "matinv4x4sf_inv.h"

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
void matinv4x4nf(void * restrict pScr, float32_t* restrict _A, int L)
{
    float32_t *X;
    int16_t *permIx;
    int L0;
    NASSERT_ALIGN(_A, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr, (2 * BBE_SIMD_WIDTH));
    if (L<=0) return;
    L0=L&~(BBE_SIMD_WIDTH/2-1);

    /* allocate data in the scratch */
    {
        X=(float32_t*)pScr;
        permIx=(int16_t*)(X+16*(BBE_SIMD_WIDTH/2));
        pScr=permIx+L;
        pScr=(void*)((((uintptr_t)pScr)+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1));
        {
            /* clean X to prevent ferret warnings */
            int k;
            xb_vecN_2xf32 *pX=(xb_vecN_2xf32 *)X;
            for (k=0; k<16; k++) BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pX,2*BBE_SIMD_WIDTH);
        }
    }
    matinv4x4nf_findPerm(permIx,_A,L);
    if (L0)
    {
        matinv4x4f_permute(_A,_A,permIx,matinv4x4_fwd_perm_tbl_bbe32,L0,e4x4_block);
        matinv4x4f_rbs4x4_inplace(_A,L0,e4x4_block);
        matinv4x4f_inv (pScr, _A,L0,e4x4_block);
        matinv4x4f_rsb4x4_inplace(_A,L0,e4x4_block);
        matinv4x4f_permute(_A,_A,permIx,matinv4x4_bkw_perm_tbl_bbe32,L0,e4x4_block);
    }
    if (L!=L0)
    {   /* process remainder */
        _A+=16*L0;
        permIx+=L0;
        L-=L0;
        matinv4x4f_permute(X,_A,permIx,matinv4x4_fwd_perm_tbl_bbe32,L,e4x4_block);
        matinv4x4f_rbs4x4_inplace(X,(L+BBE_SIMD_WIDTH/2-1)&~(BBE_SIMD_WIDTH/2-1),e4x4_block);
        matinv4x4f_inv (pScr, X,BBE_SIMD_WIDTH/2,e4x4_block);
        matinv4x4f_rsb4x4_inplace(X,(L+BBE_SIMD_WIDTH/2-1)&~(BBE_SIMD_WIDTH/2-1),e4x4_block);
        matinv4x4f_permute(_A,X,permIx,matinv4x4_bkw_perm_tbl_bbe32,L,e4x4_block);
    }
}

/* Return the scratch area size, in bytes. */
size_t matinv4x4nf_getScratchSize ( int L )
{
    size_t sz=0;
    int L0;
    L0=L&~(BBE_SIMD_WIDTH/2-1);
    if (L<=0) return 0;
    sz+=16*(BBE_SIMD_WIDTH/2)*sizeof(float32_t); /* for X part */
    sz+=L*sizeof(int16_t);  /* permIx */
    sz=(sz+2*BBE_SIMD_WIDTH-1)&~(2*BBE_SIMD_WIDTH-1);
    sz+=matinv4x4f_inv_getScratchSize(XT_MAX(BBE_SIMD_WIDTH/2,L0)); /* for matinv4x4f_inv */
    return sz;
}

#else
DISCARD_FUN(void, matinv4x4nf,(void * restrict pScr, float32_t* restrict A, int L))
size_t matinv4x4nf_getScratchSize ( int L )
{
  return (0);
}
#endif
