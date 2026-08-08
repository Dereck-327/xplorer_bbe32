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
#ifndef MATINV4X4_COMMON_H__
#define MATINV4X4_COMMON_H__
/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
#include "matinvf_common.h"

/*-------------------------------------------------------------------------
    find permutation for L matrices written in block/stream order
    Input:
    X[L*16]  input matrices
    L        number of matrices
    Output:
    permIx[L] permutation indices
-------------------------------------------------------------------------*/
void  matinv4x4nf_findPerm(int16_t *permIx,const float32_t *X, int L);
void cmatinv4x4nf_findPerm(int16_t *permIx,const complex_float *X, int L);
void  matinv4x4sf_findPerm(int16_t *permIx,const float32_t *X, int L);
void cmatinv4x4sf_findPerm(int16_t *permIx,const complex_float *X, int L);

/*-------------------------------------------------------------------------
    find permutation for L matrices written in stream order
    and return permuted matrices
    Input:
    X[L*16]  input matrices
    L        number of matrices
    Output:
    permIx[L] permutation indices
-------------------------------------------------------------------------*/
void  matinv4x4sf_findAndPerm (int16_t *permIx,float32_t *X, int L);
void  cmatinv4x4sf_findAndPerm(int16_t *permIx,complex_float *X, int L);

/*-------------------------------------------------------------------------
    backward permutation with conversion to block/stream format
    Input:
    X[16*L']   4x4 matrices in the halfstream/stream order
    permIx[L]  permutation indices
    L          number of matrices
    L' smallest bigger integer which is multiple of BBE_SIMD_WIDTH/2 for
    real or BBE_SIMD_WIDTH/4 for complex data
    Output:
    Y[16*L]    4x4 matrices in the block/stream order

    NOTE: X might be damaged 
    Restrictions:
    all matrices should be aligned
    L>0
-------------------------------------------------------------------------*/
void matinv4x4nf_permbkw(float32_t * restrict Y,
                         float32_t * restrict X,
                         const int16_t *permIx,
                         int L);
void cmatinv4x4nf_permbkw(complex_float * restrict Y,
                          complex_float * restrict X,
                          const int16_t *permIx,
                          int L);
void matinv4x4sf_permbkw(
                         float32_t * restrict X,
                         const int16_t *permIx,
                         int L);
void cmatinv4x4sf_permbkw(
                          complex_float * restrict X,
                          const int16_t *permIx,
                          int L);
/*-------------------------------------------------------------------------
    forward permutation with conversion to halfstream/stream format
    Input:
    X[16*L]    input 4x4 matrices in the block/stream order
    permIx[L]  permutation indices
    L          number of matrices
    L' smallest bigger integer which is multiple of BBE_SIMD_WIDTH/2 for
    real or BBE_SIMD_WIDTH/4 for complex data
    Output:
    Y[16*L']    4x4 matrices in halfstream/stream order

    Restrictions:
    all matrices should be aligned
    L>0
-------------------------------------------------------------------------*/
void matinv4x4nf_permfwd(float32_t * restrict Y,
                         const float32_t * restrict X,
                         const int16_t *permIx,
                         int L);
void cmatinv4x4nf_permfwd(complex_float * restrict Y,
                         const complex_float * restrict X,
                         const int16_t *permIx,
                         int L);
void matinv4x4sf_permfwd(
                               float32_t * restrict X,
                         const int16_t *permIx,
                         int L);
void cmatinv4x4sf_permfwd(
                               complex_float * restrict X,
                         const int16_t *permIx,
                         int L);

/*-------------------------------------------------
    inplace stream-to-block/block-to-stream conversion:
    Input:
    X[16*L] - L matrices 4x4 in stream order
    Output:
    X       - L matrices written in strange format
              elements of 1-st matrix is written to the 
              place of x00,x01 of original input etc.
              so, 8 matrices replaces 
-------------------------------------------------*/
void matinv4x4f_rsb4x4_inplace(float32_t *X,int L,eLayout layout);
void matinv4x4f_rbs4x4_inplace(float32_t *X,int L,eLayout layout);
void cmatinv4x4f_csb4x4_inplace(complex_float *X,int L,eLayout layout);
void cmatinv4x4f_cbs4x4_inplace(complex_float *X,int L,eLayout layout);

/*-------------------------------------------------------------------------
    inplace forward/backward permutation
    Input/output:
    X[16*L]  - matrices in the intermediate format (after conversion by
               cmatinv4x4sf_csb4x4_inplace)
    Input:
    permIx[L] - permutation indexes
    permTbl[36*32] permutation table (matinv4x4_bkw_perm_tbl_bbe32/
                   matinv4x4_fwd_perm_tbl_bbe32)
-------------------------------------------------------------------------*/
void cmatinv4x4f_permute(
                               complex_float * restrict Y,
                         const complex_float * restrict X,
                         const int16_t *permIx,
                         const int16_t *permTbl,
                         int L,eLayout layout);
void matinv4x4f_permute(
                               float32_t * restrict Y,
                         const float32_t * restrict X,
                         const int16_t *permIx,
                         const int16_t *permTbl,
                         int L,
                         eLayout layout);

/*------------------------------------------
vectorized search algorithm for BBE32EP
(truncated variant - look for 8 possible
combinations of possible 9 for given upper 
left position
Input:
X[16*L]   - L matrices
permIx[L] - upper left positions
Output:
permIx[L] - permutation index from 
matinv4x4_fwd_perm_tbl[]
------------------------------------------*/
void matinv4x4f_truncatedSearch (int16_t* permIx, const float32_t *X, int L, eLayout layout);
void cmatinv4x4f_truncatedSearch(int16_t* permIx, const complex_float *X, int L, eLayout layout);

/*------------------------------------------
find position of upper left corner
Look for an element of maximum absolute value.
Input:
X[16*L]  -L matrices - interleaved in intermediate format
Output:
pos[L]   - index of position of element with 
maximum absolute value
------------------------------------------*/
void cmatinv4x4f_findUL(int16_t *pos,const complex_float *X, int L, eLayout layout);
#endif
