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
    Direct inversion of 3x3 floating point matrices 
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#ifndef MATINV3X3_COMMON_H__
#define MATINV3X3_COMMON_H__
/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
#include "matinvf_common.h"

/*-------------------------------------------------------------------------
    find permutation for L matrices written in block order
    Input:
    X[L*SX]  input matrices
    L        number of matrices
    SX=16 for real block data 
    SX=12 for complex block data 
    SX=9  for stream order
    Output:
    permIx[L] permutation indices
-------------------------------------------------------------------------*/
void  matinv3x3nf_findPerm(int16_t *permIx,const float32_t *X, int L);
void cmatinv3x3nf_findPerm(int16_t *permIx,const complex_float *X, int L);
void  matinv3x3sf_findPerm(int16_t *permIx,const float32_t *X, int L);
void cmatinv3x3sf_findPerm(int16_t *permIx,const complex_float *X, int L);
/*-------------------------------------------------------------------------
    find permutation for L matrices written in stream order
    and return permuted matrices
    Input:
    X[L*9]   input matrices
    L        number of matrices
    Output:
    permIx[L] permutation indices
-------------------------------------------------------------------------*/
void  matinv3x3sf_findAndPerm (int16_t *permIx,float32_t *X, int L);
void  cmatinv3x3sf_findAndPerm(int16_t *permIx,complex_float *X, int L);

/*-------------------------------------------------------------------------
    backward permutation with conversion to halfstream format
    Input:
    X[SX*L']   3x3 matrices in halfstream format
    permIx[L]  permutation indices
    L          number of matrices
    L' smallest bigger integer which is multiple of BBE_SIMD_WIDTH/2 for
    real or BBE_SIMD_WIDTH/4 for complex data
    Output:
    Y[SX*L]    3x3 matrices in block format
    SX=16 for real block data 
    SX=12 for complex block data 
    SX=9  for stream order

    NOTE: X might be damaged 
    Restrictions:
    all matrices should be aligned
    L>0
-------------------------------------------------------------------------*/
void matinv3x3nf_permbkw(float32_t * restrict Y,
                         float32_t * restrict X,
                         const int16_t *permIx,
                         int L);
void cmatinv3x3nf_permbkw(complex_float * restrict Y,
                         complex_float * restrict X,
                         const int16_t *permIx,
                         int L);
void matinv3x3sf_permbkw(
                         float32_t * restrict X,
                         const int16_t *permIx,
                         int L);
void cmatinv3x3sf_permbkw(
                         complex_float * restrict X,
                         const int16_t *permIx,
                         int L);

/*-------------------------------------------------------------------------
    forward permutation with conversion to halfstream format
    Input:
    X[SX*L]    input 3x3 matrices
    permIx[L]  permutation indices
    L          number of matrices
    L' smallest bigger integer which is multiple of BBE_SIMD_WIDTH/2 for
    real or BBE_SIMD_WIDTH/4 for complex data
    Output:
    Y[SX*L']   3x3 matrices in halfstream format
    SX=16 for real block data 
    SX=12 for complex block data 
    SX=9  for stream order

    Restrictions:
    all matrices should be aligned
    L>0
-------------------------------------------------------------------------*/
void matinv3x3nf_permfwd(      float32_t * restrict Y,
                         const float32_t * restrict X,
                         const int16_t *permIx,
                         int L);
void cmatinv3x3nf_permfwd(     complex_float * restrict Y,
                         const complex_float * restrict X,
                         const int16_t *permIx,
                         int L);
void matinv3x3sf_permfwd(      
                               float32_t * restrict X,
                         const int16_t *permIx,
                         int L);
void cmatinv3x3sf_permfwd(     
                               complex_float * restrict X,
                         const int16_t *permIx,
                         int L);

/*-------------------------------------------------
    inplace stream-to-block/block-to-stream conversion:
    Input:
    X[9*L] - L matrices 3x3 in stream order 
             (last element is not processed!)
    Output:
    X       - L matrices written in strange format
              elements of 1-st matrix is written to the 
              place of x00,x01 of original input etc.
              so, 8 matrices replaces 8 original elements
              last element is left unchanged
-------------------------------------------------*/
void matinv3x3f_rsb3x3_inplace(float32_t *X,int L);
void matinv3x3f_rbs3x3_inplace(float32_t *X,int L);
void cmatinv3x3f_csb3x3_inplace(complex_float *X,int L,eLayout layout);
void cmatinv3x3f_cbs3x3_inplace(complex_float *X,int L,eLayout layout);

/*-------------------------------------------------------------------------
    inplace forward/backward permutation
    Input/output:
    X[9*L]  - matrices in the intermediate format (after conversion by
              matinv3x3sf_rsb3x3_inplace)
    Input:
    permIx[L] - permutation indexes
    permTbl[8*32] permutation table (matinv3x3_bkw_perm_tbl_bbe32/
                  matinv3x3_fwd_perm_tbl_bbe32)
-------------------------------------------------------------------------*/
void matinv3x3f_permute(
                               float32_t * restrict Y,
                         const float32_t * restrict X,
                         const int16_t *permIx,
                         const int16_t *permTbl,
                         int L, eLayout layout);
void cmatinv3x3f_permute(
                               complex_float * restrict Y,
                         const complex_float * restrict X,
                         const int16_t *permIx,
                         const int16_t *permTbl,
                         int L, eLayout layout);

/*-----------------------------------------------
    Search for permutation in block/stream ordered 
    arrays
    Input:
    X[]  - 3x3 block/stream ordered data
    Output:
    permIx[L] - permutation index
-----------------------------------------------*/
void matinv3x3f_search (int16_t *permIx,const float32_t *X, int L, eLayout layout);
void cmatinv3x3f_search(int16_t *permIx,const complex_float *X, int L, eLayout layout);

#endif
