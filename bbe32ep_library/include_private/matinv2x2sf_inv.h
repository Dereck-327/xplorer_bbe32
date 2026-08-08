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
/* inversion without permutation */

/* Portable data types. */
#include "NatureDSP_types.h"

#ifndef MATINV2X2SF_INV_H__
#define MATINV2X2SF_INV_H__

#include "matinvf_common.h"
/*-------------------------------------------------------------------------
Direct Matrix Inversion For Real Matrices without permutation
Input/output come in stream format

Input/output:
W  - BBE_SIMD_WIDTH/2 for real or BBE_SIMD_WIDTH/4 for complex matrices
X[L*2x2]  matrices 3x3 in the stream order
L         number of matrices

Temporary
pScr            scratch memory
 
Restrictions:
all matrices have to be aligned
L  should be a multiple of W
-------------------------------------------------------------------------*/
void matinv2x2f_inv (void * restrict pScr, 
                     float32_t* restrict X, 
                     int L, eLayout layout);
void cmatinv2x2f_inv(void * restrict pScr, 
                      complex_float* restrict X, 
                      int L, eLayout layout);

/* Return the scratch area size, in bytes. */
size_t matinv2x2f_inv_getScratchSize  ( int L );
size_t cmatinv2x2f_inv_getScratchSize ( int L );

#endif
