/* ------------------------------------------------------------------------ */
/* Copyright (c) 2017 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs (“Cadence    */
/* Libraries”) are the copyrighted works of Cadence Design Systems Inc.     */
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
#ifndef MATINV3X3HFS_H__
#define MATINV3X3HFS_H__
/* half-stream DMI without permutations */

/* Portable data types. */
#include "NatureDSP_types.h"

/*-------------------------------------------------------------------------
Direct Matrix Inversion For Real Matrices without permutation
Input/output come in half-stream format

Matlab formulas:
A=inv(A);
T=C*A;
D=D-T*B;
D=inv(D);
C=-D*T;
T=A*B;
A=A-T*C;
B=-T*D;

Input/output:
W  - BBE_SIMD_WIDTH/2 for real or BBE_SIMD_WIDTH/4 for complex matrices
X[L/W][SX][W]  matrices 3x3 in halfstream format
L               number of matrices

SX=16 for real and =12 for complex matrices

Temporary
pScr            scratch memory
 
Restrictions:
all matrices have to be aligned
L  should be a multiple of W
-------------------------------------------------------------------------*/
void matinv3x3hfs(void * restrict pScr, 
                  float32_t* restrict X,
                  int L);
void cmatinv3x3hfs(void * restrict pScr, 
                  complex_float* restrict X,
                  int L);

/* Return the scratch area size, in bytes. */
size_t  matinv3x3hfs_getScratchSize ( int L );
size_t cmatinv3x3hfs_getScratchSize ( int L );

#endif
