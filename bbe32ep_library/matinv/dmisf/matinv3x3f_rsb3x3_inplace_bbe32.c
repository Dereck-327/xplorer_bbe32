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

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
#include "matinv3x3_common.h"
#if HAVE_VFPU

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
void matinv3x3f_rbs3x3_inplace(float32_t *X,int L)
{
    /* bingo! for bbe32 matinv3x3sf_rbs3x3_inplace 
       and matinv3x3sf_rsb3x3_inplace are fully equivalent 
       no need to duplicate the code
    */
    matinv3x3f_rsb3x3_inplace(X,L);
}

#endif
