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
  NatureDSP Signal Processing Library. IIR part
    Lattice Block Real IIR, floating point
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/
#ifndef LATCF_COMMON_H
#define LATCF_COMMON_H

/* Portable data types. */
#include "NatureDSP_types.h"

struct tag_latcf_t;
typedef void (latcf_processFxns)(struct tag_latcf_t *handle, complex_float * restrict r,const complex_float *x, int N );

/* Filter instance structure. */
typedef struct tag_latcf_t
{
  int                 magic;   // Instance pointer validation number
  int                 M;       // Filter order
  complex_float     * delLine; // Lattice delay line
  float32_t         * coef;    // Filter coefficients
  float32_t           gain;    // Total gain
  latcf_processFxns * fxns;    // processing function
} latcf_t, *latcf_ptr_t;

latcf_processFxns latcf_process2;
latcf_processFxns latcf_process4;
latcf_processFxns latcf_process6;
latcf_processFxns latcf_process8;
latcf_processFxns latcf_processX;

#endif // LATCF_COMMON_H
