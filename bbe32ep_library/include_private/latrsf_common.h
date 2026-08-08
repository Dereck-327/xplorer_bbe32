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
    Lattice Block Real IIR, floating point, streaming version
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/
#ifndef LATRSF_COMMON_H
#define LATRSF_COMMON_H

/* Portable data types. */
#include "NatureDSP_types.h"

struct tag_latrsf_t;
typedef void (latrsf_processFxns)(struct tag_latrsf_t *handle, float32_t * restrict r,const float32_t *x, int N );

/* Filter instance structure. */
typedef struct tag_latrsf_t
{
  int                 magic;   // Instance pointer validation number
  int                 M;       // Filter order
  int                 L;       // Number of streams
  float32_t         * delLine; // Lattice delay line
  float32_t         * coef;    // Filter coefficients
  float32_t           gain;    // Total gain
  latrsf_processFxns * fxns;   // processing function
} latrsf_t, *latrsf_ptr_t;

latrsf_processFxns latrsf_process1;
latrsf_processFxns latrsf_process2;
latrsf_processFxns latrsf_process3;
latrsf_processFxns latrsf_process4;
latrsf_processFxns latrsf_process5;
latrsf_processFxns latrsf_process6;
latrsf_processFxns latrsf_process7;
latrsf_processFxns latrsf_process8;
latrsf_processFxns latrsf_processX;

#endif // LATRSF_COMMON_H
