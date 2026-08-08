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
 * NatureDSP_Baseband Library API
 * IIR Filters
 * Annotations
 */
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_iir.h"
#include "common.h"

ANNOTATE_FUN(bqriir_process   ,"Biquad IIR (real 16-bit data)");
ANNOTATE_FUN(bqriirf_process  ,"Biquad IIR (real floating-point data)");
ANNOTATE_FUN(bqciir_process   ,"Biquad IIR (complex 16-bit data)");
ANNOTATE_FUN(bqciirf_process  ,"Biquad IIR (complex floating-point data)");
ANNOTATE_FUN(bqriirs_process  ,"Biquad IIR (streaming format, real 16-bit data)");
ANNOTATE_FUN(bqriirsf_process ,"Biquad IIR (streaming format, real floating-point data)");
ANNOTATE_FUN(latr_process     ,"Lattice IIR (real 16-bit data)");
ANNOTATE_FUN(latrf_process    ,"Lattice IIR (real floating-point data)");
ANNOTATE_FUN(latc_process     ,"Lattice IIR (complex 16-bit data)");
ANNOTATE_FUN(latcf_process    ,"Lattice IIR (complex floating-point data)");
ANNOTATE_FUN(latrs_process    ,"Lattice IIR (streaming format, real 16-bit data)");
ANNOTATE_FUN(latrsf_process   ,"Lattice IIR (streaming format, real floating-point data)");
