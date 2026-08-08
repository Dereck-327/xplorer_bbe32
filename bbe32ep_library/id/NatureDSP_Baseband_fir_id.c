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
 * FIR Filters and Related Functions
 * Annotations
 */

/* Portable data types. */
#include "NatureDSP_types.h"
/* NatureDSP_Baseband Library FIR Filters and Related Functions. */
#include "NatureDSP_Baseband_fir.h"
/* Common utility declarations. */
#include "common.h"

ANNOTATE_FUN(bkfir_process,    "Block FIR filter (real 16-bit data)");
ANNOTATE_FUN(cxfir_process,    "Block FIR filter (complex 16-bit data)");
ANNOTATE_FUN(rcfir_process,    "Block FIR filter (complex 16-bit samples, real 16-bit coefficients)");
ANNOTATE_FUN(srfir_process,    "Block FIR filter with symmetric impulse response (real 16-bit data)");
ANNOTATE_FUN(srcfir_process,   "Block FIR filter with symmetric impulse response (complex 16-bit samples, real 16-bit coefficients)");
ANNOTATE_FUN(firdec_process,   "Decimating block FIR filter (complex 16-bit samples, real 16-bit coefficients)");
ANNOTATE_FUN(firinterp_process,"Interpolating block FIR filter (complex 16-bit samples, real 16-bit coefficients)");
ANNOTATE_FUN(fir_convol_circ,  "Circular convolution (complex 16-bit data)");
ANNOTATE_FUN(fir_convol_lin,   "Linear convolution (complex 16-bit data)");
ANNOTATE_FUN(fir_xcorr_circ,   "Circular cross-correlation (complex 16-bit data)");
ANNOTATE_FUN(fir_xcorr_lin,    "Linear cross-correlation (complex 16-bit data)");
ANNOTATE_FUN(fir_acorr,        "Autocorrelation (complex 16-bit data)");
ANNOTATE_FUN(despread4,        "Direct-Sequence Spread Spectrum modulation despreading");
ANNOTATE_FUN(despread8,        "Direct-Sequence Spread Spectrum modulation despreading");
ANNOTATE_FUN(fir_blms,         "Blockwise adaptive LMS filter (complex 16-bit data)");
ANNOTATE_FUN(fir_blms4,        "Blockwise adaptive LMS filter (complex 16-bit data)");
ANNOTATE_FUN(fir_blms8,        "Blockwise adaptive LMS filter (complex 16-bit data)");
ANNOTATE_FUN(fir_blms16,       "Blockwise adaptive LMS filter (complex 16-bit data)");
ANNOTATE_FUN(fir_blms32,       "Blockwise adaptive LMS filter (complex 16-bit data)");

ANNOTATE_FUN(bkfirf_process,     "Block FIR filter (real floating-point data)");
ANNOTATE_FUN(cxfirf_process,     "Block FIR filter (complex floating-point data)");
ANNOTATE_FUN(rcfirf_process,     "Block FIR filter (complex floating-point samples, real floating-point coefficients)");
ANNOTATE_FUN(firdecf_process,    "Decimating block FIR filter (complex floating-point samples, real floating-point coefficients)");
ANNOTATE_FUN(firinterpf_process, "Interpolating block FIR filter (complex floating-point samples, real floating-point coefficients)");
ANNOTATE_FUN(fir_convolf_circ,   "Circular convolution (complex floating-point data)");
ANNOTATE_FUN(fir_convolf_lin,    "Linear convolution (complex floating-point data)");
ANNOTATE_FUN(fir_rconvolf_circ,  "Circular convolution (real floating-point data)");
ANNOTATE_FUN(fir_rconvolf_lin,   "Linear convolution (real floating-point data)");
ANNOTATE_FUN(fir_xcorrf_circ,    "Circular cross-correlation (complex floating-point data)");
ANNOTATE_FUN(fir_xcorrf_lin,     "Linear cross-correlation (complex floating-point data)");
ANNOTATE_FUN(fir_rxcorrf_circ,   "Circular cross-correlation (real floating-point data)");
ANNOTATE_FUN(fir_rxcorrf_lin,    "Linear cross-correlation (real floating-point data)");
ANNOTATE_FUN(fir_acorrf,         "Autocorrelation (complex floating-point data)");
ANNOTATE_FUN(fir_racorrf,        "Autocorrelation (real floating-point data)");
