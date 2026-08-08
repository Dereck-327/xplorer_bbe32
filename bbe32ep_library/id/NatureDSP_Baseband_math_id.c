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
 * Math Functions
 * Annotations
 */

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_math.h"
#include "common.h"

#ifdef __cplusplus
#define ANNOTATE_FUN1(fun,text) \
  extern "C" const char* ANNOTATE_ATTR ANNOTATE_FUN_REF(fun) = (text)
#else
#define ANNOTATE_FUN1(fun,text) \
  const char* ANNOTATE_ATTR ANNOTATE_FUN_REF(fun)[] = (text)
#endif

ANNOTATE_FUN(vrecip16   , "Vector Reciprocal ( 16-bit data )                   ");
ANNOTATE_FUN(vrecipf    , "Vector Reciprocal ( floating-point data )           ");
ANNOTATE_FUN(vfastrecipf, "Fast Vector Reciprocal ( floating-point data )      ");
ANNOTATE_FUN(srecip16   , "Scalar Reciprocal ( 16-bit data )                   ");
ANNOTATE_FUN(srecipf    , "Scalar Reciprocal ( floating-point data )           ");
ANNOTATE_FUN(vrsqrt     , "Vector Reciprocal Square Root ( 32-bit data )       ");
ANNOTATE_FUN(vrsqrtf    , "Vector Reciprocal Square Root ( floating-point data )");
ANNOTATE_FUN(vfastrsqrtf, "Fast Vector Reciprocal Square Root ( floating-point data )");
ANNOTATE_FUN(srsqrt     , "Scalar Reciprocal Square Root ( 32-bit data )       ");
ANNOTATE_FUN(srsqrtf    , "Scalar Reciprocal Square Root ( floating-point data )");
ANNOTATE_FUN(vsqrtf     , "Vector Square Root ( floating-point data )          ");
ANNOTATE_FUN(vfastsqrtf , "Fast Vector Square Root ( floating-point data )     ");
ANNOTATE_FUN(ssqrtf     , "Scalar Square Root ( floating-point data )          ");
ANNOTATE_FUN(vfastrecip16, "Vector Reciprocal ( pseudo-floating point data )   ");
ANNOTATE_FUN(vfastrsqrt , "Vector Inverse Square Root ( pseudo-floating point data )");
ANNOTATE_FUN(vdivide    , "Vector Division ( 16-bit data )                     ");
ANNOTATE_FUN(vdividef   , "Vector Division ( floating-point data )             ");
ANNOTATE_FUN(sdivide    , "Scalar Division ( 16-bit data )                     ");
ANNOTATE_FUN(sdividef   , "Scalar Division ( floating-point data )             ");
ANNOTATE_FUN(vfmodf     , "Vector Remainder ( floating-point data )            ");
ANNOTATE_FUN(sfmodf     , "Scalar Remainder ( floating-point data )            ");
ANNOTATE_FUN(vlog2      , "Vector Base-2 Logarithm ( 32-bit input data )       ");
ANNOTATE_FUN(vlog10     , "Vector Base-10 Logarithm ( 32-bit input data )      ");
ANNOTATE_FUN(vlogn      , "Vector Natural Logarithm ( 32-bit input data )      ");
ANNOTATE_FUN(vlog2f     , "Vector Base-2 Logarithm ( floating-point data )     ");
ANNOTATE_FUN(vlog10f    , "Vector Base-10 Logarithm ( floating-point data )    ");
ANNOTATE_FUN(vlognf     , "Vector Natural Logarithm ( floating-point data )    ");
ANNOTATE_FUN(slog2      , "Scalar Base-2 Logarithm ( 32-bit input data )       ");
ANNOTATE_FUN(slog10     , "Scalar Base-10 Logarithm ( 32-bit input data )      ");
ANNOTATE_FUN(slogn      , "Scalar Natural Logarithm ( 32-bit input data )      ");
ANNOTATE_FUN(slog2f     , "Scalar Base-2 Logarithm ( floating-point data )     ");
ANNOTATE_FUN(slog10f    , "Scalar Base-10 Logarithm ( floating-point data )    ");
ANNOTATE_FUN(slognf     , "Scalar Natural Logarithm ( floating-point data )    ");
ANNOTATE_FUN(valog10f   , "Vector Base-10 Antilogarithm ( floating-point data )");
ANNOTATE_FUN(valognf    , "Vector Natural Antilogarithm ( floating-point data )");
ANNOTATE_FUN(salog10f   , "Scalar Base-10 Antilogarithm ( floating-point data )");
ANNOTATE_FUN(salognf    , "Scalar Natural Antilogarithm ( floating-point data )");
ANNOTATE_FUN(vldexpf    , "Vector Modify the Exponent of a Floating-Point Number      ");
ANNOTATE_FUN(sldexpf    , "Scalar Modify the Exponent of a Floating-Point Number      ");
ANNOTATE_FUN(vpowf      , "Vector Rise To a Power ( floating-point data )      ");
ANNOTATE_FUN(spowf      , "Scalar Rise To a Power ( floating-point data )      ");
ANNOTATE_FUN(vsine      , "Vector Sine ( 16-bit data )                         ");
ANNOTATE_FUN(vcos       , "Vector Cosine ( 16-bit data )                       ");
ANNOTATE_FUN(vsinef     , "Vector Sine ( floating-point data )                 ");
ANNOTATE_FUN(vcosf      , "Vector Cosine ( floating-point data )               ");
ANNOTATE_FUN(vfastsinef , "Fast Vector Sine ( floating-point data )            ");
ANNOTATE_FUN(vfastcosf  , "Fast Vector Cosine ( floating-point data )          ");
ANNOTATE_FUN(ssine      , "Scalar Sine ( 16-bit data )                         ");
ANNOTATE_FUN(scos       , "Scalar Cosine ( 16-bit data )                       ");
ANNOTATE_FUN(ssinef     , "Scalar Sine ( floating-point data )                 ");
ANNOTATE_FUN(scosf      , "Scalar Cosine ( floating-point data )               ");
ANNOTATE_FUN(vtan        , "Vector Tangent ( 16-bit data )                     ");
ANNOTATE_FUN(vtanf       , "Vector Tangent ( floating-point data )             ");
ANNOTATE_FUN(vfasttanf   , "Fast Vector Tangent ( floating-point data )        ");
ANNOTATE_FUN(stan        , "Scalar Tangent ( 16-bit data )                     ");
ANNOTATE_FUN(stanf       , "Scalar Tangent ( floating-point data )             ");
ANNOTATE_FUN(vcotf       , "Vector Cotangent ( floating-point data )           ");
ANNOTATE_FUN(vfastcotf   , "Fast Vector Cotangent ( floating-point data )      ");
ANNOTATE_FUN(scotf       , "Scalar Cotangent ( floating-point data )           ");
ANNOTATE_FUN(vasinf      , "Vector Arcsine ( floating-point data )             ");
ANNOTATE_FUN(vfastasinf  , "Fast Vector Arcsine ( floating-point data )        ");
ANNOTATE_FUN(sasinf      , "Scalar Arcsine ( floating-point data )             ");
ANNOTATE_FUN(vacosf      , "Vector Arccosine ( floating-point data )           ");
ANNOTATE_FUN(vfastacosf  , "Fast Vector Arccosine ( floating-point data )      ");
ANNOTATE_FUN(sacosf      , "Scalar Arccosine ( floating-point data )           ");
ANNOTATE_FUN(vatan2_16   , "Vector Full Arctangent ( 16-bit data )             ");
ANNOTATE_FUN(satan2_16   , "Scalar Full Arctangent ( 16-bit data )             ");
ANNOTATE_FUN(vatan2f     , "Vector Full Arctangent ( floating-point data )     ");
ANNOTATE_FUN(vfastatan2f , "Fast Vector Full Arctangent ( floating-point data )");
ANNOTATE_FUN(satan2f     , "Scalar Full Arctangent ( floating-point data )     ");
ANNOTATE_FUN(vatan16     , "Vector Arctangent ( 16-bit data )                  ");
ANNOTATE_FUN(vatanf      , "Vector Arctangent ( floating-point data )          ");
ANNOTATE_FUN(vfastatanf  , "Fast Vector Arctangent ( floating-point data )     ");
ANNOTATE_FUN(satan16     , "Scalar Arctangent ( 16-bit data )                  ");
ANNOTATE_FUN(satanf      , "Scalar Arctangent ( floating-point data )          ");
ANNOTATE_FUN(vtanhf      , "Vector Hyperbolic Tangent ( floating-point data )     ");
ANNOTATE_FUN(vfasttanhf  , "Fast Vector Hyperbolic Tangent ( floating-point data )");
ANNOTATE_FUN(stanhf      , "Scalar Hyperbolic Tangent ( floating-point data )     ");
ANNOTATE_FUN(vsinhf      , "Vector Hyperbolic Sine ( floating-point data )        ");
ANNOTATE_FUN(vfastsinhf  , "Fast Vector Hyperbolic Sine ( floating-point data )   ");
ANNOTATE_FUN(ssinhf      , "Scalar Hyperbolic Sine ( floating-point data )        ");
ANNOTATE_FUN(vcoshf      , "Vector Hyperbolic Cosine ( floating-point data )      ");
ANNOTATE_FUN(vfastcoshf  , "Fast Vector Hyperbolic Cosine ( floating-point data ) ");
ANNOTATE_FUN(scoshf      , "Scalar Hyperbolic Cosine ( floating-point data )      ");
ANNOTATE_FUN(vint2float  , "Integer to Floating Value Vector Conversion                  ");
ANNOTATE_FUN(sint2float  , "Integer to Floating Value Scalar Conversion                  ");
ANNOTATE_FUN(vfloat2int  , "Integer to Floating Value Vector Conversion                  ");
ANNOTATE_FUN(sfloat2int  , "Integer to Floating Value Scalar Conversion                  ");
ANNOTATE_FUN(fix2hf      , "Fixed to Half-Precision Floating Point Vector Conversion     ");
ANNOTATE_FUN(hf2fix      , "Half-Precision Floating Point to Fixed Point Vector Conversion");
ANNOTATE_FUN(vfloorf     , "Vector Floating-point Floor                                  ");
ANNOTATE_FUN(sfloorf     , "Scalar Floating-point Floor                                  ");
ANNOTATE_FUN(vceilf      , "Vector Floating-point Ceil                                   ");
ANNOTATE_FUN(sceilf      , "Scalar Floating-point Ceil                                   ");
ANNOTATE_FUN(vabsf       , "Vector Absolute Value ( floating point data )         ");
ANNOTATE_FUN(sabsf       , "Scalar Absolute Value ( floating point data )         ");
ANNOTATE_FUN(vcopysignf  , "Vector Copy Sign ( floating point data )              ");
ANNOTATE_FUN(scopysignf  , "Scalar Copy Sign ( floating point data )              ");
ANNOTATE_FUN(vclipf      , "Vector Clipping ( floating point data )               ");
ANNOTATE_FUN(sclipf      , "Scalar Clipping ( floating point data )               ");
ANNOTATE_FUN(vavgf       , "Vector Average of Two Arguments ( floating point data )     ");
ANNOTATE_FUN(vfastavgf   , "Fast Vector Average of Two Arguments ( floating point data )");
ANNOTATE_FUN(savgf       , "Scalar Average of Two Arguments ( floating point data )     ");
ANNOTATE_FUN(vcountones16, "Vector Count One Bits in a Word ( 16-bit data )        ");
ANNOTATE_FUN(vcountones32, "Vector Count One Bits in a Word ( 32-bit data )        ");
ANNOTATE_FUN(countones16 , "Scalar Count One Bits in a Word ( 16-bit data )        ");
ANNOTATE_FUN(countones32 , "Scalar Count One Bits in a Word ( 32-bit data )        ");



