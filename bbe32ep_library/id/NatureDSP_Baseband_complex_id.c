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
 * Complex Math Functions
 * Annotations
 */

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_complex.h"
#include "common.h"

#ifdef __cplusplus
#define ANNOTATE_FUN1(fun,text) \
  extern "C" const char* ANNOTATE_ATTR ANNOTATE_FUN_REF(fun) = (text)
#else
#define ANNOTATE_FUN1(fun,text) \
  const char* ANNOTATE_ATTR ANNOTATE_FUN_REF(fun)[] = (text)
#endif

ANNOTATE_FUN(vdividecx      , "Vector Division ( complex 16-bit data )                                              ");
ANNOTATE_FUN(vdividecxf     , "Vector Division ( complex floating-point data )                                      ");
ANNOTATE_FUN(vfastdividecxf , "Fast Vector Division ( complex floating-point data )                                 ");
ANNOTATE_FUN(cdivf          , "Scalar Division ( complex floating-point data )                                      ");
ANNOTATE_FUN(vpolar         , "Polar to Cartesian Vector Conversion ( 16-bit data )                                 ");
ANNOTATE_FUN(vpolarf        , "Polar to Cartesian Vector Conversion ( floating-point data )                         ");
ANNOTATE_FUN(vfastpolarf    , "Polar to Cartesian Vector Conversion ( floating-point data )                         ");
ANNOTATE_FUN(spolarf        , "Polar to Cartesian Scalar Conversion ( floating-point data )                         ");
ANNOTATE_FUN(vcartesianf    , "Cartesian to Polar Vector Conversion ( floating-point data )                         ");
ANNOTATE_FUN(vfastcartesianf, "Cartesian to Polar Vector Conversion ( floating-point data )                         ");
ANNOTATE_FUN(scartesianf    , "Cartesian to Polar Scalar Conversion ( floating-point data )                         ");
ANNOTATE_FUN(vcabs          , "Vector Magnitude ( complex 16-bit data )                                             ");
ANNOTATE_FUN(vcabsf         , "Vector Magnitude ( complex floating-point data )                                     ");
ANNOTATE_FUN(vfastcabsf     , "Vector Magnitude ( complex floating-point data )                                     ");
ANNOTATE_FUN(scabsf         , "Scalar Magnitude ( complex floating-point data )                                     ");
ANNOTATE_FUN(vconjf         , "Vector Conjugate ( complex floating-point data )                                     ");
ANNOTATE_FUN(sconjf         , "Scalar Conjugate ( complex floating-point data )                                     ");
ANNOTATE_FUN(vargf          , "Vector Phase Angle ( complex floating-point data )                                   ");
ANNOTATE_FUN(vfastargf      , "Fast Vector Phase Angle ( complex floating-point data )                              ");
ANNOTATE_FUN(sargf          , "Scalar Phase Angle ( complex floating-point data )                                   ");
ANNOTATE_FUN(cbexp          , "Vector Complex Exponential ( 16-bit data )                                           ");
ANNOTATE_FUN(cbexpf         , "Vector Complex Exponential ( floating-point data )                                   ");
ANNOTATE_FUN(cbfastexpf     , "Fast Vector Complex Exponential ( floating-point data )                              ");
ANNOTATE_FUN(sexpf          , "Scalar Complex Exponential ( floating-point data )                                   ");
ANNOTATE_FUN(cbexplin       , "Vector Complex Exponential with Linearly Evolving Phase ( 16-bit data )              ");
ANNOTATE_FUN(cbexplin_fast  , "Fast Vector Complex Exponential with Linearly Evolving Phase ( 16-bit data )         ");
ANNOTATE_FUN(cbexplinf      , "Vector Complex Exponential with Linearly Evolving Phase ( floating-point data )      ");
ANNOTATE_FUN(cbfastexplinf  , "Fast Vector Complex Exponential with Linearly Evolving Phase ( floating-point data ) ");
ANNOTATE_FUN(cexplinf       , "Scalar Complex Exponential with Linearly Evolving Phase ( floating-point data )      ");
ANNOTATE_FUN(vslope         , "Correction of Vector Slope ( complex 16-bit data )                                   ");

