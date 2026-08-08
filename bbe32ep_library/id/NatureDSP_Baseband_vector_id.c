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
 * Vector Operations
 * Annotations
 */

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_vector.h"
#include "common.h"

#ifdef __cplusplus
#define ANNOTATE_FUN1(fun,text) \
  extern "C" const char* ANNOTATE_ATTR ANNOTATE_FUN_REF(fun) = (text)
#else
#define ANNOTATE_FUN1(fun,text) \
  const char* ANNOTATE_ATTR ANNOTATE_FUN_REF(fun)[] = (text)
#endif

ANNOTATE_FUN(rvdot        , "Vector Product ( real 16-bit data )                          ");
ANNOTATE_FUN(rvdotf       , "Vector Product ( real floating-point data )                  ");
ANNOTATE_FUN(cvdot        , "Vector Product ( complex 16-bit data )                       ");
ANNOTATE_FUN(cvdotf       , "Vector Product ( complex floating-point data )               ");
ANNOTATE_FUN(rvadd        , "Vector Sum ( real 16-bit data )                              ");
ANNOTATE_FUN(rvaddf       , "Vector Sum ( real floating-point data )                      ");
ANNOTATE_FUN(cvadd        , "Vector Sum ( complex 16-bit data )                           ");
ANNOTATE_FUN(cvaddf       , "Vector Sum ( complex floating-point data )                   ");
ANNOTATE_FUN(vpower       , "Sum of Squares of a Vector ( 16-bit data )                   ");
ANNOTATE_FUN(vpowerf      , "Sum of Squares of a Vector ( floating-point data )           ");
ANNOTATE_FUN(vmag         , "Square Root of Sum of Squares ( 16-bit data )                ");
ANNOTATE_FUN(vmagf        , "Square Root of Sum of Squares ( floating-point data )        ");
ANNOTATE_FUN(vnorm        , "Vector Normalization ( 16-bit data )                         ");
ANNOTATE_FUN(vbexp        , "Common Exponent ( 16-bit data )                              ");
ANNOTATE_FUN(vbexp_fast   , "Fast Common Exponent ( 16-bit data )                         ");
ANNOTATE_FUN(vbexpf       , "Common Exponent ( floating-point data )                      ");
ANNOTATE_FUN(sbexp        , "Exponent ( 32-bit data )                                     ");
ANNOTATE_FUN(sbexpf       , "Exponent ( floating-point data )                             ");
ANNOTATE_FUN(rvmean       , "Mean of Vector Elements ( real 16-bit data )                 ");
ANNOTATE_FUN(rvmeanf      , "Mean of Vector Elements ( real floating-point data )         ");
ANNOTATE_FUN(cvmean       , "Mean of Vector Elements ( complex 16-bit data )              ");
ANNOTATE_FUN(cvmeanf      , "Mean of Vector Elements ( complex floating-point data )      ");
ANNOTATE_FUN(vmax         , "Dual Peak Search of Maximum Values ( 16-bit data )           ");
ANNOTATE_FUN(vmaxf        , "Dual Peak Search of Maximum Values ( floating-point data )   ");
ANNOTATE_FUN(vmin         , "Dual Peak Search of Minimum Values ( 16-bit data )           ");
ANNOTATE_FUN(vminf        , "Dual Peak Search of Minimum Values ( floating-point data )   ");
ANNOTATE_FUN(vmax8        , "Search for 8 Top Values ( 16-bit data )                      ");
ANNOTATE_FUN(vmax8f       , "Search for 8 Top Values ( floating-point data )              ");
ANNOTATE_FUN(vthreshold   , "Find Values above Threshold ( 16-bit data )                  ");
ANNOTATE_FUN(vthresholdf  , "Find Values above Threshold ( floating-point data )          ");
ANNOTATE_FUN(rmelements2  , "Move the Elements at Given Indices ( real 16-bit data )      ");
ANNOTATE_FUN(rmelements3  , "Move the Elements at Given Indices ( real 16-bit data )      ");
ANNOTATE_FUN(rmelements4  , "Move the Elements at Given Indices ( real 16-bit data )      ");
ANNOTATE_FUN(rmelements6  , "Move the Elements at Given Indices ( real 16-bit data )      ");
ANNOTATE_FUN(rmelements12 , "Move the Elements at Given Indices ( real 16-bit data )      ");
ANNOTATE_FUN(rmelements24 , "Move the Elements at Given Indices ( real 16-bit data )      ");
ANNOTATE_FUN(cmelements2  , "Move the Elements at Given Indices ( complex 16-bit data )   ");
ANNOTATE_FUN(cmelements3  , "Move the Elements at Given Indices ( complex 16-bit data )   ");
ANNOTATE_FUN(cmelements4  , "Move the Elements at Given Indices ( complex 16-bit data )   ");
ANNOTATE_FUN(cmelements6  , "Move the Elements at Given Indices ( complex 16-bit data )   ");
ANNOTATE_FUN(cmelements12 , "Move the Elements at Given Indices ( complex 16-bit data )   ");
ANNOTATE_FUN(cmelements24 , "Move the Elements at Given Indices ( complex 16-bit data )   ");
ANNOTATE_FUN(rrelements2  , "Remove the Elements at Given Indices ( real 16-bit data )    ");
ANNOTATE_FUN(rrelements3  , "Remove the Elements at Given Indices ( real 16-bit data )    ");
ANNOTATE_FUN(rrelements4  , "Remove the Elements at Given Indices ( real 16-bit data )    ");
ANNOTATE_FUN(rrelements6  , "Remove the Elements at Given Indices ( real 16-bit data )    ");
ANNOTATE_FUN(rrelements12 , "Remove the Elements at Given Indices ( real 16-bit data )    ");
ANNOTATE_FUN(rrelements24 , "Remove the Elements at Given Indices ( real 16-bit data )    ");
ANNOTATE_FUN(crelements2  , "Remove the Elements at Given Indices ( complex 16-bit data ) ");
ANNOTATE_FUN(crelements3  , "Remove the Elements at Given Indices ( complex 16-bit data ) ");
ANNOTATE_FUN(crelements4  , "Remove the Elements at Given Indices ( complex 16-bit data ) ");
ANNOTATE_FUN(crelements6  , "Remove the Elements at Given Indices ( complex 16-bit data ) ");
ANNOTATE_FUN(crelements12 , "Remove the Elements at Given Indices ( complex 16-bit data ) ");
ANNOTATE_FUN(crelements24 , "Remove the Elements at Given Indices ( complex 16-bit data ) ");
ANNOTATE_FUN(aunwrap      , "Angle Unwrapping ( 16-bit data )                             ");
ANNOTATE_FUN(aunwrapf     , "Angle Unwrapping ( floating-point data )                     ");


