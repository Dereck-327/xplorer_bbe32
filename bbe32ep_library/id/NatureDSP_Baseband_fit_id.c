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
 * Fitting and Interpolation Routines
 * Annotations
 */
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_fit.h"
#include "common.h"

ANNOTATE_FUN(pfit_grid1  ,"Polynomial fitting, preprocessing (real 16-bit data)");
ANNOTATE_FUN(pfit_grid2  ,"Polynomial fitting, preprocessing (real 16-bit data)");
ANNOTATE_FUN(pfit_grid3  ,"Polynomial fitting, preprocessing (real 16-bit data)");
ANNOTATE_FUN(pfit_grid4  ,"Polynomial fitting, preprocessing (real 16-bit data)");
ANNOTATE_FUN(pfit_grid5  ,"Polynomial fitting, preprocessing (real 16-bit data)");
ANNOTATE_FUN(pfit_grid6  ,"Polynomial fitting, preprocessing (real 16-bit data)");
ANNOTATE_FUN(pfitf_grid1 ,"Polynomial fitting, preprocessing (real floating-point data)");
ANNOTATE_FUN(pfitf_grid2 ,"Polynomial fitting, preprocessing (real floating-point data)");
ANNOTATE_FUN(pfitf_grid3 ,"Polynomial fitting, preprocessing (real floating-point data)");
ANNOTATE_FUN(pfitf_grid4 ,"Polynomial fitting, preprocessing (real floating-point data)");
ANNOTATE_FUN(pfitf_grid5 ,"Polynomial fitting, preprocessing (real floating-point data)");
ANNOTATE_FUN(pfitf_grid6 ,"Polynomial fitting, preprocessing (real floating-point data)");

ANNOTATE_FUN(pfit_process1 ,"Polynomial fitting, calculate coefficients (real 16-bit data)");
ANNOTATE_FUN(pfit_process2 ,"Polynomial fitting, calculate coefficients (real 16-bit data)");
ANNOTATE_FUN(pfit_process3 ,"Polynomial fitting, calculate coefficients (real 16-bit data)");
ANNOTATE_FUN(pfit_process4 ,"Polynomial fitting, calculate coefficients (real 16-bit data)");
ANNOTATE_FUN(pfit_process5 ,"Polynomial fitting, calculate coefficients (real 16-bit data)");
ANNOTATE_FUN(pfit_process6 ,"Polynomial fitting, calculate coefficients (real 16-bit data)");
ANNOTATE_FUN(pfitf_process1,"Polynomial fitting, calculate coefficients (real floating-point data)");
ANNOTATE_FUN(pfitf_process2,"Polynomial fitting, calculate coefficients (real floating-point data)");
ANNOTATE_FUN(pfitf_process3,"Polynomial fitting, calculate coefficients (real floating-point data)");
ANNOTATE_FUN(pfitf_process4,"Polynomial fitting, calculate coefficients (real floating-point data)");
ANNOTATE_FUN(pfitf_process5,"Polynomial fitting, calculate coefficients (real floating-point data)");
ANNOTATE_FUN(pfitf_process6,"Polynomial fitting, calculate coefficients (real floating-point data)");

ANNOTATE_FUN(pfit_eval1  ,"Polynomial fitting, interpolation (real 16-bit data)");
ANNOTATE_FUN(pfit_eval2  ,"Polynomial fitting, interpolation (real 16-bit data)");
ANNOTATE_FUN(pfit_eval3  ,"Polynomial fitting, interpolation (real 16-bit data)");
ANNOTATE_FUN(pfit_eval4  ,"Polynomial fitting, interpolation (real 16-bit data)");
ANNOTATE_FUN(pfit_eval5  ,"Polynomial fitting, interpolation (real 16-bit data)");
ANNOTATE_FUN(pfit_eval6  ,"Polynomial fitting, interpolation (real 16-bit data)");
ANNOTATE_FUN(pfitf_eval1 ,"Polynomial fitting, interpolation (real floating-point data)");
ANNOTATE_FUN(pfitf_eval2 ,"Polynomial fitting, interpolation (real floating-point data)");
ANNOTATE_FUN(pfitf_eval3 ,"Polynomial fitting, interpolation (real floating-point data)");
ANNOTATE_FUN(pfitf_eval4 ,"Polynomial fitting, interpolation (real floating-point data)");
ANNOTATE_FUN(pfitf_eval5 ,"Polynomial fitting, interpolation (real floating-point data)");
ANNOTATE_FUN(pfitf_eval6 ,"Polynomial fitting, interpolation (real floating-point data)");

ANNOTATE_FUN(prootf   ,"Find Roots of a Polynomial (complex floating point data)");
ANNOTATE_FUN(rprootf  ,"Find Roots of a Polynomial (real floating point data)");
ANNOTATE_FUN(vinterp_neib ,"Neighbor interpolation (real 16-bit data)");
ANNOTATE_FUN(vinterp_cneib,"Neighbor interpolation (complex 16-bit data)");
ANNOTATE_FUN(vinterp_lin  ,"Linear interpolation (real 16-bit data)");
ANNOTATE_FUN(vinterp_clin ,"Linear interpolation (complex 16-bit data)");
ANNOTATE_FUN(interp_2tap  ,"2-tap interpolation (real 16-bit data)");
ANNOTATE_FUN(interp_4tap  ,"4-tap interpolation (real 16-bit data)");
