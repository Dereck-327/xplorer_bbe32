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
    tables for log(x) approximation, Q31
*/

#include "NatureDSP_types.h"
#include "log32b_tbl.h"
#include "common.h"

/* A=round(pow2(2.^((1-(2:17))/16),31)); 
   First Cody coefficient is exact 1Q31 and can be omitted
*/
const int32_t ALIGN(32) Alog2_Q31_tbl[]=
    {
    2056437387,1969251188,1885761398,1805811301,
    1729250827,1655936265,1585730000,1518500250,
    1454120821,1392470869,1333434672,1276901417,
    1222764986,1170923762,1121280436,1073741824};
const int32_t Plog2_Q31_tbl[]={2796203,6557};
