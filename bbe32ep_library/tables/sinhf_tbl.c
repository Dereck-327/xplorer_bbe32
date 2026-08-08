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
    tables for sinhf(x) approximation
*/

#include "NatureDSP_types.h"
#include "sinhf_tbl.h"
#include "common.h"

/* Maximum magnitude of input value such that sinhf does not overflow. */
const union ufloat32uint32 sinhf_maxarg = { 0x42b2d4fc }; /* 89.41598510742188 */

/* polynomial coefficients for sinh(x) in range +/-1. 
  
   s=pow2(2,-16); x=(s:s:1); x=[-x(end:-1:1) x]; y=sinh(x)./x;
   p=polyfit(x,y,6);
   p=p(1:2:end); p(end)=[];  
*/

const union ufloat32uint32 ALIGN(32) polysinhf_tbl[4]=
{
        {0x395582a1}, /*2.0361926104e-004*/
        {0x3c087be2}, /*8.3303174043e-003*/
        {0x3e2aaad0}, /*1.6666721614e-001*/
        {0x0}         /* Padding to allow for vector loads */
};
