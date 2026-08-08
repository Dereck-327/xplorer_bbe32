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
#include "NatureDSP_types.h"
#include "common.h"
#include "cbexpTbl.h"

/* table for cbexp/vslope functions 
    - first 16 entries: sin(pi*x) for x=0...2*pi in 16 steps
    - next 8 entries:   sin(pi*(x+1/128)) for x=0...pi/8 in 8 steps
    - next 8 entries:   cos(pi*(x+1/128)) for x=0...pi/8 in 8 steps
*/
const int16_t ALIGN(32) cbexpTbl[32]=
{
    0,12540,23170,30274,32767,30274,23170,12540,0,-12540,-23170,-30274,-32768,-30274,-23170,-12540,
    804,2411,4011,5602,7180,8740,10279,11793,
    32758,32679,32522,32286,31972,31581,31114,30572
};
