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
    tables for tanhf() approximation
*/

#include "NatureDSP_types.h"
#include "tanhf_tbl.h"

#if TANHF_ALG==0
/* polynomial approximation of tanh(x) in range [-log(3)/2...-log(3)/2]
    only odd coefficients are non zero
    s=pow2(2,-16);
    x=[s:s:log(3)/2+0.008]; x=[-x(end:-1:1) x];
    y=tanh(x); z=tanh(x)./x;
    p=polyfit(x,z,8);
    p=p(1:2:end); p(end)=[];
*/
const union ufloat32uint32 ALIGN(32) polytanhf_tbl[]=
{
    {0x3c86a7d1},/* 1.6437442973e-002*/
    {0xbd57b3ab},/*-5.2661579102e-002*/
    {0x3e086615},/* 1.3320191205e-001*/
    {0xbeaaaa0f} /*-3.3332869411e-001*/
};
#elif TANHF_ALG==1
/* polynomial approximation of tanh(x) in range [-log(3)/2...-log(3)/2]
    only odd coefficients are non zero
    s=pow2(2,-16);
    x=[s:s:log(3)/2]; x=[-x(end:-1:1) x];
    y=tanh(x); z=tanh(x)./x;
    p=polyfit(x,z,10);
    p=p(1:2:end); p(end)=[];
*/
const union ufloat32uint32 ALIGN(32) polytanhf_tbl[]=
{
    {0xbbcfaa0c},/* -6.3374100424e-003*/
    {0x3cacf027},/*- 2.1110607141e-002*/
    {0xbd5c9c1f},/* -5.3859826154e-002*/
    {0x3e0886af},/*  1.3332627271e-001*/
    {0xbeaaaaa5} /*--3.3333316535e-001*/
};


#else
#error wrong TANHF_ALG
#endif

const union ufloat32uint32 halfln3={0x3f0c9f54} ; /* log(3)/2 - tanh(log(3)/2)==0.5 */
