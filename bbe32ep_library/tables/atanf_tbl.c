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
    tables for atanf(x) approximation
*/
/* Portable data types. */
#include "NatureDSP_types.h"
#include "common.h"
#include "atanf_tbl.h"
#include "common.h"

const union ufloat32uint32 ALIGN(8) atanftbl1[8] =
{
  { 0x3dbc14c0 },/* 9.183645248413086e-002 */
  { 0xbe30c39c },/*-1.726211905479431e-001 */
  { 0x3b2791e4 },/* 2.556913532316685e-003 */
  { 0x3e4dac9d },/* 2.008537799119949e-001 */
  { 0xb97d9a57 },/*-2.418545627733693e-004 */
  { 0xbeaaa7b5 },/*-3.333107531070709e-001 */
  { 0xb54f34c8 },/*-7.719031600572635e-007 */
  { 0x31cf3fa2 } /* 6.031727117772334e-009 */
};

const union ufloat32uint32 ALIGN(8) atanftbl2[8] =
{
  { 0xbcccc037 },/*-2.499399892985821e-002 */
  { 0x3e217c35 },/* 1.577003747224808e-001 */
  { 0xbecf4163 },/*-4.047957360744476e-001 */
  { 0x3ef7b762 },/* 4.838209748268127e-001 */
  { 0xbdf35059 },/*-1.188055947422981e-001 */
  { 0xbe9b8b75 },/*-3.037983477115631e-001 */
  { 0xbb80ed5c },/*-3.934545442461968e-003 */
  { 0x3956fc52 } /* 2.050262701231986e-004 */
};

/*  Matlab code for coefficients:
    order=7;
    c0=0;   c1=0.5;   c2=1;
    x1=(pow2(1,-20):pow2(1,-20):c1);% x1=[-x1(end:-1:1) x1];
    x2=(c1:pow2(1,-20):1);
   
    z1=atan(x1)./x1-1;
    z2=atan(x2)./x2-1;
    p1=polyfit(x1-c0,z1,order);
    p2=polyfit(x2-c2,z2,order);
*/
const union ufloat32uint32 ALIGN(32) atanftbl1a[8] =
{
    {0x3dbc1591}, /*  9.1838009655e-002*/
    {0xbe30c450}, /* -1.7262387276e-001*/
    {0x3b27b0ac}, /*  2.5587482378e-003*/
    {0x3e4dac72}, /*  2.0085313916e-001*/
    {0xb97d7a54}, /* -2.4173530983e-004*/
    {0xbeaaa7b6}, /* -3.3331078291e-001*/
    {0xb54f127b}, /* -7.7140401800e-007*/
    {0x31cf0483}, /*  6.0250058276e-009*/
};
const union ufloat32uint32 ALIGN(32) atanftbl2a[8]=
{
    {0xbcccc051}, /* -2.4994047359e-002*/
    {0xbc8d6020}, /* -1.7257750034e-002*/
    {0x3c876efd}, /*  1.6532415524e-002*/
    {0xbd4a83ca}, /* -4.9442090094e-002*/
    {0x3d43896c}, /*  4.7738477588e-002*/
    {0x3d10efe0}, /*  3.5385012627e-002*/
    {0xbe921fc2}, /* -2.8539854288e-001*/
    {0xbe5bc096}, /* -2.1460184455e-001*/    
};

/*  Matlab code for coefficients:

    order=9;
    x1=(pow2(1,-18):pow2(1,-18):1);
    z1=atan(x1)./x1-1;

    p1=polyfit(x1,z1,order);
*/
const union ufloat32uint32 ALIGN(32) atanftbl_10ord[10]=
{
    {0x3c9cfd46}, /*  1.9163738868559731e-002 */
    {0xbdf59de5}, /* -1.1993006791495249e-001 */
    {0x3e9a31ab}, /*  3.0116017501554360e-001 */
    {0xbeb1ace1}, /* -3.4702209574263121e-001 */
    {0x3daae6a0}, /*  8.3447695083062212e-002 */
    {0x3e37564d}, /*  1.7904014634945500e-001 */
    {0x3b4bd271}, /*  3.1100773808500366e-003 */
    {0xbeaacafc}, /* -3.3357989361779528e-001 */    
    {0x370d895a}, /*  8.4362349722692591e-006 */
    {0xb3995f3e}  /* -7.1419421484951209e-008 */    
};

