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
#include "sqrtTbl.h"

/* 3-rd order polynomials for sqrt(x*3/4+0.25) for x=0...1 with no centering
   Matlab code:
   P=3;
   X=[];Y=[]; Z=[];
   for n=1:16
       x=(n-1)/16+(0:pow2(1,-16):1/16);
       z=sqrt(x*3/4+0.25);
       x0=(n-1.0)/16;
       p=polyfit(x-x0,z,P);
       y=polyval(p,x-x0);
       X=[X x];
       Y=[Y y];
       Z=[Z z];
   end
*/
const int16_t ALIGN(32) sqrtTbl[4*16]=
{
 22178, 14918, 10596, 7844, 5997, 4707, 3774, 3081, 2553, 2144, 1822, 1563, 1353, 1180, 1037,  917,
-18225,-14126,-11360,-9390,-7930,-6812,-5934,-5229,-4653,-4176,-3775,-3434,-3142,-2889,-2668,-2474,
 24573, 22551, 20958,19660,18577,17656,16859,16161,15543,14991,14494,14043,13632,13255,12908,12587,
 16384, 17854, 19212,20480,21674,22806,23884,24915,25905,26859,27780,28672,29537,30377,31194,31991
};
