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
#include "atanTbl.h"

/* table for arctangent functions 
   
   tablated coefficients for 2-nd order min/max polynomials 
   around 16 points for atan(x)/(pi/4) in range -1...1
   coefficient formats are q17,16 and q15 for subsequent 
   orders
*/

/*
        Matlab code for generation of min/max coefficients
P=2;
N=16;
X=[];Y=[]; Z=[];
for n=1:N
    x=1*((n-1)/N+(0:pow2(1,-16):1/N));
    z=atan(x)/(pi/4);
    x0=1*(n-0.5)/N;
    p=polyfit(x-x0,z,P);
    y=polyval(p,x-x0);
    X=[X x];
    Y=[Y y];
    Z=[Z z];
    q=zeros(1,P+1); q=q+15;
    q(P+1)=15;
    for k=P:-1:1
        q(k)=q(k+1)-1;
    end
end

*/
const int16_t ALIGN(32) atanTbl[3*16]=
{
        -325,-960,-1552,-2077,-2518,-2866,-3120,-3285,-3369,-3384,-3343,-3259,-3142,-3003,-2850,-2689,   /* Q17 */
        20836,20675,20360,19905,19329,18654,17904,17102,16269,15423,14582,13756,12955,12187,11455,10762, /* Q16 */
        1303,3900,6467,8985,11439,13814,16100,18288,20374,22355,24230,26001,27670,29241,30718,32106      /* Q15 */
};
