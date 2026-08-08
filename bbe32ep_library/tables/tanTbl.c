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
#include "tanTbl.h"

/* table for tangent functions
   
   tablated coefficients for 3-nd order min/max polynomials 
   around 16 points for sin(x*pi/8) in range 0...2
   coefficient formats are q18, q17,16 and q15 for subsequent 
   orders
*/

/*
P=3;
N=16;
X=[];Y=[]; Z=[];
coef=zeros(N,P+1);
for n=1:N
    x=2*((n-1)/N+(0:pow2(1,-16):1/N));
    z=tan(x*pi/8);
    x0=2*(n-0.5)/N;
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
    coef(n,:)=round(pow2(p,q));
    fprintf(1,'%d ',round(pow2(p,q)));fprintf(1,'\n');
end
figure(1); plot(X,Y-Z); grid on;
figure(2); plot(X,Y,X,Z); grid on;
fprintf(1,'%d,',coef); fprintf(1,'\n');

*/
const int16_t ALIGN(32) tanTbl[4*16]=
{83,85,88,93,100,110,122,137,157,182,214,255,308,376,467,587,
31,94,158,226,298,377,463,560,670,797,946,1121,1332,1587,1902,2295,
6438,6469,6532,6628,6758,6927,7136,7391,7698,8064,8499,9014,9625,10353,11221,12266,
804,2417,4042,5686,7358,9068,10825,12640,14525,16494,18563,20750,23078,25572,28266,31198};
