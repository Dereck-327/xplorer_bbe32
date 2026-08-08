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
#include "sineTbl.h"

/* table for sine/cosine functions 
   
   tablated coefficients for 3-nd order min/max polynomials 
   around 16 points for sin(x*pi/2) in range 0...2
   coefficient formats are q18, q17,16 and q15 for subsequent 
   orders
*/

/*
        Matlab code for generation of min/max coefficients
        P=3 or 2;
        N=16;
        X=[];Y=[]; Z=[];
        for n=1:N
            x=2*((n-1)/N+(0:pow2(1,-16):1/N));
            z=sin(x*pi/2);
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
        end

*/
const int16_t ALIGN(32) sineTbl[5*16]=
{
    -2632,-2531,-2332,-2044,-1678,-1247,-768,-259,259,768,1247,1678,2044,2332,2531,2632,                 /* Q18 */
    -990,-2932,-4761,-6407,-7807,-8907,-9665,-10051,-10051,-9665,-8907,-7807,-6407,-4761,-2932,-990,     /* Q17 */
    25612,24628,22697,19894,16327,12132,7471,2523,-2523,-7471,-12132,-16327,-19894,-22697,-24628,-25612, /* Q16 */
    25587,24604,22675,19875,16311,12120,7464,2520,-2520,-7464,-12120,-16311,-19875,-22675,-24604,-25587, /* Q16 for P=2 */
    3212,9512,15447,20788,25330,28899,31357,32610,32610,31357,28899,25330,20788,15447,9512,3212          /* Q15 */
};
