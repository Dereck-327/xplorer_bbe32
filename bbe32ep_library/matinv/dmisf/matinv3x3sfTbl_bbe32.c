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
    NatureDSP_Baseband library. Direct Matrix Inversion
    Direct inversion of 3x3 floating point matrices 
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
#include "matinv3x3Tbl.h"
#if HAVE_VFPU

const int16_t ALIGN(32) matinv3x3_fwd_perm_tbl_bbe32[32*MATINV3X3_PERM_TBL_SIZE] = 
{
    0*2,0*2+1,1*2,1*2+1,2*2,2*2+1,3*2,3*2+1,4*2,4*2+1,5*2,5*2+1,6*2,6*2+1,7*2,7*2+1,8*2,8*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0*2,0*2+1,2*2,2*2+1,1*2,1*2+1,3*2,3*2+1,5*2,5*2+1,4*2,4*2+1,6*2,6*2+1,8*2,8*2+1,7*2,7*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1*2,1*2+1,2*2,2*2+1,0*2,0*2+1,4*2,4*2+1,5*2,5*2+1,3*2,3*2+1,7*2,7*2+1,8*2,8*2+1,6*2,6*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0*2,0*2+1,1*2,1*2+1,2*2,2*2+1,6*2,6*2+1,7*2,7*2+1,8*2,8*2+1,3*2,3*2+1,4*2,4*2+1,5*2,5*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0*2,0*2+1,2*2,2*2+1,1*2,1*2+1,6*2,6*2+1,8*2,8*2+1,7*2,7*2+1,3*2,3*2+1,5*2,5*2+1,4*2,4*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1*2,1*2+1,2*2,2*2+1,0*2,0*2+1,7*2,7*2+1,8*2,8*2+1,6*2,6*2+1,4*2,4*2+1,5*2,5*2+1,3*2,3*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    3*2,3*2+1,4*2,4*2+1,5*2,5*2+1,6*2,6*2+1,7*2,7*2+1,8*2,8*2+1,0*2,0*2+1,1*2,1*2+1,2*2,2*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    4*2,4*2+1,5*2,5*2+1,3*2,3*2+1,7*2,7*2+1,8*2,8*2+1,6*2,6*2+1,1*2,1*2+1,2*2,2*2+1,0*2,0*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

const int16_t ALIGN(32) matinv3x3_bkw_perm_tbl_bbe32[MATINV3X3_PERM_TBL_SIZE*32] ={
0*2,0*2+1,1*2,1*2+1,2*2,2*2+1,3*2,3*2+1,4*2,4*2+1,5*2,5*2+1,6*2,6*2+1,7*2,7*2+1,8*2,8*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0*2,0*2+1,1*2,1*2+1,2*2,2*2+1,6*2,6*2+1,7*2,7*2+1,8*2,8*2+1,3*2,3*2+1,4*2,4*2+1,5*2,5*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
6*2,6*2+1,7*2,7*2+1,8*2,8*2+1,0*2,0*2+1,1*2,1*2+1,2*2,2*2+1,3*2,3*2+1,4*2,4*2+1,5*2,5*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0*2,0*2+1,2*2,2*2+1,1*2,1*2+1,3*2,3*2+1,5*2,5*2+1,4*2,4*2+1,6*2,6*2+1,8*2,8*2+1,7*2,7*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0*2,0*2+1,2*2,2*2+1,1*2,1*2+1,6*2,6*2+1,8*2,8*2+1,7*2,7*2+1,3*2,3*2+1,5*2,5*2+1,4*2,4*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
6*2,6*2+1,8*2,8*2+1,7*2,7*2+1,0*2,0*2+1,2*2,2*2+1,1*2,1*2+1,3*2,3*2+1,5*2,5*2+1,4*2,4*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
2*2,2*2+1,0*2,0*2+1,1*2,1*2+1,5*2,5*2+1,3*2,3*2+1,4*2,4*2+1,8*2,8*2+1,6*2,6*2+1,7*2,7*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
8*2,8*2+1,6*2,6*2+1,7*2,7*2+1,2*2,2*2+1,0*2,0*2+1,1*2,1*2+1,5*2,5*2+1,3*2,3*2+1,4*2,4*2+1,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

/*
    Special search table derived from original matinv3x3_fwd_perm_tbl[] using 
    following Matlab code:
------------------------------------------------------
matinv3x3_fwd_perm_tbl=[
    0,1,2,3,4,5,6,7,8;...
    0,2,1,3,5,4,6,8,7;...
    1,2,0,4,5,3,7,8,6;...
    0,1,2,6,7,8,3,4,5;...
    0,2,1,6,8,7,3,5,4;...
    1,2,0,7,8,6,4,5,3;...
    3,4,5,6,7,8,0,1,2;...
    4,5,3,7,8,6,1,2,0];
a=zeros(1,8);
b=zeros(1,8);
c=zeros(1,8);
d=zeros(1,8);
for ix=1:8
    perm = matinv3x3_fwd_perm_tbl(ix,:);
    a(ix)=perm(1);
    b(ix)=perm(2);
    c(ix)=perm(4);
    d(ix)=perm(5);
end
% convert to N_2 way select
a=[a*2;a*2+1]; a=reshape(a,1,16);
b=[b*2;b*2+1]; b=reshape(b,1,16);
c=[c*2;c*2+1]; c=reshape(c,1,16);
d=[d*2;d*2+1]; d=reshape(d,1,16);
fprintf(1,'%d,',a);     fprintf(1,'\n');
fprintf(1,'%d,',b);     fprintf(1,'\n');
fprintf(1,'%d,',c);     fprintf(1,'\n');
fprintf(1,'%d,',d);     fprintf(1,'\n');
------------------------------------------------------
*/
const int16_t ALIGN(32) matinv3x3sf_searchTbl[16*5]=
{
0,1,0,1,2,3,0,1,0,1,2,3,6,7,8,9,
2,3,4,5,4,5,2,3,4,5,4,5,8,9,10,11,
6,7,6,7,8,9,12,13,12,13,14,15,12,13,14,15,
8,9,10,11,10,11,14,15,16,17,16,17,14,15,16,17,
//
0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0
};

#endif
