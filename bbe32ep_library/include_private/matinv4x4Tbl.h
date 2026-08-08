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
    Direct inversion of 4x4 floating point matrices 
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
/* permutation tables for 4x4 inversion */
#ifndef MATINV4X4TBL_H__
#define MATINV4X4TBL_H__

/* Portable data types. */
#include "NatureDSP_types.h"

#define MATINV4X4_PERM_TBL_SIZE 36

extern const int16_t  matinv4x4_fwd_perm_tbl[MATINV4X4_PERM_TBL_SIZE*16];
extern const int16_t  matinv4x4_bkw_perm_tbl[MATINV4X4_PERM_TBL_SIZE*16];

/* search tables for ported variants */
extern const int16_t matinv4x4_fwd_perm_tbl_bbe32[MATINV4X4_PERM_TBL_SIZE*32];
extern const int16_t matinv4x4_bkw_perm_tbl_bbe32[MATINV4X4_PERM_TBL_SIZE*32];

/*
    Special search table derived from original matinv4x4_fwd_perm_tbl[] using 
    following Matlab code:
------------------------------------------------------
M = [0,1,2;...
     0,3,4;...
     1,3,5;...
     2,4,5 ]; 
a=zeros(1,9);
b=zeros(1,9);
c=zeros(1,9);
d=zeros(1,9);
for maxIx=(0:15)
    row = bitshift((maxIx),-2);
    col = bitand ( (maxIx), 3 );
    p=0;
    for m=0:2
        for n=0:2
            p=p+1;
            ix(p) = 6*M(row+1,m+1) + M(col+1,n+1);
            perm = matinv4x4_fwd_perm_tbl(ix(p)+1,:);
            a(p)=perm(1);
            b(p)=perm(2);
            c(p)=perm(5);
            d(p)=perm(6);
        end
    end
    % convert to N_2 way select and remove last 9-th element
     a=a(1:8); a=[a*2;a*2+1]; a=reshape(a,1,16);
     b=b(1:8); b=[b*2;b*2+1]; b=reshape(b,1,16);
     c=c(1:8); c=[c*2;c*2+1]; c=reshape(c,1,16);
     d=d(1:8); d=[d*2;d*2+1]; d=reshape(d,1,16);
     ix=ix(1:8); ix=[ix;ix*0]; ix=reshape(ix,1,16);
     fprintf(1,'%d,',a);     fprintf(1,'\n'); 
     fprintf(1,'%d,',b);     fprintf(1,'\n'); 
     fprintf(1,'%d,',c);     fprintf(1,'\n'); 
     fprintf(1,'%d,',d);     fprintf(1,'\n'); 
     fprintf(1,'%d,',ix);    fprintf(1,'\n'); 
end
------------------------------------------------------
*/
extern const int16_t ALIGN(32) matinv4x4sf_searchTbl[16*16*5];
#endif
