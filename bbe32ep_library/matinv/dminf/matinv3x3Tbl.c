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
/* permutation tables for 3x3 inversion */

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
#include "matinv3x3Tbl.h"
#if HAVE_VFPU
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

/*------------------------------------------------------------------------------
  MATLAB code for tables generation:

  %
  % Build forward and inverse permutation tsbles for the reference C code
  %

  % 3 reorderings of 3 elements such that the first two elements from 3 rows
  % contain all "2 of 3" combinations: (1,2), (1,3), (2,3).
  ep = [1,2,3;1,3,2;2,3,1];
  I3 = eye(3);
  A = reshape(0:8,3,3)';
  P = zeros(3,3,3);
  for k=1:3, P(:,:,k) = I3(:,ep(k,:)); end;
  fwd_perm_tbl = zeros(9,9);
  bkw_perm_tbl = zeros(9,9);
  for m=1:3
    for n=1:3
      fwd_perm_tbl((m-1)*3+n,:) = reshape((P(:,:,m)'*A*P(:,:,n))',1,9);
      bkw_perm_tbl((m-1)*3+n,:) = reshape((P(:,:,n)*A*P(:,:,m)')',1,9);
    end
  end

  % Throw away one permutation to have a nice number of 8 variants.
  fwd_perm_tbl(8,:) = [];
  bkw_perm_tbl(8,:) = [];

  fprintf('#define PERM_TBL_SIZE    8\n');
  fprintf('\n');
  fprintf('static const int fwd_perm_tbl[PERM_TBL_SIZE*9] =\n{\n');
  for m=1:8
    fprintf('  ');
    fprintf('%d,',fwd_perm_tbl(m,:));
    fprintf('\n');
  end
  fprintf('};\n');
  fprintf('\n');
  fprintf('static const int bkw_perm_tbl[PERM_TBL_SIZE*9] =\n{\n');
  for m=1:8
    fprintf('  ');
    fprintf('%d,',bkw_perm_tbl(m,:));
    fprintf('\n');
  end
  fprintf('};\n\n');

------------------------------------------------------------------------------*/
const int16_t ALIGN(32) matinv3x3_fwd_perm_tbl[9*MATINV3X3_PERM_TBL_SIZE] = 
{
    0,1,2,3,4,5,6,7,8,
    0,2,1,3,5,4,6,8,7,
    1,2,0,4,5,3,7,8,6,
    0,1,2,6,7,8,3,4,5,
    0,2,1,6,8,7,3,5,4,
    1,2,0,7,8,6,4,5,3,
    3,4,5,6,7,8,0,1,2,
    //3,5,4,6,8,7,0,2,1,
    4,5,3,7,8,6,1,2,0
};

const int16_t ALIGN(32) matinv3x3_bkw_perm_tbl[9*MATINV3X3_PERM_TBL_SIZE] =
{
0,1,2,3,4,5,6,7,8,
0,1,2,6,7,8,3,4,5,
6,7,8,0,1,2,3,4,5,
0,2,1,3,5,4,6,8,7,
0,2,1,6,8,7,3,5,4,
6,8,7,0,2,1,3,5,4,
2,0,1,5,3,4,8,6,7,
//2,0,1,8,6,7,5,3,4,
8,6,7,2,0,1,5,3,4
};
#endif
