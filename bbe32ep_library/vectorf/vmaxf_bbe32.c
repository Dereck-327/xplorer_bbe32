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
  NatureDSP_Baseband library. Vector Operations
    Dual Peak Search
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Vector Operations. */
#include "NatureDSP_Baseband_vector.h"
/* Infinities for single precision routines */
#include "inff_tbl.h"

/*-------------------------------------------------------------------------
Dual Peak Search 

Description: These functions retrieve the maximum (minimum) and next-to-
maximum (next-to-minimum) values of vector elements. They output both the 
peak values and their indices. 

Representation:
vmax,vmin    16-bit signed fixed-point format
             Special values for vmax and vmin (see note 1) are -32768 and
             32767, respectively
vmaxf,vminf  IEEE-754 Std. single precision floating-point format
             Special values for vmaxf and vminf are -HUGE_VALF and +HUGE_VALF,
             respectively

Notes:
1. Each kind of dual-peak search function reserves a special value (see
   above) to maintain internal invariants during the search process. However,
   it is still legal for the input vector to contain special values among
   other data. If this is the case, then the following peculiarities should
   be considered:
     A) Input vector elements equal to the special value are ignored.
     B) If the total number K of input vector elements distinct from the special
        value is less than 2 (i.e. K==0 or K==1), then the last 2-K entries of output
        vector idx[2] are assigned zero, and the last 2-K elements of output vector
        m[2] are assigned the special value.
2. If the peak value is encountered more than once in the input data vector, then
   it will be reported twice in m[2], and idx[2] will contain indices of the first
   two occurencies of the peak value.
3. For floating-point functions NAN values in vectors are ignored. If, however,
   the input vector contains NANs only, then the entries of output vector idx[2] are
   assigned zero, and the elements of output vector m[2] are assigned the special value.
Parameters:
Input:
x[N]     Input data vector
N        Length of input data vector
Output:
m[2]     2 peak values is descending (vmax) or ascending (vmin) order
idx[2]   Indices of 2 peak elements; optional

Restrictions:
x        Aligned on 32-byte boundary
x,m,idx  Must not overlap
N        Multiple of 16 (vmax, vmin) or 8 (maxf, xminf)
-------------------------------------------------------------------------*/
#if !HAVE_VFPU
DISCARD_FUN(void,vmaxf,( float32_t * restrict m, 
             int16_t   * restrict idx, 
      const  float32_t * restrict x, int N ))
#else
void vmaxf ( float32_t * restrict m, 
             int16_t   * restrict idx, 
      const  float32_t * restrict x, int N )
{
  int n;
  float32_t a, b;
  xb_vecN_2xf32 x0, x1;
  xb_vecN_2xf32 max0, max1, maxx, z1;
  vboolN_2 b0, b1;
  vboolN bx0, bx1, bxx;
  xb_vecNx16 id0, id1, tmp, id, id01, zm, _8;
  xb_int16 q0, q1;
  const xb_vecN_2xf32  * restrict pX = (const xb_vecN_2xf32  *)x;
  NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);
  if (N <= 0) { m[0] = m[1] = minusInff.f; idx[0] = idx[1] = 0; return; }
  id = BBE_SEQNX16();
  id = BBE_SHFLNX16I(id, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
  _8 = 8; _8 = BBE_REPNX16(_8, 0);
  id0 = id1 = 0;
  max0 = minusInff.f; max1 = minusInff.f;
  zm = (N - 1);              /* replicate (N-1) */
  for (n = 0; n<(N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
  {
    BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
    b0 = BBE_OGTN_2XF32(x0, max0);
    tmp = BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(max0),BBE_MOVNX16_FROMN_2XF32(x0), BBE_MOVN_FROMN_2(b0));
    x1 = BBE_MOVN_2XF32_FROMNX16(tmp);
    maxx = BBE_MAXN_2XF32(x0, max0);

    bx0 = BBE_JOINB(b0, b0);
    tmp = BBE_MOVNX16T(id0, id, (bx0));
    id0 = BBE_MOVNX16T(id, id0, (bx0));

    b1 = BBE_OGTN_2XF32(x1, max1);
    max1 = BBE_MAXN_2XF32(x1, max1);
    bx1 = BBE_JOINB(b1, b1);
    id1 = BBE_MOVNX16T(tmp, id1, bx1);

    max0 = maxx;
    id = BBE_ADDNX16(id, _8);
  }

  a = BBE_RMAXNUMN_2XF32(max0); /* first minimum value */
  x0 = a;
  b0 = BBE_OEQN_2XF32(x0, max0);
  bx0 = BBE_JOINB(b0, b0);
  id = BBE_MOVNX16T(id0, zm, bx0);
  q0 = BBE_RMINNX16(id);        /* first index */

  tmp = BBE_MOVNX16_FROM16(q0);
  tmp = BBE_REPNX16(tmp, 0);
  bxx = BBE_EQNX16(tmp, id0);

  BBE_EXTRACTB16(b0, b1, bxx);
  z1 = minusInff.f;
  tmp = BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(z1), BBE_MOVNX16_FROMN_2XF32(max0), BBE_MOVN_FROMN_2(b1));
  x1 = BBE_MOVN_2XF32_FROMNX16(tmp);
  b1 = BBE_OEQN_2XF32(x1, max1);            /* flags: "min values is equivalent" */
  BBE_BMINNX16(bx1, id01, id0, id1);        /* flags: "min index of min value"   */

  bx0 = BBE_ANDB(BBE_JOINB(b1, b1), bx1);   /* flags: "choise min value with min index from id0 (else from id1)" */

  b1 = BBE_OGTN_2XF32(x1, max1);
  max1 = BBE_MAXN_2XF32( x1, max1);
  bx1 = BBE_ORB(BBE_JOINB(b1, b1), bx0);    /* use all flags */
  max1 = BBE_MAXN_2XF32(x1, max1);
  b = BBE_RMAXNUMN_2XF32(max1);             /* second maximum value */
  x1 = (b);
  b0 = BBE_OEQN_2XF32(x1, max1);
  bx0 = BBE_JOINB(b0, b0);
  id1 = BBE_MOVNX16T(id0, id1, bx1);
  id1 = BBE_MOVNX16T(id1, zm, bx0);
  q1 = BBE_RMINNX16(id1);                   /* second index */
  m[0] = a; m[1]=b;
  idx[0] = q0; idx[1]=q1;
  return; 

} /* vmaxf() */
#endif
