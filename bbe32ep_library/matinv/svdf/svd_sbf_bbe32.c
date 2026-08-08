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
    NatureDSP_Baseband library. Singular Value Decomposition
    Special stream-to-block data order conversions
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* SVD common declarations */
#include "svd_common.h"

#if HAVE_VFPU

#define sz_f32   sizeof(float32_t)
#define sz_f32c  sizeof(complex_float)

#if 0
/* Calculate the number of data elements occupied by a matrix/vector stored
 * in block order. elemNum is the number of payload data elements, each of
 * elemSize bytes. */
static int getSpace( size_t elemSize, int elemNum )
{
  int stp, wid;
  /* Base-2 log of SIMD width for elemSize-byte elements. */
  wid = LOG2_BBE_SIMD_WIDTH - ( 29 - XT_NSA(elemSize) );
  /* At least one element must fit into a SIMD vector! */
  NASSERT(wid>=0);
  /* Select the allocation step from the number of elements: the next
    * power of two, not greater than the SIMD vector size. */
  stp = 30 - XT_NSA(elemNum);
  if (stp>wid) stp = wid;
  /* Allocation size is the storage size rounded up to the next
    * multiple of allocation step. */ 
  return ((1+((elemNum-1)>>stp))<<stp);
}
#endif

/*
 * Convert matrices from stream to block order. Note that the number of
 * matrices to be converted and the total number of matrices in the
 * input stream are specified through separate input arguments: L and
 * stride, respectively.
 * Input:
 *   M,N             Matrix size
 *   L               Number of matrices to be converted
 *   stride          Number of matrices in the input stream
 *   x[M*N][stride]  Input matrices in stream order
 * Output:
 *   y[L][SY]        Output matrices in block order
 * Restrictions:
 *   y,x             Must not overlap and must be aligned on 
 *                   2*BBE_SIMD_WIDTH-byte boundary
 *   L<=stride       Number of converted matrices cannot exceed the number
 *                   of matrices in the input stream
 *   M,N             The product M*N must be a multiple of 4
 * where SY = denotes the number of data entries to store M*N elements in
 * block order with proper alignment.
 */

/* Convert complex data from stream order to block order.
 * Restrictions: 
 *   stride  Must be a multiple of BBE_SIMD_WIDTH/4 */
void svd_csbmxnxsf( complex_float * restrict y, 
              const complex_float *          x,
              int M, int N, int L, int stride )
{
#if 1
  xb_vecN_2xf32 * restrict Y0;
  xb_vecN_2xf32 * restrict Y1;
  xb_vecN_2xf32 * restrict Y2;
  xb_vecN_2xf32 * restrict Y3;
  const xb_vecN_2xf32 * X;
  xb_vecN_2xf32 a0,a1,a2,a3;
  xb_vecN_2xf32 b0,b1,b2,b3;
  vboolN_2 p1,p2,p3;
  int n,k,MN,S;

  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH);
  NASSERT(L<=stride);
  NASSERT(0==((M*N)%4));
  NASSERT(0==(stride%(BBE_SIMD_WIDTH/4)));

  MN = M*N;
  S = (MN+BBE_SIMD_WIDTH/4-1)/(BBE_SIMD_WIDTH/4)*(BBE_SIMD_WIDTH/4);
  Y0 = (xb_vecN_2xf32*)((uintptr_t)y + 0*S*sz_f32c);
  Y1 = (xb_vecN_2xf32*)((uintptr_t)y + 1*S*sz_f32c);
  Y2 = (xb_vecN_2xf32*)((uintptr_t)y + 2*S*sz_f32c);
  Y3 = (xb_vecN_2xf32*)((uintptr_t)y + 3*S*sz_f32c);
  for ( k=0; k<L; k+=BBE_SIMD_WIDTH/4 ) {
    p1 = BBE_MOVN_2_FROMN(BBE_LTNX16(BBE_MOVVA16(k+1),BBE_MOVVA16(L)));
    p2 = BBE_MOVN_2_FROMN(BBE_LTNX16(BBE_MOVVA16(k+2),BBE_MOVVA16(L)));
    p3 = BBE_MOVN_2_FROMN(BBE_LTNX16(BBE_MOVVA16(k+3),BBE_MOVVA16(L)));
    X = (xb_vecN_2xf32*)((uintptr_t)x + k*sz_f32c);
    for ( n=0; n<MN/4; n++ ) {
      BBE_LVN_2XF32_XP(a0, X, stride*sz_f32c);
      BBE_LVN_2XF32_XP(a1, X, stride*sz_f32c);
      BBE_LVN_2XF32_XP(a2, X, stride*sz_f32c);
      BBE_LVN_2XF32_XP(a3, X, stride*sz_f32c);

      BBE_DSELN_2XF32I(b1, b0, a2, a0, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELN_2XF32I(b3, b2, a3, a1, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELN_2XF32I(a1, a0, b2, b0, BBE_DSELI_INTERLEAVE_4);
      BBE_DSELN_2XF32I(a3, a2, b3, b1, BBE_DSELI_INTERLEAVE_4);

      BBE_SVN_2XF32_IP (a0, Y0, BBE_SIMD_WIDTH/4*sz_f32c    );
      BBE_SVN_2XF32T_IP(a1, Y1, BBE_SIMD_WIDTH/4*sz_f32c, p1);
      BBE_SVN_2XF32T_IP(a2, Y2, BBE_SIMD_WIDTH/4*sz_f32c, p2);
      BBE_SVN_2XF32T_IP(a3, Y3, BBE_SIMD_WIDTH/4*sz_f32c, p3);
    } /* n */

    Y0 = (xb_vecN_2xf32*)((uintptr_t)Y0 + 3*S*sz_f32c);
    Y1 = (xb_vecN_2xf32*)((uintptr_t)Y1 + 3*S*sz_f32c);
    Y2 = (xb_vecN_2xf32*)((uintptr_t)Y2 + 3*S*sz_f32c);
    Y3 = (xb_vecN_2xf32*)((uintptr_t)Y3 + 3*S*sz_f32c);
  } /* k */
#else
  int k,n,S;
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH);
  NASSERT(L<=stride);
  NASSERT(0==((M*N)%4));
  NASSERT(0==(stride%(BBE_SIMD_WIDTH/4)));
  S = getSpace(sz_f32c, M*N);
  for ( k=0; k<L; k++ ) {
    for ( n=0; n<M*N; n++ ) {
      y[k*S+n] = x[n*stride+k];
    }
  }
#endif
} /* svd_csbmxnxsf() */

/* Convert real data from stream order to block order.
 * Restrictions: 
 *   stride  Must be a multiple of BBE_SIMD_WIDTH/2 */
void svd_rsbmxnxsf( float32_t * restrict y, 
              const float32_t * x, 
              int M, int N, int L, int stride )
{
#if 1
  const xb_vecN_2xf32 * X;
  xb_vecN_2xf32 * restrict Y;
  valign va;
  vboolN_2 p1,p2,p3,p4,p5,p6,p7;
  xb_vecN_2xf32 a0,a1,a2,a3,a4,a5,a6,a7;
  xb_vecN_2xf32 b0,b1,b2,b3,b4,b5,b6,b7;
  int k,n,MN,P,S;

  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH);
  NASSERT(L<=stride);
  NASSERT(0==((M*N)%4));
  NASSERT(0==(stride%(BBE_SIMD_WIDTH/2)));

  MN = M*N;
  if (MN>4) {
    S = (M*N+(BBE_SIMD_WIDTH/2-1))/(BBE_SIMD_WIDTH/2)*(BBE_SIMD_WIDTH/2);
    for ( k=0; k<L; k+=(BBE_SIMD_WIDTH/2) ) {
      p1 = BBE_MOVN_2_FROMN(BBE_LTNX16(BBE_MOVVA16(k+1),BBE_MOVVA16(L)));
      p2 = BBE_MOVN_2_FROMN(BBE_LTNX16(BBE_MOVVA16(k+2),BBE_MOVVA16(L)));
      p3 = BBE_MOVN_2_FROMN(BBE_LTNX16(BBE_MOVVA16(k+3),BBE_MOVVA16(L)));
      p4 = BBE_MOVN_2_FROMN(BBE_LTNX16(BBE_MOVVA16(k+4),BBE_MOVVA16(L)));
      p5 = BBE_MOVN_2_FROMN(BBE_LTNX16(BBE_MOVVA16(k+5),BBE_MOVVA16(L)));
      p6 = BBE_MOVN_2_FROMN(BBE_LTNX16(BBE_MOVVA16(k+6),BBE_MOVVA16(L)));
      p7 = BBE_MOVN_2_FROMN(BBE_LTNX16(BBE_MOVVA16(k+7),BBE_MOVVA16(L)));
      X = (xb_vecN_2xf32*)((uintptr_t)x + k*sz_f32);
      Y = (xb_vecN_2xf32*)((uintptr_t)y + k*S*sz_f32);
      for ( n=0; n<MN/(BBE_SIMD_WIDTH/2); n++ ) {
        BBE_LVN_2XF32_XP(a0, X, stride*sz_f32);
        BBE_LVN_2XF32_XP(a1, X, stride*sz_f32);
        BBE_LVN_2XF32_XP(a2, X, stride*sz_f32);
        BBE_LVN_2XF32_XP(a3, X, stride*sz_f32);
        BBE_LVN_2XF32_XP(a4, X, stride*sz_f32);
        BBE_LVN_2XF32_XP(a5, X, stride*sz_f32);
        BBE_LVN_2XF32_XP(a6, X, stride*sz_f32);
        BBE_LVN_2XF32_XP(a7, X, stride*sz_f32);

        BBE_DSELN_2XF32I(b1, b0, a4, a0, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(b3, b2, a5, a1, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(b5, b4, a6, a2, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(b7, b6, a7, a3, BBE_DSELI_INTERLEAVE_2);

        BBE_DSELN_2XF32I(a1, a0, b4, b0, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(a3, a2, b5, b1, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(a5, a4, b6, b2, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(a7, a6, b7, b3, BBE_DSELI_INTERLEAVE_2);

        BBE_DSELN_2XF32I(b1, b0, a4, a0, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(b3, b2, a5, a1, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(b5, b4, a6, a2, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(b7, b6, a7, a3, BBE_DSELI_INTERLEAVE_2);

        BBE_SVN_2XF32_XP (b0, Y, S*sz_f32    );
        BBE_SVN_2XF32T_XP(b1, Y, S*sz_f32, p1);
        BBE_SVN_2XF32T_XP(b2, Y, S*sz_f32, p2);
        BBE_SVN_2XF32T_XP(b3, Y, S*sz_f32, p3);
        BBE_SVN_2XF32T_XP(b4, Y, S*sz_f32, p4);
        BBE_SVN_2XF32T_XP(b5, Y, S*sz_f32, p5);
        BBE_SVN_2XF32T_XP(b6, Y, S*sz_f32, p6);
        BBE_SVN_2XF32T_XP(b7, Y, (BBE_SIMD_WIDTH/2-7*S)*sz_f32, p7);
      } /* n */
      if (MN&(BBE_SIMD_WIDTH/4)) {
        BBE_LVN_2XF32_XP(a0, X, stride*sz_f32);
        BBE_LVN_2XF32_XP(a1, X, stride*sz_f32);
        BBE_LVN_2XF32_XP(a2, X, stride*sz_f32);
        BBE_LVN_2XF32_XP(a3, X, stride*sz_f32);

        BBE_DSELN_2XF32I(b1, b0, a2, a0, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(b3, b2, a3, a1, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(a1, a0, b2, b0, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(a3, a2, b3, b1, BBE_DSELI_INTERLEAVE_2);
        b0 = a0; b1 = BBE_SELN_2XF32I(a0, a0, BBE_SELI_ROTATE_RIGHT_8);
        b2 = a1; b3 = BBE_SELN_2XF32I(a1, a1, BBE_SELI_ROTATE_RIGHT_8);
        b4 = a2; b5 = BBE_SELN_2XF32I(a2, a2, BBE_SELI_ROTATE_RIGHT_8);
        b6 = a3; b7 = BBE_SELN_2XF32I(a3, a3, BBE_SELI_ROTATE_RIGHT_8);

        BBE_SVN_2XF32_XP (b0, Y, S*sz_f32    );
        BBE_SVN_2XF32T_XP(b1, Y, S*sz_f32, p1);
        BBE_SVN_2XF32T_XP(b2, Y, S*sz_f32, p2);
        BBE_SVN_2XF32T_XP(b3, Y, S*sz_f32, p3);
        BBE_SVN_2XF32T_XP(b4, Y, S*sz_f32, p4);
        BBE_SVN_2XF32T_XP(b5, Y, S*sz_f32, p5);
        BBE_SVN_2XF32T_XP(b6, Y, S*sz_f32, p6);
        BBE_SVN_2XF32T_XP(b7, Y, (BBE_SIMD_WIDTH/2-7*S)*sz_f32, p7);
      } /* MN */
    } /* k */
  } else { /* MN */
    Y = (xb_vecN_2xf32*)y;
    X = (xb_vecN_2xf32*)x;
    for ( k=0; k<L/(BBE_SIMD_WIDTH/2); k++ ) {
      BBE_LVN_2XF32_XP(a0, X, stride*sz_f32);
      BBE_LVN_2XF32_XP(a1, X, stride*sz_f32);
      BBE_LVN_2XF32_XP(a2, X, stride*sz_f32);
      BBE_LVN_2XF32_XP(a3, X, (BBE_SIMD_WIDTH/2-3*stride)*sz_f32);
      BBE_DSELN_2XF32I(b1, b0, a2, a0, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELN_2XF32I(b3, b2, a3, a1, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELN_2XF32I(a1, a0, b2, b0, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELN_2XF32I(a3, a2, b3, b1, BBE_DSELI_INTERLEAVE_2);
      BBE_SVN_2XF32_IP(a0, Y, BBE_SIMD_WIDTH/2*sz_f32);
      BBE_SVN_2XF32_IP(a1, Y, BBE_SIMD_WIDTH/2*sz_f32);
      BBE_SVN_2XF32_IP(a2, Y, BBE_SIMD_WIDTH/2*sz_f32);
      BBE_SVN_2XF32_IP(a3, Y, BBE_SIMD_WIDTH/2*sz_f32);
    } /* k */
    if (L&(BBE_SIMD_WIDTH/2-1)) {
      BBE_LVN_2XF32_XP(a0, X, stride*sz_f32);
      BBE_LVN_2XF32_XP(a1, X, stride*sz_f32);
      BBE_LVN_2XF32_XP(a2, X, stride*sz_f32);
      BBE_LVN_2XF32_XP(a3, X, stride*sz_f32);
      BBE_DSELN_2XF32I(b1, b0, a2, a0, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELN_2XF32I(b3, b2, a3, a1, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELN_2XF32I(a1, a0, b2, b0, BBE_DSELI_INTERLEAVE_2);
      BBE_DSELN_2XF32I(a3, a2, b3, b1, BBE_DSELI_INTERLEAVE_2);
      va = BBE_ZALIGN(); P = (L&(BBE_SIMD_WIDTH/2-1));
      BBE_SAVN_2XF32_XP(a0, va, Y, (4*P-0*(BBE_SIMD_WIDTH/2))*sz_f32);
      BBE_SAVN_2XF32_XP(a1, va, Y, (4*P-1*(BBE_SIMD_WIDTH/2))*sz_f32);
      BBE_SAVN_2XF32_XP(a2, va, Y, (4*P-2*(BBE_SIMD_WIDTH/2))*sz_f32);
      BBE_SAVN_2XF32_XP(a3, va, Y, (4*P-3*(BBE_SIMD_WIDTH/2))*sz_f32);
      BBE_SAVN_2XF32POS_FP(va, Y);
    } /* L */
  } /* MN */
#else
  int k,n,S;
  NASSERT_ALIGN(y, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2*BBE_SIMD_WIDTH);
  NASSERT(L<=stride);
  NASSERT(0==((M*N)%4));
  NASSERT(0==(stride%(BBE_SIMD_WIDTH/2)));
  S = getSpace(sz_f32, M*N);
  for ( k=0; k<L; k++ ) {
    for ( n=0; n<M*N; n++ ) {
      y[k*S+n] = x[n*stride+k];
    }
  }
#endif
} /* svd_rsbmxnxsf() */

#endif /* HAVE_VFPU */
