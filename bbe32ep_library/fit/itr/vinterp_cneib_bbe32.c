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
  NatureDSP_Baseband library. Fitting and Interpolation Routines
    Vector Interpolation
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fit.h"

/*---------------------------------------------------------------------------
Find values of the target function at points xi[N] by locating the
nearest control point.

Representation: 16-bit signed fixed-point

Parameters:
Input:
  C          1 for real data (vinterp_neib), 2 for complex (vinterp_cneib)
  xi[N]      Desired points
  y[C*M]     Values of the target function sampled at control points
  x[M]       Control points
Output:
  yi[C*N]    Values of the target function at points xi[N]
Temporary:
  scr[]      Scratch memory area. To determine the scratch area size required by
             a function vinterp_<fun>, use the respective helper function 
             vinterp_<fun>_getScratchSize(M,N)

Restrictions:
  scr,x[M],y[C*M],xi[N],yi[C*N]  Must not overlap
  x[M],xi[N]                     Must be sorted in the ascending order
  scr,x[M],xi[N],yi[C*N]         Must be aligned on 32-byte boundary 
  y[C*M]                         Must be aligned on (C*2)-byte boundary
  M>=1                           The minimum number of control points is 1
  N                              Must be a multiple of 16
---------------------------------------------------------------------------*/

void vinterp_cneib ( void *    restrict scr,
                     int16_t * restrict yi,
               const int16_t *          xi,
               const int16_t *          y,
               const int16_t *          x,
               int M, int N )
{
  int16_t * restrict a_xm;

  xb_vecNx16 * restrict XM_wr;
  xb_vecNx16 * restrict YI_wr;

  int m, n;

  NASSERT_ALIGN32(scr);
  NASSERT_ALIGN32(xi);
  NASSERT_ALIGN32(yi);
  NASSERT_ALIGN32(x);

  NASSERT_ALIGN4(y);

  NASSERT(!(N & (BBE_SIMD_WIDTH - 1)));

  NASSERT(M >= 1);

  if (N<=0) return;

  //
  // Partition the scratch memory.
  //

  {
    void * ptr = scr;

    a_xm = (int16_t*)ptr;
    ptr = a_xm + M + (-M & (BBE_SIMD_WIDTH - 1));

    ASSERT((int8_t*)ptr - (int8_t*)scr <= (int)vinterp_cneib_getScratchSize(M,N));
  }

  //----------------------------------------------------------------------------
  // Compute the middle point for each segment between control points.

  {
    const xb_vecNx16 * X_rd;
    valign X_va;

    xb_vecNx16 x0, x1, x2, xm0;
    xb_vecNx40 w0;

    vsaN vsa0;

    X_rd = (const xb_vecNx16*)x;
    XM_wr = (xb_vecNx16*)a_xm;

    vsa0 = BBE_MOVVSA32(1);

    X_va = BBE_LAVNX16_PP(X_rd);
    BBE_LAVNX16_XP(x0, X_va, X_rd, (uint8_t*)(x+M) - (uint8_t*)X_rd);

    for (m = 0; m<M - 1; m += BBE_SIMD_WIDTH)
    {
      BBE_LAVNX16_XP(x1, X_va, X_rd, (uint8_t*)(x+M) - (uint8_t*)X_rd);

      x2 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_RIGHT_1);

      w0 = BBE_ADDWNX16(x2, x0);
      x0 = x1;

      xm0 = BBE_PACKVNX40(w0, vsa0);

      BBE_SVNX16_IP(xm0, XM_wr, +2 * BBE_SIMD_WIDTH);
    }

    __Pragma("no_reorder");

    // Append a pseudo point. It is necessary for the right-most control point
    // to absorb all the desired point lying to the right of it.
    a_xm[M - 1] = 0x7fff;
  }

  //------------------------------------------------------------------------------
  // C reference code of the search algorithm:
  //
  //  for ( k=m=n=0; k<N+M-1; k++ )
  //  {
  //    yi[2*n+0] = y[2*m+0];
  //    yi[2*n+1] = y[2*m+1];
  //
  //    if ( (int32_t)xi[n] <= a_xm[m] )
  //    {
  //      n = MIN( n+1, N-1 );
  //    }
  //    else
  //    {
  //      m = MIN( m+1, M-1 );
  //    }
  //  }
  //

  if (N > BBE_SIMD_WIDTH)
  {
    //----------------------------------------------------------------------------
    // Whenever the total number of desired points is at least twice the SIMD
    // vector width, use a two-channel search procedure which performs two
    // instances of the same search algorithm independently and simultaneously
    // for even and odd vectors of desired points.
    //

    const xb_vecNx16 * XI_rd;
    const int16_t    * XM_rd;
    const int16_t    * Y_rd;

    xb_vecNx16 t0, t1;
    xb_vecNx16 xm0, xi0, y0;

    vboolN vb0, vb1;
    vboolN vb_mask0, vb_mask1;

    uint32_t ix_m_re, ix_m_cx;            // Control point real/complex indices 
    uint32_t ix_n_re, ix_n_cx0, ix_n_cx1; // Desired point real/complex indices

    int N0, N1;

    //                     4     3    2    1    0
    // ixN_<...> layout: n_cx1 n_cx0 n_re m_cx m_re
    xb_vecNx16 ix0_v, ix1_v;      // Index vectors
    xb_vecNx16 ix0_lim, ix1_lim;  // Index limits
    xb_vecNx16 ix_upd0, ix_upd1;  // Index update terms

    static const int16_t ALIGN(32) cst[][BBE_SIMD_WIDTH] = {
      // m_re m_cx       n_re                n_cx0                 n_cx1
      { 2, 4, 0, 0, 0 }, // ix_upd0
      { -2, -4, 2 * 2 * BBE_SIMD_WIDTH, 4 * 4 * BBE_SIMD_WIDTH / 2, 4 * 4 * BBE_SIMD_WIDTH / 2 }, // ix_upd1
      { 0, 0, 0, 0 * 4 * BBE_SIMD_WIDTH / 2, 1 * 4 * BBE_SIMD_WIDTH / 2 }, // ix0_v
      { 0, 0, 2 * BBE_SIMD_WIDTH, 2 * 4 * BBE_SIMD_WIDTH / 2, 3 * 4 * BBE_SIMD_WIDTH / 2 }  // ix1_v
    };

    // Upper limits for control points' real/complex data.
    ix_m_re = 2 * (M - 1);
    ix_m_cx = 4 * (M - 1);

    //
    // Setup the even search instance.
    //

    // The number of desired points the even search instance has to process.
    N0 = (((N / BBE_SIMD_WIDTH) + 1) >> 1)*BBE_SIMD_WIDTH;

    // Upper limit for desired points' real data.
    ix_n_re = 2 * (2 * N0 - 2 * BBE_SIMD_WIDTH);

    // Upper limit for desired points' complex data.
    ix_n_cx0 = 4 * (2 * N0 - 2 * BBE_SIMD_WIDTH) + 0 * 4 * BBE_SIMD_WIDTH / 2;
    ix_n_cx1 = 4 * (2 * N0 - 2 * BBE_SIMD_WIDTH) + 1 * 4 * BBE_SIMD_WIDTH / 2;

    ix0_lim = 0;

    // Make a tuple of all index limits.
    t0 = BBE_MOVVA16(ix_n_cx1); ix0_lim = BBE_SELNX16I(ix0_lim, t0, BBE_SELI_ROTATE_LEFT_1);
    t0 = BBE_MOVVA16(ix_n_cx0); ix0_lim = BBE_SELNX16I(ix0_lim, t0, BBE_SELI_ROTATE_LEFT_1);
    t0 = BBE_MOVVA16(ix_n_re); ix0_lim = BBE_SELNX16I(ix0_lim, t0, BBE_SELI_ROTATE_LEFT_1);
    t0 = BBE_MOVVA16(ix_m_cx); ix0_lim = BBE_SELNX16I(ix0_lim, t0, BBE_SELI_ROTATE_LEFT_1);
    t0 = BBE_MOVVA16(ix_m_re); ix0_lim = BBE_SELNX16I(ix0_lim, t0, BBE_SELI_ROTATE_LEFT_1);

    // Load initial index values.
    ix0_v = BBE_LVNX16_I((const xb_vecNx16*)cst, 2 * 2 * BBE_SIMD_WIDTH);

    //
    // Setup the odd search instance.
    //

    // The number of desired points the odd search instance has to process.
    N1 = (((N / BBE_SIMD_WIDTH) + 0) >> 1)*BBE_SIMD_WIDTH;

    // Upper limit for desired points' real data.
    ix_n_re = 2 * (2 * N1 - 1 * BBE_SIMD_WIDTH);

    // Upper limit for desired points' complex data.
    ix_n_cx0 = 4 * (2 * N1 - 1 * BBE_SIMD_WIDTH) + 0 * 4 * BBE_SIMD_WIDTH / 2;
    ix_n_cx1 = 4 * (2 * N1 - 1 * BBE_SIMD_WIDTH) + 1 * 4 * BBE_SIMD_WIDTH / 2;

    ix1_lim = 0;

    // Make a tuple of all index limits.
    t0 = BBE_MOVVA16(ix_n_cx1); ix1_lim = BBE_SELNX16I(ix1_lim, t0, BBE_SELI_ROTATE_LEFT_1);
    t0 = BBE_MOVVA16(ix_n_cx0); ix1_lim = BBE_SELNX16I(ix1_lim, t0, BBE_SELI_ROTATE_LEFT_1);
    t0 = BBE_MOVVA16(ix_n_re); ix1_lim = BBE_SELNX16I(ix1_lim, t0, BBE_SELI_ROTATE_LEFT_1);
    t0 = BBE_MOVVA16(ix_m_cx); ix1_lim = BBE_SELNX16I(ix1_lim, t0, BBE_SELI_ROTATE_LEFT_1);
    t0 = BBE_MOVVA16(ix_m_re); ix1_lim = BBE_SELNX16I(ix1_lim, t0, BBE_SELI_ROTATE_LEFT_1);

    // Load initial index values.
    ix1_v = BBE_LVNX16_I((const xb_vecNx16*)cst, 3 * 2 * BBE_SIMD_WIDTH);

    //
    // Do the search.
    //

    t0 = 0;

    vb_mask0 = vb_mask1 = BBE_EQNX16(t0, t0);

    ix_upd0 = BBE_LVNX16_I((const xb_vecNx16*)cst, 0 * 2 * BBE_SIMD_WIDTH);
    ix_upd1 = BBE_LVNX16_I((const xb_vecNx16*)cst, 1 * 2 * BBE_SIMD_WIDTH);

    XM_rd = (const int16_t   *)a_xm;
    Y_rd = (const int16_t   *)y;
    XI_rd = (const xb_vecNx16*)xi;
    YI_wr = (xb_vecNx16*)yi;

    __Pragma("concurrent");
    for (n = 0; n<(N / BBE_SIMD_WIDTH + 1) / 2 + M - 1; n++)
    {
      //------------------------------------------------------------------------
      // Even search instance.

      // Extract individual indices from the index tuple.
      ix_m_re = BBE_EXTRNX16(ix0_v, 0);
      ix_n_re = BBE_EXTRNX16(ix0_v, 2);
      ix_n_cx0 = BBE_EXTRNX16(ix0_v, 3);
      ix_n_cx1 = BBE_EXTRNX16(ix0_v, 4);
#if 0
      // Better schedule when compiled for 27june
      ix_m_cx = BBE_EXTRNX16(ix0_v, 1);
#else
      // Makes happy 23oct!
      ix_m_cx = ix_m_re * 2;
#endif
      // Load the current vector of desired point positions.
      xi0 = BBE_LVNX16_X(XI_rd, ix_n_re);

      // Load the current control point position and value.
      xm0 = BBE_LSNX16_X(XM_rd, ix_m_re); xm0 = BBE_REPNX16(xm0, 0);
      y0 = BBE_LPNX16_X(Y_rd, ix_m_cx); y0 = BBE_REPNX16C(y0, 0);

      BBE_EXTRACTB(vb1, vb0, vb_mask0);

      // Assign desired points the control point value. Desired points assigned
      // on previous iterations are protected of overwriting by the mask.
      BBE_SVNX16T_X(y0, YI_wr, ix_n_cx0, vb0);
      BBE_SVNX16T_X(y0, YI_wr, ix_n_cx1, vb1);

      //
      // Look if all desired points of the current vector has been attached to
      // respective control points. That is, check if the last desired point
      // matches the current control point.
      //

      t0 = BBE_REPNX16(xi0, BBE_SIMD_WIDTH - 1);

      vb0 = BBE_LENX16(t0, xm0);

      //
      // Conditionally update the index tuple: if all the desired points of the
      // current vector has been attached, then we go to the next vector,
      // otherwise we step to the next control point.
      //

      t0 = BBE_REPNX16(ix0_v, 2);

      ix0_v = BBE_ADDNX16(ix0_v, ix_upd0);
      BBE_ADDNX16T(ix0_v, ix0_v, ix_upd1, vb0);

      // Saturate updated indices!
      ix0_v = BBE_MINUNX16(ix0_v, ix0_lim);

      t1 = BBE_REPNX16(ix0_v, 2);

      // Detect if we have actually switched to the next desired points vector.
      vb1 = BBE_NEQNX16(t1, t0);

      //
      // Update the storage mask for the current vector of desired points.
      //

      // Set non-zero for those desired points that still haven't been attached
      // to a control point.
      vb_mask0 = BBE_LTNX16(xm0, xi0);
      // If we have just switched to the next vector, then all points are to be
      // attached, so all flags must be non-zero.
      vb_mask0 = BBE_ORB(vb_mask0, vb1);

      //------------------------------------------------------------------------
      // Odd search instance. Exactly the same operations as of the even one.

      ix_m_re = BBE_EXTRNX16(ix1_v, 0);
      ix_m_cx = BBE_EXTRNX16(ix1_v, 1);
      ix_n_re = BBE_EXTRNX16(ix1_v, 2);
      ix_n_cx0 = BBE_EXTRNX16(ix1_v, 3);
      ix_n_cx1 = BBE_EXTRNX16(ix1_v, 4);

      xi0 = BBE_LVNX16_X(XI_rd, ix_n_re);

      xm0 = BBE_LSNX16_X(XM_rd, ix_m_re); xm0 = BBE_REPNX16(xm0, 0);
      y0 = BBE_LPNX16_X(Y_rd, ix_m_cx); y0 = BBE_REPNX16C(y0, 0);

      BBE_EXTRACTB(vb1, vb0, vb_mask1);

      BBE_SVNX16T_X(y0, YI_wr, ix_n_cx0, vb0);
      BBE_SVNX16T_X(y0, YI_wr, ix_n_cx1, vb1);

      t0 = BBE_REPNX16(xi0, BBE_SIMD_WIDTH - 1);

      vb0 = BBE_LENX16(t0, xm0);

      t0 = BBE_REPNX16(ix1_v, 2);

      ix1_v = BBE_ADDNX16(ix1_v, ix_upd0);
      BBE_ADDNX16T(ix1_v, ix1_v, ix_upd1, vb0);

      ix1_v = BBE_MINUNX16(ix1_v, ix1_lim);

      t1 = BBE_REPNX16(ix1_v, 2);

      vb1 = BBE_NEQNX16(t1, t0);

      vb_mask1 = BBE_LTNX16(xm0, xi0);
      vb_mask1 = BBE_ORB(vb_mask1, vb1);
    }
  }
  else
  {
    const int16_t * XM_rd;
    const int16_t * Y_rd;

    xb_vecNx16 xm0, y0;
    xb_vecNx16 xi0, yi1;

    vboolN vb_mask, vb_cx0, vb_cx1;

    XM_rd = (const int16_t   *)a_xm;
    Y_rd = (const int16_t   *)y;
    YI_wr = (xb_vecNx16*)yi;

    xi0 = BBE_LVNX16_I((const xb_vecNx16*)xi, 0);

    yi1 = 0;

    vb_mask = BBE_EQNX16(yi1, yi1);

    __Pragma("loop_count min=1");
    for (m = 0; m<M; m++)
    {
      BBE_LSNX16_IP(xm0, XM_rd, +2);
      BBE_LPNX16_IP(y0, Y_rd, +4);

      xm0 = BBE_REPNX16(xm0, 0);
      y0 = BBE_REPNX16C(y0, 0);

      BBE_EXTRACTB(vb_cx1, vb_cx0, vb_mask);

      BBE_SVNX16T_I(y0, YI_wr, 0, vb_cx0);

      yi1 = BBE_MOVNX16T(y0, yi1, vb_cx1);

      vb_mask = BBE_LTNX16(xm0, xi0);
    }

    BBE_SVNX16_I(yi1, YI_wr, 4 * BBE_SIMD_WIDTH / 2);
  }

} /* vinterp_cneib() */

/* Return the scratch area size, in bytes. */
size_t vinterp_cneib_getScratchSize ( int M, int N )
{
  return (2 * (((M)+15)&(~15)));

} /* vinterp_cneib_getScratchSize() */
