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

#if !( HAVE_RECIP && HAVE_VSAMATH)
DISCARD_FUN(size_t, vinterp_clin_getScratchSize, (int M, int N))
DISCARD_FUN(void,vinterp_clin,( void *    restrict scr,
                   int16_t * restrict yi,
             const int16_t *          xi,
             const int16_t *          y,
             const int16_t *          x,
             int M, int N ))
#else
/*---------------------------------------------------------------------------
Find values of the target function at points xi[N] by linear
interpolation between the two closest control points.

Representation: 16-bit signed fixed-point

Parameters:
Input:
  C          1 for real data (vinterp_lin), 2 for complex (vinterp_clin)
  scr[]      Scratch area, VINTERP_[C]LIN_SCRATCH(M,N) bytes
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
  x[M]                           Must be distinct
  scr,x[M],y[C*M],xi[N],yi[C*N]  Must be aligned on 32-byte boundary 
  M>=2                           The minimum number of control points is 2
  N                              Must be a multiple of 16

Note:
  Interpolation result is not defined for those desired points that do
  not belong to the closed interval between the left-most and right-most
  control points: xi<x[0] || xi>x[M-1].
---------------------------------------------------------------------------*/
inline_ void * align_ptr(void * ptr, size_t align)
{
  return (void*)(((uintptr_t)ptr + align - 1) & ~(align - 1));

} // align_ptr()

void vinterp_clin ( void *    restrict scr,
                    int16_t * restrict yi,
              const int16_t *          xi,
              const int16_t *          y,
              const int16_t *          x,
              int M, int N )
{
  // Temporary arrays for interpolation parameters for M control points.
  int16_t * a_dy_m;     // [2*M] First derivative, 16-bit, normalized, Q(29-dy_exp)
  int16_t * a_dy_exp_m; // [2*M] First derivative exponent

  // Temporary arrays for control points data replicated for N desired points.
  int16_t * a_x0;       // [N]   Left-hand control point position, Q15
  int16_t * a_y0;       // [2*N] Left-hand control point value, Q15
  int16_t * a_dy_n;     // [2*N] First derivative, Q(29-dy_exp)
  int16_t * a_dy_exp_n; // [2*N] First derivative exponent

  xb_vecNx16 * restrict DY_M_wr;
  xb_vecNx16 * restrict DY_EXP_M_wr;
  xb_vecNx16 * restrict X0_wr;
  xb_vecNx16 * restrict Y0_wr;
  xb_vecNx16 * restrict DY_N_wr;
  xb_vecNx16 * restrict DY_EXP_N_wr;
  xb_vecNx16 * restrict YI_wr;

  int k, m, n;

  NASSERT(!(N & (BBE_SIMD_WIDTH - 1)));

  NASSERT(M >= 2);

  NASSERT_ALIGN32(scr);
  NASSERT_ALIGN32(yi);
  NASSERT_ALIGN32(xi);
  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  if (N<=0) return;
  //
  // Partition the scratch memory area.
  //

  {
    void * ptr = scr;

    a_dy_m = (int16_t*)align_ptr(ptr, 32);
    ptr = a_dy_m + 2 * M + 2 * BBE_SIMD_WIDTH;
    a_dy_exp_m = (int16_t*)align_ptr(ptr, 32);
    ptr = a_dy_exp_m + 2 * M + 2 * BBE_SIMD_WIDTH;
    a_x0 = (int16_t*)align_ptr(ptr, 32);
    ptr = a_x0 + N;
    a_y0 = (int16_t*)align_ptr(ptr, 32);
    ptr = a_y0 + 2 * N;
    a_dy_n = (int16_t*)align_ptr(ptr, 32);
    ptr = a_dy_n + 2 * N;
    a_dy_exp_n = (int16_t*)align_ptr(ptr, 32);
    ptr = a_dy_exp_n + 2 * N;

    ASSERT( (int8_t*)ptr - (int8_t*)scr <= (int)vinterp_clin_getScratchSize(M,N) );
  }

  //----------------------------------------------------------------------------
  // For each segment between the control points compute the first derivative of
  // the target function: dy(m) = (y(m+1)-y(m))/(x(m+1)-x(m)).
  {
    const xb_vecNx16 * Y_rd;
    const xb_vecNx16 * X_rd;
    valign X_va, Y_va;

    xb_vecNx16 x0, x1;
    xb_vecNx16 r_dx0_lo, r_dx0_hi;
    xb_vecNx16 r_dx1_lo, r_dx1_hi;
    xb_vecNx16 dx0_exp, dx1_exp;
    xb_vecNx16 y00, y01, y02;
    xb_vecNx16 y10, y11;
    xb_vecNx16 dy0, dy1;
    xb_vecNx16 dy0_exp, dy1_exp;
    xb_vecNx16 t0, t1;
    xb_vecNx16 c0, c1;
    xb_vecNx40 w0, w1, w2, w3, w4, w5;

    vsaN vsa0, vsa1, vsa2, vsa3;

    vsa1 = BBE_MOVVSA32(16);

    Y_rd = (const xb_vecNx16*)y;
    X_rd = (const xb_vecNx16*)x;
    DY_M_wr = (xb_vecNx16*)a_dy_m;
    DY_EXP_M_wr = (xb_vecNx16*)a_dy_exp_m;

    X_va = BBE_LAVNX16_PP(X_rd);
    Y_va = BBE_LAVNX16_PP(Y_rd);

    BBE_LAVNX16_XP(x0, X_va, X_rd, (uint8_t*)(x+M) - (uint8_t*)X_rd);
    BBE_LAVNX16_XP(y00, Y_va, Y_rd, (uint8_t*)(y+2*M) - (uint8_t*)Y_rd);

    //__Pragma("ymemory( X_rd )");
    //__Pragma("ymemory( Y_rd )");
    for (m = 0; m<(M - 1); m += BBE_SIMD_WIDTH)
    {
      // Q15
      BBE_LAVNX16_XP(x1, X_va, X_rd, (uint8_t*)(x+M) - (uint8_t*)X_rd);

      // CQ15
      BBE_LAVNX16_XP(y01, Y_va, Y_rd, (uint8_t*)(y+2*M) - (uint8_t*)Y_rd);
      BBE_LAVNX16_XP(y02, Y_va, Y_rd, (uint8_t*)(y+2*M) - (uint8_t*)Y_rd);

      //
      // Calculate distances between control points and find their reciprocals.
      //

      t0 = BBE_SELNX16I(x1, x0, BBE_SELI_ROTATE_RIGHT_1);

      w0 = BBE_UNPKSNX16(t0);

      c1 = BBE_MOVVINX16(BBE_MOVVI_INT16_1);

      // Q1.15 <- Q15 - Q15
      BBE_MULSNX16(w0, x0, c1);
      x0 = x1;

      vsa0 = BBE_NSANX40(w0);

      // Q(15+dx_exp) <- Q15 + dx_exp
      w0 = BBE_SLLNX40(w0, vsa0);

      // w1: Q(61-d_exp) <- 2^76/Q(15+dx_exp)
      // t0: Q(62-2*dx_exp) <- 2^92/Q(15+dx_exp)^2, 16-bit unsigned
      // t1: Q(dx_exp-1) <- Q(15+dx_exp) - 16
      BBE_RECIPLUNX40_0(w1, t1, t0, w0);
      BBE_RECIPLUNX40_1(w1, t1, t0, w0);

      // Q(61-dx_exp) <- Q(62-2*dx_exp)*Q(dx_exp-1)
      BBE_MULUSANX16(w1, t0, t1);

      // Q(54-dx_exp) <- Q(61-dx_exp) - 7 w/o rounding, w/o saturation
      w1 = BBE_SRAINX40(w1, 7);

      // 32-bit signed, positive, normalized!
      t0 = BBE_MOVVWL(w1);
      t1 = BBE_MOVVWH(w1);

      // Separate high and low parts of all 32-bit words.
      BBE_DSELNX16I(t1, t0, t1, t0, BBE_DSELI_DEINTERLEAVE_1);

      // Duplicate low parts of 1/dx, Q(54-dx_exp)
      BBE_DSELNX16I(r_dx1_lo, r_dx0_lo, t0, t0, BBE_DSELI_INTERLEAVE_1);
      // Duplicate high parts of 1/dx,  Q(38-dx_exp)
      BBE_DSELNX16I(r_dx1_hi, r_dx0_hi, t1, t1, BBE_DSELI_INTERLEAVE_1);

      //
      // Calculate r_dx*(y(m+1)-y(m)) by 32x16-bit multiplications
      //

      y10 = BBE_SELNX16I(y01, y00, BBE_SELI_ROTATE_RIGHT_2);
      y11 = BBE_SELNX16I(y02, y01, BBE_SELI_ROTATE_RIGHT_2);

      // CQ(53-dx_exp) <- Q(38-dx_exp)*CQ15
      w0 = BBE_MULNX16(r_dx0_hi, y10);
      w1 = BBE_MULNX16(r_dx1_hi, y11);

      // CQ(53-dx_exp) <- CQ(53-dx_exp) - Q(38-dx_exp)*CQ15
      BBE_MULSNX16(w0, r_dx0_hi, y00);
      BBE_MULSNX16(w1, r_dx1_hi, y01);

      // CQ(69-dx_exp) <- Q(54-dx_exp)*CQ15
      w2 = BBE_MULUSRNX16(r_dx0_lo, y10, vsa1);
      w3 = BBE_MULUSRNX16(r_dx1_lo, y11, vsa1);

      // CQ(69-dx_exp) <- Q(54-dx_exp)*CQ15
      w4 = BBE_MULUSNX16(r_dx0_lo, y00);
      w5 = BBE_MULUSNX16(r_dx1_lo, y01);
      y00 = y02;

      // CQ(69-dx_exp)
      w2 = BBE_SUBNX40(w2, w4);
      w3 = BBE_SUBNX40(w3, w5);

      // CQ(53-dx_exp) <- CQ(69-dx_exp) - 16 w/ rounding
      w2 = BBE_SRAINX40(w2, 16);
      w3 = BBE_SRAINX40(w3, 16);

      // CQ(53-dx_exp)
      w0 = BBE_ADDNX40(w0, w2);
      w1 = BBE_ADDNX40(w1, w3);

      //
      // Normalize derivatives and save them as 16-bit normalized values.
      //

      vsa2 = BBE_NSANX40(w0);
      vsa3 = BBE_NSANX40(w1);

      dy0_exp = BBE_MOVVVS(vsa2);
      dy1_exp = BBE_MOVVVS(vsa3);

      vsa2 = BBE_SUBSAVSN(24, vsa2);
      vsa3 = BBE_SUBSAVSN(24, vsa3);

      w0 = BBE_RNDADJNX40(w0, vsa2);
      w1 = BBE_RNDADJNX40(w1, vsa3);

      // CQ(29+dy_exp-dx_exp) <- CQ(53-dx_exp) + dy_exp - 24
      dy0 = BBE_PACKVNX40(w0, vsa2);
      dy1 = BBE_PACKVNX40(w1, vsa3);

      BBE_SVNX16_IP(dy0, DY_M_wr, +4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(dy1, DY_M_wr, +4 * BBE_SIMD_WIDTH / 2);

      t0 = BBE_MOVVVS(vsa0);

      BBE_DSELNX16I(dx1_exp, dx0_exp, t0, t0, BBE_DSELI_INTERLEAVE_1);

      dy0_exp = BBE_SUBNX16(dy0_exp, dx0_exp);
      dy1_exp = BBE_SUBNX16(dy1_exp, dx1_exp);

      BBE_SVNX16_IP(dy0_exp, DY_EXP_M_wr, +2 * BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(dy1_exp, DY_EXP_M_wr, +2 * BBE_SIMD_WIDTH);
    }

    c0 = 0;

    __Pragma("no_reorder");

    BBE_SPNX16_X(c0, a_dy_m, 4 * (M - 1));
    BBE_SPNX16_X(c0, a_dy_exp_m, 4 * (M - 1));
  }
  __Pragma("no_reorder");

  //----------------------------------------------------------------------------
  // For each desired point, find its position in the control points grid and
  // copy interpolation parameters to point's position in temporary arrays.
  //
  // C reference code of the search algorithms:
  //
  //  for ( k=0, m=M-1, n=N-1; k<N+M-1; k++ )
  //  {
  //    if ( xi[n] >= x[m] )
  //    {
  //      a_x0[n] = x[m];
  //
  //      a_y0[2*n+0] = y[2*m+0];
  //      a_y0[2*n+1] = y[2*m+1];
  //
  //      a_dy_n[2*n+0] = a_dy_m[2*m+0];
  //      a_dy_n[2*n+1] = a_dy_m[2*m+1];
  //
  //      a_dy_exp_n[2*n+0] = a_dy_exp_m[2*m+0];
  //      a_dy_exp_n[2*n+1] = a_dy_exp_m[2*m+1];
  //
  //      n = MAX( 0, n-1 );
  //    }
  //    else
  //    {
  //      m = MAX( 0, m-1 );
  //    }
  //  }
  //

  {
    const xb_vecNx16 * XI_rd;
    const int16_t    * Y_rd;
    const int16_t    * X_rd;
    const int16_t    * DY_M_rd;
    const int16_t    * DY_EXP_M_rd;

    xb_vecNx16 x0, y0, dy, dy_exp;
    xb_vecNx16 xi0;
    xb_vecNx16 c0, t0, t1;

    uint32_t ix_m_re, ix_m_cx;            // Control point real/complex indices
    uint32_t ix_n_re, ix_n_cx0, ix_n_cx1; // Desired point real/complex indices

    //                     4        3        2       1       0
    // ix_<...> layout: ix_n_cx1 ix_n_cx0 ix_n_re ix_m_cx ix_m_re
    xb_vecNx16 ix_v, ix_upd0, ix_upd1;

    vboolN vb_cx0, vb_cx1;
    vboolN vb_mask, vb0, vb1;

    static const int16_t ALIGN(32) cst[][BBE_SIMD_WIDTH] = {
      // ix_m_re ix_m_cx       ix_n_re              ix_n_cx0               ix_n_cx1
      { -2, -4, 0, 0, 0 },
      { +2, +4, -2 * BBE_SIMD_WIDTH, -2 * 4 * BBE_SIMD_WIDTH / 2, -2 * 4 * BBE_SIMD_WIDTH / 2 }
    };

    //
    // Setup pointers for control points' data.
    //

    X_rd = (const int16_t*)x;
    Y_rd = (const int16_t*)y;
    DY_M_rd = (const int16_t*)a_dy_m;
    DY_EXP_M_rd = (const int16_t*)a_dy_exp_m;

    //
    // Setup pointers from desired points' data.
    //

    XI_rd = (const xb_vecNx16*)xi;
    X0_wr = (xb_vecNx16*)a_x0;
    Y0_wr = (xb_vecNx16*)a_y0;
    DY_N_wr = (xb_vecNx16*)a_dy_n;
    DY_EXP_N_wr = (xb_vecNx16*)a_dy_exp_n;

    //
    // Setup data indices. The search array starts from the last points and
    // moves backwards.
    //

    ix_m_re = (M - 1) * 2;
    ix_m_cx = (M - 1) * 4;

    ix_n_re = (N - BBE_SIMD_WIDTH) * 2;

    ix_n_cx0 = (N - 2 * BBE_SIMD_WIDTH / 2) * 4;
    ix_n_cx1 = (N - 1 * BBE_SIMD_WIDTH / 2) * 4;

    ix_v = 0;

    // Load data indices to a vector variable.
    t0 = BBE_MOVVA16(ix_n_cx1); ix_v = BBE_SELNX16I(ix_v, t0, BBE_SELI_ROTATE_LEFT_1);
    t0 = BBE_MOVVA16(ix_n_cx0); ix_v = BBE_SELNX16I(ix_v, t0, BBE_SELI_ROTATE_LEFT_1);
    t0 = BBE_MOVVA16(ix_n_re); ix_v = BBE_SELNX16I(ix_v, t0, BBE_SELI_ROTATE_LEFT_1);
    t0 = BBE_MOVVA16(ix_m_cx); ix_v = BBE_SELNX16I(ix_v, t0, BBE_SELI_ROTATE_LEFT_1);
    t0 = BBE_MOVVA16(ix_m_re); ix_v = BBE_SELNX16I(ix_v, t0, BBE_SELI_ROTATE_LEFT_1);

    ix_upd0 = BBE_LVNX16_I((const xb_vecNx16*)cst, 0 * 2 * BBE_SIMD_WIDTH);
    ix_upd1 = BBE_LVNX16_I((const xb_vecNx16*)cst, 1 * 2 * BBE_SIMD_WIDTH);

    t0 = 0;

    vb_mask = BBE_EQNX16(t0, t0);

    for (k = 0; k<N / BBE_SIMD_WIDTH + M - 1; k++)
    {
      //
      // Extract data indices from the vector.
      //

      ix_m_re = BBE_EXTRNX16(ix_v, 0);
      ix_m_cx = BBE_EXTRNX16(ix_v, 1);
      ix_n_re = BBE_EXTRNX16(ix_v, 2);
      ix_n_cx0 = BBE_EXTRNX16(ix_v, 3);
      ix_n_cx1 = BBE_EXTRNX16(ix_v, 4);

      //
      // Load control point interpolation data.
      //

      x0 = BBE_LSNX16_X(X_rd, ix_m_re);
      y0 = BBE_LPNX16_X(Y_rd, ix_m_cx);
      dy = BBE_LPNX16_X(DY_M_rd, ix_m_cx);
      dy_exp = BBE_LPNX16_X(DY_EXP_M_rd, ix_m_cx);

      x0 = BBE_REPNX16(x0, 0);
      y0 = BBE_REPNX16C(y0, 0);
      dy = BBE_REPNX16C(dy, 0);
      dy_exp = BBE_REPNX16C(dy_exp, 0);

      //
      // Replicate control point data for the current vector of desired points.
      // Those desired points that have already been attached on previous
      // iterations are protected of overwriting with a mask.
      //

      BBE_EXTRACTB(vb_cx1, vb_cx0, vb_mask);

      BBE_SVNX16T_X(x0, X0_wr, ix_n_re, vb_mask);

      BBE_SVNX16T_X(y0, Y0_wr, ix_n_cx0, vb_cx0);
      BBE_SVNX16T_X(y0, Y0_wr, ix_n_cx1, vb_cx1);
      BBE_SVNX16T_X(dy, DY_N_wr, ix_n_cx0, vb_cx0);
      BBE_SVNX16T_X(dy, DY_N_wr, ix_n_cx1, vb_cx1);
      BBE_SVNX16T_X(dy_exp, DY_EXP_N_wr, ix_n_cx0, vb_cx0);
      BBE_SVNX16T_X(dy_exp, DY_EXP_N_wr, ix_n_cx1, vb_cx1);

      //
      // If the least significant desired point of the current vector has been
      // assigned to the current control point, then we step back in the desired
      // points vector. Otherwise we step back in the control points array.
      //

      // Load the current vector of desired points.
      xi0 = BBE_LVNX16_X(XI_rd, ix_n_re);

      t0 = BBE_REPNX16(xi0, 0);

      // Check if the least significant desired point resides to the right of
      // the current control point (if so, then all the points of the current
      // vector are done).
      vb0 = BBE_LENX16(x0, t0);

      t0 = BBE_REPNX16(ix_v, 2);

      ix_v = BBE_ADDNX16(ix_v, ix_upd0);

      // Update either the control point, or the desired point vector index.
      BBE_ADDNX16T(ix_v, ix_v, ix_upd1, vb0);

      c0 = 0;
      // Prevent indices from falling to negative values.
      ix_v = BBE_MAXNX16(ix_v, c0);

      t1 = BBE_REPNX16(ix_v, 2);

      // Check if the desired point vector index has actually changed.
      vb1 = BBE_NEQNX16(t0, t1);

      //
      // Update the assignment mask for the next iteration.
      //

      // Set ones for all desired points of the current vector residing to the
      // left of the current control point.
      vb_mask = BBE_LTNX16(xi0, x0);
      // If the desired points vector index has just changed, then all new
      // points are to be attached, so all flags are forced to non-zero.
      vb_mask = BBE_ORB(vb_mask, vb1);
    }
  }

  __Pragma("no_reorder");

  //----------------------------------------------------------------------------
  // Perform linear interpolation for each desired point.
  {
    const xb_vecNx16 * XI_rd;
    const xb_vecNx16 * X0_rd;
    const xb_vecNx16 * Y0_rd;
    const xb_vecNx16 * DY_N_rd;
    const xb_vecNx16 * DY_EXP_N_rd;

    xb_vecNx16 xi0, xi1, yi0, yi1;
    xb_vecNx16 x0, x1, y0, y1;
    xb_vecNx16 dy0, dy1, dy0_exp, dy1_exp;
    xb_vecNx40 w0, w1;
    xb_vecNx16 c1, c29, c38;

    vsaN vsa0, vsa1;

    YI_wr = (xb_vecNx16*)yi;
    XI_rd = (const xb_vecNx16*)xi;
    X0_rd = (const xb_vecNx16*)a_x0;
    Y0_rd = (const xb_vecNx16*)a_y0;
    DY_N_rd = (const xb_vecNx16*)a_dy_n;
    DY_EXP_N_rd = (const xb_vecNx16*)a_dy_exp_n;

    for (n = 0; n<N / BBE_SIMD_WIDTH; n++)
    {
      // Q15
      BBE_LVNX16_IP(xi0, XI_rd, +2 * BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(x0, X0_rd, +2 * BBE_SIMD_WIDTH);

      // Duplicate each entry.
      BBE_DSELNX16I(xi1, xi0, xi0, xi0, BBE_DSELI_INTERLEAVE_1);
      BBE_DSELNX16I(x1, x0, x0, x0, BBE_DSELI_INTERLEAVE_1);

      // CQ15
      BBE_LVNX16_IP(y0, Y0_rd, +4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(y1, Y0_rd, +4 * BBE_SIMD_WIDTH / 2);

      // CQ(29+dy_exp)
      BBE_LVNX16_IP(dy0, DY_N_rd, +4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(dy1, DY_N_rd, +4 * BBE_SIMD_WIDTH / 2);

      BBE_LVNX16_IP(dy0_exp, DY_EXP_N_rd, +4 * BBE_SIMD_WIDTH / 2);
      BBE_LVNX16_IP(dy1_exp, DY_EXP_N_rd, +4 * BBE_SIMD_WIDTH / 2);

      c29 = BBE_MOVVINT16(29);

      dy0_exp = BBE_ADDNX16(c29, dy0_exp);
      dy1_exp = BBE_ADDNX16(c29, dy1_exp);

      c38 = BBE_MOVVINT16(38);

      dy0_exp = BBE_MINNX16(c38, dy0_exp);
      dy1_exp = BBE_MINNX16(c38, dy1_exp);

      vsa0 = BBE_MOVVSV(dy0_exp, 0);
      vsa1 = BBE_MOVVSV(dy1_exp, 0);

      // CQ(44+dy_exp) <- CQ(29+dy_exp)*Q15
      w0 = BBE_MULRNX16(dy0, xi0, vsa0);
      w1 = BBE_MULRNX16(dy1, xi1, vsa1);

      BBE_MULSNX16(w0, dy0, x0);
      BBE_MULSNX16(w1, dy1, x1);

      // CQ1.15 <- CQ(44+dy_exp) - 29 - dy_exp w/ rounding
      w0 = BBE_SRANX40(w0, vsa0);
      w1 = BBE_SRANX40(w1, vsa1);

      c1 = BBE_MOVVINX16(BBE_MOVVI_INT16_1);

      BBE_MULANX16(w0, y0, c1);
      BBE_MULANX16(w1, y1, c1);

      // CQ15 <- sat16( CQ1.15 )
      yi0 = BBE_PACKSNX40(w0);
      yi1 = BBE_PACKSNX40(w1);

      BBE_SVNX16_IP(yi0, YI_wr, +4 * BBE_SIMD_WIDTH / 2);
      BBE_SVNX16_IP(yi1, YI_wr, +4 * BBE_SIMD_WIDTH / 2);
    }
  }

} /* vinterp_clin() */

/* Return the scratch area size, in bytes. */
size_t vinterp_clin_getScratchSize ( int M, int N )
{
  size_t size;
  size = (
  /* a_dy_m[2*M]     */  4 * (M+16 + (-(M)& 15)) + 
  /* a_dy_exp_m[2*M] */  4 * (M+16 + (-(M)& 15)) + 
  /* a_x0[N]         */  2 * (N) + 
  /* a_y0[2*N]       */  4 * (N) + 
  /* a_dy_n[2*N]     */  4 * (N) + 
  /* a_dy_exp_n[2*N] */  4 * (N) );
  return (size);
} /* vinterp_clin_getScratchSize() */
#endif
