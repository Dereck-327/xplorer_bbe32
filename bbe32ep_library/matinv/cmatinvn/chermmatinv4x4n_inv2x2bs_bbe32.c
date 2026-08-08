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
    Blockwise inversion of 4x4 Hermitian matrices
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Common declarations for blockwise inversion of 4x4 complex matrices. */
#include "cmatinvn_common.h"

#if HAVE_HERMMATINV4X4

/*-------------------------------------------------------------------
 In-place matrix inversion for 2x2 complex matrices stored in
 block-streaming format.

 Internal operation for l-th matrix:
   CQ(30-Q-eX[l]) <- (1<<30)/( CQ(Q) + eX[l] )
 where eX[l] is a normalization shift amount computed internally
 for each input matrix.

 Fixed point position for input data Q is actually used only to compute
 the fixed point position for output data Qi[L]!

  Temporary:
    pScr[(L/16)*16*4]  Scratch memory area
  Input:
    X[L/8][4][8*2]     2x2 complex matrices in block-streaming format, CQ(Q)
    Q                  Fixed point position for input matrices
  Output:
    X[L/8][4][8*2]     Inverse 2x2 matrices, CQ(Qi[L])
    Qi[L]              Fixed point position for each inverse matrix
  Restrictions:
    X,Qi,pScr          Must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
    L                  Must be a multiple of 8
-------------------------------------------------------------------*/

/* Special version for Hermitian matrices. */
void chermmatinv4x4n_inv2x2bs( int16_t * restrict pScr,
                               int16_t * restrict X,
                               int16_t * restrict Qi,
                               int L, int Q )
{
  const xb_vecNx16 * X0_rd;
  const xb_vecNx16 * X1_rd;
        xb_vecNx16 * X_wr;
  const xb_vecNx16 * QI_rd;
        xb_vecNx16 * QI_wr;
  const xb_vecNx16 * P_rd;
        xb_vecNx16 * P_wr;

  xb_vecNx16 x00, x01, x02, x03;
  xb_vecNx16 x10, x11, x12, x13;
  xb_vecNx40 w00, w01, w03;
  xb_vecNx40 w10, w11, w13;

  xb_vecNx16  p0, p1;
  xb_vecNx40  w0, w1;

  xb_vecNx16 d_exp;
  xb_vecNx16 x_exp0, x_exp1;

  xb_vecNx16 t0, t1;
  xb_vecNx16 c15, c38q;

  vsaN vsa0, vsa1, vsa_c23;

  vboolN vb0;

  int l, _L;

  NASSERT_ALIGN( pScr, BBE_SIMD_WIDTH*2 );
  NASSERT_ALIGN( X   , BBE_SIMD_WIDTH*2 );
  NASSERT_ALIGN( Qi  , BBE_SIMD_WIDTH*2 );

  NASSERT( 0 == ( L % (BBE_SIMD_WIDTH/2) ) );

  // L truncated to a multiple of BBE_SIMD_WIDTH.
  _L = ( L & ~(BBE_SIMD_WIDTH-1) );

  c38q = BBE_MOVVA16( 30+8 - Q );

  vsa_c23 = BBE_MOVVSA32( 23 );

  //----------------------------------------------------------------------------
  // Compute determinants and their reciprocals for _L matrices, and store
  // results in the scratch memory.
  //  d = X(1,1)*X(2,2)-X(1,2)*X(2,1);
  //  p = d'/(d*d');
  //

  X0_rd  = (xb_vecNx16*)( (uintptr_t)X + 0*4*BBE_SIMD_WIDTH/2 );
  X1_rd  = (xb_vecNx16*)( (uintptr_t)X + 4*4*BBE_SIMD_WIDTH/2 );
  QI_wr  = (xb_vecNx16*)Qi;
  P_wr   = (xb_vecNx16*)pScr;

  for ( l=0; l<_L/BBE_SIMD_WIDTH; l++ )
  {
    // Load BBE_SIMD_WIDTH 2x2 Hermitian matrices for the next iteration.
    // Formal conversion: CQ15 <- CQ(Q)/2^(15-Q) + 15 - Q
    BBE_LVNX16_IP( x00, X0_rd, +1*4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x01, X0_rd, +2*4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x03, X0_rd, +5*4*BBE_SIMD_WIDTH/2 );

    BBE_LVNX16_IP( x10, X1_rd, +1*4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x11, X1_rd, +2*4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x13, X1_rd, +5*4*BBE_SIMD_WIDTH/2 );

    // Extract real parts of elements (0,0),(1,1).
    x00 = BBE_SELNX16I( x10, x00, BBE_SELI_EXTRACT_1_OF_2_OFF_0 );
    x03 = BBE_SELNX16I( x13, x03, BBE_SELI_EXTRACT_1_OF_2_OFF_0 );

    // Squared magnitude of off-diagonal elements (0,1),(1,0).
    // Q1.30 <- |CQ15|^2
    w0 = BBE_MAGINX16C( x11, x01 );

    // Deinterleave squared magnitudes.
    w0 = BBE_SHFLNX40I( w0, BBE_W_SHFLI_DITLV_1 );

    // Negated determinant (for a Hermitian matrix it's a real number).
    // Q2.30 <- Q1.30 - Q15*Q15
    BBE_MULSNX16( w0, x00, x03 );

    //
    // Take abs. values and remember the original sign.
    //

    w1 = BBE_ZERONX40();

    vb0 = BBE_GTNX40( w0, w1 );

    w0 = BBE_ABSNX40( w0 );

    //
    // 40-bit normalization
    //

    vsa0 = BBE_NSANX40( w0 );

    w0 = BBE_SLLNX40( w0, vsa0 );

    //
    // Compute the reciprocal value: rd = 1/d
    //

    // w1: Q(46-d_exp) <- 2^76/Q(30+d_exp)
    // t0: Q(32-2*d_exp) <- 2^92/Q(30+d_exp)^2; 16-bit unsigned
    // t1: Q(14+d_exp) <- Q(30+d_exp) - 16
    BBE_RECIPLUNX40_0( w1, t1, t0, w0 );
    BBE_RECIPLUNX40_1( w1, t1, t0, w0 );

    // Q(46-d_exp) <- Q(46-d_exp) + Q(32-2*d_exp)*Q(14+d_exp)
    BBE_MULUSANX16( w1, t0, t1 );

    w1 = BBE_RNDSADJNX40( w1, vsa_c23 );

    // Q(23-d_exp) <- Q(46-d_exp) - 23 w/ rounding and saturation
    t0 = BBE_PACKVNX40( w1, vsa_c23 );

    //
    // Restore the sign and duplicate reciprocal determinants.
    //

    BBE_NEGSNX16T( t0, t0, vb0 );

    p0 = BBE_SHFLNX16I( t0, BBE_SHFLI_DOUBLE_1_LO );
    p1 = BBE_SHFLNX16I( t0, BBE_SHFLI_DOUBLE_1_HI );

    BBE_SVNX16_IP( p0, P_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( p1, P_wr, +4*BBE_SIMD_WIDTH/2 );

    //
    // Compute exponent: 30 + 8 - d_exp
    //

    d_exp = BBE_MOVVVS( vsa0 );

    t0 = BBE_SUBNX16( c38q, d_exp );

    BBE_SVNX16_IP( t0, QI_wr, +2*BBE_SIMD_WIDTH );
  }

  //----------------------------------------------------------------------------
  // Compute inverse matrices for _L input matrices.
  //  Y = p*[X(2,2),-X(1,2);-X(2,1),X(1,1)]; 
  //

  c15 = BBE_MOVVINT16( 15 );

  X0_rd = (xb_vecNx16*)X;
  X_wr  = (xb_vecNx16*)X;
  QI_rd = (xb_vecNx16*)( (uintptr_t)QI_wr - 2*_L );
  QI_wr = (xb_vecNx16*)Qi;
  P_rd  = (xb_vecNx16*)( (uintptr_t)P_wr - 4*_L );

  __Pragma( "concurrent" );
  for ( l=0; l<_L/BBE_SIMD_WIDTH; l++ )
  {
    // Load BBE_SIMD_WIDTH 2x2 Hermitian matrices.
    // Formal conversion: CQ15 <- CQ(Q)/2^(15-Q) + 15 - Q
    x01 = BBE_LVNX16_I( X0_rd, 1*4*BBE_SIMD_WIDTH/2 );
    x03 = BBE_LVNX16_I( X0_rd, 3*4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x00, X0_rd, 4*4*BBE_SIMD_WIDTH/2 );

    x11 = BBE_LVNX16_I( X0_rd, 1*4*BBE_SIMD_WIDTH/2 );
    x13 = BBE_LVNX16_I( X0_rd, 3*4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x10, X0_rd, 4*4*BBE_SIMD_WIDTH/2 );

    // Load BBE_SIMD_WIDTH reciprocal determinants
    BBE_LVNX16_IP( p0, P_rd, 4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( p1, P_rd, 4*BBE_SIMD_WIDTH/2 );

    //
    // X = sat16([X(2,2),-X(1,2);-X(2,1),X(1,1)]);
    //

    {
      xb_vecNx16 x00_ = x00, x01_ = x01, x03_ = x03;
      xb_vecNx16 x10_ = x10, x11_ = x11, x13_ = x13;

      x00 = x03_;                 x10 = x13_;
      x01 = BBE_NEGSNX16( x01_ ); x11 = BBE_NEGSNX16( x11_ );
      x03 = x00_;                 x13 = x10_;
    }

    BBE_BMAXABSNX16( vb0, t0, x01, x00 );
    BBE_BMAXABSNX16( vb0, t0,  t0, x03 );
    vsa0 = BBE_NSANX16C( t0 );
    x_exp0 = BBE_MOVVVS( vsa0 );
    vsa0 = BBE_SUBSAVSN( 15, x_exp0 );

    BBE_BMAXABSNX16( vb0, t0, x11, x10 );
    BBE_BMAXABSNX16( vb0, t0,  t0, x13 );
    vsa1 = BBE_NSANX16C( t0 );
    x_exp1 = BBE_MOVVVS( vsa1 );
    vsa1 = BBE_SUBSAVSN( 15, x_exp1 );

    // CQ(38-d_exp) <- CQ15*Q(23-d_exp)
    w00 = BBE_MULRNX16( x00, p0, vsa0 ); w10 = BBE_MULRNX16( x10, p1, vsa1 );
    w01 = BBE_MULRNX16( x01, p0, vsa0 ); w11 = BBE_MULRNX16( x11, p1, vsa1 );
    w03 = BBE_MULRNX16( x03, p0, vsa0 ); w13 = BBE_MULRNX16( x13, p1, vsa1 );

    // CQ(23-d_exp+x_exp) <- CQ(38-d_exp) - ( 15 - x_exp ) w/ rounding
    x00 = BBE_PACKVNX40( w00, vsa0 ); x10 = BBE_PACKVNX40( w10, vsa1 );
    x01 = BBE_PACKVNX40( w01, vsa0 ); x11 = BBE_PACKVNX40( w11, vsa1 );
    x03 = BBE_PACKVNX40( w03, vsa0 ); x13 = BBE_PACKVNX40( w13, vsa1 );

    x02 = BBE_CONJSNX16C( x01 ); x12 = BBE_CONJSNX16C( x11 );

    // Formal conversion:
    // CQ(38-Q-d_exp+x_exp) <- CQ(23-d_exp+x_exp)/2^(15-Q) + 15 - Q
    BBE_SVNX16_IP( x00, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x01, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x02, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x03, X_wr, +4*BBE_SIMD_WIDTH/2 );

    BBE_SVNX16_IP( x10, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x11, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x12, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x13, X_wr, +4*BBE_SIMD_WIDTH/2 );

    BBE_LVNX16_IP( t0, QI_rd, +2*BBE_SIMD_WIDTH );

    t1 = BBE_SELNX16I( x_exp1, x_exp0, BBE_SELI_EXTRACT_1_OF_2_OFF_0 );
    t0 = BBE_ADDNX16( t0, t1 );

    BBE_SVNX16_IP( t0, QI_wr, +2*BBE_SIMD_WIDTH );
  }

  //----------------------------------------------------------------------------
  // Process the last BBE_SIMD_WIDTH/2 matrices if L is not a multiple
  // of BBE_SIMD_WIDTH:
  //  d = X(1,1)*X(2,2)-X(1,2)*X(2,1);
  //  Y = d'/(d*d')*[X(2,2),-X(1,2);-X(2,1),X(1,1)]; 
  //

  if ( _L < L )
  {
    // Load BBE_SIMD_WIDTH 2x2 Hermitian matrices.
    // Formal conversion: CQ15 <- CQ(Q)/2^(15-Q) + 15 - Q
    x00 = BBE_LVNX16_I( X0_rd, 0*4*BBE_SIMD_WIDTH/2 );
    x01 = BBE_LVNX16_I( X0_rd, 1*4*BBE_SIMD_WIDTH/2 );
    x03 = BBE_LVNX16_I( X0_rd, 3*4*BBE_SIMD_WIDTH/2 );

    x10 = BBE_LVNX16_I( X0_rd, 4*4*BBE_SIMD_WIDTH/2 );
    x11 = BBE_LVNX16_I( X0_rd, 5*4*BBE_SIMD_WIDTH/2 );
    x13 = BBE_LVNX16_I( X0_rd, 7*4*BBE_SIMD_WIDTH/2 );

    // Duplicate real parts of elements (0,0),(1,1).
    t0 = BBE_SHFLNX16I( x00, BBE_SHFLI_DUPLICATE_1_EVEN );
    t1 = BBE_SHFLNX16I( x03, BBE_SHFLI_DUPLICATE_1_EVEN );

    // Squared magnitude of off-diagonal elements (0,1),(1,0), duplicated.
    // Q1.30 <- |CQ15|^2
    w0 = BBE_MAGINX16C( x01, x01 );

    // Negated determinants (for a Hermitian matrix it's a real number).
    // Q2.30 <- Q1.30 - Q15*Q15
    BBE_MULSNX16( w0, t0, t1 );

    //
    // Take abs. values and remember the original sign.
    //

    w1 = 0;

    vb0 = BBE_GTNX40( w0, w1 );

    w0 = BBE_ABSNX40( w0 );

    //
    // 40-bit normalization
    //

    vsa0 = BBE_NSANX40( w0 );

    w0 = BBE_SLLNX40( w0, vsa0 );

    d_exp = BBE_MOVVVS( vsa0 );

    //
    // Compute the reciprocal value: rd = 1/d
    //

    // w1: Q(46-d_exp) <- 2^76/Q(30+d_exp)
    // t0: Q(32-2*d_exp) <- 2^92/Q(30+d_exp)^2; 16-bit unsigned
    // t1: Q(14+d_exp) <- Q(30+d_exp) - 16
    BBE_RECIPLUNX40_0( w1, t1, t0, w0 );
    BBE_RECIPLUNX40_1( w1, t1, t0, w0 );

    // Q(46-d_exp) <- Q(46-d_exp) + Q(32-2*d_exp)*Q(14+d_exp)
    BBE_MULUSANX16( w1, t0, t1 );

    w1 = BBE_RNDSADJNX40( w1, vsa_c23 );

    // Q(23-d_exp) <- Q(46-d_exp) - 23; w/ rounding and saturation
    p0 = BBE_PACKVNX40( w1, vsa_c23 );

    // Restore the sign.
    BBE_NEGSNX16T( p0, p0, vb0 );

    //
    // X = sat16([X(2,2),-X(1,2);-X(2,1),X(1,1)]);
    //

    {
      xb_vecNx16 x00_ = x00, x01_ = x01, x03_ = x03;

      x00 = x03_;                
      x01 = BBE_NEGSNX16( x01_ );
      x03 = x00_;                
    }

    BBE_BMAXABSNX16( vb0, t0, x01, x00 );
    BBE_BMAXABSNX16( vb0, t0,  t0, x03 );
    vsa0 = BBE_NSANX16C( t0 );
    x_exp0 = BBE_MOVVVS( vsa0 );
    vsa0 = BBE_SUBSAVSN( 15, x_exp0 );

    // CQ(38-d_exp) <- CQ15*Q(23-d_exp)
    w00 = BBE_MULRNX16( x00, p0, vsa0 );
    w01 = BBE_MULRNX16( x01, p0, vsa0 );
    w03 = BBE_MULRNX16( x03, p0, vsa0 );

    // CQ(23-d_exp+x_exp) <- CQ(38-d_exp) - ( 15 - x_exp ) w/ rounding
    x00 = BBE_PACKVNX40( w00, vsa0 );
    x01 = BBE_PACKVNX40( w01, vsa0 );
    x03 = BBE_PACKVNX40( w03, vsa0 );

    x02 = BBE_CONJSNX16C( x01 );

    // Formal conversion:
    // CQ(38-Q-d_exp+x_exp) <- CQ(23-d_exp+x_exp)/2^(15-Q) + 15 - Q
    BBE_SVNX16_IP( x00, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x01, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x02, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x03, X_wr, +4*BBE_SIMD_WIDTH/2 );

    t0 = BBE_MOVVA16( 30+8 - Q );   
    t0 = BBE_SUBNX16( t0, d_exp ); 
    t0 = BBE_ADDNX16( t0, x_exp0 ); 

    t0 = BBE_SELNX16I( t0, t0, BBE_SELI_EXTRACT_1_OF_2_OFF_0 );
    t1 = BBE_SELNX16I( t0, t0, BBE_SELI_ROTATE_RIGHT_4       );

    BBE_SV4X16_I( t0, QI_wr, 0*4*2 );
    BBE_SV4X16_I( t1, QI_wr, 1*4*2 );
  }

} /* chermmatinv4x4n_inv2x2bs() */

#endif
