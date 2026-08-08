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
    Direct Matrix Inversion
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Common utility declarations. */
#include "cmatinvn_common.h"

#if HAVE_MATINV4X4

/*-------------------------------------------------------------------
 In-place matrix inversion for 2x2 complex matrices stored in
 block-streaming format. This function is similar to inv2x2bs(),
 but it allows to specify the fixed point position for each
 input matrix, and it doesn't normalize input matrices.

 Internal operation for the l-th matrix:
   CQ(30-Q[l]-eX[l]) <- (1<<30)/( CQ(Q[l]) + eX[l] )
 where eX[l] is a normalization shift amount computed internally
 for each inverse matrix.

 Fixed point position for input data is actually used only to compute
 the fixed point position for output data!

  Temporary:
    pScr[(L/16)*16*4]  Scratch memory area
  Input:
    X[L/8][4][8*2]     2x2 complex matrices in block-streaming format, CQ(Q)
    Q[L]               Fixed point position for each input matrix
  Output:
    X[L/8][4][8*2]     Inverse 2x2 matrices, CQ(Qi[L])
    Q[L]               Fixed point position for each inverse matrix
  Restrictions:
    X,Q,pScr           Must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
    L                  Must be a multiple of 8
-------------------------------------------------------------------*/
void cmatinv4x4n_inv2x2bs_q( int16_t * restrict pScr,
                             int16_t * restrict X,
                             int16_t * restrict Q,
                             int L )
{
#if 0
  //
  // MATLAB code for a single matrix operation, Q = 15:
  //  % Q30 <- Q15*Q15
  //  d = X(1,1)*X(2,2)-X(1,2)*X(2,1);
  //  d_exp = nsa40(d)-7;
  //  % Q(2.13+d_exp) <- Q30 + d_exp - 17 w/o rounding
  //  d = sat16(floor(d*2^(d_exp-17)));
  //  % Q(5.11+2*d_exp) <- Q(2.13+d_exp)^2 + Q(2.13+d_exp)^2 - 15 w/ rounding
  //  d2 = floor(d*d'*2^-15+0.5); if d2>=2^16, d2 = d2-2^16; end;
  //  % Q(18-2*d_exp) <- Q29/Q(5.11+2*d_exp)
  //  rd2 = sat16u(floor(2^29/d2));
  //  % Q(15-d_exp) <- Q(2.13+d_exp)*Q(18-2*d_exp) - 16 w/ rounding
  //  p = sat16(floor(d*rd2*2^-16+0.5+0.5j));
  //  % Q(15-d_exp+X_exp) <- Q(15-d_exp)*Q15 + X_exp - 15 w/ rounding
  //  X = sat16([X(2,2),-X(1,2);-X(2,1),X(1,1)]);
  //  X = sat16(floor(p'*X*2^-15+0.5+0.5j));
  //  Qi = 15-d_exp;
  //
  const xb_vecNx16 *          X_rd;
        xb_vecNx16 * restrict X_wr;
  const xb_vecNx16 *          Q_rd;
        xb_vecNx16 * restrict Q_wr;
  const xb_vecNx16 *          P_rd;
        xb_vecNx16 * restrict P_wr;

  xb_vecNx16 x00, x01, x02, x03;
  xb_vecNx16 x10, x11, x12, x13;
  xb_vecNx40 w00, w01, w02, w03;
  xb_vecNx40 w10, w11, w12, w13;

  xb_vecNx16  d0, d1;
  xb_vecNx16  d2, rd2, rd2_0, rd2_1;
  xb_vecNx16U d2u, rd2u, rem;
  xb_vecNx16  p0, p1;
  xb_vecNx40  w0, w1;

  xb_vecNx16 d_exp0, d_exp1;

  xb_vecNx16 t0, t1;
  xb_vecNx16 c37;
  xb_vecNx40 c1q29;

  vsaN vsa0, vsa1, vsa_c15, vsa_c16, vsa_c24;

  int l;

  NASSERT_ALIGN32( X );
  NASSERT_ALIGN32( Q );

  NASSERT( !(L&7) );

  c37 = BBE_MOVVINT16( 37 );

  c1q29 = BBE_MOVWINX40( BBE_MOVWI_Q9_30_1 );
  c1q29 = BBE_SRAINX40( c1q29, 1 );
 
  vsa_c16 = BBE_MOVVSA32( 16 );
  vsa_c24 = BBE_MOVVSA32( 24 );

  //----------------------------------------------------------------------------
  // Compute the determinant and its reciprocal for L&~15 matrices, and store
  // results in the scratch memory.
  //  d = X(1,1)*X(2,2)-X(1,2)*X(2,1);
  //  p = d'/(d*d');
  //

  X_rd = (const xb_vecNx16*)X;
  Q_rd = (const xb_vecNx16*)Q;
  Q_wr = (      xb_vecNx16*)Q;
  P_wr = (      xb_vecNx16*)pScr;

  for ( l=0; l<L/16; l++ )
  {
    // Formal conversion: CQ15 <- CQ(Q)/2^(15-Q) + 15 - Q
    BBE_LVNX16_IP( x00, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x01, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x02, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x03, X_rd, +4*BBE_SIMD_WIDTH/2 );

    BBE_LVNX16_IP( x10, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x11, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x12, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x13, X_rd, +4*BBE_SIMD_WIDTH/2 );

    // CQ2.30 <- CQ15*CQ15
    w0 = BBE_MULNX16C( x00, x03 );
    w1 = BBE_MULNX16C( x10, x13 );
    // CQ2.30 <- CQ2.30 - CQ15*CQ15
    BBE_MULSNX16C( w0, x01, x02 );
    BBE_MULSNX16C( w1, x11, x12 );

    //  d_exp = nsa40(d)-7;
    //  Q(2.13+d_exp) <- Q30 - [(17+7) - (d_exp+7)] w/o rounding
    vsa0   = BBE_NSANX40C( w0 );             vsa1   = BBE_NSANX40C( w1 );  
	d_exp0 = BBE_MOVVVS( vsa0 );			 d_exp1 = BBE_MOVVVS( vsa1 );
    w0     = BBE_SLANX40( w0, vsa0 );        w1     = BBE_SLANX40( w1, vsa1 );
    d0     = BBE_PACKVNX40( w0, vsa_c24 );   d1     = BBE_PACKVNX40( w1, vsa_c24 );
	

    // Q(5.11+2*d_exp) <- Q(2.13+d_exp)^2 + Q(2.13+d_exp)^2 - 15 w/ rounding
    // 16-bit unsigned; 4.0 <= d2 < 32.0
    d2 = BBE_MAGINX16CPACKQU(d1, d0);

    // 32-by-16 unsigned division
    // Q(18-2*d_exp) <- Q29/Q(5.11+2*d_exp)
    d2u = xb_vecNx16_rtor_xb_vecNx16U( d2 );
    BBE_DIVNX32U( rd2u, rem, c1q29, d2u );
    rd2 = xb_vecNx16U_rtor_xb_vecNx16( rd2u );

    rd2_0 = BBE_SHFLNX16I( rd2, BBE_SHFLI_DUPLICATE_1_EVEN );
    rd2_1 = BBE_SHFLNX16I( rd2, BBE_SHFLI_DUPLICATE_1_ODD  );

    // 16-by-16 bit signed/unsigned multiply
    // CQ(31-d_exp) <- CQ(2.13+d_exp)*Q(18-2*d_exp)
    w0 = BBE_MULUSRNX16( rd2_0, d0, vsa_c16 );
    w1 = BBE_MULUSRNX16( rd2_1, d1, vsa_c16 );
    // CQ(15-d_exp) <- CQ(31-d_exp) - 16 w/ rounding
    p0 = BBE_PACKVNX40( w0, vsa_c16 );
    p1 = BBE_PACKVNX40( w1, vsa_c16 );

    BBE_SVNX16_IP( p0, P_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( p1, P_wr, +4*BBE_SIMD_WIDTH/2 );

    
	t0 = BBE_SELNX16I( d_exp1, d_exp0, BBE_SELI_EXTRACT_1_OF_2_OFF_0 );

    BBE_LVNX16_IP( t1, Q_rd, +2*BBE_SIMD_WIDTH );

    // t0 <- 30+7 - (d_exp+7) = 30 - d_exp
    t0 = BBE_SUBNX16( c37, t0 );
    // t0 <- 30 - d_exp - Q
    t0 = BBE_SUBNX16( t0, t1 );

    BBE_SVNX16_IP( t0, Q_wr, +2*BBE_SIMD_WIDTH );
  }

  //----------------------------------------------------------------------------
  // Compute inverse matrices for L&~15 input matrices.
  //  Y = p*[X(2,2),-X(1,2);-X(2,1),X(1,1)]; 
  //
#ifdef COMPILER_XTENSA // !!!!
  #pragma no_reorder 
#endif
  vsa_c15 = BBE_MOVVSA32( 15 );

  X_rd = (const xb_vecNx16*)X;
  X_wr = (      xb_vecNx16*)X;
  P_rd = (const xb_vecNx16*)( (uintptr_t)P_wr  - 4*(L/16)*16 );

  for ( l=0; l<L/16; l++ )
  {
    // Formal conversion: CQ15 <- CQ(Q)/2^(15-Q) + 15 - Q
    BBE_LVNX16_IP( x00, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x01, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x02, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x03, X_rd, +4*BBE_SIMD_WIDTH/2 );

    BBE_LVNX16_IP( x10, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x11, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x12, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x13, X_rd, +4*BBE_SIMD_WIDTH/2 );

    BBE_LVNX16_IP( p0, P_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( p1, P_rd, +4*BBE_SIMD_WIDTH/2 );

    //
    // X = sat16([X(2,2),-X(1,2);-X(2,1),X(1,1)]);
    //

    {
      xb_vecNx16 x00_ = x00, x01_ = x01, x02_ = x02, x03_ = x03;
      xb_vecNx16 x10_ = x10, x11_ = x11, x12_ = x12, x13_ = x13;

      x00 = x03_;                 x10 = x13_;
      x01 = BBE_NEGSNX16( x01_ ); x11 = BBE_NEGSNX16( x11_ );
      x02 = BBE_NEGSNX16( x02_ ); x12 = BBE_NEGSNX16( x12_ );
      x03 = x00_;                 x13 = x10_;
    }

    // CQ(30-d_exp) <- CQ15*CQ(15-d_exp)
    w00 = BBE_MULRNX16J( x00, p0, vsa_c15 );
    w10 = BBE_MULRNX16J( x10, p1, vsa_c15 );
    w01 = BBE_MULRNX16J( x01, p0, vsa_c15 );
    w11 = BBE_MULRNX16J( x11, p1, vsa_c15 );
    w02 = BBE_MULRNX16J( x02, p0, vsa_c15 );
    w12 = BBE_MULRNX16J( x12, p1, vsa_c15 );
    w03 = BBE_MULRNX16J( x03, p0, vsa_c15 );
    w13 = BBE_MULRNX16J( x13, p1, vsa_c15 );

    // CQ(15-d_exp) <- CQ(30-d_exp) - 15 w/ rounding
    x00 = BBE_PACKVNX40( w00, vsa_c15 );
    x10 = BBE_PACKVNX40( w10, vsa_c15 );
    x01 = BBE_PACKVNX40( w01, vsa_c15 );
    x11 = BBE_PACKVNX40( w11, vsa_c15 );
    x02 = BBE_PACKVNX40( w02, vsa_c15 );
    x12 = BBE_PACKVNX40( w12, vsa_c15 );
    x03 = BBE_PACKVNX40( w03, vsa_c15 );
    x13 = BBE_PACKVNX40( w13, vsa_c15 );

    // Formal conversion:
    // CQ(30-Q-d_exp) <- CQ(15-d_exp)/2^(15-Q) + 30 - Q
    BBE_SVNX16_IP( x00, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x01, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x02, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x03, X_wr, +4*BBE_SIMD_WIDTH/2 );

    BBE_SVNX16_IP( x10, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x11, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x12, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x13, X_wr, +4*BBE_SIMD_WIDTH/2 );
  }

  //----------------------------------------------------------------------------
  // Process the last 8 matrices if L is not a multiple of 16:
  //  d = X(1,1)*X(2,2)-X(1,2)*X(2,1);
  //  Y = d'/(d*d')*[X(2,2),-X(1,2);-X(2,1),X(1,1)]; 
  //

  if ( L&8 )
  {
    // Formal conversion: CQ15 <- CQ(Q)/2^(15-Q) + 15 - Q
    BBE_LVNX16_IP( x00, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x01, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x02, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x03, X_rd, +4*BBE_SIMD_WIDTH/2 );

    // CQ2.30 <- CQ15*CQ15
    w0 = BBE_MULNX16C( x00, x03 );
    // CQ2.30 <- CQ2.30 - CQ15*CQ15
    BBE_MULSNX16C( w0, x01, x02 );

    //  d_exp = nsa40(d)-7;
    //  Q(2.13+d_exp) <- Q30 - [(17+7) - (d_exp+7)] w/o rounding
    vsa0   = BBE_NSANX40C( w0 );
    w0     = BBE_SLLNX40( w0, vsa0 );     
    d0     = BBE_PACKVNX40( w0, vsa_c24 );

    // Q(5.11+2*d_exp) <- Q(2.13+d_exp)^2 + Q(2.13+d_exp)^2 - 15 w/ rounding
    // 16-bit unsigned; 4.0 <= d2 < 32.0
    d2 = BBE_MAGINX16CPACKQU(d0, d0);

    // 32-by-16 unsigned division
    // Q(18-2*d_exp) <- Q29/Q(5.11+2*d_exp)
    d2u = xb_vecNx16_rtor_xb_vecNx16U( d2 );
    	
	// divide half of vector (even elements only) 
	BBE_DIVNX32U_4STEP0_0( c1q29, d2u );
	BBE_DIVNX16U_4STEP_0( d2u );
	BBE_DIVNX16U_4STEP_0( d2u );
	rd2u = BBE_DIVNX16U_4STEPN_0( d2u );

    rd2 = xb_vecNx16U_rtor_xb_vecNx16( rd2u );

    rd2_0 = BBE_SHFLNX16I( rd2, BBE_SHFLI_DUPLICATE_1_EVEN );

    // 16-by-16 bit signed/unsigned multiply
    // CQ(31-d_exp) <- CQ(2.13+d_exp)*Q(18-2*d_exp)
    w0 = BBE_MULUSNX16( rd2_0, d0 );
    // CQ(15-d_exp) <- CQ(31-d_exp) - 16 w/ rounding
    w0 = BBE_RNDSADJNX40( w0, vsa_c16 );
    p0 = BBE_PACKVNX40( w0, vsa_c16 );

    //
    // X = sat16([X(2,2),-X(1,2);-X(2,1),X(1,1)]);
    //

    {
      xb_vecNx16 x00_ = x00, x01_ = x01, x02_ = x02, x03_ = x03;

      x00 = x03_;                
      x01 = BBE_NEGSNX16( x01_ );
      x02 = BBE_NEGSNX16( x02_ );
      x03 = x00_;                
    }

    // CQ(30-d_exp) <- CQ15*CQ(15-d_exp)
    w00 = BBE_MULRNX16J( x00, p0, vsa_c15 );
    w01 = BBE_MULRNX16J( x01, p0, vsa_c15 );
    w02 = BBE_MULRNX16J( x02, p0, vsa_c15 );
    w03 = BBE_MULRNX16J( x03, p0, vsa_c15 );

    // CQ(15-d_exp) <- CQ(30-d_exp) - 15 w/ rounding
    x00 = BBE_PACKVNX40( w00, vsa_c15 );
    x01 = BBE_PACKVNX40( w01, vsa_c15 );
    x02 = BBE_PACKVNX40( w02, vsa_c15 );
    x03 = BBE_PACKVNX40( w03, vsa_c15 );

    // Formal conversion:
    // CQ(30-Q-d_exp) <- CQ(15-d_exp)/2^(15-Q) + 30 - Q
    BBE_SVNX16_IP( x00, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x01, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x02, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x03, X_wr, +4*BBE_SIMD_WIDTH/2 );

	d_exp0 = BBE_MOVVVS( vsa0 );  
    t0 = BBE_SELNX16I( d_exp0, d_exp0, BBE_SELI_EXTRACT_1_OF_2_OFF_0 );

    BBE_LVNX16_IP( t1, Q_rd, +2*BBE_SIMD_WIDTH );

    // t0 <- 30+7 - (d_exp+7) = 30 - d_exp
    t0 = BBE_SUBNX16( c37, t0 );
    // t0 <- 30 - d_exp - Q
    t0 = BBE_SUBNX16( t0, t1 );

    t1 = BBE_SELNX16I( t0, t0, BBE_SELI_ROTATE_RIGHT_4 );

    BBE_SV4X16_I( t0, Q_wr, 0*4*2 );
    BBE_SV4X16_I( t1, Q_wr, 1*4*2 );
  }
#else
  //
  // MATLAB code for a single matrix operation, Q = 15:
  //  % Q30 <- Q15*Q15
  //  d = X(1,1)*X(2,2)-X(1,2)*X(2,1);
  //  d_exp = nsa40(d)-7;
  //  % Q(2.13+d_exp) <- Q30 + d_exp - 17 w/o rounding
  //  d = sat16(floor(d*2^(d_exp-17)));
  //  % Q(5.11+2*d_exp) <- Q(2.13+d_exp)^2 + Q(2.13+d_exp)^2 - 15 w/ rounding
  //  d2 = floor(d*d'*2^-15+0.5); if d2>=2^16, d2 = d2-2^16; end;
  //  % Q(18-2*d_exp) <- Q29/Q(5.11+2*d_exp)
  //  rd2 = sat16u(floor(2^29/d2));
  //  % Q(15-d_exp) <- Q(2.13+d_exp)*Q(18-2*d_exp) - 16 w/ rounding
  //  p = sat16(floor(d*rd2*2^-16+0.5+0.5j));
  //  % Q(15-d_exp+X_exp) <- Q(15-d_exp)*Q15 + X_exp - 15 w/ rounding
  //  X = sat16([X(2,2),-X(1,2);-X(2,1),X(1,1)]);
  //  X = sat16(floor(p'*X*2^-15+0.5+0.5j));
  //  Qi = 15-d_exp;
  //

  const xb_vecNx16 *          X_rd;
        xb_vecNx16 * restrict X_wr;
  const xb_vecNx16 *          Q_rd;
        xb_vecNx16 * restrict Q_wr;
  const xb_vecNx16 *          P_rd;
        xb_vecNx16 * restrict P_wr;

  xb_vecNx16 x00, x01, x02, x03;
  xb_vecNx16 x10, x11, x12, x13;

  xb_vecNx16  d0, d1;
  xb_vecNx16  d2, rd2, rd2_0, rd2_1;
  xb_vecNx16U d2u, rd2u, rem;
  xb_vecNx16  p0, p1;
  xb_vecNx40  w0, w1;

  xb_vecNx16 d_exp0, d_exp1;

  xb_vecNx16 t0, t1;
  xb_vecNx16 c37;
  xb_vecNx40 c1q29;

  vsaN vsa0, vsa1, vsa_c16, vsa_c24;
  valign aQ_wr=BBE_ZALIGN();

  int l;

  NASSERT_ALIGN( X   , (2*BBE_SIMD_WIDTH) );
  NASSERT_ALIGN( pScr, (2*BBE_SIMD_WIDTH) );
  NASSERT_ALIGN( Q   , (2*BBE_SIMD_WIDTH) );

  NASSERT( L%(BBE_SIMD_WIDTH/2)==0 );

  c37 = BBE_MOVVINT16( 37 );

  c1q29 = BBE_MOVWINX40( BBE_MOVWI_Q9_30_1 );
  c1q29 = BBE_SRAINX40( c1q29, 1 );
 
  vsa_c16 = BBE_MOVVSA32( 16 );
  vsa_c24 = BBE_MOVVSA32( 24 );

  //----------------------------------------------------------------------------
  // Compute the determinant and its reciprocal for L&~(BBE_SIMD_WIDTH-1) 
  // matrices, and store results in the scratch memory.
  //  d = X(1,1)*X(2,2)-X(1,2)*X(2,1);
  //  p = d'/(d*d');
  //

  X_rd = (const xb_vecNx16*)X;
  Q_rd = (const xb_vecNx16*)Q;
  Q_wr = (      xb_vecNx16*)Q;
  P_wr = (      xb_vecNx16*)pScr;

  for ( l=0; l<(L>>LOG2_BBE_SIMD_WIDTH); l++ )
  {
    // Formal conversion: CQ15 <- CQ(Q)/2^(15-Q) + 15 - Q
    BBE_LVNX16_IP( x00, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x01, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x02, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x03, X_rd, +4*BBE_SIMD_WIDTH/2 );

    BBE_LVNX16_IP( x10, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x11, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x12, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x13, X_rd, +4*BBE_SIMD_WIDTH/2 );

    // CQ2.30 <- CQ15*CQ15
    w0 = BBE_MULNX16C( x00, x03 );
    w1 = BBE_MULNX16C( x10, x13 );
    // CQ2.30 <- CQ2.30 - CQ15*CQ15
    BBE_MULSNX16C( w0, x01, x02 );
    BBE_MULSNX16C( w1, x11, x12 );

    //  d_exp = nsa40(d)-7;
    //  Q(2.13+d_exp) <- Q30 - [(17+7) - (d_exp+7)] w/o rounding
    vsa0   = BBE_NSANX40C( w0 );             vsa1   = BBE_NSANX40C( w1 );  
    d_exp0 = BBE_MOVVVS( vsa0 );             d_exp1 = BBE_MOVVVS( vsa1 );
    w0     = BBE_SLANX40( w0, vsa0 );        w1     = BBE_SLANX40( w1, vsa1 );
    d0     = BBE_PACKVNX40( w0, vsa_c24 );   d1     = BBE_PACKVNX40( w1, vsa_c24 );

    // Q(5.11+2*d_exp) <- Q(2.13+d_exp)^2 + Q(2.13+d_exp)^2 - 15 w/ rounding
    // 16-bit unsigned; 4.0 <= d2 < 32.0
    d2 = BBE_MAGINX16CPACKQU(d1, d0);

    // 32-by-16 unsigned division
    // Q(18-2*d_exp) <- Q29/Q(5.11+2*d_exp)
    d2u = xb_vecNx16_rtor_xb_vecNx16U( d2 );
    BBE_DIVNX32U( rd2u, rem, c1q29, d2u );
    rd2 = xb_vecNx16U_rtor_xb_vecNx16( rd2u );

    rd2_0 = BBE_SHFLNX16I( rd2, BBE_SHFLI_DUPLICATE_1_EVEN );
    rd2_1 = BBE_SHFLNX16I( rd2, BBE_SHFLI_DUPLICATE_1_ODD  );

    // 16-by-16 bit signed/unsigned multiply
    // CQ(31-d_exp) <- CQ(2.13+d_exp)*Q(18-2*d_exp)
    w0 = BBE_MULUSRNX16( rd2_0, d0, vsa_c16 );
    w1 = BBE_MULUSRNX16( rd2_1, d1, vsa_c16 );
    // CQ(15-d_exp) <- CQ(31-d_exp) - 16 w/ rounding
    p0 = BBE_PACKVNX40( w0, vsa_c16 );
    p1 = BBE_PACKVNX40( w1, vsa_c16 );

    BBE_SVNX16_IP( p0, P_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( p1, P_wr, +4*BBE_SIMD_WIDTH/2 );
    
    t0 = BBE_SELNX16I( d_exp1, d_exp0, BBE_SELI_EXTRACT_1_OF_2_OFF_0 );

    BBE_LVNX16_IP( t1, Q_rd, +2*BBE_SIMD_WIDTH );

    // t0 <- 30+7 - (d_exp+7) = 30 - d_exp
    t0 = BBE_SUBNX16( c37, t0 );
    // t0 <- 30 - d_exp - Q
    t0 = BBE_SUBNX16( t0, t1 );

    BBE_SVNX16_IP( t0, Q_wr, +2*BBE_SIMD_WIDTH );
  }

  __Pragma( "no_reorder" );

  //----------------------------------------------------------------------------
  // Compute inverse matrices for L&~(BBE_SIMD_WIDTH-1) input matrices.
  //  Y = p*[X(2,2),-X(1,2);-X(2,1),X(1,1)]; 
  //

  X_rd = (const xb_vecNx16*)X;
  X_wr = (      xb_vecNx16*)X;
  P_rd = (const xb_vecNx16*)( (uintptr_t)P_wr  - 4*(L/BBE_SIMD_WIDTH)*BBE_SIMD_WIDTH );

  for ( l=0; l<(L>>LOG2_BBE_SIMD_WIDTH); l++ )
  {
    // Formal conversion: CQ15 <- CQ(Q)/2^(15-Q) + 15 - Q
    BBE_LVNX16_IP( x00, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x01, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x02, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x03, X_rd, +4*BBE_SIMD_WIDTH/2 );

    BBE_LVNX16_IP( x10, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x11, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x12, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x13, X_rd, +4*BBE_SIMD_WIDTH/2 );

    BBE_LVNX16_IP( p0, P_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( p1, P_rd, +4*BBE_SIMD_WIDTH/2 );

    //
    // X = sat16([X(2,2),-X(1,2);-X(2,1),X(1,1)]);
    //

    {
      xb_vecNx16 x00_ = x00, x01_ = x01, x02_ = x02, x03_ = x03;
      xb_vecNx16 x10_ = x10, x11_ = x11, x12_ = x12, x13_ = x13;

      x00 = x03_;                 x10 = x13_;
      x01 = BBE_NEGSNX16( x01_ ); x11 = BBE_NEGSNX16( x11_ );
      x02 = BBE_NEGSNX16( x02_ ); x12 = BBE_NEGSNX16( x12_ );
      x03 = x00_;                 x13 = x10_;
    }

    // CQ(15-d_exp) <- CQ15*CQ(15-d_exp) - 15 w/ rounding
    x00 = BBE_MULNX16JPACKQ( x00, p0 );
    x01 = BBE_MULNX16JPACKQ( x01, p0 );
    x02 = BBE_MULNX16JPACKQ( x02, p0 );
    x03 = BBE_MULNX16JPACKQ( x03, p0 );
    x10 = BBE_MULNX16JPACKQ( x10, p1 );
    x11 = BBE_MULNX16JPACKQ( x11, p1 );
    x12 = BBE_MULNX16JPACKQ( x12, p1 );
    x13 = BBE_MULNX16JPACKQ( x13, p1 );

    // Formal conversion:
    // CQ(30-Q-d_exp) <- CQ(15-d_exp)/2^(15-Q) + 30 - Q
    BBE_SVNX16_IP( x00, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x01, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x02, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x03, X_wr, +4*BBE_SIMD_WIDTH/2 );

    BBE_SVNX16_IP( x10, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x11, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x12, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x13, X_wr, +4*BBE_SIMD_WIDTH/2 );
  }

  //----------------------------------------------------------------------------
  // Process the last BBE_SIMD_WIDTH/2 matrices if L is not a multiple of 
  // BBE_SIMD_WIDTH:
  //  d = X(1,1)*X(2,2)-X(1,2)*X(2,1);
  //  Y = d'/(d*d')*[X(2,2),-X(1,2);-X(2,1),X(1,1)]; 
  //

  if ( L&(BBE_SIMD_WIDTH/2) )
  {
    // Formal conversion: CQ15 <- CQ(Q)/2^(15-Q) + 15 - Q
    BBE_LVNX16_IP( x00, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x01, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x02, X_rd, +4*BBE_SIMD_WIDTH/2 );
    BBE_LVNX16_IP( x03, X_rd, +4*BBE_SIMD_WIDTH/2 );

    // CQ2.30 <- CQ15*CQ15
    w0 = BBE_MULNX16C( x00, x03 );
    // CQ2.30 <- CQ2.30 - CQ15*CQ15
    BBE_MULSNX16C( w0, x01, x02 );

    //  d_exp = nsa40(d)-7;
    //  Q(2.13+d_exp) <- Q30 - [(17+7) - (d_exp+7)] w/o rounding
    vsa0   = BBE_NSANX40C( w0 );
    w0     = BBE_SLLNX40( w0, vsa0 );     
    d0     = BBE_PACKVNX40( w0, vsa_c24 );

    // Q(5.11+2*d_exp) <- Q(2.13+d_exp)^2 + Q(2.13+d_exp)^2 - 15 w/ rounding
    // 16-bit unsigned; 4.0 <= d2 < 32.0
    d2 = BBE_MAGINX16CPACKQU(d0, d0);

    // 32-by-16 unsigned division
    // Q(18-2*d_exp) <- Q29/Q(5.11+2*d_exp)
    d2u = xb_vecNx16_rtor_xb_vecNx16U( d2 );
        
    // divide half of vector (even elements only) 
    BBE_DIVNX32U_4STEP0_0( c1q29, d2u );
    BBE_DIVNX16U_4STEP_0( d2u );
    BBE_DIVNX16U_4STEP_0( d2u );
    rd2u = BBE_DIVNX16U_4STEPN_0( d2u );

    rd2 = xb_vecNx16U_rtor_xb_vecNx16( rd2u );

    rd2_0 = BBE_SHFLNX16I( rd2, BBE_SHFLI_DUPLICATE_1_EVEN );

    // 16-by-16 bit signed/unsigned multiply
    // CQ(31-d_exp) <- CQ(2.13+d_exp)*Q(18-2*d_exp)
    w0 = BBE_MULUSNX16( rd2_0, d0 );
    // CQ(15-d_exp) <- CQ(31-d_exp) - 16 w/ rounding
    w0 = BBE_RNDSADJNX40( w0, vsa_c16 );
    p0 = BBE_PACKVNX40( w0, vsa_c16 );

    //
    // X = sat16([X(2,2),-X(1,2);-X(2,1),X(1,1)]);
    //

    {
      xb_vecNx16 x00_ = x00, x01_ = x01, x02_ = x02, x03_ = x03;

      x00 = x03_;                
      x01 = BBE_NEGSNX16( x01_ );
      x02 = BBE_NEGSNX16( x02_ );
      x03 = x00_;                
    }

    // CQ(15-d_exp) <- CQ15*CQ(15-d_exp) - 15 w/ rounding
    x00 = BBE_MULNX16JPACKQ( x00, p0 );
    x01 = BBE_MULNX16JPACKQ( x01, p0 );
    x02 = BBE_MULNX16JPACKQ( x02, p0 );
    x03 = BBE_MULNX16JPACKQ( x03, p0 );

    // Formal conversion:
    // CQ(30-Q-d_exp) <- CQ(15-d_exp)/2^(15-Q) + 30 - Q
    BBE_SVNX16_IP( x00, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x01, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x02, X_wr, +4*BBE_SIMD_WIDTH/2 );
    BBE_SVNX16_IP( x03, X_wr, +4*BBE_SIMD_WIDTH/2 );

    d_exp0 = BBE_MOVVVS( vsa0 );  
    t0 = BBE_SELNX16I( d_exp0, d_exp0, BBE_SELI_EXTRACT_1_OF_2_OFF_0 );

    BBE_LVNX16_IP( t1, Q_rd, +2*BBE_SIMD_WIDTH );

    // t0 <- 30+7 - (d_exp+7) = 30 - d_exp
    t0 = BBE_SUBNX16( c37, t0 );
    // t0 <- 30 - d_exp - Q
    t0 = BBE_SUBNX16( t0, t1 );

    BBE_SAVNX16_XP(t0,aQ_wr,Q_wr,BBE_SIMD_WIDTH);
    BBE_SAPOS_FP(aQ_wr,Q_wr);
  }
#endif
} // cmatinv4x4n_inv2x2bs_q()

#endif
