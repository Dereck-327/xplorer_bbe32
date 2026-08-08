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

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Common utility declarations. */
#include "cmatinvn_common.h"

#if HAVE_HERMMATINV4X4

/*-------------------------------------------------------------------
 Normalize input 4x4 matrices and search for a suboptimal permutation 
 for each matrix, so that the 2x2 subblock at the upper left of a
 permuted matrix has a large determinant.
  Input:
    X[L][16*2]             Input 4x4 complex matrices, CQ(qX)
    cmatinv4x4n_searchPermutation:
        search_tbl[16][2*8]    For each of 16 element locations, enumerates
                               8 variants of the upper left 2x2 subblock
                               selection patterns that relocate the element
                               of interest to that subblock.
    chermmatinv4x4n_searchPermutation:
        search_tbl[8]          Enumerates upper left 2x2 subblock
                               selection patterns for 6 variants of
                               symmetric permutation
  Output:
    X[L][16*2]             Normalized 4x4 complex matrices, CQ(qX+eX)
    eX[L]                  Normalization shift amount for each matrix
    perm_ix[L]             Permutation indices
  Restrictions:
    X  Must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
-------------------------------------------------------------------*/
void chermmatinv4x4n_searchPermutation( int16_t * restrict X,
                                        int16_t * restrict eX,
                                        int16_t * restrict perm_ix,
                                  const int16_t *          search_tbl,
                                        int L )
{
  const xb_vecNx16 *          X_rd;
        xb_vecNx16 * restrict X_wr;
        int16_t    * restrict IX_wr;

  int l;

  NASSERT_ALIGN( X, (2*BBE_SIMD_WIDTH) );

  {
    xb_vecNx16 x00, x01, x10, x11, t0, t1;
    xb_vecNx16 x0r, x1r, x0i, x1i;
    xb_vecNx16 a00r, a01r, a01i, a11r;
    xb_vecNx16 tbl_ix, ix, t, m0, m1;
    xb_vecNx40 w0;
    xb_int16 nsa0, nsa1;
    vsaN sa0, sa1;
    vboolN b0, b1;
    int dum;

    vselN sel00, sel01, sel11, sel_ix;

    vsaN _17 = BBE_MOVVSA32( 17 );

    tbl_ix = BBE_SLLINX16( BBE_SEQNX16(), 5 );
    /* Search table contains selection patterns for 6 variants of the
     * upper left 2x2 subblock. */
    m0 = BBE_LVNX16_I( (xb_vecNx16*)search_tbl, 0 );
    m1 = BBE_MOVVINT16(0x1F);/* mask to extract permutation patterns */

    t = BBE_ANDNX16(m0, m1); m0 = BBE_SRLINX16(m0, 5);
    sel00 = BBE_MOVVSV( t, 0 ); /* 6 variants of element (0,0) */
    t = BBE_ANDNX16(m0, m1);
    sel01 = BBE_MOVVSV( t, 0 ); /* 6 variants of element (0,1) */
    sel11 = BBE_MOVVSV( m0, 5 ); /* 6 variants of element (1,1) */

    X_rd  = (xb_vecNx16*)X;
    X_wr  = (xb_vecNx16*)X;
    IX_wr = perm_ix;
    
    for ( l=0; l<L/2; l++ )
    {
      /* Load 2 4x4 complex matrices, CQ(qX) */
      BBE_LVNX16_IP( x00, X_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( x01, X_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( x10, X_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( x11, X_rd, BBE_SIMD_WIDTH/2*4 );

      /* Deinterleave re/im parts. */
      BBE_DSELNX16I( x0i, x0r, x01, x00, BBE_DSELI_DEINTERLEAVE_1 );
      BBE_DSELNX16I( x1i, x1r, x11, x10, BBE_DSELI_DEINTERLEAVE_1 );

      /* Reduced max absolute values. */
      BBE_BMAXABSNX16( b0, t0, x01, x00 );
      BBE_BMAXABSNX16( b1, t1, x11, x10 );
      nsa0 = BBE_RMAXNX16( t0 );
      nsa1 = BBE_RMAXNX16( t1 );

      /* Determine the normalization shift amount. */
      sa0 = BBE_NSANX16( BBE_REPNX16( BBE_MOVNX16_FROM16(nsa0), 0 ) );
      sa1 = BBE_NSANX16( BBE_REPNX16( BBE_MOVNX16_FROM16(nsa1), 0 ) );
      nsa0 = BBE_MOVVVS(sa0);
      nsa1 = BBE_MOVVVS(sa1);

      /* Save the exponent. */
      BBE_SSNX16_IP( nsa0, eX, 2 );
      BBE_SSNX16_IP( nsa1, eX, 2 );

      /* Normalization with formal conversion to Q15:
       * CQ15 <- [CQ(qX)+eX]/2^(15-qX-eX) + 15 - qX - eX */
      x00 = BBE_SLANX16( x0r, sa0 );
      x01 = BBE_SLANX16( x0i, sa0 );
      x10 = BBE_SLANX16( x1r, sa1 );
      x11 = BBE_SLANX16( x1i, sa1 );

      /* Save deinterleaved and normalized data */
      BBE_SVNX16_IP( x00, X_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( x01, X_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( x10, X_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( x11, X_wr, BBE_SIMD_WIDTH/2*4 );

      /* Retrieve elements of the subblock. Diagonal elements (0,0) and (1,1) 
       * are pure real, off-diagonal element (1,0) is conjugate of (1,0). */
      a00r = BBE_SELNX16( x1r, x0r, sel00 );
      a01r = BBE_SELNX16( x1r, x0r, sel01 );
      a01i = BBE_SELNX16( x1i, x0i, sel01 );
      a11r = BBE_SELNX16( x1r, x0r, sel11 );

      /* Calculate determinant for each of 6 variants.
       * Q2.30 <- Q15*Q15 - Q15*Q15 - Q15*Q15 */
      w0 = BBE_MULRNX16( a00r, a11r, _17 );
      BBE_MULSNX16 ( w0, a01r, a01r );
      BBE_MULSNX16 ( w0, a01i, a01i );
      /* Q2.13 <- Q2.30 - 17 */
      t = BBE_PACKVNX40( w0, _17 );

      /*
       * Take absolute value and find the maximum.
       * 't' contains results for 2 matrices, searching is performed separately.
       */
      t = BBE_ABSSNX16(t);
      m0 = BBE_SELNX16I( t, t, BBE_SELI_EXTRACT_LO_HALVES );
      m1 = BBE_SELNX16I( t, t, BBE_SELI_EXTRACT_HI_HALVES );


      /* Retrieve the permutation index that maximizes the absolute value
       * of the determinant. */
      /* 1st matrix */
      BBE_RBMAXNX16( b0, m0, m0 );
      BBE_SQZN( sel_ix, dum, b0 );
      ix = BBE_SHFLNX16( tbl_ix, sel_ix ); 
      BBE_SSNX16_IP( ix, IX_wr, 2 );
      /* 2nd matrix */
      BBE_RBMAXNX16( b0, m1, m1 );
      BBE_SQZN( sel_ix, dum, b0 );
      ix = BBE_SHFLNX16( tbl_ix, sel_ix ); 
      BBE_SSNX16_IP( ix, IX_wr, 2 );
    }
    if (L&1)
    {
      /* Load 4x4 complex matrix, CQ(qX) */
      BBE_LVNX16_IP( x00, X_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( x01, X_rd, BBE_SIMD_WIDTH/2*4 );

      /* Deinterleave re/im parts. */
      BBE_DSELNX16I( x0i, x0r, x01, x00, BBE_DSELI_DEINTERLEAVE_1 );

      /* Reduced max absolute values. */
      BBE_BMAXABSNX16( b0, t0, x01, x00 );
      nsa0 = BBE_RMAXNX16( t0 );

      /* Determine the normalization shift amount. */
      sa0 = BBE_NSANX16( BBE_REPNX16( BBE_MOVNX16_FROM16(nsa0), 0 ) );
      nsa0 = BBE_MOVVVS(sa0);

      /* Save the exponent. */
      BBE_SSNX16_IP( nsa0, eX, 2 );

      /* Normalization with formal conversion to Q15:
       * CQ15 <- [CQ(qX)+eX]/2^(15-qX-eX) + 15 - qX - eX */
      x00 = BBE_SLANX16( x0r, sa0 );
      x01 = BBE_SLANX16( x0i, sa0 );

      /* Save deinterleaved and normalized data */
      BBE_SVNX16_IP( x00, X_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( x01, X_wr, BBE_SIMD_WIDTH/2*4 );

      /* Retrieve elements of the subblock. Diagonal elements (0,0) and (1,1) 
       * are pure real, off-diagonal element (1,0) is conjugate of (1,0). */
      a00r = BBE_SHFLNX16( x0r, sel00 );
      a01r = BBE_SHFLNX16( x0r, sel01 );
      a01i = BBE_SHFLNX16( x0i, sel01 );
      a11r = BBE_SHFLNX16( x0r, sel11 );

      /* Calculate determinant for each of 6 variants.
       * Q2.30 <- Q15*Q15 - Q15*Q15 - Q15*Q15 */
      w0 = BBE_MULRNX16( a00r, a11r, _17 );
      BBE_MULSNX16 ( w0, a01r, a01r );
      BBE_MULSNX16 ( w0, a01i, a01i );
      /* Q2.13 <- Q2.30 - 17 */
      t = BBE_PACKVNX40( w0, _17 );

      /*
       * Take absolute value and find the maximum.
       * 't' contains results for 2 matrices, searching is performed separately.
       */
      t = BBE_ABSSNX16(t);

      /* Retrieve the permutation index that maximizes the absolute value
       * of the determinant. */
      BBE_RBMAXNX16( b0, t, t );
      BBE_SQZN( sel_ix, dum, b0 );
      ix = BBE_SHFLNX16( BBE_SLLINX16( BBE_SEQNX16(), 5 ), sel_ix ); 
      BBE_SSNX16_IP( ix, IX_wr, 2 );
    }
    (void)dum;
  }

} /* chermmatinv4x4n_searchPermutation() */

#endif /* HAVE_VSAMATH && HAVE_DIV */
