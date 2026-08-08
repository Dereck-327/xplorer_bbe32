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
void cmatinv4x4n_searchPermutation( int16_t * restrict X,
                                    int16_t * restrict eX,
                                    int16_t * restrict perm_ix,
                              const int16_t *          search_tbl,
                                    int L )
{
  const xb_vecNx16 *          X_rd;
        xb_vecNx16 * restrict X_wr;
  const int16_t    * restrict IX_rd;
        int16_t    * restrict IX_wr;
  const xb_vecNx16 *          STBL;

  int l;

  NASSERT_ALIGN( X, (2*BBE_SIMD_WIDTH) );

  //--------------------------------------------------------------------------
  // For each input matrix: normalize data and find the element of greatest 
  // magnitude.

  {
    xb_vecNx16 x0, x1, z0, z1, p0;
    xb_vecNx16 sel_ix, ix, out;
    xb_vecNx40 w0;
    vsaN sa0, sa1, _16;
    xb_int16 nsa;
    vboolN b0;

    static const int16_t ALIGN(32) sel_ix_i[16] = { 0*32,8*32,1*32,9*32,2*32,10*32,3*32,11*32,4*32,12*32,5*32,13*32,6*32,14*32,7*32,15*32 };

    sel_ix = BBE_LVNX16_I( (xb_vecNx16*)sel_ix_i, 0 );
    out = BBE_MOVQINT16(1);

    X_rd  = (xb_vecNx16*)X;
    X_wr  = (xb_vecNx16*)X;
    IX_wr = perm_ix;

    _16 = BBE_MOVVSA32(16);

    __Pragma("loop_count min=1");
    for ( l=0; l<L; l++ )
    {
      /* Load 4x4 complex matrix, CQ(qX) */
      BBE_LVNX16_IP( x0, X_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( x1, X_rd, BBE_SIMD_WIDTH/2*4 );

      /*
       * Find element of maximal magnitude.
       */

      /* Q2.30 <- Q15*Q15 + Q15*Q15 */
      w0 = BBE_MAGIRNX16C( x1, x0, _16 );
      /* Q1.14 <- Q2.30 - 16 w/ rounding */
      p0 = BBE_PACKVNX40( w0, _16 );

      BBE_RBMAXNX16( b0, p0, p0 );

      /* Determine the index of the maximal element.
       * Use it as byte offset in the search_tbl. */
      ix = BBE_MOVNX16T( sel_ix, out, b0 );
      ix = BBE_RMINNX16( ix );

      BBE_SSNX16_IP( ix, IX_wr, 2 );

      /*
       * Normalize data.
       */

      sa0 = BBE_NSANX16( x0 );
      sa1 = BBE_NSANX16( x1 );

      sa0 = BBE_MINVSN( sa1, sa0 );
      nsa = BBE_RMINNX16( BBE_MOVVVS(sa0) );

      sa0 = BBE_MOVVSV( BBE_REPNX16( BBE_MOVNX16_FROM16(nsa), 0 ), 0 );

      /* Normalization with formal conversion to Q15:
       * CQ15 <- [CQ(qX)+eX]/2^(15-qX-eX) + 15 - qX - eX */
      x0 = BBE_SLLNX16( x0, sa0 );
      x1 = BBE_SLLNX16( x1, sa0 );

      /* Save the exponent. */
      BBE_SSNX16_IP( nsa, eX, 2 );

      /* Deinterleave re/im parts. */
      BBE_DSELNX16I( z1, z0, x1, x0, BBE_DSELI_DEINTERLEAVE_1 );

      /* Save deinterleaved and normalized data */
      BBE_SVNX16_IP( z0, X_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z1, X_wr, BBE_SIMD_WIDTH/2*4 );
    }

  }

  __Pragma( "no_reorder" );
  //--------------------------------------------------------------------------
  // Try 8 permutations that move the element of maximum amplitude to
  // the upper left 2x2 subblock. For each permutation compute the
  // determinant of upper left 2x2 subblock, and select the permutation
  // that provides the maximum absolute value.
  //

  {
    xb_vecNx16 x0r, x1r, x0i, x1i;
    xb_vecNx16 a0r, a0i, a1r, a1i, a2r, a2i, a3r, a3i;
    xb_vecNx16 s0, d0, d1, mag, ix, zero;
    xb_vecNx16 mask0;
    xb_vecNx40 w0, w1;
    vboolN b0;
    vselN sel0, sel1, sel2, sel3, sel_ix, c16;
    int dum, off0, off1;

    vsaN _17 = BBE_MOVVSA32(17);
    zero = BBE_ZERONX16();
    mask0 = BBE_MOVVA16(0x0F0F);
    s0 = BBE_MOVVINT16(0x10);
    s0 = BBE_SELNX16I( s0, zero, BBE_SELI_EXTRACT_HI_HALVES );
    c16 = BBE_MOVVSV(s0, 0);

    X_rd  = (xb_vecNx16*)X;
    X_wr  = (xb_vecNx16*)X;
    STBL  = (xb_vecNx16*)search_tbl;
    IX_wr = perm_ix;
    IX_rd = IX_wr;

    for ( l=0; l<(L/2); l++ )
    {
      /* Load 2 4x4 complex matrices, CQ15 */
      BBE_LVNX16_IP( x0r, X_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( x0i, X_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( x1r, X_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( x1i, X_rd, BBE_SIMD_WIDTH/2*4 );

      /*
       * Load selection patterns to check 8 variants of upper left subblock.
       */

      off0 = *IX_rd++;
      off1 = *IX_rd++;

      d0 = BBE_LVNX16_X( STBL, off0 );
      d1 = BBE_LVNX16_X( STBL, off1 );
      s0 = BBE_SELNX16I( d1, d0, BBE_SELI_EXTRACT_LO_HALVES );
      d0 = BBE_ANDNX16(s0, mask0); s0 = BBE_SRLINX16(s0, 4);
      d1 = BBE_ANDNX16(s0, mask0);
      sel0 = BBE_MOVVSV( d0, 0 );
      sel1 = BBE_MOVVSV( d1, 0 );
      sel2 = BBE_MOVVSV( d0, 8 );
      sel3 = BBE_MOVVSV( d1, 8 );
      /* add offset for correct permutation of 2nd matrix */
      sel0 = BBE_ADDSVSN(sel0, c16);
      sel1 = BBE_ADDSVSN(sel1, c16);
      sel2 = BBE_ADDSVSN(sel2, c16);
      sel3 = BBE_ADDSVSN(sel3, c16);

      /*
       * Form 8 variants of upper left subblock.
       */

      a0r = BBE_SELNX16( x1r, x0r, sel0 );
      a0i = BBE_SELNX16( x1i, x0i, sel0 );
      a1r = BBE_SELNX16( x1r, x0r, sel1 );
      a1i = BBE_SELNX16( x1i, x0i, sel1 );
      a2r = BBE_SELNX16( x1r, x0r, sel2 );
      a2i = BBE_SELNX16( x1i, x0i, sel2 );
      a3r = BBE_SELNX16( x1r, x0r, sel3 );
      a3i = BBE_SELNX16( x1i, x0i, sel3 );

      /*
       * Compute squared magnitudes of 8 determinants.
       */

      /* CQ2.30 <- CQ15*CQ15 - CQ15*CQ15 */
      w0 = BBE_MULRNX16( a0r, a3r, _17 );
      BBE_MULSNX16 ( w0, a0i, a3i );
      BBE_MULSNX16 ( w0, a1r, a2r );
      BBE_MULANX16 ( w0, a1i, a2i );

      w1 = BBE_MULRNX16( a0r, a3i, _17 );
      BBE_MULANX16 ( w1, a0i, a3r );
      BBE_MULSNX16 ( w1, a1r, a2i );
      BBE_MULSNX16 ( w1, a1i, a2r );

      /* CQ2.13 <- CQ2.30 - 17 */
      d0 = BBE_PACKVNX40( w0, _17 );
      d1 = BBE_PACKVNX40( w1, _17 );

      /* Q5.26 <- Q2.13*Q2.13 + Q2.13*Q2.13 */
      w0 = BBE_MULNX16( d0, d0 );
      BBE_MULANX16( w0, d1, d1 );
      /* 2x2 determinant does not exceed 4 in its magnitude, so it's safe to use
       * signed Q4.11 for the squared magnitude. */
      /* Q4.11 <- sature16( Q5.26 - 15 ) */
      mag = BBE_PACKQNX40( w0 );

      /*
       * Find the greatest magnitude and save the respective permutation index.
       * 'mag' contains results for 2 matrices, searching is performed separately.
       */
      d0 = BBE_SELNX16I( mag, mag, BBE_SELI_EXTRACT_LO_HALVES );
      d1 = BBE_SELNX16I( mag, mag, BBE_SELI_EXTRACT_HI_HALVES );
      
      /* 1st matrix */
      s0 = BBE_LVNX16_X( STBL, off0 );
      ix = BBE_SELNX16I( s0, s0, BBE_SELI_EXTRACT_HI_HALVES );
      BBE_RBMAXNX16( b0, d0, d0 );
      BBE_SQZN( sel_ix, dum, b0 );
      ix = BBE_SHFLNX16( ix, sel_ix );
      BBE_SSNX16_IP( ix, IX_wr, 2 );
      /* 2nd matrix */
      s0 = BBE_LVNX16_X( STBL, off1 );
      ix = BBE_SELNX16I( s0, s0, BBE_SELI_EXTRACT_HI_HALVES );
      BBE_RBMAXNX16( b0, d1, d1 );
      BBE_SQZN( sel_ix, dum, b0 );
      ix = BBE_SHFLNX16( ix, sel_ix );
      BBE_SSNX16_IP( ix, IX_wr, 2 );
    }
    if (L&1)
    {
      /* Load 4x4 complex matrix, CQ15 */
      BBE_LVNX16_IP( x0r, X_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( x0i, X_rd, BBE_SIMD_WIDTH/2*4 );

      /*
       * Load selection patterns to check 8 variants of upper left subblock.
       */
      
      off0 = *IX_rd;

      s0 = BBE_LVNX16_X( STBL, off0 );
      ix = BBE_SELNX16I( zero, s0, BBE_SELI_EXTRACT_HI_HALVES );
      s0 = BBE_SELNX16I( zero, s0, BBE_SELI_EXTRACT_LO_HALVES );

      {
          xb_vecNx16 mask, t0, t1;
          mask = BBE_MOVVA16(0x0F0F);
          t0 = BBE_ANDNX16(s0, mask); s0 = BBE_SRLINX16(s0, 4);
          t1 = BBE_ANDNX16(s0, mask);
          sel0 = BBE_MOVVSV( t0, 0 );
          sel1 = BBE_MOVVSV( t1, 0 );
          sel2 = BBE_MOVVSV( t0, 8 );
          sel3 = BBE_MOVVSV( t1, 8 );
      }

      /*
       * Form 8 variants of upper left subblock.
       */

      a0r = BBE_SHFLNX16( x0r, sel0 );
      a0i = BBE_SHFLNX16( x0i, sel0 );
      a1r = BBE_SHFLNX16( x0r, sel1 );
      a1i = BBE_SHFLNX16( x0i, sel1 );
      a2r = BBE_SHFLNX16( x0r, sel2 );
      a2i = BBE_SHFLNX16( x0i, sel2 );
      a3r = BBE_SHFLNX16( x0r, sel3 );
      a3i = BBE_SHFLNX16( x0i, sel3 );

      /*
       * Compute squared magnitudes of 8 determinants.
       */

      /* CQ2.30 <- CQ15*CQ15 - CQ15*CQ15 */
      w0 = BBE_MULRNX16( a0r, a3r, _17 );
      BBE_MULSNX16 ( w0, a0i, a3i );
      BBE_MULSNX16 ( w0, a1r, a2r );
      BBE_MULANX16 ( w0, a1i, a2i );

      w1 = BBE_MULRNX16( a0r, a3i, _17 );
      BBE_MULANX16 ( w1, a0i, a3r );
      BBE_MULSNX16 ( w1, a1r, a2i );
      BBE_MULSNX16 ( w1, a1i, a2r );

      /* CQ2.13 <- CQ2.30 - 17 */
      d0 = BBE_PACKVNX40( w0, _17 );
      d1 = BBE_PACKVNX40( w1, _17 );

      /* Q5.26 <- Q2.13*Q2.13 + Q2.13*Q2.13 */
      w0 = BBE_MULNX16( d0, d0 );
      BBE_MULANX16( w0, d1, d1 );
      /* 2x2 determinant does not exceed 4 in its magnitude, so it's safe to use
       * signed Q4.11 for the squared magnitude. */
      /* Q4.11 <- sature16( Q5.26 - 15 ) */
      mag = BBE_PACKQNX40( w0 );

      /*
       * Find the greatest magnitude and save the respective permutation index.
       */

      BBE_RBMAXNX16( b0, mag, mag );
      BBE_SQZN( sel_ix, dum, b0 );
      ix = BBE_SHFLNX16( ix, sel_ix );

      BBE_SSNX16_IP( ix, IX_wr, 2 );
    }
    (void)dum;
  }

} /* cmatinv4x4n_searchPermutation() */

#endif
