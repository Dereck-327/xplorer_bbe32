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
 Apply permutation to each 4x4 complex matrix, then convert the whole
 bulk of matrices to block-streaming format: each group of 8 consecutive
 matrices is converted to streaming format [4*4][8*2], and that groups
 are stored sequentally, i.e. we arrive to format [L/8][4*4][8*2]. Finally,
 2x2 subblocks are separated into 4 arrays A, B, C and D: [L/8][2*2][8*2].
  Input:
    X[L][16*2]              4x4 complex matrices in block format
    perm_ix[L]              Permutation index for each matrix, 0..34
    pattern_tbl[35][16]     35-entry permutation pattern table
  Output:
    A[(L+7)/8][4][8*2]      Upper left 2x2 subblock in block-streaming format
    B[(L+7)/8][4][8*2]      Upper right 2x2 subblock in block-streaming format
    C[(L+7)/8][4][8*2]      Lower left 2x2 subblock in block-streaming format
    D[(L+7)/8][4][8*2]      Lower right 2x2 subblock in block-streaming format
  Note:
    L may be not a multiple of 8, in which case the last block is complemented
    with the respective volume of dummy data.
  Restrictions:
    A,B,C,D,X must be aligned on 2*BBE_SIMD_WIDTH-byte boundary
-------------------------------------------------------------------*/
void cmatinv4x4n_permute_b2bs( int16_t * restrict A,
                               int16_t * restrict B,
                               int16_t * restrict C,
                               int16_t * restrict D,
                               int16_t * restrict X,
                         const int16_t * restrict perm_ix,
                         const int16_t * restrict pattern_tbl,
                         int L )
{
  const xb_vecNx16 *          X0_rd;
  const xb_vecNx16 *          X1_rd;
        xb_vecNx16 * restrict X_wr;
  const xb_vecNx16 *          PTBL;
        xb_vecNx16 * restrict A_wr;
        xb_vecNx16 * restrict B_wr;
        xb_vecNx16 * restrict C_wr;
        xb_vecNx16 * restrict D_wr;

  int l;

  NASSERT_ALIGN( A, BBE_SIMD_WIDTH*2 );
  NASSERT_ALIGN( B, BBE_SIMD_WIDTH*2 );
  NASSERT_ALIGN( C, BBE_SIMD_WIDTH*2 );
  NASSERT_ALIGN( D, BBE_SIMD_WIDTH*2 );
  NASSERT_ALIGN( X, BBE_SIMD_WIDTH*2 );

  NASSERT_ALIGN( pattern_tbl, BBE_SIMD_WIDTH*2 );
  NASSERT(L > 0);

  //--------------------------------------------------------------------------
  // Perform forward permutation for each input matrix.

  {
    xb_vecNx16 x0r, x0i;
    xb_vecNx16 y0r, y0i;
    xb_vecNx16 z0, z1;
    xb_vecNx16 p;
    vselN s0;

    X0_rd = (xb_vecNx16*)X;
    X_wr  = (xb_vecNx16*)X;
    PTBL  = (xb_vecNx16*)pattern_tbl;

    __Pragma( "loop_count min=1" );
    for ( l=0; l<L; l++ )
    {
      /* Load the permutation pattern */
      p = BBE_LVNX16_X( PTBL, perm_ix[l] );

      s0 = BBE_MOVVSV( p, 0 );

      /* Load input matrix; re/im parts are deinterleaved. */
      BBE_LVNX16_IP( x0r, X0_rd, BBE_SIMD_WIDTH*2 );
      BBE_LVNX16_IP( x0i, X0_rd, BBE_SIMD_WIDTH*2 );

      /* Apply the permutation. */
      y0r = BBE_SHFLNX16( x0r, s0 );
      y0i = BBE_SHFLNX16( x0i, s0 );

      /* Interleave re/im components. */
      BBE_DSELNX16I( z1, z0, y0i, y0r, BBE_DSELI_INTERLEAVE_1 );

      BBE_SVNX16_IP( z0, X_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z1, X_wr, BBE_SIMD_WIDTH/2*4 );
    }

  }

  __Pragma( "no_reorder" );
  //--------------------------------------------------------------------------
  // Break each 4x4 matrix into 4 2x2 subblocks A,B,C,D, and convert them to
  // block-stremaing format.

  {
    int nb;
    valign X_va;

    xb_vecNx16 x0, x1, x2, x3, x4, x5, x6, x7;
    xb_vecNx16 y00, y01, y10, y11, z00, z01, z10, z11;

    A_wr = (xb_vecNx16*)A;
    B_wr = (xb_vecNx16*)B;
    C_wr = (xb_vecNx16*)C;
    D_wr = (xb_vecNx16*)D;
    
    X0_rd = (xb_vecNx16*)( (uintptr_t)X );
    X1_rd = (xb_vecNx16*)( (uintptr_t)X + 8*4 );

    for ( l=0; l<L/(BBE_SIMD_WIDTH/2); l++ )
    {
      /* Load upper half (2x4) of SIMD_WIDTH/2 matrices */
      BBE_LVNX16_IP( x0, X0_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x1, X0_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x2, X0_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x3, X0_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x4, X0_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x5, X0_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x6, X0_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x7, X0_rd, BBE_SIMD_WIDTH*4 );

      /* interleave */
      BBE_DSELNX16I(x1, x0, x1, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x3, x2, x3, x2, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x5, x4, x5, x4, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x6, x7, x6, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x2, x0, x2, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x3, x1, x3, x1, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x6, x4, x6, x4, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x5, x7, x5, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(y10, y00, x4, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(y11, y01, x5, x1, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(z10, z00, x6, x2, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(z11, z01, x7, x3, BBE_DSELI_DEINTERLEAVE_2);

      /* Save upper left 2x2 subblock A */
      BBE_SVNX16_IP( y00, A_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( y01, A_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( y10, A_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( y11, A_wr, BBE_SIMD_WIDTH/2*4 );

      /* Save upper right 2x2 subblock B */
      BBE_SVNX16_IP( z00, B_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z01, B_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z10, B_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z11, B_wr, BBE_SIMD_WIDTH/2*4 );
    }
    for ( l=0; l<L/(BBE_SIMD_WIDTH/2); l++ )
    {
      /* Load lower half (2x4) of SIMD_WIDTH/2 matrices */
      BBE_LVNX16_IP( x0, X1_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x1, X1_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x2, X1_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x3, X1_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x4, X1_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x5, X1_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x6, X1_rd, BBE_SIMD_WIDTH*4 );
      BBE_LVNX16_IP( x7, X1_rd, BBE_SIMD_WIDTH*4 );

      /* interleave */
      BBE_DSELNX16I(x1, x0, x1, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x3, x2, x3, x2, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x5, x4, x5, x4, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x6, x7, x6, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x2, x0, x2, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x3, x1, x3, x1, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x6, x4, x6, x4, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x5, x7, x5, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(y10, y00, x4, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(y11, y01, x5, x1, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(z10, z00, x6, x2, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(z11, z01, x7, x3, BBE_DSELI_DEINTERLEAVE_2);

      /* Save lower left 2x2 subblock C */
      BBE_SVNX16_IP( y00, C_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( y01, C_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( y10, C_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( y11, C_wr, BBE_SIMD_WIDTH/2*4 );

      /* Save lower right 2x2 subblock D */
      BBE_SVNX16_IP( z00, D_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z01, D_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z10, D_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z11, D_wr, BBE_SIMD_WIDTH/2*4 );
    }
    if ( L & (BBE_SIMD_WIDTH/2-1) )
    {
      xb_vecNx16 x8, x9, xA, xB, xC, xD, xE, xF;
      /* Number of bytes left in the input array. */
      nb = ( L & (BBE_SIMD_WIDTH/2-1) )*16*4;

      X_va = BBE_LAVNX16_PP(X0_rd);

      /* Load 1st 4x4 matrix */
      BBE_LAVNX16_XP( x0, X_va, X0_rd, nb );
      BBE_LAVNX16_XP( x1, X_va, X0_rd, nb );
      nb -= 2*BBE_SIMD_WIDTH/2*4;
      /* Conditionally load next remainder 4x4 matrices */
      BBE_LAVNX16_XP( x2, X_va, X0_rd, nb );
      BBE_LAVNX16_XP( x3, X_va, X0_rd, nb );
      nb -= 2*BBE_SIMD_WIDTH/2*4;
      BBE_LAVNX16_XP( x4, X_va, X0_rd, nb );
      BBE_LAVNX16_XP( x5, X_va, X0_rd, nb );
      nb -= 2*BBE_SIMD_WIDTH/2*4;
      BBE_LAVNX16_XP( x6, X_va, X0_rd, nb );
      BBE_LAVNX16_XP( x7, X_va, X0_rd, nb );
      nb -= 2*BBE_SIMD_WIDTH/2*4;
      BBE_LAVNX16_XP( x8, X_va, X0_rd, nb );
      BBE_LAVNX16_XP( x9, X_va, X0_rd, nb );
      nb -= 2*BBE_SIMD_WIDTH/2*4;
      BBE_LAVNX16_XP( xA, X_va, X0_rd, nb );
      BBE_LAVNX16_XP( xB, X_va, X0_rd, nb );
      nb -= 2*BBE_SIMD_WIDTH/2*4;
      BBE_LAVNX16_XP( xC, X_va, X0_rd, nb );
      BBE_LAVNX16_XP( xD, X_va, X0_rd, nb );
      nb -= 2*BBE_SIMD_WIDTH/2*4;

      xE = xF = BBE_ZERONX16();

      /* interleave upper half */
      BBE_DSELNX16I(x2, x0, x2, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x6, x4, x6, x4, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xA, x8, xA, x8, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xE, xC, xE, xC, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x4, x0, x4, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x6, x2, x6, x2, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xC, x8, xC, x8, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xE, xA, xE, xA, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(y10, y00, x8, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(y11, y01, xA, x2, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(z10, z00, xC, x4, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(z11, z01, xE, x6, BBE_DSELI_DEINTERLEAVE_2);

      /* Save upper left 2x2 subblock A */
      BBE_SVNX16_IP( y00, A_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( y01, A_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( y10, A_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( y11, A_wr, BBE_SIMD_WIDTH/2*4 );

      /* Save upper right 2x2 subblock B */
      BBE_SVNX16_IP( z00, B_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z01, B_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z10, B_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z11, B_wr, BBE_SIMD_WIDTH/2*4 );

      /* interleave lower half */
      BBE_DSELNX16I(x3, x1, x3, x1, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x5, x7, x5, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xB, x9, xB, x9, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xF, xD, xF, xD, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x5, x1, x5, x1, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x3, x7, x3, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xD, x9, xD, x9, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xF, xB, xF, xB, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(y10, y00, x9, x1, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(y11, y01, xB, x3, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(z10, z00, xD, x5, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(z11, z01, xF, x7, BBE_DSELI_DEINTERLEAVE_2);

      /* Save lower left 2x2 subblock C */
      BBE_SVNX16_IP( y00, C_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( y01, C_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( y10, C_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( y11, C_wr, BBE_SIMD_WIDTH/2*4 );

      /* Save lower right 2x2 subblock D */
      BBE_SVNX16_IP( z00, D_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z01, D_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z10, D_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z11, D_wr, BBE_SIMD_WIDTH/2*4 );
    }
  }

} /* cmatinv4x4n_permute_b2bs() */

#endif
