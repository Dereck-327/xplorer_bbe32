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
 For each group of 8 matrices, merge 2x2 subblocks A, B, C and D into
 4x4 matrices, then convert 8 4x4 matrices from streaming format to block
 format. After that, permute each 4x4 matrix according to its permutation
 index.
  Input:
    A[(L+7)/8][4][8*2]      Upper left 2x2 subblock in block-streaming format
    B[(L+7)/8][4][8*2]      Upper right 2x2 subblock in block-streaming format
    C[(L+7)/8][4][8*2]      Lower left 2x2 subblock in block-streaming format
    D[(L+7)/8][4][8*2]      Lower right 2x2 subblock in block-streaming format
    perm_ix[L]              Permutation index for each matrix, 0..34
    pattern_tbl[35][16]     35-entry permutation pattern table
  Output:
    Y[L][16*2]              4x4 complex matrices in block format
  Note:
    L may be not a multiple of 8.
  Restrictions:
    Y,A,B,C,D  Must be aligned on BBE_SIMD_WIDTH*2-byte boundary
-------------------------------------------------------------------*/
void cmatinv4x4n_permute_bs2b( int16_t * restrict Y,
                         const int16_t * restrict A,
                         const int16_t * restrict B,
                         const int16_t * restrict C,
                         const int16_t * restrict D,
                         const int16_t * restrict perm_ix,
                         const int16_t * restrict pattern_tbl,
                         int L )
{
  const xb_vecNx16 *          Y_rd;
        xb_vecNx16 * restrict Y0_wr;
        xb_vecNx16 * restrict Y1_wr;
  const xb_vecNx16 *          A_rd;
  const xb_vecNx16 *          B_rd;
  const xb_vecNx16 *          C_rd;
  const xb_vecNx16 *          D_rd;
  const xb_vecNx16 *          PTBL;

  int l;

  NASSERT_ALIGN( Y, (BBE_SIMD_WIDTH*2) );
  NASSERT_ALIGN( A, (BBE_SIMD_WIDTH*2) );
  NASSERT_ALIGN( B, (BBE_SIMD_WIDTH*2) );
  NASSERT_ALIGN( C, (BBE_SIMD_WIDTH*2) );
  NASSERT_ALIGN( D, (BBE_SIMD_WIDTH*2) );

  NASSERT(L > 0);

  //--------------------------------------------------------------------------
  // Combine each set of 2x2 subblocks A,B,C,D into a 4x4 matrix. Then convert
  // 4x4 matrices to block format.

  {
    int nb;
    valign Y_va;
    
    xb_vecNx16 x0, x1, x2, x3, x4, x5, x6, x7;
    xb_vecNx16 y00, y01, y10, y11, z00, z01, z10, z11;
    
    
    Y0_wr = (xb_vecNx16*)( (uintptr_t)Y );
    Y1_wr = (xb_vecNx16*)( (uintptr_t)Y + 8*4 );
    A_rd = (xb_vecNx16*)A;
    B_rd = (xb_vecNx16*)B;
    C_rd = (xb_vecNx16*)C;
    D_rd = (xb_vecNx16*)D;

    __Pragma( "ymemory( A_rd )" );
    __Pragma( "ymemory( B_rd )" );
    __Pragma( "ymemory( C_rd )" );
    __Pragma( "ymemory( D_rd )" );
    for ( l=0; l<L/(BBE_SIMD_WIDTH/2); l++ )
    {
      /* Load upper left 2x2 subblock A */
      BBE_LVNX16_IP( y00, A_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( y01, A_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( y10, A_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( y11, A_rd, BBE_SIMD_WIDTH/2*4 );

      /* Load upper right 2x2 subblock B */
      BBE_LVNX16_IP( z00, B_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( z01, B_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( z10, B_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( z11, B_rd, BBE_SIMD_WIDTH/2*4 );

      /* interleave */
      BBE_DSELNX16I(x1, x0, y01, y00, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x3, x2, z01, z00, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x5, x4, y11, y10, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x6, z11, z10, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x2, x0, x2, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x3, x1, x3, x1, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x6, x4, x6, x4, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x5, x7, x5, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x4, x0, x4, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x5, x1, x5, x1, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x6, x2, x6, x2, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x3, x7, x3, BBE_DSELI_DEINTERLEAVE_2);

      /* Save 8 upper halves of SIMD_WIDTH/2 matrices */
      BBE_SVNX16_IP( x0, Y0_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x1, Y0_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x2, Y0_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x3, Y0_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x4, Y0_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x5, Y0_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x6, Y0_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x7, Y0_wr, BBE_SIMD_WIDTH*4 );
    }
    for ( l=0; l<L/(BBE_SIMD_WIDTH/2); l++ )
    {
      /* Load lower left 2x2 subblock C */
      BBE_LVNX16_IP( y00, C_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( y01, C_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( y10, C_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( y11, C_rd, BBE_SIMD_WIDTH/2*4 );

      /* Load lower right 2x2 subblock D */
      BBE_LVNX16_IP( z00, D_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( z01, D_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( z10, D_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( z11, D_rd, BBE_SIMD_WIDTH/2*4 );

      /* interleave */
      BBE_DSELNX16I(x1, x0, y01, y00, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x3, x2, z01, z00, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x5, x4, y11, y10, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x6, z11, z10, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x2, x0, x2, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x3, x1, x3, x1, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x6, x4, x6, x4, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x5, x7, x5, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x4, x0, x4, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x5, x1, x5, x1, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x6, x2, x6, x2, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x3, x7, x3, BBE_DSELI_DEINTERLEAVE_2);

      /* Save 8 lower halves of SIMD_WIDTH/2 matrices */
      BBE_SVNX16_IP( x0, Y1_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x1, Y1_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x2, Y1_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x3, Y1_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x4, Y1_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x5, Y1_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x6, Y1_wr, BBE_SIMD_WIDTH*4 );
      BBE_SVNX16_IP( x7, Y1_wr, BBE_SIMD_WIDTH*4 );
    }

    if ( L & (BBE_SIMD_WIDTH/2-1) )
    {
      xb_vecNx16 x8, x9, xA, xB, xC, xD, xE, xF;
      /* Load upper left 2x2 subblock A */
      BBE_LVNX16_IP( y00, A_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( y01, A_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( y10, A_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( y11, A_rd, BBE_SIMD_WIDTH/2*4 );

      /* Load upper right 2x2 subblock B */
      BBE_LVNX16_IP( z00, B_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( z01, B_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( z10, B_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( z11, B_rd, BBE_SIMD_WIDTH/2*4 );

      /* interleave upper half */
      BBE_DSELNX16I(x2, x0, y01, y00, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x6, x4, z01, z00, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xA, x8, y11, y10, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xE, xC, z11, z10, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x4, x0, x4, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x6, x2, x6, x2, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xC, x8, xC, x8, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xE, xA, xE, xA, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x8, x0, x8, x0, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xA, x2, xA, x2, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xC, x4, xC, x4, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xE, x6, xE, x6, BBE_DSELI_DEINTERLEAVE_2);

      /* Load lower left 2x2 subblock C */
      BBE_LVNX16_IP( y00, C_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( y01, C_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( y10, C_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( y11, C_rd, BBE_SIMD_WIDTH/2*4 );

      /* Load lower right 2x2 subblock D */
      BBE_LVNX16_IP( z00, D_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( z01, D_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( z10, D_rd, BBE_SIMD_WIDTH/2*4 );
      BBE_LVNX16_IP( z11, D_rd, BBE_SIMD_WIDTH/2*4 );

      /* interleave lower half */
      BBE_DSELNX16I(x3, x1, y01, y00, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x5, z01, z00, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xB, x9, y11, y10, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xF, xD, z11, z10, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x5, x1, x5, x1, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x7, x3, x7, x3, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xD, x9, xD, x9, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xF, xB, xF, xB, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(x9, x1, x9, x1, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xB, x3, xB, x3, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xD, x5, xD, x5, BBE_DSELI_DEINTERLEAVE_2);
      BBE_DSELNX16I(xF, x7, xF, x7, BBE_DSELI_DEINTERLEAVE_2);

      /* Number of bytes to be stored to the output array. */
      nb = ( L & (BBE_SIMD_WIDTH/2-1) )*16*4;

      Y_va = BBE_ZALIGN();

      /* Save 1st 4x4 matrix */
      BBE_SAVNX16_XP( x0, Y_va, Y0_wr, nb );
      BBE_SAVNX16_XP( x1, Y_va, Y0_wr, nb );
      nb -= 2*BBE_SIMD_WIDTH/2*4;
      /* Conditionally save next remainder 4x4 matrices */
      BBE_SAVNX16_XP( x2, Y_va, Y0_wr, nb );
      BBE_SAVNX16_XP( x3, Y_va, Y0_wr, nb );
      nb -= 2*BBE_SIMD_WIDTH/2*4;
      BBE_SAVNX16_XP( x4, Y_va, Y0_wr, nb );
      BBE_SAVNX16_XP( x5, Y_va, Y0_wr, nb );
      nb -= 2*BBE_SIMD_WIDTH/2*4;
      BBE_SAVNX16_XP( x6, Y_va, Y0_wr, nb );
      BBE_SAVNX16_XP( x7, Y_va, Y0_wr, nb );
      nb -= 2*BBE_SIMD_WIDTH/2*4;
      BBE_SAVNX16_XP( x8, Y_va, Y0_wr, nb );
      BBE_SAVNX16_XP( x9, Y_va, Y0_wr, nb );
      nb -= 2*BBE_SIMD_WIDTH/2*4;
      BBE_SAVNX16_XP( xA, Y_va, Y0_wr, nb );
      BBE_SAVNX16_XP( xB, Y_va, Y0_wr, nb );
      nb -= 2*BBE_SIMD_WIDTH/2*4;
      BBE_SAVNX16_XP( xC, Y_va, Y0_wr, nb );
      BBE_SAVNX16_XP( xD, Y_va, Y0_wr, nb );

      BBE_SANX16POS_FP(Y_va, Y0_wr);
    }

  }

  __Pragma( "no_reorder" );
  //--------------------------------------------------------------------------
  // Perform backward permutation for each input matrix.

  {
    xb_vecNx16 x0r, x0i;
    xb_vecNx16 y0r, y0i;
    xb_vecNx16 z0, z1;
    xb_vecNx16 p;
    vselN s0;

    Y_rd  = (xb_vecNx16*)Y;
    Y0_wr = (xb_vecNx16*)Y;
    PTBL  = (xb_vecNx16*)pattern_tbl;

    __Pragma( "loop_count min=1" );
    for ( l=0; l<L; l++ )
    {
      /* Load the permutation pattern */
      p = BBE_LVNX16_X( PTBL, perm_ix[l] );

      s0 = BBE_MOVVSV( p, 8 );

      /* Load input matrix; re/im parts are interleaved. */
      BBE_LVNX16_IP( z0, Y_rd, BBE_SIMD_WIDTH*2 );
      BBE_LVNX16_IP( z1, Y_rd, BBE_SIMD_WIDTH*2 );
      
      /* Deinterleave re/im components. */
      BBE_DSELNX16I( x0i, x0r, z1, z0, BBE_DSELI_DEINTERLEAVE_1 );

      /* Apply the permutation. */
      y0r = BBE_SHFLNX16( x0r, s0 );
      y0i = BBE_SHFLNX16( x0i, s0 );

      /* Interleave re/im components. */
      BBE_DSELNX16I( z1, z0, y0i, y0r, BBE_DSELI_INTERLEAVE_1 );

      BBE_SVNX16_IP( z0, Y0_wr, BBE_SIMD_WIDTH/2*4 );
      BBE_SVNX16_IP( z1, Y0_wr, BBE_SIMD_WIDTH/2*4 );
    }
  }

} /* cmatinv4x4n_permute_bs2b() */

#endif
