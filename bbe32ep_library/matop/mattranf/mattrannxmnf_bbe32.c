/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
  NatureDSP_Baseband library. Matrix Operations
    Real Matrix Transpose
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"

/*----------------------------------------------------*/
/* Restricted implementations of matrix transpose     */

/* N is multiple of 4 and M is multiple of 4 */
static void mattrannxmnf_N4M4 ( float32_t * restrict y, 
                                float32_t * restrict x, 
                                int N, int M, int L );
/* N is multiple of 4 */
static void mattrannxmnf_N4 ( float32_t * restrict y, 
                              float32_t * restrict x, 
                              int N, int M, int L );
/* M is multiple of 4 */
static void mattrannxmnf_M4 ( float32_t * restrict y, 
                              float32_t * restrict x, 
                              int N, int M, int L );
/* N is multiple of 8 and M is multiple of 8 */
static void mattrannxmnf_N8M8 ( float32_t * restrict y, 
                                float32_t * restrict x, 
                                int N, int M, int L );

/*----------------------------------------------------*/

/*-------------------------------------------------------------------------
Real Matrix Transpose

Description: These functions perform transposition for each matrix from input
sequence and store results to output sequence. Both the block order and
streaming order are allowed for input/output matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Note:
The functions mattrannxnn(), mattrannxmn() and mattrannxmnf() (real matrix
transpose for the block order) may distort the input matrices sequence x[L*S].

Parameters:
Input:
x[L*S]  Sequence of input matrices
N,M     Matrix dimensions 
L       Number of matrices
Output:
y[L*S]  Sequence of output matrices

Restrictions:
x,y     Aligned on 32-byte boundary
x,y     Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, Floating-Point, MxN->NxM, S=MxN
   Restrictions:
     N,M must be multiples of 4
*/
void mattrannxmnf ( float32_t * restrict y, 
                    float32_t * restrict x, 
                    int N, int M, int L )
{
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT(N % 4 == 0);
    NASSERT(M % 4 == 0);

    if (N<=0 || M<=0 || L<=0) return;

    if ((N&4) && (M&4))
    {
        mattrannxmnf_N4M4(y, x, N, M, L);
    }
    else if (N&4)
    {
        mattrannxmnf_N4(y, x, N, M, L);
    }
    else if (M&4)
    {
        mattrannxmnf_M4(y, x, N, M, L);
    }
    else
    {
        mattrannxmnf_N8M8(y, x, N, M, L);
    }
} /* mattrannxmnf() */

/* N is multiple of 4 and M is multiple of 4 */
void mattrannxmnf_N4M4 ( float32_t * restrict y, 
                         float32_t * restrict x, 
                         int N, int M, int L )
{
    const xb_vecNx16 * restrict px0;
    const xb_vecNx16 * restrict px1;
          xb_vecNx16 * restrict py0;
          xb_vecNx16 * restrict py1;
          xb_vecNx16 * restrict py2;
          xb_vecNx16 * restrict py3;
          xb_vecNx16 * restrict py4;
          xb_vecNx16 * restrict py5;
          xb_vecNx16 * restrict py6;
          xb_vecNx16 * restrict py7;
    int n, m, l;

    xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7;
    xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
    xb_vecNx16 t1, t3, t5, t7;
    vboolN mask4;
    valign aly1, aly3, aly5, aly7;

    py7 = (xb_vecNx16 *)(y);
    mask4 = BBE_LTRNI(8);
    aly1 = aly3 = aly5 = aly7 = BBE_ZALIGN();

    for (l=0; l<L; l++)
    {
        px0 = (const xb_vecNx16 *)(x+l*M*N);
        px1 = (const xb_vecNx16 *)(x+l*M*N+N-4);

        for (n=0; n<(N>>3); n++)
        {
            py0 = py7;
            py1 = (xb_vecNx16 *)((float32_t *)py0 + M);
            py2 = (xb_vecNx16 *)((float32_t *)py1 + M);
            py3 = (xb_vecNx16 *)((float32_t *)py2 + M);
            py4 = (xb_vecNx16 *)((float32_t *)py3 + M);
            py5 = (xb_vecNx16 *)((float32_t *)py4 + M);
            py6 = (xb_vecNx16 *)((float32_t *)py5 + M);
            py7 = (xb_vecNx16 *)((float32_t *)py6 + M);

            for (m=0; m<(M>>(LOG2_BBE_SIMD_WIDTH-1)); m++)
            {
                /* Load input matrix X */
                BBE_LVNX16_XP(X0, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X2, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X4, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X6, px0, 2*N*sizeof(float32_t));

                t1 = BBE_LVNX16_I(px1, 2*BBE_SIMD_WIDTH);
                BBE_LVNX16_XP(X1, px1, 2*N*sizeof(float32_t));
                t3 = BBE_LVNX16_I(px1, 2*BBE_SIMD_WIDTH);
                BBE_LVNX16_XP(X3, px1, 2*N*sizeof(float32_t));
                t5 = BBE_LVNX16_I(px1, 2*BBE_SIMD_WIDTH);
                BBE_LVNX16_XP(X5, px1, 2*N*sizeof(float32_t));
                t7 = BBE_LVNX16_I(px1, 2*BBE_SIMD_WIDTH);
                BBE_LVNX16_XP(X7, px1, 2*N*sizeof(float32_t));
                X1 = BBE_SELNX16I(t1, X1, BBE_SELI_ROTATE_RIGHT_8);
                X3 = BBE_SELNX16I(t3, X3, BBE_SELI_ROTATE_RIGHT_8);
                X5 = BBE_SELNX16I(t5, X5, BBE_SELI_ROTATE_RIGHT_8);
                X7 = BBE_SELNX16I(t7, X7, BBE_SELI_ROTATE_RIGHT_8);

                /* Transpose */
                BBE_DSELNX16I(X1, X0, X1, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X2, X3, X2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X5, X4, X5, X4, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X7, X6, X7, X6, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X6, X4, X6, X4, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X7, X5, X7, X5, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y4, Y0, X4, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y5, Y1, X5, X1, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y6, Y2, X6, X2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y7, Y3, X7, X3, BBE_DSELI_DEINTERLEAVE_2);

                /* Save results */
                BBE_SVNX16_IP(Y0, py0, 2*BBE_SIMD_WIDTH);
                BBE_SANX16_IP(Y1, aly1, py1);
                BBE_SVNX16_IP(Y2, py2, 2*BBE_SIMD_WIDTH);
                BBE_SANX16_IP(Y3, aly3, py3);
                BBE_SVNX16_IP(Y4, py4, 2*BBE_SIMD_WIDTH);
                BBE_SANX16_IP(Y5, aly5, py5);
                BBE_SVNX16_IP(Y6, py6, 2*BBE_SIMD_WIDTH);
                BBE_SANX16_IP(Y7, aly7, py7);
            }
            /* Transpose last 4x(BBE_SIMD_WIDTH/2) values */
            {
                /* Load input matrix X */
                BBE_LVNX16_XP(X0, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X2, px0, 2*N*sizeof(float32_t));
                t1 = BBE_LVNX16_I(px1, 2*BBE_SIMD_WIDTH);
                BBE_LVNX16_XP(X1, px1, 2*N*sizeof(float32_t));
                t3 = BBE_LVNX16_I(px1, 2*BBE_SIMD_WIDTH);
                BBE_LVNX16_XP(X3, px1, 2*N*sizeof(float32_t));
                X1 = BBE_SELNX16I(t1, X1, BBE_SELI_ROTATE_RIGHT_8);
                X3 = BBE_SELNX16I(t3, X3, BBE_SELI_ROTATE_RIGHT_8);

                /* Transpose */
                BBE_DSELNX16I(X1, X0, X1, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X2, X3, X2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y4, Y0, X0, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y5, Y1, X1, X1, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y6, Y2, X2, X2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y7, Y3, X3, X3, BBE_DSELI_DEINTERLEAVE_2);

                /* Save results ( (BBE_SIMD_WIDTH/2)x4 values ) */
                BBE_SVNX16T_I(Y0, py0, 0, mask4);
                BBE_SAVNX16_XP(Y1, aly1, py1, 4*sizeof(float32_t));
                BBE_SVNX16T_I(Y2, py2, 0, mask4);
                BBE_SAVNX16_XP(Y3, aly3, py3, 4*sizeof(float32_t));
                BBE_SVNX16T_I(Y4, py4, 0, mask4);
                BBE_SAVNX16_XP(Y5, aly5, py5, 4*sizeof(float32_t));
                BBE_SVNX16T_I(Y6, py6, 0, mask4);
                BBE_SAVNX16_XP(Y7, aly7, py7, 4*sizeof(float32_t));
            }

            BBE_SANX16POS_FP(aly1, py1);
            BBE_SANX16POS_FP(aly3, py3);
            BBE_SANX16POS_FP(aly5, py5);
            BBE_SANX16POS_FP(aly7, py7);

            px0 = (xb_vecNx16 *)((float32_t *)px0 - M*N + BBE_SIMD_WIDTH/2);
            px1 = (xb_vecNx16 *)((float32_t *)px1 - M*N + BBE_SIMD_WIDTH/2);
        }
        /* Transpose last 4 columns of matrix X */
        {
            py0 = py7;
            py1 = (xb_vecNx16 *)((float32_t *)py0 + M);
            py2 = (xb_vecNx16 *)((float32_t *)py1 + M);
            py3 = (xb_vecNx16 *)((float32_t *)py2 + M);

            for (m=0; m<(M>>(LOG2_BBE_SIMD_WIDTH-1)); m++)
            {
                /* Load input matrix X */
                BBE_LVNX16_XP(X0, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X2, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X4, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X6, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X1, px1, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X3, px1, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X5, px1, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X7, px1, 2*N*sizeof(float32_t));
                X1 = BBE_SHFLNX16I(X1, BBE_SHFLI_SWAP_8);
                X3 = BBE_SHFLNX16I(X3, BBE_SHFLI_SWAP_8);
                X5 = BBE_SHFLNX16I(X5, BBE_SHFLI_SWAP_8);
                X7 = BBE_SHFLNX16I(X7, BBE_SHFLI_SWAP_8);

                /* Transpose */
                BBE_DSELNX16I(X1, X0, X1, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X2, X3, X2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X5, X4, X5, X4, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X7, X6, X7, X6, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X6, X4, X6, X4, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X7, X5, X7, X5, BBE_DSELI_DEINTERLEAVE_2);
                Y0 = BBE_SELNX16I(X4, X0, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
                Y1 = BBE_SELNX16I(X5, X1, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
                Y2 = BBE_SELNX16I(X6, X2, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
                Y3 = BBE_SELNX16I(X7, X3, BBE_SELI_EXTRACT_2_OF_4_OFF_0);

                /* Save results */
                BBE_SVNX16_IP(Y0, py0, 2*BBE_SIMD_WIDTH);
                BBE_SANX16_IP(Y1, aly1, py1);
                BBE_SVNX16_IP(Y2, py2, 2*BBE_SIMD_WIDTH);
                BBE_SANX16_IP(Y3, aly3, py3);
            }
            /* Transpose last 4x4 values */
            {
                /* Load input matrix X */
                BBE_LVNX16T_XP(X0, px0, 2*N*sizeof(float32_t), mask4);
                BBE_LVNX16T_XP(X2, px0, 2*N*sizeof(float32_t), mask4);
                BBE_LVNX16F_XP(X1, px1, 2*N*sizeof(float32_t), mask4);
                BBE_LVNX16F_XP(X3, px1, 2*N*sizeof(float32_t), mask4);
                X1 = BBE_SHFLNX16I(X1, BBE_SHFLI_SWAP_8);
                X3 = BBE_SHFLNX16I(X3, BBE_SHFLI_SWAP_8);

                /* Transpose */
                BBE_DSELNX16I(X1, X0, X1, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X2, X3, X2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_DEINTERLEAVE_2);
                Y0 = BBE_SELNX16I(X0, X0, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
                Y1 = BBE_SELNX16I(X1, X1, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
                Y2 = BBE_SELNX16I(X2, X2, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
                Y3 = BBE_SELNX16I(X3, X3, BBE_SELI_EXTRACT_2_OF_4_OFF_0);

                /* Save results ( (BBE_SIMD_WIDTH/2)x4 values ) */
                BBE_SVNX16T_I(Y0, py0, 0, mask4);
                BBE_SAVNX16_XP(Y1, aly1, py1, 4*sizeof(float32_t));
                BBE_SVNX16T_I(Y2, py2, 0, mask4);
                BBE_SAVNX16_XP(Y3, aly3, py3, 4*sizeof(float32_t));
            }

            BBE_SANX16POS_FP(aly1, py1);
            BBE_SANX16POS_FP(aly3, py3);

            py7 = py3;
        }
    }
}

/* N is multiple of 4 */
void mattrannxmnf_N4 ( float32_t * restrict y, 
                       float32_t * restrict x, 
                       int N, int M, int L )
{
    const xb_vecNx16 * restrict px0;
    const xb_vecNx16 * restrict px1;
    const xb_vecNx16 * restrict px2;
    const xb_vecNx16 * restrict px3;
    const xb_vecNx16 * restrict px4;
    const xb_vecNx16 * restrict px5;
    const xb_vecNx16 * restrict px6;
    const xb_vecNx16 * restrict px7;
          xb_vecNx16 * restrict py0;
          xb_vecNx16 * restrict py1;
    int n, m, l;

    xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7;
    xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
    vboolN mask4;
    valign alx1, alx3, alx5, alx7;

    px7 = (const xb_vecNx16 *)(x);
    mask4 = BBE_LTRNI(8);

    for (l=0; l<L; l++)
    {
        py0 = (xb_vecNx16 *)(y+l*M*N);
        py1 = (xb_vecNx16 *)(y+l*M*N+M);

        __Pragma("loop_count min=1");
        for (m=0; m<(M>>3); m++)
        {
            px0 = px7;
            px1 = (const xb_vecNx16 *)((float32_t *)px0 + N);
            px2 = (const xb_vecNx16 *)((float32_t *)px1 + N);
            px3 = (const xb_vecNx16 *)((float32_t *)px2 + N);
            px4 = (const xb_vecNx16 *)((float32_t *)px3 + N);
            px5 = (const xb_vecNx16 *)((float32_t *)px4 + N);
            px6 = (const xb_vecNx16 *)((float32_t *)px5 + N);
            px7 = (const xb_vecNx16 *)((float32_t *)px6 + N);

            alx1 = BBE_LANX16_PP(px1);
            alx3 = BBE_LANX16_PP(px3);
            alx5 = BBE_LANX16_PP(px5);
            alx7 = BBE_LANX16_PP(px7);

            for (n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-1)); n++)
            {
                /* Load input matrix X */
                BBE_LVNX16_IP(X0, px0, 2*BBE_SIMD_WIDTH);
                BBE_LANX16_IP(X1, alx1, px1);
                BBE_LVNX16_IP(X2, px2, 2*BBE_SIMD_WIDTH);
                BBE_LANX16_IP(X3, alx3, px3);
                BBE_LVNX16_IP(X4, px4, 2*BBE_SIMD_WIDTH);
                BBE_LANX16_IP(X5, alx5, px5);
                BBE_LVNX16_IP(X6, px6, 2*BBE_SIMD_WIDTH);
                BBE_LANX16_IP(X7, alx7, px7);

                /* Transpose */
                BBE_DSELNX16I(X4, X0, X4, X0, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(X5, X1, X5, X1, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(X6, X2, X6, X2, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(X7, X3, X7, X3, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(X6, X4, X6, X4, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(X7, X5, X7, X5, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(Y1, Y0, X1, X0, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(Y3, Y2, X3, X2, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(Y5, Y4, X5, X4, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(Y7, Y6, X7, X6, BBE_DSELI_INTERLEAVE_2);

                /* Save results */
                BBE_SVNX16_XP(Y0, py0, 2*M*sizeof(float32_t));
                BBE_SVNX16_XP(Y1, py1, 2*M*sizeof(float32_t));
                BBE_SVNX16_XP(Y2, py0, 2*M*sizeof(float32_t));
                BBE_SVNX16_XP(Y3, py1, 2*M*sizeof(float32_t));
                BBE_SVNX16_XP(Y4, py0, 2*M*sizeof(float32_t));
                BBE_SVNX16_XP(Y5, py1, 2*M*sizeof(float32_t));
                BBE_SVNX16_XP(Y6, py0, 2*M*sizeof(float32_t));
                BBE_SVNX16_XP(Y7, py1, 2*M*sizeof(float32_t));
            }
            /* Transpose last (BBE_SIMD_WIDTH/2)x4 values */
            {
                /* Load input matrix X */
                X0 = BBE_LVNX16T_I(px0, 0, mask4);
                BBE_LAVNX16_XP(X1, alx1, px1, 4*sizeof(float32_t));
                X2 = BBE_LVNX16T_I(px2, 0, mask4);
                BBE_LAVNX16_XP(X3, alx3, px3, 4*sizeof(float32_t));
                X4 = BBE_LVNX16T_I(px4, 0, mask4);
                BBE_LAVNX16_XP(X5, alx5, px5, 4*sizeof(float32_t));
                X6 = BBE_LVNX16T_I(px6, 0, mask4);
                BBE_LAVNX16_XP(X7, alx7, px7, 4*sizeof(float32_t));

                /* Transpose */
                BBE_DSELNX16I(X4, X0, X4, X0, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(X5, X1, X5, X1, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(X6, X2, X6, X2, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(X7, X3, X7, X3, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(Y1, Y0, X1, X0, BBE_DSELI_INTERLEAVE_2);
                BBE_DSELNX16I(Y3, Y2, X3, X2, BBE_DSELI_INTERLEAVE_2);

                /* Save results */
                BBE_SVNX16_XP(Y0, py0, 2*M*sizeof(float32_t));
                BBE_SVNX16_XP(Y1, py1, 2*M*sizeof(float32_t));
                BBE_SVNX16_XP(Y2, py0, 2*M*sizeof(float32_t));
                BBE_SVNX16_XP(Y3, py1, 2*M*sizeof(float32_t));
            }

            py0 = (xb_vecNx16 *)((float32_t *)py0 - M*N + BBE_SIMD_WIDTH/2);
            py1 = (xb_vecNx16 *)((float32_t *)py1 - M*N + BBE_SIMD_WIDTH/2);
        }
    }
}

/* M is multiple of 4 */
void mattrannxmnf_M4 ( float32_t * restrict y, 
                       float32_t * restrict x, 
                       int N, int M, int L )
{
    const xb_vecNx16 * restrict px0;
    const xb_vecNx16 * restrict px1;
          xb_vecNx16 * restrict py0;
          xb_vecNx16 * restrict py1;
          xb_vecNx16 * restrict py2;
          xb_vecNx16 * restrict py3;
          xb_vecNx16 * restrict py4;
          xb_vecNx16 * restrict py5;
          xb_vecNx16 * restrict py6;
          xb_vecNx16 * restrict py7;
    int n, m, l;

    xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7;
    xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
    vboolN mask4;
    valign aly1, aly3, aly5, aly7;

    py7 = (xb_vecNx16 *)(y);
    mask4 = BBE_LTRNI(8);
    aly1 = aly3 = aly5 = aly7 = BBE_ZALIGN();

    for (l=0; l<L; l++)
    {
        px0 = (const xb_vecNx16 *)(x+l*M*N);
        px1 = (const xb_vecNx16 *)(x+l*M*N+N);

        __Pragma("loop_count min=1");
        for (n=0; n<(N>>3); n++)
        {
            py0 = py7;
            py1 = (xb_vecNx16 *)((float32_t *)py0 + M);
            py2 = (xb_vecNx16 *)((float32_t *)py1 + M);
            py3 = (xb_vecNx16 *)((float32_t *)py2 + M);
            py4 = (xb_vecNx16 *)((float32_t *)py3 + M);
            py5 = (xb_vecNx16 *)((float32_t *)py4 + M);
            py6 = (xb_vecNx16 *)((float32_t *)py5 + M);
            py7 = (xb_vecNx16 *)((float32_t *)py6 + M);

            for (m=0; m<(M>>(LOG2_BBE_SIMD_WIDTH-1)); m++)
            {
                /* Load input matrix X */
                BBE_LVNX16_XP(X0, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X1, px1, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X2, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X3, px1, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X4, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X5, px1, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X6, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X7, px1, 2*N*sizeof(float32_t));

                /* Transpose */
                BBE_DSELNX16I(X1, X0, X1, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X2, X3, X2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X5, X4, X5, X4, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X7, X6, X7, X6, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X6, X4, X6, X4, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X7, X5, X7, X5, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y4, Y0, X4, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y5, Y1, X5, X1, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y6, Y2, X6, X2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y7, Y3, X7, X3, BBE_DSELI_DEINTERLEAVE_2);

                /* Save results */
                BBE_SVNX16_IP(Y0, py0, 2*BBE_SIMD_WIDTH);
                BBE_SANX16_IP(Y1, aly1, py1);
                BBE_SVNX16_IP(Y2, py2, 2*BBE_SIMD_WIDTH);
                BBE_SANX16_IP(Y3, aly3, py3);
                BBE_SVNX16_IP(Y4, py4, 2*BBE_SIMD_WIDTH);
                BBE_SANX16_IP(Y5, aly5, py5);
                BBE_SVNX16_IP(Y6, py6, 2*BBE_SIMD_WIDTH);
                BBE_SANX16_IP(Y7, aly7, py7);
            }
            /* Transpose last 4x(BBE_SIMD_WIDTH/2) values */
            {
                /* Load input matrix X */
                BBE_LVNX16_XP(X0, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X1, px1, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X2, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X3, px1, 2*N*sizeof(float32_t));

                /* Transpose */
                BBE_DSELNX16I(X1, X0, X1, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X2, X3, X2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y4, Y0, X0, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y5, Y1, X1, X1, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y6, Y2, X2, X2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y7, Y3, X3, X3, BBE_DSELI_DEINTERLEAVE_2);

                /* Save results ( (BBE_SIMD_WIDTH/2)x4 values ) */
                BBE_SVNX16T_I(Y0, py0, 0, mask4);
                BBE_SAVNX16_XP(Y1, aly1, py1, 4*sizeof(float32_t));
                BBE_SVNX16T_I(Y2, py2, 0, mask4);
                BBE_SAVNX16_XP(Y3, aly3, py3, 4*sizeof(float32_t));
                BBE_SVNX16T_I(Y4, py4, 0, mask4);
                BBE_SAVNX16_XP(Y5, aly5, py5, 4*sizeof(float32_t));
                BBE_SVNX16T_I(Y6, py6, 0, mask4);
                BBE_SAVNX16_XP(Y7, aly7, py7, 4*sizeof(float32_t));
            }

            BBE_SANX16POS_FP(aly1, py1);
            BBE_SANX16POS_FP(aly3, py3);
            BBE_SANX16POS_FP(aly5, py5);
            BBE_SANX16POS_FP(aly7, py7);

            px0 = (xb_vecNx16 *)((float32_t *)px0 - M*N + BBE_SIMD_WIDTH/2);
            px1 = (xb_vecNx16 *)((float32_t *)px1 - M*N + BBE_SIMD_WIDTH/2);
        }
    }
} /* mattrannxmnf_M4() */

/* N is multiple of 8 and M is multiple of 8 */
void mattrannxmnf_N8M8 ( float32_t * restrict y, 
                         float32_t * restrict x, 
                         int N, int M, int L )
{
    const xb_vecNx16 * restrict px0;
    const xb_vecNx16 * restrict px1;
          xb_vecNx16 * restrict py0;
          xb_vecNx16 * restrict py1;
          xb_vecNx16 * restrict py2;
          xb_vecNx16 * restrict py3;
          xb_vecNx16 * restrict py4;
          xb_vecNx16 * restrict py5;
          xb_vecNx16 * restrict py6;
          xb_vecNx16 * restrict py7;
    int n, m, l;

    xb_vecNx16 X0, X1, X2, X3, X4, X5, X6, X7;
    xb_vecNx16 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;

    py7 = (xb_vecNx16 *)(y);

    for (l=0; l<L; l++)
    {
        px0 = (const xb_vecNx16 *)(x+l*M*N);
        px1 = (const xb_vecNx16 *)(x+l*M*N+N);

        __Pragma("loop_count min=1");
        for (n=0; n<(N>>3); n++)
        {
            py0 = py7;
            py1 = (xb_vecNx16 *)((float32_t *)py0 + M);
            py2 = (xb_vecNx16 *)((float32_t *)py1 + M);
            py3 = (xb_vecNx16 *)((float32_t *)py2 + M);
            py4 = (xb_vecNx16 *)((float32_t *)py3 + M);
            py5 = (xb_vecNx16 *)((float32_t *)py4 + M);
            py6 = (xb_vecNx16 *)((float32_t *)py5 + M);
            py7 = (xb_vecNx16 *)((float32_t *)py6 + M);

            __Pragma("loop_count min=1");
            for (m=0; m<(M>>(LOG2_BBE_SIMD_WIDTH-1)); m++)
            {
                /* Load input matrix X */
                BBE_LVNX16_XP(X0, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X1, px1, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X2, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X3, px1, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X4, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X5, px1, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X6, px0, 2*N*sizeof(float32_t));
                BBE_LVNX16_XP(X7, px1, 2*N*sizeof(float32_t));

                /* Transpose */
                BBE_DSELNX16I(X1, X0, X1, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X2, X3, X2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X5, X4, X5, X4, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X7, X6, X7, X6, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X2, X0, X2, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X3, X1, X3, X1, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X6, X4, X6, X4, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(X7, X5, X7, X5, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y4, Y0, X4, X0, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y5, Y1, X5, X1, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y6, Y2, X6, X2, BBE_DSELI_DEINTERLEAVE_2);
                BBE_DSELNX16I(Y7, Y3, X7, X3, BBE_DSELI_DEINTERLEAVE_2);

                /* Save results */
                BBE_SVNX16_IP(Y0, py0, 2*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(Y1, py1, 2*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(Y2, py2, 2*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(Y3, py3, 2*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(Y4, py4, 2*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(Y5, py5, 2*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(Y6, py6, 2*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(Y7, py7, 2*BBE_SIMD_WIDTH);
            }

            px0 = (xb_vecNx16 *)((float32_t *)px0 - M*N + BBE_SIMD_WIDTH/2);
            px1 = (xb_vecNx16 *)((float32_t *)px1 - M*N + BBE_SIMD_WIDTH/2);
        }
    }
}/* mattrannxmnf_N8M8() */
