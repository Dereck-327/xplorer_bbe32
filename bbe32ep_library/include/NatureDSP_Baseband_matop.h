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
 * NatureDSP_Baseband Library API
 * Matrix Operations
 */

#ifndef __NATUREDSP_BASEBAND_MATOP_H
#define __NATUREDSP_BASEBAND_MATOP_H

#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
Matrix Operations
matmul,
matvmul            Real matrix-matrix/matrix-vector multiply
cmatmul,
cmatvmul           Complex matrix-matrix/matrix-vector multiply
rcmatmul
rcmatvmul          Real matrix by complex matrix/vector multiply
mattran            Real matrix transpose
cmattran           Complex matrix conjugate transpose
cmatherm           Matrix Hermitian product
<r|c>sb            Streaming to packed conversion for real and complex matrices
<r|c>bs            Packed to streaming conversion for real and complex matrices
<r|c>interleave    M-to-1 complex/real streams interleave
<r|c>deinterleave  1-to-M complex/real streams deinterleave
===========================================================================*/

/*-------------------------------------------------------------------------
Real Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of real matrices or vectors. Both the block order and streaming 
order are allowed for input/output matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand input matrices
y[L*Sy]     Sequence of right-hand input matrices
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/
/* Block Order, Fixed-Point, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
    L must be a multiple of 4
*/
void matmul2x2n ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q );

/* Block Order, Floating-Point, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
    L must be a multiple of 2
*/
void matmul2x2nf( float32_t * restrict z,
            const float32_t * restrict x,
            const float32_t * restrict y,
            int L );

/* Block Order, Fixed-Point, 3x3*3x3->3x3, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/                                                                                                                                                      
void matmul3x3n ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q );

/* Block Order, Floating-Point, 3x3*3x3->3x3, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/                                                                                                                                                      
void matmul3x3nf ( float32_t * restrict z, 
             const float32_t * restrict x, 
             const float32_t * restrict y, 
             int L );

/* Block Order, Fixed-Point, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/
void matmul4x4n ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q );

/* Block Order, Floating-Point, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/
void matmul4x4nf ( float32_t * restrict z, 
             const float32_t * restrict x, 
             const float32_t * restrict y, 
             int L );

/* Block Order, Fixed-Point, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     None
*/
void matmul8x8n ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q );

/* Block Order, Floating-Point, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     None
*/
void matmul8x8nf ( float32_t * restrict z, 
             const float32_t * restrict x, 
             const float32_t * restrict y, 
             int L );

/* Block Order, Fixed-Point, 16x16*16x16->16x16, Sx=256, Sy=256, Sz=256
   Restrictions:
     None
*/
void matmul16x16n ( int16_t * restrict z, 
              const int16_t * restrict x, 
              const int16_t * restrict y, 
              int L, int Q );

/* Block Order, Fixed-Point, NxN*NxN->NxN, Sx=NxN, Sy=NxN, Sz=NxN
   Restrictions:  
     N must be a multiple of 4
*/
void matmulnxnn ( void * pScr,
                  int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int N, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t matmulnxnn_getScratchSize ( int N, int L );

/* Block Order, Fixed-Point, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:  
     N,M must be multiples of 4
*/
void matmulnxmn ( void * pScr,
                  int16_t * restrict z,
            const int16_t * restrict x,
            const int16_t * restrict y,
            int N, int M, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t matmulnxmn_getScratchSize ( int N, int M, int L );

/* Block Order, Floating-Point, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:  
     N,M must be multiples of 4
*/
void matmulnxmnf ( void * pScr,
                   float32_t * restrict z,
             const float32_t * restrict x,
             const float32_t * restrict y,
             int N, int M, int L );

/* Return the scratch area size, in bytes. */
size_t matmulnxmnf_getScratchSize ( int N, int M, int L );

/* Block Order, Fixed-Point, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions: 
     L must be a multiple of 8 
*/
void matvmul2x2n ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q );

/* Block Order, Floating-Point, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions: 
     L must be a multiple of 4
*/
void matvmul2x2nf ( float32_t * restrict z, 
              const float32_t * restrict x, 
              const float32_t * restrict y, 
              int L );

/* Block Order, Fixed-Point, 3x3*3x1->3x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 4
*/
void matvmul3x3n ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q );

/* Block Order, Floating-Point, 3x3*3x1->3x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 2
*/
void matvmul3x3nf ( float32_t * restrict z, 
              const float32_t * restrict x, 
              const float32_t * restrict y, 
              int L );

/* Block Order, Fixed-Point, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
    L must be a multiple of 4
*/
void matvmul4x4n ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q );

/* Block Order, Floating-Point, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
    L must be a multiple of 2
*/
void matvmul4x4nf ( float32_t * restrict z, 
              const float32_t * restrict x, 
              const float32_t * restrict y, 
              int L );

/* Block Order, Fixed-Point, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
      L must be even
*/
void matvmul8x8n ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q );

/* Block Order, Floating-Point, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
      None
*/
void matvmul8x8nf ( float32_t * restrict z, 
              const float32_t * restrict x, 
              const float32_t * restrict y, 
              int L );

/* Block Order, Fixed-Point, 16x16*16x1->16x1, Sx=256, Sy=16, Sz=16
   Restrictions:
     None
*/
void matvmul16x16n ( int16_t * restrict z, 
               const int16_t * restrict x, 
               const int16_t * restrict y, 
               int L, int Q );

/* Block Order, Fixed-Point, NxN*Nx1->Nx1, Sx=NxN, Sy=(N>8)?((N+15)&~15):N, Sz=(N>8)?((N+15)&~15):N
   Restrictions:
     L must be a multiple of 4
     N must be a multiple of 4
*/
void  matvmulnxnn ( void * pScr,
                    int16_t * restrict z, 
              const int16_t * restrict x, 
              const int16_t * restrict y, 
              int N, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t matvmulnxnn_getScratchSize( int N, int L );

/* Block Order, Fixed-Point, MxN*Nx1->Mx1, Sx=MxN, Sy=(N>8)?((N+15)&~15):N, Sz=(M>8)?((M+15)&~15):M
   Restrictions:
     L must be a multiple of 4
     N, M must be multiples of 4
*/
void matvmulnxmn ( void * pScr,
                   int16_t * restrict z,
             const int16_t * restrict x,
             const int16_t * restrict y,
             int N, int M, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t matvmulnxmn_getScratchSize ( int N, int M, int L );

/* Block Order, Floating-Point, MxN*Nx1->Mx1, Sx=MxN, Sy=(N>4)?((N+7)&~7):N, Sz=(M>4)?((M+7)&~7):M
   Restrictions:
     L must be a multiple of 2
     N, M must be multiples of 4
*/
void matvmulnxmnf ( void * pScr,
                    float32_t * restrict z,
              const float32_t * restrict x,
              const float32_t * restrict y,
              int N, int M, int L );

/* Return the scratch area size, in bytes. */
size_t matvmulnxmnf_getScratchSize ( int N, int M, int L );

/* Streaming Order, Fixed-Point, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 16
*/
void matmul2x2s ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q );

/* Streaming Order, Floating-Point, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 8
*/
void matmul2x2sf ( float32_t * restrict z, 
             const float32_t * restrict x, 
             const float32_t * restrict y, 
             int L );

/* Streaming Order, Fixed-Point, 3x3*3x3->3x3, Sx=9, Sy=9, Sz=9
   Restrictions:
     L must be a multiple of 16
*/
void matmul3x3s ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q );

/* Streaming Order, Floating-Point, 3x3*3x3->3x3, Sx=9, Sy=9, Sz=9
   Restrictions:
     L must be a multiple of 8
*/
void matmul3x3sf ( float32_t * restrict z, 
             const float32_t * restrict x, 
             const float32_t * restrict y, 
             int L );

/* Streaming Order, Fixed-Point, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     L must be a multiple of 16
*/
void matmul4x4s ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q );

/* Streaming Order, Floating-Point, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     L must be a multiple of 8
*/
void matmul4x4sf ( float32_t * restrict z, 
             const float32_t * restrict x, 
             const float32_t * restrict y, 
             int L );

/* Streaming Order, Fixed-Point, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     L must be a multiple of 16
*/
void matmul8x8s ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int L, int Q );

/* Streaming Order, Floating-Point, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     L must be a multiple of 8
*/
void matmul8x8sf ( float32_t * restrict z, 
             const float32_t * restrict x, 
             const float32_t * restrict y, 
             int L );

/* Streaming Order, Fixed-Point, 16x16*16x16->16x16, Sx=256, Sy=256, Sz=256
   Restrictions:
     L must be a multiple of 16
*/
void matmul16x16s ( int16_t * restrict z, 
              const int16_t * restrict x, 
              const int16_t * restrict y, 
              int L, int Q );

/* Streaming Order, Fixed-Point, NxN*NxN->NxN, Sx=NxN, Sy=NxN, Sz=NxN
   Restrictions:
     L must be a multiple of 16
*/
void matmulnxns ( int16_t * restrict z, 
            const int16_t * restrict x, 
            const int16_t * restrict y, 
            int N, int L, int Q );

/* Streaming Order, Fixed-Point, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
     L must be a multiple of 16
*/
void matmulnxms ( int16_t * restrict z,
            const int16_t * restrict x,
            const int16_t * restrict y,
            int N, int M, int L, int Q );

/* Streaming Order, Floating-Point, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
     L must be a multiple of 8
*/
void matmulnxmsf ( float32_t * restrict z,
             const float32_t * restrict x,
             const float32_t * restrict y,
             int N, int M, int L );

/* Streaming Order, Fixed-Point, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions:
     L must be a multiple of 16
*/
void matvmul2x2s ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q );

/* Streaming Order, Floating-Point, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions:
     L must be a multiple of 8
*/
void matvmul2x2sf ( float32_t * restrict z, 
              const float32_t * restrict x, 
              const float32_t * restrict y, 
              int L );

/* Streaming Order, Fixed-Point, 3x3*3x1->3x1, Sx=9, Sy=3, Sz=3
   Restrictions:
     L must be a multiple of 16
*/
void matvmul3x3s ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q );

/* Streaming Order, Floating-Point, 3x3*3x1->3x1, Sx=9, Sy=3, Sz=3
   Restrictions:
     L must be a multiple of 8
*/
void matvmul3x3sf ( float32_t * restrict z, 
              const float32_t * restrict x, 
              const float32_t * restrict y, 
              int L );

/* Streaming Order, Fixed-Point, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 16
*/
void matvmul4x4s ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q );

/* Streaming Order, Floating-Point, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 8
*/
void matvmul4x4sf ( float32_t * restrict z, 
              const float32_t * restrict x, 
              const float32_t * restrict y, 
              int L );

/* Streaming Order, Fixed-Point, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     L must be a multiple of 16
*/
void matvmul8x8s ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int L, int Q );

/* Streaming Order, Floating-Point, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     L must be a multiple of 8
*/
void matvmul8x8sf ( float32_t * restrict z, 
              const float32_t * restrict x, 
              const float32_t * restrict y, 
              int L );

/* Streaming Order, Fixed-Point, 16x16*16x1->16x1, Sx=256, Sy=16, Sz=16
   Restrictions:
     L must be a multiple of 16
*/
void matvmul16x16s ( int16_t * restrict z, 
               const int16_t * restrict x, 
               const int16_t * restrict y, 
               int L, int Q );

/* Streaming Order, Fixed-Point, NxN*Nx1->Nx1, Sx=NxN, Sy=N, Sz=N
   Restrictions:
     L must be a multiple of 16
*/
void matvmulnxns ( int16_t * restrict z, 
             const int16_t * restrict x, 
             const int16_t * restrict y, 
             int N, int L, int Q );

/* Streaming Order, Fixed-Point, MxN*Nx1->Mx1, Sx=MxN, Sy=N, Sz=M
   Restrictions:
     L must be a multiple of 16
*/
void  matvmulnxms ( int16_t * restrict z,
              const int16_t * restrict x,
              const int16_t * restrict y,
              int N, int M, int L, int Q );

/* Streaming Order, Floating-Point, MxN*Nx1->Mx1, Sx=MxN, Sy=N, Sz=M
   Restrictions:
     L must be a multiple of 8
*/
void  matvmulnxmsf ( float32_t * restrict z,
               const float32_t * restrict x,
               const float32_t * restrict y,
               int N, int M, int L );

/*-------------------------------------------------------------------------
Complex Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of complex matrices or vectors. Both the block order and 
streaming order are allowed for input/output matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand complex matrices
y[L*Sy]     Sequence of right-hand complex matrices
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices 
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of complex result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, Fixed-Point, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
     L must be even
*/
void cmatmul2x2n ( complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int L, int Q );

/* Block Order, Floating-Point, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
     None
*/
void cmatmul2x2nf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int L );

/* Block Order, Fixed-Point, 3x3*3x3->3x3, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/
void cmatmul3x3n ( complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y,
             int L, int Q );

/* Block Order, Floating-Point, 3x3*3x3->3x3, Sx=12, Sy=12, Sz=12
   Restrictions:
     None
*/
void cmatmul3x3nf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y,
              int L );

/* Block Order, Fixed-Point, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/
void cmatmul4x4n ( complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int L, int Q );

/* Block Order, Floating-Point, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/
void cmatmul4x4nf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int L );

/* Block Order, Fixed-Point, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     None
*/
void cmatmul8x8n ( complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int L, int Q );

/* Block Order, Floating-Point, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     None
*/
void cmatmul8x8nf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int L );

/* Block Order, Fixed-Point, 16x16*16x16->16x16, Sx=256, Sy=256, Sz=256
   Restrictions:
     None
*/
void cmatmul16x16n ( complex_fract16 * restrict z, 
               const complex_fract16 * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q );

/* Block Order, Fixed-Point, NxN*NxN->NxN, Sx=NxN, Sy=NxN, Sz=NxN
   Restrictions:
     N must be a multiple of 4
*/
void cmatmulnxnn ( void * pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t cmatmulnxnn_getScratchSize ( int N, int L );

/* Block Order, Fixed-Point, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
     N,M must be multiples of 4
*/
void cmatmulnxmn ( void * pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t cmatmulnxmn_getScratchSize ( int N, int M, int L );

/* Block Order, Floating-Point, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
     N,M must be multiples of 4
*/
void cmatmulnxmnf ( void * pScr,
                    complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int N, int M, int L );

/* Return the scratch area size, in bytes. */
size_t cmatmulnxmnf_getScratchSize ( int N, int M, int L );

/* Block Order, Fixed-Point, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions:
     L must be a multiple of 4
*/
void cmatvmul2x2n ( complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Block Order, Floating-Point, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions:
     L must be even
*/
void cmatvmul2x2nf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int L );

/* Block Order, Fixed-Point, 3x3*3x1->3x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be even
*/
void cmatvmul3x3n ( complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Block Order, Floating-Point, 3x3*3x1->3x1, Sx=12, Sy=4, Sz=4
   Restrictions:
     None
*/
void cmatvmul3x3nf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int L );

/* Block Order, Fixed-Point, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be even
*/
void cmatvmul4x4n ( complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Block Order, Floating-Point, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     None
*/
void cmatvmul4x4nf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int L );

/* Block Order, Fixed-Point, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     None
*/
void cmatvmul8x8n ( complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Block Order, Floating-Point, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     None
*/
void cmatvmul8x8nf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int L );

/* Block Order, Fixed-Point, 16x16*16x1->16x1, Sx=256, Sy=16, Sz=16
   Restrictions:
     None
*/
void cmatvmul16x16n ( complex_fract16 * restrict z, 
                const complex_fract16 * restrict x, 
                const complex_fract16 * restrict y, 
                int L, int Q );

/* Block Order, Fixed-Point, NxN*Nx1->Nx1, Sx=NxN, Sy=(N>4)?((N+7)&~7):N, Sz=(N>4)?((N+7)&~7):N
   Restrictions:
     L must be even
     N must be a multiple of 4
*/
void cmatvmulnxnn ( void * pScr,
                    complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int N, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t cmatvmulnxnn_getScratchSize ( int N, int L );

/* Block Order, Fixed-Point, MxN*Nx1->Mx1, Sx=MxN, Sy=(N>4)?((N+7)&~7):N, Sz=(M>4)?((M+7)&~7):M
   Restrictions:
     L must be even
     N,M must be multiples of 4
*/
void cmatvmulnxmn ( void * pScr,
                    complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int N, int M, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t cmatvmulnxmn_getScratchSize ( int N, int M, int L );

/* Block Order, Floating-Point, MxN*Nx1->Mx1, Sx=MxN, Sy=N, Sz=M
   Restrictions:
     N,M must be multiples of 4
*/
void cmatvmulnxmnf ( void * pScr,
                     complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int N, int M, int L );

/* Return the scratch area size, in bytes. */
size_t cmatvmulnxmnf_getScratchSize ( int N, int M, int L );

/* Streaming Order, Fixed-Point, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 8
*/
void cmatmul2x2s ( complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int L, int Q );

/* Streaming Order, Floating-Point, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 4
*/
void cmatmul2x2sf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int L );

/* Streaming Order, Fixed-Point, 3x3*3x3->3x3, Sx=9, Sy=9, Sz=9
   Restrictions:
     L must be a multiple of 8
*/
void cmatmul3x3s ( complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int L, int Q );

/* Streaming Order, Floating-Point, 3x3*3x3->3x3, Sx=9, Sy=9, Sz=9
   Restrictions:
     L must be a multiple of 4
*/
void cmatmul3x3sf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int L );

/* Streaming Order, Fixed-Point, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     L must be a multiple of 8
*/
void cmatmul4x4s ( complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int L, int Q );

/* Streaming Order, Floating-Point, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     L must be a multiple of 4
*/
void cmatmul4x4sf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int L );

/* Streaming Order, Fixed-Point, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     L must be a multiple of 8
*/
void cmatmul8x8s ( complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int L, int Q );

/* Streaming Order, Floating-Point, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     L must be a multiple of 4
*/
void cmatmul8x8sf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int L );

/* Streaming Order, Fixed-Point, 16x16*16x16->16x16, Sx=256, Sy=256, Sz=256
   Restrictions:
     L must be a multiple of 8
*/
void cmatmul16x16s ( complex_fract16 * restrict z, 
               const complex_fract16 * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q );

/* Streaming Order, Fixed-Point, NxN*NxN->NxN, Sx=NxN, Sy=NxN, Sz=NxN
   Restrictions:
     L must be a multiple of 8
*/
void cmatmulnxns ( complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int L, int Q );

/* Streaming Order, Fixed-Point, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
     L must be a multiple of 8
*/
void cmatmulnxms ( complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q );

/* Streaming Order, Floating-Point, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
     L must be a multiple of 4
*/
void cmatmulnxmsf ( complex_float * restrict z, 
              const complex_float * restrict x, 
              const complex_float * restrict y, 
              int N, int M, int L );

/* Streaming Order, Fixed-Point, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions:
     L must be a multiple of 8
*/
void cmatvmul2x2s ( complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Streaming Order, Floating-Point, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions:
     L must be a multiple of 4
*/
void cmatvmul2x2sf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int L );

/* Streaming Order, Fixed-Point, 3x3*3x1->3x1, Sx=9, Sy=3, Sz=3
   Restrictions:
     L must be a multiple of 8
*/
void cmatvmul3x3s ( complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Streaming Order, Floating-Point, 3x3*3x1->3x1, Sx=9, Sy=3, Sz=3
   Restrictions:
     L must be a multiple of 4
*/
void cmatvmul3x3sf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int L );

/* Streaming Order, Fixed-Point, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 8
*/
void cmatvmul4x4s ( complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Streaming Order, Floating-Point, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 4
*/
void cmatvmul4x4sf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int L );

/* Streaming Order, Fixed-Point, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     L must be a multiple of 8
*/
void cmatvmul8x8s ( complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Streaming Order, Floating-Point, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     L must be a multiple of 4
*/
void cmatvmul8x8sf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int L );

/* Streaming Order, Fixed-Point, 16x16*16x1->16x1, Sx=256, Sy=16, Sz=16
   Restrictions:
     L must be a multiple of 8
*/
void cmatvmul16x16s ( complex_fract16 * restrict z, 
                const complex_fract16 * restrict x, 
                const complex_fract16 * restrict y, 
                int L, int Q );

/* Streaming Order, Fixed-Point, NxN*Nx1->Nx1, Sx=NxN, Sy=N, Sz=N
   Restrictions:
     L must be a multiple of 8
*/
void cmatvmulnxns ( complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int N, int L, int Q );

/* Streaming Order, Fixed-Point, MxN*Nx1->Mx1, Sx=MxN, Sy=N, Sz=M
   Restrictions:
     L must be a multiple of 8
*/
void cmatvmulnxms ( complex_fract16 * restrict z, 
              const complex_fract16 * restrict x, 
              const complex_fract16 * restrict y, 
              int N, int M, int L, int Q );

/* Streaming Order, Floating-Point, MxN*Nx1->Mx1, Sx=MxN, Sy=N, Sz=M
   Restrictions:
     L must be a multiple of 4
*/
void cmatvmulnxmsf ( complex_float * restrict z, 
               const complex_float * restrict x, 
               const complex_float * restrict y, 
               int N, int M, int L );

/*-------------------------------------------------------------------------
Real Matrix by Complex Matrix/Vector Multiply 

Description: These functions perform pairwise multiplication of left-hand
real matrices by right-hand complex matrices or vectors. Both the block order
and streaming order are allowed for input/output matrix sequences.

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand real matrices
y[L*Sy]     Sequence of right-hand complex matrices or vectors
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices 
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of complex result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 4
*/
void rcmatmul2x2n ( complex_fract16 * restrict z, 
              const int16_t         * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Block Order, 3x3*3x3->3x3, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/
void rcmatmul3x3n ( complex_fract16 * restrict z, 
              const int16_t         * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Block Order, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/
void rcmatmul4x4n ( complex_fract16 * restrict z, 
              const int16_t         * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Block Order, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     None
*/
void rcmatmul8x8n ( complex_fract16 * restrict z, 
              const int16_t         * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Block Order, 16x16*16x16->16x16, Sx=256, Sy=256, Sz=256
   Restrictions:
     None
*/
void rcmatmul16x16n ( complex_fract16 * restrict z, 
                const int16_t         * restrict x, 
                const complex_fract16 * restrict y, 
                int L, int Q );

/* Block Order, NxN*NxN->NxN, Sx=NxN, Sy=NxN, Sz=NxN
   Restrictions:
     N must be a multiple of 4
*/
void rcmatmulnxnn ( void * pScr,
                    complex_fract16 * restrict z, 
              const int16_t         * restrict x, 
              const complex_fract16 * restrict y, 
              int N, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t rcmatmulnxnn_getScratchSize ( int N, int L );

/* Block Order, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
      N,M must be multiples of 4
*/
void rcmatmulnxmn ( void * pScr,
                    complex_fract16 * restrict z, 
              const int16_t         * restrict x, 
              const complex_fract16 * restrict y, 
              int N, int M, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t rcmatmulnxmn_getScratchSize ( int N, int M, int L );

/* Block Order, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions:
     L must be a multiple of 4
*/
void rcmatvmul2x2n ( complex_fract16 * restrict z, 
               const int16_t         * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q );

/* Block Order, 3x3*3x1->3x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be even
*/
void rcmatvmul3x3n ( complex_fract16 * restrict z, 
               const int16_t         * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q );

/* Block Order, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be even
*/
void rcmatvmul4x4n ( complex_fract16 * restrict z, 
               const int16_t         * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q );

/* Block Order, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     None
*/
void rcmatvmul8x8n ( complex_fract16 * restrict z, 
               const int16_t         * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q );

/* Block Order, 16x16*16x1->16x1, Sx=256, Sy=16, Sz=16
   Restrictions:
     None
*/
void rcmatvmul16x16n ( complex_fract16 * restrict z, 
                 const int16_t         * restrict x, 
                 const complex_fract16 * restrict y, 
                 int L, int Q );

/* Block Order, NxN*Nx1->Nx1,  Sx=NxN, Sy=N, Sz=N
   Restrictions:
     N,L must be multiples of 4
*/
void rcmatvmulnxnn ( void * pScr,
                     complex_fract16 * restrict z, 
               const int16_t         * restrict x, 
               const complex_fract16 * restrict y, 
               int N, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t rcmatvmulnxnn_getScratchSize ( int N, int L );

/* Block Order, MxN*Nx1->Mx1, Sx=MxN, Sy=N, Sz=M
   Restrictions:
     N,M,L must be multiples of 4
*/
void rcmatvmulnxmn ( void * pScr,
                     complex_fract16 * restrict z, 
               const int16_t         * restrict x, 
               const complex_fract16 * restrict y, 
               int N, int M, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t rcmatvmulnxmn_getScratchSize ( int N, int M, int L );

/* Streaming Order, 2x2*2x2->2x2, Sx=4, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 16
*/
void rcmatmul2x2s ( complex_fract16 * restrict z, 
              const int16_t         * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Streaming Order, 3x3*3x3->3x3, Sx=9, Sy=9, Sz=9
   Restrictions:
     L must be a multiple of 16
*/
void rcmatmul3x3s ( complex_fract16 * restrict z, 
              const int16_t         * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Streaming Order, 4x4*4x4->4x4, Sx=16, Sy=16, Sz=16
   Restrictions:
     L must be a multiple of 16
*/
void rcmatmul4x4s ( complex_fract16 * restrict z, 
              const int16_t         * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Streaming Order, 8x8*8x8->8x8, Sx=64, Sy=64, Sz=64
   Restrictions:
     L must be a multiple of 16
*/
void rcmatmul8x8s ( complex_fract16 * restrict z, 
              const int16_t         * restrict x, 
              const complex_fract16 * restrict y, 
              int L, int Q );

/* Streaming Order, 16x16*16x16->16x16, Sx=256, Sy=256, Sz=256
   Restrictions:
     L must be a multiple of 16
*/
void rcmatmul16x16s ( complex_fract16 * restrict z, 
                const int16_t         * restrict x, 
                const complex_fract16 * restrict y, 
                int L, int Q );

/* Streaming Order, NxN*NxN->NxN, Sx=NxN, Sy=NxN, Sz=NxN
   Restrictions:
     L must be a multiple of 16
*/
void rcmatmulnxns ( complex_fract16 * restrict z, 
              const int16_t         * restrict x, 
              const complex_fract16 * restrict y, 
              int N, int L, int Q );

/* Streaming Order, MxN*NxM ->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:
     L must be a multiple of 16
*/
void rcmatmulnxms ( complex_fract16 * restrict z, 
              const int16_t         * restrict x, 
              const complex_fract16 * restrict y, 
              int N, int M, int L, int Q );

/* Streaming Order, 2x2*2x1->2x1, Sx=4, Sy=2, Sz=2
   Restrictions:
     L must be a multiple of 16
*/
void rcmatvmul2x2s ( complex_fract16 * restrict z, 
               const int16_t         * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q );

/* Streaming Order, 3x3*3x1->3x1, Sx=9, Sy=3, Sz=3
   Restrictions:
     L must be a multiple of 16
*/
void rcmatvmul3x3s ( complex_fract16 * restrict z, 
               const int16_t         * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q );

/* Streaming Order, 4x4*4x1->4x1, Sx=16, Sy=4, Sz=4
   Restrictions:
     L must be a multiple of 16
*/
void rcmatvmul4x4s ( complex_fract16 * restrict z, 
               const int16_t         * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q );

/* Streaming Order, 8x8*8x1->8x1, Sx=64, Sy=8, Sz=8
   Restrictions:
     L must be a multiple of 16
*/
void rcmatvmul8x8s ( complex_fract16 * restrict z, 
               const int16_t         * restrict x, 
               const complex_fract16 * restrict y, 
               int L, int Q );

/* Streaming Order, 16x16*16x1->16x1, Sx=256, Sy=16, Sz=16
   Restrictions:
     L must be a multiple of 16
*/
void rcmatvmul16x16s ( complex_fract16 * restrict z, 
                 const int16_t         * restrict x, 
                 const complex_fract16 * restrict y, 
                 int L, int Q );

/* Streaming Order, NxN*Nx1->Nx1, Sx=NxN, Sy=N, Sz=N
   Restrictions:
     L must be a multiple of 16
*/
void rcmatvmulnxns ( complex_fract16 * restrict z, 
               const int16_t         * restrict x, 
               const complex_fract16 * restrict y, 
               int N, int L, int Q );

/* Streaming Order, MxN*Nx1->Mx1, Sx=MxN, Sy=N, Sz=M
   Restrictions:
     L must be a multiple of 16
*/
void rcmatvmulnxms ( complex_fract16 * restrict z, 
               const int16_t         * restrict x, 
               const complex_fract16 * restrict y, 
               int N, int M, int L, int Q );

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
/* Block Order, Fixed-Point, 2x2->2x2, S=4
   Restrictions:
     L must be a multiple of 4
*/
void mattran2x2n ( int16_t * restrict y, 
             const int16_t * restrict x, 
             int L );

/* Block Order, Floating-Point, 2x2->2x2, S=4
   Restrictions:
     L must be even
*/
void mattran2x2nf ( float32_t * restrict y, 
              const float32_t * restrict x, 
              int L );

/* Block Order, Fixed-Point, 4x4->4x4, S=16
   Restrictions:
     None
 */
void mattran4x4n ( int16_t * restrict y, 
             const int16_t * restrict x, 
             int L );

/* Block Order, Floating-Point, 4x4->4x4, S=16
   Restrictions:
     None
 */
void mattran4x4nf ( float32_t * restrict y, 
              const float32_t * restrict x, 
              int L );

/* Block Order, Fixed-Point, 8x8->8x8, S=64
   Restrictions:
     None
*/
void mattran8x8n ( int16_t * restrict y, 
             const int16_t * restrict x, 
             int L );

/* Block Order, Floating-Point, 8x8->8x8, S=64
   Restrictions:
     None
*/
void mattran8x8nf ( float32_t * restrict y, 
              const float32_t * restrict x, 
              int L );

/* Block Order, Fixed-Point, 16x16->16x16, S=256
   Restrictions:
     None
*/
void mattran16x16n ( int16_t * restrict y, 
               const int16_t * restrict x, 
               int L );

/* Block Order, Fixed-Point, NxN->NxN, S=NxN
   Restrictions:
     N must be a multiple of 4
*/
void mattrannxnn ( int16_t * restrict y, 
                   int16_t * restrict x, 
                   int N, int L );

/* Block Order, Fixed-Point, MxN->NxM, S=MxN
   Restrictions:
     N,M must be multiples of 4
*/
void mattrannxmn ( int16_t * restrict y, 
                   int16_t * restrict x, 
                   int N, int M, int L );

/* Block Order, Floating-Point, MxN->NxM, S=MxN
   Restrictions:
     N,M must be multiples of 4
*/
void mattrannxmnf ( float32_t * restrict y, 
                    float32_t * restrict x, 
                    int N, int M, int L );

/* Streaming Order, Fixed-Point, 2x2->2x2, S=4
   Restrictions:
     L must be a multiple of 16
*/
void mattran2x2s ( int16_t * restrict y, 
             const int16_t * restrict x, 
             int L );

/* Streaming Order, Floating-Point, 2x2->2x2, S=4
   Restrictions:
     L must be a multiple of 8
*/
void mattran2x2sf ( float32_t * restrict y, 
              const float32_t * restrict x, 
              int L );

/* Streaming Order, Fixed-Point, 4x4->4x4, S=16
   Restrictions:
     L must be a multiple of 16
*/
void mattran4x4s ( int16_t * restrict y, 
             const int16_t * restrict x, 
             int L );

/* Streaming Order, Floating-Point, 4x4->4x4, S=16
   Restrictions:
     L must be a multiple of 8
*/
void mattran4x4sf ( float32_t * restrict y, 
              const float32_t * restrict x, 
              int L );

/* Streaming Order, Fixed-Point, 8x8->8x8, S=64
   Restrictions:
     L must be a multiple of 16
*/
void mattran8x8s ( int16_t * restrict y, 
             const int16_t * restrict x, 
             int L );

/* Streaming Order, Floating-Point, 8x8->8x8, S=64
   Restrictions:
     L must be a multiple of 8
*/
void mattran8x8sf ( float32_t * restrict y, 
              const float32_t * restrict x, 
              int L );

/* Streaming Order, Fixed-Point, 16x16->16x16, S=256
   Restrictions:
     L must be a multiple of 16
*/
void mattran16x16s ( int16_t * restrict y, 
               const int16_t * restrict x, 
               int L );

/* Streaming Order, Fixed-Point, NxN->NxN, S=NxN
   Restrictions:
     L must be a multiple of 16
*/
void mattrannxns ( int16_t * restrict y, 
             const int16_t * restrict x, 
             int N, int L );

/* Streaming Order, Fixed-Point, MxN->NxM, S=MxN
   Restrictions:
     L must be a multiple of 16
*/
void mattrannxms ( int16_t * restrict y, 
             const int16_t * restrict x, 
             int N, int M, int L );

/* Streaming Order, Floating-Point, MxN->NxM, S=MxN
   Restrictions:
     L must be a multiple of 8
*/
void mattrannxmsf ( float32_t * restrict y, 
              const float32_t * restrict x, 
              int N, int M, int L );

/*-------------------------------------------------------------------------
Complex Matrix Conjugate Transpose

Description: These functions perform transposition and then take the complex
conjugate for each matrix from input sequence. Results are stored to output
sequence. Both the block order and streaming order are allowed for input/output
matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Note:
1. Complex conjugation of fixed-point data may involve 16-bit saturation of
   imaginary components
2. The functions cmattrannxnn(), cmattrannxmn() and cmattrannxmnf() (conjugate
   transpose for the block order) may distort the input matrices sequence x[L*S].

Parameters:
Input:
x[L*S]  Sequence of input matrices.
N,M     Matrix dimensions 
L       Number of matrices
Output:
y[L*S]  Sequence of output matrices

Restrictions:
x,y     Aligned on 32-byte boundary
x,y     Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, Fixed-Point, 2x2->2x2, S=4
   Restrictions:
    L must be even
*/
void cmattran2x2n ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L );

/* Block Order, Floating-Point, 2x2->2x2, S=4
   Restrictions:
    None
*/
void cmattran2x2nf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L );

/* Block Order, Fixed-Point, 4x4->4x4, S=16
   Restrictions:
     None
*/
void cmattran4x4n ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L );

/* Block Order, Floating-Point, 4x4->4x4, S=16
   Restrictions:
     None
*/
void cmattran4x4nf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L );

/* Block Order, Fixed-Point, 8x8->8x8, S=64
   Restrictions:
     None
*/
void cmattran8x8n ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L );

/* Block Order, Floating-Point, 8x8->8x8, S=64
   Restrictions:
     None
*/
void cmattran8x8nf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L );

/* Block Order, Fixed-Point, 16x16->16x16, S=256
   Restrictions:
     None
*/
void cmattran16x16n ( complex_fract16 * restrict y, 
                const complex_fract16 * restrict x, 
                int L );

/* Block Order, Fixed-Point, NxN->NxN, S=NxN
   Restrictions:
     N must be a multiple of 4
*/
void cmattrannxnn ( complex_fract16 * restrict y, 
                    complex_fract16 * restrict x, 
                    int N, int L );

/* Block Order, Fixed-Point, MxN->NxM, S=MxN
   Restrictions:
     N,M must be multiples of 4
*/
void cmattrannxmn ( complex_fract16 * restrict y, 
                    complex_fract16 * restrict x, 
                    int N, int M, int L );

/* Block Order, Floating-Point, MxN->NxM, S=MxN
   Restrictions:
     N,M must be multiples of 4
*/
void cmattrannxmnf ( complex_float * restrict y, 
                     complex_float * restrict x, 
                     int N, int M, int L );

/* Streaming Order, Fixed-Point, 2x2->2x2, S=4
   Restrictions:
     L must be a multiple of 8
*/
void cmattran2x2s ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L );

/* Streaming Order, Floating-Point, 2x2->2x2, S=4
   Restrictions:
     L must be a multiple of 4
*/
void cmattran2x2sf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L );

/* Streaming Order, Fixed-Point, 4x4->4x4, S=16
   Restrictions:
     L must be a multiple of 8
*/
void cmattran4x4s ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L );

/* Streaming Order, Floating-Point, 4x4->4x4, S=16
   Restrictions:
     L must be a multiple of 4
*/
void cmattran4x4sf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L );

/* Streaming Order, Fixed-Point, 8x8->8x8, S=64
   Restrictions:
     L must be a multiple of 8
*/
void cmattran8x8s ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L );

/* Streaming Order, Floating-Point, 8x8->8x8, S=64
   Restrictions:
     L must be a multiple of 4
*/
void cmattran8x8sf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L );

/* Streaming Order, Fixed-Point, 16x16->16x16, S=256
   Restrictions:
     L must be a multiple of 8
*/
void cmattran16x16s ( complex_fract16 * restrict y, 
                const complex_fract16 * restrict x, 
                int L );

/* Streaming Order, Fixed-Point, NxN->NxN, S=NxN
   Restrictions:
     L must be a multiple of 8
*/
void cmattrannxns ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int N, int L );

/* Streaming Order, Fixed-Point, MxN->NxM, S=MxN
   Restrictions:
     L must be a multiple of 8
*/
void cmattrannxms ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int N, int M, int L );

/* Streaming Order, Floating-Point, MxN->NxM, S=MxN
   Restrictions:
     L must be a multiple of 4
*/
void cmattrannxmsf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int N, int M, int L );

/*-------------------------------------------------------------------------
Matrix Hermitian Product

Description: These functions left multiply each complex input matrix by its
conjugate transpose. The result is a Hermitian (or self-adjoint) matrix. Both
the block order and streaming order are allowed for input/output matrix
sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr      Scratch memory area. To determine the scratch area size required by
          a function <fun>, use the respective helper function 
          <fun>_getScratchSize()
Input:
x[L*Sx]   Complex input matrices
M         Matrix dimension 
N         Matrix dimension (columnar for MxN)
L         Number of matrices 
Q         Position of fractional point in matrix representation, 0..16
Output:
y[L*Sy]   Complex output matrices

Restrictions:
pScr,x,y  Aligned on 32-byte boundary
pScr,x,y  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/
/* Block Order, Fixed-Point, 2x2*2x2->2x2, Sx=4, Sy=4
   Restrictions:
     L must be even
*/
void cmatherm2x2n ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L, int Q );

/* Block Order, Floating-Point, 2x2*2x2->2x2, Sx=4, Sy=4
   Restrictions:
     None
*/
void cmatherm2x2nf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L );

/* Block Order, Fixed-Point, 4x4*4x4->4x4, Sx=16, Sy=16
   Restrictions:
     None
*/
void cmatherm4x4n ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L, int Q );

/* Block Order, Floating-Point, 4x4*4x4->4x4, Sx=16, Sy=16
   Restrictions:
     None
*/
void cmatherm4x4nf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L );

/* Block Order, Fixed-Point, 8x8*8x8->8x8, Sx=64, Sy=64
   Restrictions:
     None
*/
void cmatherm8x8n ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L, int Q );

/* Block Order, Floating-Point, 8x8*8x8->8x8, Sx=64, Sy=64
   Restrictions:
     None
*/
void cmatherm8x8nf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L );

/* Block Order, Fixed-Point, 16x16*16x16->16x16, Sx=256, Sy=256
   Restrictions:
     None
*/
void cmatherm16x16n ( complex_fract16 * restrict y, 
                const complex_fract16 * restrict x, 
                int L, int Q );

/* Block Order, Fixed-Point, NxN*NxN->NxN, Sx=NxN, Sy=NxN
   Restrictions:
     N must be a multiple of 4
*/
void cmathermnxnn ( void * pScr,
                    complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int N, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t cmathermnxnn_getScratchSize ( int N, int L );

/* Block Order, Fixed-Point, NxM*MxN->NxN, Sx=MxN, Sy=NxN
   Restrictions:
     M,N must be multiples of 4
*/
void cmathermnxmn ( void * pScr,
                    complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int N, int M, int L, int Q );

/* Return the scratch area size, in bytes. */
size_t cmathermnxmn_getScratchSize ( int N, int M, int L );

/* Block Order, Floating-Point, NxM*MxN->NxN, Sx=MxN, Sy=NxN
   Restrictions:
     M,N must be multiples of 4
*/
void cmathermnxmnf ( void * pScr,
                     complex_float * restrict y, 
               const complex_float * restrict x, 
               int N, int M, int L );

/* Return the scratch area size, in bytes. */
size_t cmathermnxmnf_getScratchSize ( int N, int M, int L );

/* Streaming Order, Fixed-Point, 2x2*2x2->2x2, Sx=4, Sy=4
   Restrictions:
     L must be a multiple of 8
*/
void cmatherm2x2s ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L, int Q );

/* Streaming Order, Floating-Point, 2x2*2x2->2x2, Sx=4, Sy=4
   Restrictions:
     L must be a multiple of 4
*/
void cmatherm2x2sf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L );

/* Streaming Order, Fixed-Point, 4x4*4x4->4x4, Sx=16, Sy=16
   Restrictions:
     L must be a multiple of 8
*/
void cmatherm4x4s ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L, int Q );

/* Streaming Order, Floating-Point, 4x4*4x4->4x4, Sx=16, Sy=16
   Restrictions:
     L must be a multiple of 4
*/
void cmatherm4x4sf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L );

/* Streaming Order, Fixed-Point, 8x8*8x8->8x8, Sx=64, Sy=64
   Restrictions:
     L must be a multiple of 8
*/
void cmatherm8x8s ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int L, int Q );

/* Streaming Order, Floating-Point, 8x8*8x8->8x8, Sx=64, Sy=64
   Restrictions:
     L must be a multiple of 4
*/
void cmatherm8x8sf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int L );

/* Streaming Order, Fixed-Point, 16x16*16x16->16x16, Sx=256, Sy=256
   Restrictions:
     L must be a multiple of 8
*/
void cmatherm16x16s ( complex_fract16 * restrict y, 
                const complex_fract16 * restrict x, 
                int L, int Q );

/* Streaming Order, Fixed-Point, NxN*NxN->NxN, Sx=NxN, Sy=NxN
   Restrictions:
     L must be a multiple of 8
*/
void cmathermnxns ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int N, int L, int Q );

/* Streaming Order, Fixed-Point, NxM*MxN->NxN, Sx=MxN, Sy=NxN
   Restrictions:
     L must be a multiple of 8
*/
void cmathermnxms ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x, 
              int N, int M, int L, int Q );

/* Streaming Order, Floating-Point, NxM*MxN->NxN, Sx=MxN, Sy=NxN
   Restrictions:
     L must be a multiple of 4
*/
void cmathermnxmsf ( complex_float * restrict y, 
               const complex_float * restrict x, 
               int N, int M, int L );

/*-------------------------------------------------------------------------
Streaming to Packed Conversion for Real and Complex Matrices

Description: convert a sequence of L MxN matrices from streaming to packed
(block) order. Use rsb*() functions for real data, and csb*() functions - for
complex data.

Representation:
<r|c>sb<size>    16-bit fixed-point data
<r|c>sb<size>f   IEEE-754 Std single precision floating-point data

Storage size SY denotes the number of data elements required to store an
MxN matrix Y in block order. If matrix size M*N is less than the SIMD vector
size for appropriate data type, then the storage size equals M*N rounded up
to the next power of two, otherwise SY equals M*N rounded up to the next
multiple of the SIMD vector size.

SIMD vector size:
  - for real fixed-point data 2*BBE_SIMD_WIDTH/sizeof(int16_t) == 16
  - for complex fixed-point data 2*BBE_SIMD_WIDTH/sizeof(complex_fract16) == 8
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4

Parameters:
Input:
x[M*N][L]  Input sequence of MxN matrices, streaming order
M          Number of rows in a matrix
N          Number of columns in a matrix
L          Number of matrices
Output:
y[L][Sy]   Output sequence of matrices, block order.

Restrictions:
x,y        Aligned on 32-byte boundary
x,y        Must not overlap
L          Must be a multiple of:
             16 for real fixed-point data 
              8 for complex fixed-point data and real floating-point data
              4 for complex floating-point data
-------------------------------------------------------------------------*/
/* Real fixed-point data, 16 bits per data element */
void rsb2x2 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=2, N=2, SY= 4 */
void rsb3x3 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=3, N=3, SY=16 */
void rsb4x4 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=4, N=4, SY=16 */
void rsb8x8 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=8, N=8, SY=64 */
void rsb2x1 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=2, N=1, SY= 2 */
void rsb4x1 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=4, N=1, SY= 4 */
void rsb8x1 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=8, N=1, SY= 8 */
void rsbmxn ( int16_t * restrict y, const int16_t * restrict x, int M, int N, int L ); /* SY=<see the description> */

/* Real floating-point data, 32 bits per data element */
void rsb2x2f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=2, N=2, SY= 4 */
void rsb3x3f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=3, N=3, SY=16 */
void rsb4x4f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=4, N=4, SY=16 */
void rsb8x8f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=8, N=8, SY=64 */
void rsb2x1f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=2, N=1, SY= 2 */
void rsb4x1f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=4, N=1, SY= 4 */
void rsb8x1f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=8, N=1, SY= 8 */
void rsbmxnf ( float32_t * restrict y, const float32_t * restrict x, int M, int N, int L ); /* SY=<see the description> */

/* Complex fixed-point data, 32 bits per data element */
void csb2x2 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=2, N=2, SY= 4 */
void csb3x3 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=3, N=3, SY=16 */
void csb4x4 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=4, N=4, SY=16 */
void csb8x8 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=8, N=8, SY=64 */
void csb2x1 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=2, N=1, SY= 2 */
void csb4x1 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=4, N=1, SY= 4 */
void csb8x1 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=8, N=1, SY= 8 */
void csbmxn ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int M, int N, int L ); /* SY=<see the description> */

/* Complex floating-point data, 64 bits per data element */
void csb2x2f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=2, N=2, SY= 4 */
void csb3x3f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=3, N=3, SY=12 */
void csb4x4f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=4, N=4, SY=16 */
void csb8x8f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=8, N=8, SY=64 */
void csb2x1f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=2, N=1, SY= 2 */
void csb4x1f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=4, N=1, SY= 4 */
void csb8x1f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=8, N=1, SY= 8 */
void csbmxnf ( complex_float * restrict y, const complex_float * restrict x, int M, int N, int L ); /* SY=<see the description> */

/*-------------------------------------------------------------------------
Packed to Streaming Conversion for Real and Complex Matrices

Description: convert a sequence of L MxN matrices from packed (block) order
to streaming order. Use rbs*() functions for real data, and cbs*() functions -
for complex data.

Representation:
<r|c>bs<size>    16-bit fixed-point data
<r|c>bs<size>f   IEEE-754 Std single precision floating-point data

Storage size SX denotes the number of data elements required to store an
MxN matrix X in block order. If matrix size M*N is less than the SIMD vector
size for appropriate data type, then the storage size equals M*N rounded up
to the next power of two, otherwise SX equals M*N rounded up to the next
multiple of the SIMD vector size.

SIMD vector size:
  - for real fixed-point data 2*BBE_SIMD_WIDTH/sizeof(int16_t) == 16
  - for complex fixed-point data 2*BBE_SIMD_WIDTH/sizeof(complex_fract16) == 8
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4

Parameters:
Input:
x[L][SX]   Input sequence of MxN matrices, block order
M          Number of rows in a matrix
N          Number of columns in a matrix
L          Number of matrices
Output:
y[M*N][L]  Output sequence of matrices, streaming order

Restrictions:
x,y      Aligned on 32-byte boundary
x,y      Must not overlap
L        Must be a multiple of:
           16 for real fixed-point data 
            8 for complex fixed-point data and real floating-point data
            4 for complex floating-point data
-------------------------------------------------------------------------*/
/* Real fixed-point data, 16 bits per data element */
void rbs2x2 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=2, N=2, SX= 4 */
void rbs3x3 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=3, N=3, SX=16 */
void rbs4x4 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=4, N=4, SX=16 */
void rbs8x8 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=8, N=8, SX=64 */
void rbs2x1 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=2, N=1, SX= 2 */
void rbs4x1 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=4, N=1, SX= 4 */
void rbs8x1 ( int16_t * restrict y, const int16_t * restrict x, int L ); /* M=8, N=1, SX= 8 */
void rbsmxn ( int16_t * restrict y, const int16_t * restrict x, int M, int N, int L ); /* SX=<see the description> */

/* Real floating-point data, 32 bits per data element */
void rbs2x2f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=2, N=2, SX= 4, Sy= 4 */
void rbs3x3f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=3, N=3, SX=16, Sy= 9 */
void rbs4x4f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=4, N=4, SX=16, Sy=16 */
void rbs8x8f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=8, N=8, SX=64, Sy=64 */
void rbs2x1f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=2, N=1, SX= 2, Sy= 2 */
void rbs4x1f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=4, N=1, SX= 4, Sy= 4 */
void rbs8x1f ( float32_t * restrict y, const float32_t * restrict x, int L ); /* M=8, N=1, SX= 8, Sy= 8 */
void rbsmxnf ( float32_t * restrict y, const float32_t * restrict x, int M, int N, int L ); /* SX=<see the description> */

/* Complex fixed-point data, 32 bits per data element */
void cbs2x2 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=2, N=2, SX= 4 */
void cbs3x3 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=3, N=3, SX=16 */
void cbs4x4 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=4, N=4, SX=16 */
void cbs8x8 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=8, N=8, SX=64 */
void cbs2x1 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=2, N=1, SX= 2 */
void cbs4x1 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=4, N=1, SX= 4 */
void cbs8x1 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int L ); /* M=8, N=1, SX= 8 */
void cbsmxn ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int M, int N, int L ); /* SX=<see the description> */

/* Complex floating-point data, 64 bits per data element */
void cbs2x2f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=2, N=2, SX= 4 */
void cbs3x3f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=3, N=3, SX=12 */
void cbs4x4f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=4, N=4, SX=16 */
void cbs8x8f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=8, N=8, SX=64 */
void cbs2x1f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=2, N=1, SX= 2 */
void cbs4x1f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=4, N=1, SX= 4 */
void cbs8x1f ( complex_float * restrict y, const complex_float * restrict x, int L ); /* M=8, N=1, SX= 8 */
void cbsmxnf ( complex_float * restrict y, const complex_float * restrict x, int M, int N, int L ); /* SX=<see the description> */

/*-------------------------------------------------------------------------
M-to-1 complex/real streams interleave

Description: interleave element by element M=2,3,4 or 8 streams into a
single stream. Use rinterleave<M>() functions for real data, and
cinterleave<M>() functions - for complex data.

Representation:
<r|c>interleave<M>   16-bit fixed-point data
<r|c>interleave<M>f  IEEE-754 Std single precision floating-point data

Parameters:
Input:
M                    Number of streams
N                    Number of elements per each input stream
x0[N],...,x<M-1>[N]  Input data streams
x[M]                 M pointers to input data streams
Output:
y[M*N]               Interleaved data streams

Restrictions:
x0,...,x<M-1>,
x[0..M-1],y          Must not overlap and must be aligned on 32-byte boundary
N                    Must be a multiple of:
                       16 for real fixed-point data 
                        8 for complex fixed-point data and real floating-point data
                        4 for complex floating-point data
-------------------------------------------------------------------------*/
/* Complex-valued fixed-point functions */
void cinterleave2 ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x0,
              const complex_fract16 * restrict x1, int N );
void cinterleave3 ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x0, 
              const complex_fract16 * restrict x1,
              const complex_fract16 * restrict x2, int N );
void cinterleave4 ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x0,
              const complex_fract16 * restrict x1,
              const complex_fract16 * restrict x2,
              const complex_fract16 * restrict x3, int N );
void cinterleave8 ( complex_fract16 * restrict y, 
              const complex_fract16 * restrict x0,
              const complex_fract16 * restrict x1,
              const complex_fract16 * restrict x2,
              const complex_fract16 * restrict x3,
              const complex_fract16 * restrict x4,
              const complex_fract16 * restrict x5,
              const complex_fract16 * restrict x6,
              const complex_fract16 * restrict x7, int N );

/* Complex-valued floating-point functions */
void cinterleave2f ( complex_float * restrict y, 
               const complex_float * restrict x0,
               const complex_float * restrict x1, int N );
void cinterleave3f ( complex_float * restrict y, 
               const complex_float * restrict x0, 
               const complex_float * restrict x1,
               const complex_float * restrict x2, int N );
void cinterleave4f ( complex_float * restrict y, 
               const complex_float * restrict x0,
               const complex_float * restrict x1,
               const complex_float * restrict x2,
               const complex_float * restrict x3, int N );
void cinterleavemf ( complex_float * restrict y,
                     complex_float * restrict x[], int M, int N );

/* Real-valued fixed-point functions */
void rinterleave2 ( int16_t * restrict y, 
              const int16_t * restrict x0,
              const int16_t * restrict x1, int N );
void rinterleave3 ( int16_t * restrict y, 
              const int16_t * restrict x0,
              const int16_t * restrict x1,
              const int16_t * restrict x2, int N );
void rinterleave4 ( int16_t * restrict y, 
              const int16_t * restrict x0, 
              const int16_t * restrict x1,
              const int16_t * restrict x2, 
              const int16_t * restrict x3, int N ); 
void rinterleave8 ( int16_t * restrict y, 
              const int16_t * restrict x0,
              const int16_t * restrict x1,
              const int16_t * restrict x2,
              const int16_t * restrict x3,
              const int16_t * restrict x4,
              const int16_t * restrict x5,
              const int16_t * restrict x6,
              const int16_t * restrict x7, int N ); 

/* Real-valued floating-point functions */
void rinterleave2f ( float32_t * restrict y, 
               const float32_t * restrict x0,
               const float32_t * restrict x1, int N );
void rinterleave3f ( float32_t * restrict y, 
               const float32_t * restrict x0, 
               const float32_t * restrict x1,
               const float32_t * restrict x2, int N );
void rinterleave4f ( float32_t * restrict y, 
               const float32_t * restrict x0,
               const float32_t * restrict x1,
               const float32_t * restrict x2,
               const float32_t * restrict x3, int N );
void rinterleavemf ( float32_t * restrict y,
                     float32_t * restrict x[], int M, int N );

/*-------------------------------------------------------------------------
1-to-M complex/real streams deinterleave

Description: decompose the input data stream into M=2,3,4 or 8 output 
streams, element by element. Use rdeinterleave<M>() functions for real
data, and cdeinterleave<M>() functions - for complex data.

Representation:
<r|c>deinterleave<M>   16-bit fixed-point data
<r|c>deinterleave<M>f  IEEE-754 Std single precision floating-point data

Parameters:
Input:
M                    Number of streams
N                    Number of elements per each output stream
x[M*N]               Input data stream
Output:
y0[N],...,y<M-1>[N]  Deinterleaved data streams
y[M]                 M pointers to output data streams

Restrictions:
x, y0,...,y<M-1>,
y[0..M-1]            Must not overlap and must be aligned on 32-byte boundary
N                    Must be a multiple of:
                       16 for real fixed-point data 
                        8 for complex fixed-point data and real floating-point data
                        4 for complex floating-point data
-------------------------------------------------------------------------*/
/* Complex-valued fixed-point functions */
void cdeinterleave2 ( complex_fract16 * restrict y0, 
                      complex_fract16 * restrict y1, 
                const complex_fract16 * restrict x, int N );
void cdeinterleave3 ( complex_fract16 * restrict y0,
                      complex_fract16 * restrict y1,
                      complex_fract16 * restrict y2, 
                const complex_fract16 * restrict x, int N );
void cdeinterleave4 ( complex_fract16 * restrict y0,
                      complex_fract16 * restrict y1,
                      complex_fract16 * restrict y2,
                      complex_fract16 * restrict y3, 
                const complex_fract16 * restrict x, int N );
void cdeinterleave8 ( complex_fract16 * restrict y0,
                      complex_fract16 * restrict y1,
                      complex_fract16 * restrict y2,
                      complex_fract16 * restrict y3,
                      complex_fract16 * restrict y4,
                      complex_fract16 * restrict y5,
                      complex_fract16 * restrict y6,
                      complex_fract16 * restrict y7, 
                const complex_fract16 * restrict x, int N );

/* Complex-valued floating-point functions */
void cdeinterleave2f ( complex_float * restrict y0, 
                       complex_float * restrict y1, 
                 const complex_float * restrict x, int N );
void cdeinterleave3f ( complex_float * restrict y0,
                       complex_float * restrict y1,
                       complex_float * restrict y2, 
                 const complex_float * restrict x, int N );
void cdeinterleave4f ( complex_float * restrict y0,
                       complex_float * restrict y1,
                       complex_float * restrict y2,
                       complex_float * restrict y3, 
                 const complex_float * restrict x, int N );
void cdeinterleavemf ( complex_float * restrict y[],
                 const complex_float * restrict x, int M, int N );

/* Real-valued fixed-point functions */
void rdeinterleave2 ( int16_t * restrict y0, 
                      int16_t * restrict y1, 
                const int16_t * restrict x, int N );
void rdeinterleave3 ( int16_t * restrict y0,
                      int16_t * restrict y1,
                      int16_t * restrict y2, 
                const int16_t * restrict x, int N );
void rdeinterleave4 ( int16_t * restrict y0,
                      int16_t * restrict y1,
                      int16_t * restrict y2,
                      int16_t * restrict y3, 
                const int16_t * restrict x, int N );
void rdeinterleave8 ( int16_t * restrict y0,
                      int16_t * restrict y1,
                      int16_t * restrict y2,
                      int16_t * restrict y3,
                      int16_t * restrict y4,
                      int16_t * restrict y5,
                      int16_t * restrict y6,
                      int16_t * restrict y7, 
                const int16_t * restrict x, int N );

/* Real-valued floating-point functions */
void rdeinterleave2f ( float32_t * restrict y0, 
                       float32_t * restrict y1, 
                 const float32_t * restrict x, int N );
void rdeinterleave3f ( float32_t * restrict y0,
                       float32_t * restrict y1,
                       float32_t * restrict y2, 
                 const float32_t * restrict x, int N );
void rdeinterleave4f ( float32_t * restrict y0,
                       float32_t * restrict y1,
                       float32_t * restrict y2,
                       float32_t * restrict y3, 
                 const float32_t * restrict x, int N );
void rdeinterleavemf ( float32_t * restrict y[], 
                 const float32_t * restrict x, int M, int N );

#ifdef __cplusplus
};
#endif

#endif /* __NATUREDSP_BASEBAND_MATOP_H */
