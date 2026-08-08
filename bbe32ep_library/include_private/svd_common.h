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
    NatureDSP_Baseband library. Singular Value Decomposition
    Common declarations
    IntegrIT, 2006-2017
*/

#ifndef __SVD_COMMON_H
#define __SVD_COMMON_H

/* Portable data types. */
#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#if HAVE_VFPU

/*
 * Reduce MxN matrices with M>=N to upper bidiagonal form by left- and
 * right-hand Householder transformations.
 * Original matrix A, upper bidiagonal form B, left- and right-hand 
 * transformation matrices U and V relate to each other through the 
 * following identities: A == U*B*V', U'*A*V == B.
 * Matrix U is comprised of orthonormal columns, matrix V is orthogonal
*  (unitary for complex data).
 * Input:
 *   M,N         Matrix dimensions
 *   L           Number of matrices
 *   needU       Set to non-zero value if left-hand transformation matrices
 *               U are required (applies to stream order functions)
 * Input/Output:
 *   A[M*N]xL    Original matrices (in); left-hand transformation matrices
 *               U (out, optional, stream order functions). Even when 
 *               not utilized to keep the resulting matrix U, A is still used
 *               as an intermediate storage, so input data are damaged in any
 *               case.
 * Output:
 *   D[N]xL      Main diagonal of reduced matrices
 *   F[N-1]xL    First superdiagonal of reduced matrices
 *   U[M*N]xL    Left-hand transformation matrix (block order functions). Set
 *               to zero if not needed.
 *   V[N*N]xL    Orthogonal matrix of accumulated right-hand transformations,
 *               optional. Set to zero if not needed.
 * Restrictiions:
 *   D,F,U,V,A   Must not overlap and must be aligned on 2*BBE_SIMD_WIDTH-byte
 *               boundary
 *   N>1         Input matrices must have at least 2 columns
 *   M>=N        Number of columns must not exceed the number of rows
 *   M,N,L       Subject to additional limitations exposed by a particular
 *               function
 */

/* Complex-valued, Stream Order
 * Note that output matrices U replace input matrices A.
 * Restrictions: 
 *   L  Must be a multiple of BBE_SIMD_WIDTH/4 */
void svd_bidiag_mxnsf (
                complex_float * restrict D, /* D[N][L]   */
                complex_float * restrict F, /* F[N-1][L] */
                complex_float * restrict V, /* V[N*N][L] */
                complex_float * restrict A, /* A[M*N][L] */
                int M, int N, int L, int needU );

/* Real-valued, Stream Order.
 * Note that output matrices U replace input matrices A.
 * Restrictions: 
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2 */
void rsvd_bidiag_mxnsf (
                float32_t * restrict D, /* D[N][L]   */
                float32_t * restrict F, /* F[N-1][L] */
                float32_t * restrict V, /* V[N*N][L] */
                float32_t * restrict A, /* A[M*N][L] */
                int M, int N, int L, int needU );

/*
 * Golub-Reinsch SVD for real or complex upper bidiagonal matrices stored 
 * in stream order. Functions compute the SVD for L matrices, each of M rows
 * and N columns.
 *
 * These functions implement the Golub-Reinsch SVD algorithm as stated in:
 * [1] "Matrix Computations" by G.H. Golub and C.F. Van Loan, 4-th Edition
 * [2] "Singular Value Decomposition and Least Squares Solutions" by 
 *     G.H. Golub abd C. Reinsch, published in "Handbook for Automatie 
 *     Computation", Vol.II, Contribution I/10.
 *
 * Temporary:
 *   pScr            Scratch area.  Required size (in bytes) is defined by 
 *                   functions [r]grsvdsf_getScratchSize(M,N,L)
 * Input:
 *   M,N             Matrix dimensions
 *   L               Number of matrices
 *   F[N-1][L]       First superdiagonal of input matrices. Note that data
 *                   will be distorted by intermediate results.
 * Input/Output:
 *   D[N][L]         Main diagonal of input matrices (in); singular
 *                   values in descending order, or NaNs if failed to
 *                   converge (out, applies to the real-valued variant).
 *   U[M*N][L]       (OPTIONAL) left-hand transformation matrix of the 
 *                   bidiagonalization transform (in); matrix of left-
 *                   singular orthonormal vectors (out).
 *   V[N*N][L]       (OPTIONAL) right-hand transformation matrix of the
 *                   bidiagonalization transform (in); orthogonal (unitary)
 *                   matrix of right-singular vectors (out).
 * Output:
 *   s[N][L]         Singular values in descending order, or NaNs if failed to
 *                   converge (applies to the complex-vlaued variant)
 * Restrictiions:
 *   N>1             Input matrices must have at least 2 columns
 *   M>=N            Number of columns must not exceed the number of rows
 *   L               Must be a multiple of BBE_SIMD_WIDTH/2
 *   pScr,s,D,F,U,V  Must not overlap and must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf ( void * pScr,
               float32_t     * restrict s,
               complex_float * restrict D,
               complex_float * restrict F,
               complex_float * restrict U,
               complex_float * restrict V,
               int M, int N, int L );

void rgrsvdsf( void * pScr,
               float32_t * restrict D,
               float32_t * restrict F,
               float32_t * restrict U,
               float32_t * restrict V,
               int M, int N, int L );

size_t grsvdsf_getScratchSize ( int M, int N, int L );
size_t rgrsvdsf_getScratchSize( int M, int N, int L );

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * Summarize index data and fetch data for Wilkinson's shift computation.
 * Input:
 *   itsLim       Iterations count limit
 *   N            Number of columns in an input matrix
 *   L            Number of matrices
 *   D[N][L]      Main diagonal of input matrices
 *   F[N-1][L]    First superdiagonal of input matrices
 *   a_l[L]       Left index of working subblocks
 *   a_k[L]       Right index of working subblocks
 * Input/Output:
 *   a_its[L]     Iteration counters
 * Output:
 *   p_l_lo       Lowest left index over input matrices
 *   p_k_up       Uppermost right index over input matrices
 *   a_m[N][L/(BBE_SIMD_WIDTH/4)] (complex variant)
 *   a_m[N][L/(BBE_SIMD_WIDTH/2)] (real variant)
 *                Boolean labels of Working subblocks
 *   a_a[L]       Fetched data: D[a_k[0..L-1]-1] 
 *   a_b[L]       Fetched data: D[a_k[0..L-1]] 
 *   a_x[L]       Fetched data: D[a_l[0..L-1]] 
 *   a_y[L]       Fetched data: F[a_l[0..L-1]] 
 *   a_c[L]       Fetched data: F[a_k[0..L-1]-1] 
 *   a_s[L]       Fetched data: F[a_k[0..L-1]-2] 
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_its,rng_l,rng_k,a_m,a_a,a_b,a_x,a_y,a_c,a_s,a_l,a_k,D,F
 *      Must not overlap
 *   a_its,a_a,a_b,a_x,a_y,a_c,a_s,a_l,a_k,D,F
 *      Must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf_gks_fetch( 
                    int16_t       * restrict p_l_lo,
                    int16_t       * restrict p_k_up,
                    int16_t       * restrict a_its,
                    vboolN_4      * restrict a_m,
                    complex_float * restrict a_a,
                    complex_float * restrict a_b,
                    complex_float * restrict a_x,
                    complex_float * restrict a_y,
                    complex_float * restrict a_c,
                    complex_float * restrict a_s,
              const int16_t       *          a_l,
              const int16_t       *          a_k,
              const complex_float *          D,
              const complex_float *          F,
              int itsLim, int N, int L );

void rgrsvdsf_gks_fetch( 
                    int16_t   * restrict p_l_lo,
                    int16_t   * restrict p_k_up,
                    int16_t   * restrict a_its,
                    vboolN_2  * restrict a_m,
                    float32_t * restrict a_a,
                    float32_t * restrict a_b,
                    float32_t * restrict a_x,
                    float32_t * restrict a_y,
                    float32_t * restrict a_c,
                    float32_t * restrict a_s,
              const int16_t   *          a_l,
              const int16_t   *          a_k,
              const float32_t *          D,
              const float32_t *          F,
              int itsLim, int N, int L );

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * Compute Wilkinson's shift.
 * Input:
 *   L            Number of matrices
 *   D[N][L]      Main diagonal of input matrices
 *   F[N-1][L]    First superdiagonal of input matrices
 *   a_c[L]       Fetched data: F[a_k[0..L-1]-1].  May be reused
 *                as a temporal storage of intermediate results.
 *   a_s[L]       Fetched data: F[a_k[0..L-1]-2]. May be reused
 *                as a temporal storage of intermediate results.
 *   a_x[L]       Fetched data: D[a_l[0..L-1]] 
 *   a_y[L]       Fetched data: F[a_l[0..L-1]] 
 * Input/Output:
 *   a_a[L]       In:  fetched data: D[a_k[0..L-1]-1] 
 *                Out: element (0,0) of B'*B-mu*I
 *   a_b[L]       In:  fetched data: D[a_k[0..L-1]] 
 *                Out: element (0,1) of B'*B-mu*I
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_a,a_b,a_x,a_y,a_c,a_s
 *      Must not overlap and must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf_gks_wilkShift(
                    complex_float * restrict a_a,
                    complex_float * restrict a_b,
                    complex_float * restrict a_c,
                    complex_float * restrict a_s,
              const complex_float *          a_x,
              const complex_float *          a_y,
              int L );

void rgrsvdsf_gks_wilkShift(
                    float32_t * restrict a_a,
                    float32_t * restrict a_b,
                    float32_t * restrict a_c,
                    float32_t * restrict a_s,
              const float32_t *          a_x,
              const float32_t *          a_y,
              int L );

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * Compute Givens's rotating matrix.
 * For real data:
 *   G(a,b) <- [c,s;-s,c]: [a,b]*[c,s;-s,c] == [*,0], c^2+s^2 == 1
 * For complex data:
 *   G(a,b) <- [c,conj(s);-s,conj(c)]: [a,b]*G(a,b) == [*,0], c*conj(c)+s*conj(s) == 1
 * Temporary:
 *   a_x[L]       Scratch array of L entries
 * Input:
 *   L            Number of matrices
 *   a_a[L]       a-values
 *   a_b[L]       b-values
 * Output:
 *   a_c[L]       G(a,b) cosine values
 *   a_s[L]       G(a,b) sine values
 *  Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_a,a_b,a_c,a_s
 *      Must not overlap and must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf_gks_givens( 
                    complex_float * restrict a_x,
                    complex_float * restrict a_c,
                    complex_float * restrict a_s,
              const complex_float *          a_a, 
              const complex_float *          a_b, 
              int L );

void rgrsvdsf_gks_givens( 
                    float32_t * restrict a_c,
                    float32_t * restrict a_s,
              const float32_t *          a_a, 
              const float32_t *          a_b, 
              int L );

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * 1st updating step of SVD step iteration: B <- B*G(a,b).
 * Input:
 *   n            Current position
 *   L            Number of matrices
 *   a_c[L]       G(a,b) cosine values
 *   a_s[L]       G(a,b) sine values
 *   a_l[L]       Left index of working subblocks
 *   a_m[N][L/(BBE_SIMD_WIDTH/4)] (complex variant)
 *   a_m[N][L/(BBE_SIMD_WIDTH/2)] (real variant)
 *                Boolean labels of Working subblocks
 * Input/Output:
 *   a_a[L]       In: B(n-1,n) Out: B(n,n)  
 *   a_b[L]       In: B(n-1,n+1) Out: B(n+1,n)
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_a,a_b,a_c,a_s,a_l,a_m
 *      Must not overlap
 *   a_a,a_b,a_c,a_s,a_l
 *      Must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf_gks_step1(
                    complex_float * restrict D,
                    complex_float * restrict F,
                    complex_float * restrict a_a,
                    complex_float * restrict a_b,
              const complex_float *          a_c,
              const complex_float *          a_s,
              const int16_t       *          a_l,
              const vboolN_4      *          a_m,
              int n, int L );

void rgrsvdsf_gks_step1(
                    float32_t * restrict D,
                    float32_t * restrict F,
                    float32_t * restrict a_a,
                    float32_t * restrict a_b,
              const float32_t *          a_c,
              const float32_t *          a_s,
              const int16_t   *          a_l,
              const vboolN_2  *          a_m,
              int n, int L );

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * 2nd updating step of SVD step iteration: B <- G(a,b)'*B.
 * Input:
 *   n            Current position
 *   L            Number of matrices
 *   a_c[L]       G(a,b) cosine values
 *   a_s[L]       G(a,b) sine values
 *   a_k[L]       Right index of working subblocks
 *   a_m[N][L/(BBE_SIMD_WIDTH/4)] (complex variant)
 *   a_m[N][L/(BBE_SIMD_WIDTH/2)] (real variant)
 *                Boolean labels of Working subblocks
 * Input/Output:
 *   a_a[L]       In: B(n,n) Out: B(n,n+1)  
 *   a_b[L]       In: B(n+1,n) Out: B(n,n+2)
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_a,a_b,a_c,a_s,a_k,a_m
 *      Must not overlap
 *   a_a,a_b,a_c,a_s,a_k
 *      Must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf_gks_step2(
                    complex_float * restrict D,  
                    complex_float * restrict F,
                    complex_float * restrict a_a,
                    complex_float * restrict a_b,
              const complex_float *          a_c,
              const complex_float *          a_s,
              const int16_t       *          a_k,
              const vboolN_4      *          a_m,
              int n, int L );

void rgrsvdsf_gks_step2(
                    float32_t * restrict D,  
                    float32_t * restrict F,
                    float32_t * restrict a_a,
                    float32_t * restrict a_b,
              const float32_t *          a_c,
              const float32_t *          a_s,
              const int16_t   *          a_k,
              const vboolN_2  *          a_m,
              int n, int L );

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * Accumulate transformations: U <- U*conj(G(a,b)) or V <- V*G(a,b),
 * where U (MxN) and V (NxN) are matrices comprised of left- and right
 * singular vectors. For real data, a single function is used for both
 * U and V updates. For complex data there are two separate functions
 * for matrices U and V.
 * Input:
 *   n            Current position
 *   M,N          Matrix dimensions
 *   L            Number of matrices
 *   a_c[L]       G(a,b) cosine values
 *   a_s[L]       G(a,b) sine values
 *   a_m[N][L/(BBE_SIMD_WIDTH/4)] (complex variant)
 *   a_m[N][L/(BBE_SIMD_WIDTH/2)] (real variant)
 *                Boolean labels of Working subblocks
 * Input/Output:
 *   W[M*N][L]    Matrix of left-singular (U) or right-singular (V) 
 *                orthonormal vectors
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   W,a_c,a_s,a_m
 *      Must not overlap
 *   W,a_c,a_s
 *      Must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void grsvdsf_gks_accumLeft(
                    complex_float * restrict W,
              const complex_float *          a_c,
              const complex_float *          a_s,
              const vboolN_4      *          a_m,
              int n, int M, int N, int L );

void grsvdsf_gks_accumRight(
                    complex_float * restrict W,
              const complex_float *          a_c,
              const complex_float *          a_s,
              const vboolN_4      *          a_m,
              int n, int N, int L );

void rgrsvdsf_gks_accum(
                    float32_t * restrict W,
              const float32_t *          a_c,
              const float32_t *          a_s,
              const vboolN_2  *          a_m,
              int n, int M, int N, int L );

/*
 * Convert matrices from stream to block order. Note that the number of
 * matrices to be converted and the total number of matrices in the
 * input stream are specified through separate input arguments: L and
 * stride, respectively.
 * Input:
 *   M,N             Matrix size
 *   L               Number of matrices to be converted
 *   stride          Number of matrices in the input stream
 *   x[M*N][stride]  Input matrices in stream order
 * Output:
 *   y[L][SY]        Output matrices in block order
 * Restrictions:
 *   y,x             Must not overlap and must be aligned on 
 *                   2*BBE_SIMD_WIDTH-byte boundary
 *   L<=stride       Number of converted matrices cannot exceed the number
 *                   of matrices in the input stream
 *   M,N             The product M*N must be a multiple of 4
 * where SY = denotes the number of data entries to store M*N elements in
 * block order with proper alignment.
 */

/* Convert complex data from stream order to block order.
 * Restrictions: 
 *   stride  Must be a multiple of BBE_SIMD_WIDTH/4 */
void svd_csbmxnxsf( complex_float * restrict y, 
              const complex_float *          x,
              int M, int N, int L, int stride );

/* Convert real data from stream order to block order.
 * Restrictions: 
 *   stride  Must be a multiple of BBE_SIMD_WIDTH/2 */
void svd_rsbmxnxsf( float32_t * restrict y, 
              const float32_t * x, 
              int M, int N, int L, int stride );

#endif /* HAVE_VFPU */

#ifdef __cplusplus
};
#endif

#endif /* __SVD_COMMON_H */
