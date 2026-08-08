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
    NatureDSP_Baseband library. Eigenvalues and eigenvectors
    Common declarations
    IntegrIT, 2006-2017
*/

#ifndef __EIGEN_COMMON_H
#define __EIGEN_COMMON_H

/* Portable data types. */
#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Reduce square input matrices to upper-Hessenberg form by a unitary 
 * similarity transformation: H <- P'*A*P. Resulting matrices H replace
 * input matrices A.
 * Input:
 *   N    Matrix size
 *   L    Number of matrices
 * Input/Output:
 *   A[]  General square matrices (in), reduced matrices (out). 
 *        Values below the first subdiagonal of resulting 
 *        matrices are not defined
 * Output:
 *   P[]  Unitary/orthogonal trasformation matrix, This argument is 
 *        OPTIONAL, set it to zero if matrix P is not required
 * Restrictions:
 *   A,P  Must not overlap and must be aligned on 2*BBE_SIMD_WIDTH-byte
 *        boundary
 *   N>2  Minimum maxtrix size is 3x3
 *   N,L  Subject to additional limitations exposed by a particular function
 */
#if 0
/* Complex-valued, Packed Order
 * Restrictions: 
 *   N  Must be a positive multiple of 4 */
void eigen_hess_nxnnf ( 
                complex_float * restrict P, /* P[L][S(sz_f32c,N*N)] */
                complex_float * restrict A, /* A[L][S(sz_f32c,N*N)] */
                int N, int L );

/* Real-valued, Packed Order
 * Restrictions: 
 *   N  Must be a positive multiple of 4 */
void reigen_hess_nxnnf (
                float32_t * restrict P,     /* P[L][S(sz_f32,N*N)] */
                float32_t * restrict A,     /* A[L][S(sz_f32,N*N)] */
                int N, int L );
#endif
/* Complex-valued, Stream Order
 * Restrictions: 
 *   L  Must be a multiple of BBE_SIMD_WIDTH/4 */
void eigen_hess_nxnsf ( 
                complex_float * restrict P, /* P[N*N][L] */
                complex_float * restrict A, /* A[N*N][L] */
                int N, int L );

/* Real-valued, Stream Order
 * Restrictions: 
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2 */
void reigen_hess_nxnsf ( 
                float32_t * restrict P,     /* P[N*N][L] */
                float32_t * restrict A,     /* A[N*N][L] */
                int N, int L );

/*
 * Convert square upper-Hessenberg matrices from full stream to compact
 * block order. Compact storage does not keep zero elements below the first 
 * subdiagonal. Note that there are N*(N+3)/2-1 non-zero elsements an NxN
 * upper-Hessenberg matrix.
 * Input:
 *   N          Matrix size
 *   L          Number of matrices
 *   x[N*N][L]  Input matrices in full stream order. Implementation may perform
 *              in-place transformations of input matrices, so INPUT DATA MAY APPEAR
 *              DAMAGED after the call.
 * Output:
 *   y[L][SY]   Output matrices in compact block order,
 * Restrictions:
 *   y,x  Must not overlap and must be aligned on 2*BBE_SIMD_WIDTH-byte
 *        boundary
 *   L    Must be a multiple of BBE_SIMD_WIDTH/2 for real data, or
 *        a multiple of BBE_SIMD_WIDTH/4 for complex data
 * where SY = denotes the number of data elements needed to store all 
 * non-zerp elements of an NxN upper-Hessenberg matrix in block order 
 * with proper alignment.
 */

void eigen_s2hn_nxnf ( complex_float * restrict y, /* y[L][S(sz_f32c, N*(N+3)/2-1)] */
                       complex_float * restrict x, /* x[N*N][L]                     */
                 int N, int L );

void reigen_s2hn_nxnf ( float32_t * restrict y,     /* y[L][S(sz_f32, N*(N+3)/2-1)] */
                        float32_t * restrict x,     /* x[N*N][L]                    */
                  int N, int L );

/*
 * Extract the upper triangular part of a Hessenberg matrix stored in
 * compact block order and convert it to compact stream order. 
 * Compactness impiles that elements below the first subdiagonal
 * (upper-Hessenberg) or below the main diagonal (upper triangular
 * matrix) are not stored in memory.
 * Number of payload elements in an NxN upper-Hessenberg matrix: N*(N+3)/2-1.
 * Number of elements in the upper triangular part of an NxN matrix, including
 * the main diagonal: N*(N+1)/2.
 * Input:
 *   N                  Matrix size
 *   x[L][SX]           Upper-Hessenberg matrices in block order
 * Output:
 *   y[N*(N+3)/2-1][L]  Upper triangular part of input matrices in stream order.
 *                      After the function completes, results are contained in
 *                      the first N*(N+1)/2*L elements of the array.
 * Restrictions:
 *   y,x  Must not overlap and must be aligned on 2*BBE_SIMD_WIDTH-byte
 *        boundary
 *   L    Must be a multiple of BBE_SIMD_WIDTH/4
 * where SX = S(sizeof(complex_float),N*(N+3)/2-1) denotes the number of
 * data elements needed to store all non-zerp elements of an NxN 
 * upper-Hessenberg matrix in block order with proper alignment.
 */

void eigen_hn2ts_nxnf ( complex_float * restrict y,
                  const complex_float * restrict x,
                  int N, int L );

/*
 * Francis QR algorithm for complex upper-Hessenberg matrix. 
 * For NxN matrix A, the function computes the full set of N eigenvalues e[N]
 * from its upper-Hessenberg form H. It also outputs a triangular form T and an
 * updated unitary matrix P (optional) such that A == P*T*P'. 
 * Input:
 *   N               Matrix size
 *   low,upp         Lower/upper boundaries of block diagonal structure. 
 *                   Set to 0,N-1 if QR should be applied to full matrix H.
 * Input/Output:
 *   H[N*(N+3)/2-1]  Upper-Hessenberg form (in), triangular form (out). Zero
 *                   elements below the first subdiagonal are not stored.
 *                   Values on the first subdiagonal of triangular output
 *                   form are not specified.
 *   P[N*N]          Transformation matrix reducing original matrix to 
 *                   upper-Hessenberg form (in), or to triangular form T (out).
 *                   This argument is optional, set to zero if transformation
 *                   matrix is not required
 * Output:
 *   e[N]            Eigenvalues of original matrix. Complex eigenvalues
 *                   come in conjugate pairs. If the QR algorithm fails to 
 *                   converge, real and imaginary components of ev[0..N-1]
 *                   are set to NaN.
 * Return Value:
 *   Non-zero if successfull, or zero if the QR algorithm failed to converge.
 * Restrictions:
 *   0<=low<=upp<N
 */

int eigen_hqr_f ( complex_float * restrict e, 
                  complex_float * restrict H, 
                  complex_float * restrict P, 
                  int low, int upp, int N );

/*
 * Francis QR algorithm for real upper-Hessenberg matrix. 
 * For NxN matrix A, the function computes the full set of N eigenvalues e[N]
 * from its upper-Hessenberg form H. It also outputs a quasi-triangular form T 
 * and an updated orthogonal matrix P (optional) such that A == P*T*P'. 
 * Input:
 *   N               Matrix size
 *   low,upp         Lower/upper boundaries of block diagonal structure. 
 *                   Set to 0,N-1 if QR should be applied to full matrix H.
 * Input/Output:
 *   H[N*(N+3)/2-1]  Upper-Hessenberg form (in), quasi-triangular form with
 *                   2x2 blocks on the main diagonal (out). Zero elements
 *                   below the first subdiagonal are not stored.
 *   P[N*N]          Transformation matrix reducing original matrix to 
 *                   upper-Hessenberg form (in), or to quasi-triangular
 *                   form T (out). This argument is optional, set to zero
 *                   if transformation matrix is not required
 * Output:
 *   e[N]            Eigenvalues of original matrix. Complex eigenvalues
 *                   come in conjugate pairs. If the QR algorithm fails to 
 *                   converge, real and imaginary components of ev[0..N-1]
 *                   are set to NaN.
 * Return Value:
 *   Non-zero if successfull, or zero if the QR algorithm failed to converge.
 * Restrictions:
 *   0<=low<=upp<N
 */

int reigen_hqr_f ( complex_float * restrict e, 
                   float32_t     * restrict H, 
                   float32_t     * restrict P, 
                   int low, int upp, int N );

/* 
 * Determine eigenvectors of a triangular complex form T by in-place 
 * backsubstitution process. 
 * Input:
 *   N               Matrix size
 * Input/Output:
 *   T[N*(N+3)/2-1]  NxN triangular form T (in); N column eigenvectors (out).
 *                   Zeros below the first subdiagonal are not stored. Values
 *                   on the subdiagonal are not specified for both input and
 *                   output.
 */

void eigen_bksubst_f ( complex_float * restrict T, int N );

/* 
 * Determine eigenvectors of a quasi-triangular real form T by in-place 
 * backsubstitution process. Real part of a complex eigenvector is stored
 * in the first vector of the conjugate pair, and the imaginary part is
 * stored in the second vector. 
 * Input:
 *   N               Matrix size
 *   e[N]            Eigenvalues, either real or in conjugate pairs
 * Input/Output:
 *   T[N*(N+3)/2-1]  NxN quasi-triangular form T with 2x2 blocks on the main
 *                   diagonal (in); N column eigenvectors (out). Zeros below
 *                   the first subdiagonal are not stored
 */

void reigen_bksubst_f ( float32_t     * restrict T, 
                  const complex_float * restrict e,
                  int N );

/*
 * Left-multiply real square upper-Hessenberg matrices H by orthogonal 
 * matrices P and store resulting matrices to U: U <- P*H. 
 * Note:
 *   TBD If balancing is implemented, the complexity of this functions may
 *       be reduced by taking into account the block diagonal structure of
 *       input matrices. See the MATLAB reference.
 * Input:
 *   N      Matrix size
 *   L      Number of matrices
 *   P[]    Orthogonal transformation matrices
 *   H[]    Upper-Hessenberg matrices. Zero elements below the first 
 *          subdiagonal are not stored
 * Output
 *   U[]    Resulting matrices
 * Restrictions:
 *   U,P,H  Must not overlap and must be aligned on 2*BBE_SIMD_WIDTH-byte
 *          boundary
 *   Variant functions may impose additional restrictions
 */

/* Real Data, Stream Order
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2 */
void reigen_hmulp_nxnsf ( float32_t * restrict U, /* U[N*N][L]         */
                    const float32_t * restrict P, /* P[N*N][L]         */
                    const float32_t * restrict H, /* H[N*(N+3)/2-1][L] */
                    int N, int L );
#if 0
/* Real Data, Packed Order
 * Restrictions:
 *   N  Must be a positive multiple of 4 */
void reigen_hmulp_nxnnf ( float32_t * restrict U, /* U[L][S(N*N)]         */
                    const float32_t * restrict P, /* P[L][S(N*N)]         */
                    const float32_t * restrict H, /* H[L][S(N*(N+3)/2-1)] */
                    int N, int L );
#endif
/*
 * Left-multiply complex square upper triangular matrices T by unitary
 * matrices P, then scale each column of a product matrix so that its
 * L2 norm is 1, and store results to output matrices: V <- P*T*S, where
 * S is a diagonal matrix with scaling factors on the main diagonal.
 * Note:
 *   TBD If balancing is implemented, the complexity of this functions may
 *       be reduced by taking into account the block diagonal structure of
 *       input matrices. See the MATLAB reference.
 * Input:
 *   N      Matrix size
 *   L      Number of matrices
 *   P[]    Unitary transformation matrices
 *   T[]    Upper triangular matrices. See function definitions for info
 *          on the storage format
 * Output
 *   V[]    Resulting matrices
 * Restrictions:
 *   V,P,H  Must not overlap and must be aligned on 2*BBE_SIMD_WIDTH-byte
 *          boundary
 *   Variant functions may impose additional restrictions
 */

/* Complex Data, Stream Order.
 * Elements of T below the MAIN DIAGONAL are not stored
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/4 */
void eigen_tmulp_nxnsf( complex_float * restrict V, /* V[N*N][L]       */
                  const complex_float * restrict P, /* P[N*N][L]       */
                  const complex_float * restrict T, /* T[N*(N+1)/2][L] */
                  int N, int L );
#if 0
/* Complex Data, Packed Order
 * Elements of T below the FIRST SUBDIAGONAL are not stored.
 * Values on the first subdiagonal of matrix T do not matter.
 * Restrictions:
 *   N  Must be a positive multiple of 4 */
void eigen_tmulp_nxnnf( complex_float * restrict V, /* V[N*N][L]         */
                  const complex_float * restrict P, /* P[N*N][L]         */
                  const complex_float * restrict T, /* T[N*(N+3)/2-1][L] */
                  int N, int L );
#endif
/*
 * Combine eigenvectors of real matrices from re/im components. The function
 * exploits the special order of conjugate eigenvalues/eigenvectors by rhqr()
 * function. Eigenvectors are also rescaled so that L2 norm of a vector is 1.
 * Input:
 *   N          Matrix size
 *   L          Number of matrices
 *   U[N*N][L]  Real matrices of re/im components of eigenvectors
 *   e[N][L]    Complex eigenvalues
 * Output:
 *   V[N*N][L]  Complex matrices of eigenvectors
 * Restrictions:
 *   V,U,e      Must not overlap and must be aligned on 2*BBE_SIMD_WIDTH-byte
 *              boundary
 *   Variant functions may impose additional restrictions
 */

/* RealData, Stream Order
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2 */
void reigen_evcomb_nxnsf ( complex_float * restrict V, /* V[N*N][L] */
                     const float32_t     * restrict U, /* U[N*N][L] */
                     const complex_float * restrict e, /* e[N][L]   */
                     int N, int L );
#if 0
/* Real Data, Packed Order
 * Restrictions:
 *   N  Must be a positive multiple of 4 */
void reigen_evcomb_nxnnf ( complex_float * restrict V, /* V[L][S(N*N)] */
                     const float32_t     * restrict U, /* U[L][S(N*N)] */
                     const complex_float * restrict e, /* e[L][S(N)]   */
                     int N, int L );
#endif
#ifdef __cplusplus
};
#endif

#endif /* __EIGEN_COMMON_H */
