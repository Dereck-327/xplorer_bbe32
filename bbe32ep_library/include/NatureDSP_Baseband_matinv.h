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
 * Matrix Decomposition and Inversion Functions
 */

#ifndef __NATUREDSP_BASEBAND_MATINV_H
#define __NATUREDSP_BASEBAND_MATINV_H

#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
Direct Matrix Inversion
cmatinv             Direct matrix inversion for complex matrices
matinv              Direct matrix inversion for real matrices
chermmatinv         Direct matrix inversion for Hermitian matrices

Matrix Inversion By Gauss-Jordan algortihm
[c]matinvgjnxnn     Gauss-Jordan algortihm for real/complex block ordered matrices
[c]matinvgjnxns     Gauss-Jordan algortihm for real/complex stream ordered matrices

LU Decomposition
[c]lunxnn           LU decomposition for real/complex block ordered matrices
[c]lunxns           LU decomposition for real/complex stream ordered matrices

Matrix determinant
[c]detnxnn          Determinant for real/complex block ordered matrices
[c]detnxns          Determinant for real/complex stream ordered matrices

Cholesky decomposition for complex stream ordered fixed-point matrices
cholmxns            Cholesky decomposition
cholfwdmxnxps       Cholesky forward recursion
cholbkwnxps         Cholesky backward recursion

Cholesky decomposition for complex block ordered fixed-point matrices
cholmxnn            Cholesky decomposition
cholfwdmxnxpn       Cholesky forward recursion
cholbkwnxpn         Cholesky backward recursion

Cholesky decomposition for real/complex stream ordered floating-point matrices
[r]cholmxnsf        Cholesky decomposition
[r]cholfwdmxnxpsf   Cholesky forward recursion
[r]cholbkwnxpsf     Cholesky backward recursion
[r]cholmmsemxnxpsf  Cholesky MMSE solver

Cholesky decomposition for real/complex block ordered floating-point matrices
[r]cholmxnnf        Cholesky decomposition
[r]cholfwdmxnxpnf   Cholesky forward recursion
[r]cholbkwnxpnf     Cholesky backward recursion
[r]cholmmsemxnxpsnf Cholesky MMSE solver

Banded Cholesky decomposition for complex block ordered fixed-point matrices
bcholwxnn           Banded Cholesky decomposition
bcholfwdwxnxpn      Banded Cholesky forward recursion
bcholbkwwxnxpn      Banded Cholesky backward recursion

QR decomposition for real/complex stream ordered fixed-point matrices
[c]qr2x2s           QR decomposition for 2x2 matrices
[c]qr_build_rmxns   Similar to [c]qrs, but returns Householder vectors instead of factor Q
[c]qr_calc_qbmxnxps Apply Householder reflections to right hand vector B
[c]qr_bkwnxps       QR back substitution

QR decomposition for complex block ordered fixed-point matrices
cqr_build_rmxnn     QR decomposition
cqr_calc_qbmxnxpn   Apply Householder reflections to right hand vector B
cqr_bkwmxnxpn       QR back substitution

QR decomposition for real/complex stream ordered floating-point matrices
[c]qr_build_rmxnsf     QR decomposition
[c]qr_calc_qbmxnxpsf   Apply Householder reflections to right hand vector B
[c]qr_bkwmxnxpsf       QR back substitution

QR decomposition for real/complex block ordered floating-point matrices
[c]qr_build_rmxnnf     QR decomposition
[c]qr_calc_qbmxnxpnf   Apply Householder reflections to right hand vector B
[c]qr_bkwmxnxpnf       QR back substitution

Singular Value Decomposition (SVD)
[r]svdmxnn          Thin SVD for real/complex block ordered matrices
[r]svdmxns          Thin SVD for real/complex stream ordered matrices

Eigenvalues and eigenvectors
[ê]eigennxnn        Eigenvalues and eigenvectors of real/complex block ordered matrices
[ê]eigennxns        Eigenvalues and eigenvectors of real/complex stream ordered matrices
===========================================================================*/

/*===========================================================================
Direct Matrix Inversion
cmatinv   Direct matrix inversion for complex matrices
matinv    Direct matrix inversion for real matrices
chermmatinv  Direct matrix inversion for Hermitian matrices
===========================================================================*/

/*-------------------------------------------------------------------------
Direct Matrix Inversion For Complex Matrices

Description: perform in-place inversion of 2x2 and 4x4 complex matrices. 
2x2 matrices are inverted by Cramer's rule. For 3x3, 4x4 matrices we employ 
the blockwise inversion algorithm encompassed with a suboptimal row/column
permutation that gains better conditioning of the block structure.

Data format and order options:
  Suffix   Data Order                 Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Notes:
1. In general, accuracy of a matrix inversion algorithm implementation is a
   function of input matrix condition number. Thus it is user's responsibility
   to qualify the reliability of numeric results. Refer to NatureDSP Baseband 
   Library Reference for details. 
2. For blockwise inversion of 4x4 fixed-point matrices, it is reasonable to
   limit the dynamic range of input data by 11..13 significant bits. This
   measure reduces the possibility of an overflow at internal computations.

Parameters:
Temporary:
pScr      Scratch memory area. To determine the scratch area size required by
          a function <fun>, use the respective helper function 
          <fun>_getScratchSize()
Input:
L         Number of matrices
qA        Number of fractional bits for fixed-point input/output data
Input/Output:
A[L][SA]  Sequence of L NxN complex input/result matrices. SA is the number
          of data elements occupied by a single NxN matrix in a block 
          (stream) ordered sequence, see function specifications.

Restrictions:
pScr,A    Aligned on 32-byte boundary

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/
/* Block Order, Fixed-Point, 2x2, SA=4
   Restrictions: 
     L must be even
*/
void cmatinv2x2n ( void * restrict pScr, complex_fract16 * restrict A, int qA, int L );

/* Return the scratch area size, in bytes. */
size_t cmatinv2x2n_getScratchSize ( int L );

/* Block Order, Fixed-Point, 4x4, SA=16
   Restrictions:
     None
*/
void cmatinv4x4n ( void * restrict pScr, complex_fract16 * restrict A, int qA, int L );

/* Return the scratch area size, in bytes. */
size_t cmatinv4x4n_getScratchSize ( int L );

/* Block Order, Floating-Point, 2x2, SA=4
   Restrictions: 
     None
*/
void cmatinv2x2nf ( void * restrict pScr, complex_float * restrict A, int L );

/* Return the scratch area size, in bytes. */
size_t cmatinv2x2nf_getScratchSize ( int L );

/* Block Order, Floating-Point, 3x3, SA=12
   Restrictions: 
     None
*/
void cmatinv3x3nf ( void * restrict pScr, complex_float * restrict A, int L );

/* Return the scratch area size, in bytes. */
size_t cmatinv3x3nf_getScratchSize ( int L );

/* Block Order, Floating-Point, 4x4, SA=16
   Restrictions: 
     None
*/
void cmatinv4x4nf ( void * restrict pScr, complex_float * restrict A, int L );

/* Return the scratch area size, in bytes. */
size_t cmatinv4x4nf_getScratchSize ( int L );

/* Stream Order, Floating-Point, 2x2, SA=2
   Restrictions: 
     L must be a multiple of 4
*/
void cmatinv2x2sf ( void * restrict pScr, complex_float * restrict A, int L );

/* Return the scratch area size, in bytes. */
size_t cmatinv2x2sf_getScratchSize ( int L );

/* Stream Order, Floating-Point, 3x3, SA=9
   Restrictions: 
     L must be a multiple of 4
*/
void cmatinv3x3sf ( void * restrict pScr, complex_float * restrict A, int L );

/* Return the scratch area size, in bytes. */
size_t cmatinv3x3sf_getScratchSize ( int L );

/* Stream Order, Floating-Point, 4x4, SA=16
   Restrictions: 
     L must be a multiple of 4
*/
void cmatinv4x4sf ( void * restrict pScr, complex_float * restrict A, int L );

/* Return the scratch area size, in bytes. */
size_t cmatinv4x4sf_getScratchSize ( int L );

/*-------------------------------------------------------------------------
Direct Matrix Inversion For Real Matrices

Description: perform in-place inversion of 2x2 and 4x4 real matrices. 
2x2 matrices are inverted by Cramer's rule. For 3x3, 4x4 matrices we employ 
the blockwise inversion algorithm encompassed with a suboptimal row/column
permutation that gains better conditioning of the block structure.

Data format and order options:
  Suffix   Data Order                Data Format   
    nf       Block     IEEE-754 Std single precision floating-point
    sf       Stream    IEEE-754 Std single precision floating-point

Note:
In general, accuracy of a matrix inversion algorithm implementation is a
function of input matrix condition number. Thus it is user's responsibility
to qualify the reliability of numeric results. Refer to NatureDSP Baseband 
Library Reference for details. 

Parameters:
Temporary:
pScr      Scratch memory area. To determine the scratch area size required by
          a function <fun>, use the respective helper function 
          <fun>_getScratchSize()
Input:
L         Number of matrices
Input/Output:
A[L][SA]  Sequence of L NxN input/result matrices. SA is the number
          of data elements occupied by a single NxN matrix in a block 
          (stream) ordered sequence, see function specifications.

Restrictions:
pScr,A    Aligned on 32-byte boundary

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/
/* Block Order, Floating-Point, 2x2, SA=4
   Restrictions: 
     L must be even
*/
void matinv2x2nf ( void * restrict pScr, float32_t * restrict A, int L );

/* Return the scratch area size, in bytes. */
size_t matinv2x2nf_getScratchSize ( int L );

/* Block Order, Floating-Point, 3x3, SA=16
   Restrictions: 
     None
*/
void matinv3x3nf ( void * restrict pScr, float32_t * restrict A, int L );

/* Return the scratch area size, in bytes. */
size_t matinv3x3nf_getScratchSize ( int L );

/* Block Order, Floating-Point, 4x4, SA=16
   Restrictions: 
     None
*/
void matinv4x4nf ( void * restrict pScr, float32_t * restrict A, int L );

/* Return the scratch area size, in bytes. */
size_t matinv4x4nf_getScratchSize ( int L );

/* Stream Order, Floating-Point, 2x2, SA=4
   Restrictions: 
     L must be a multiple of 8
*/
void matinv2x2sf ( void * restrict pScr, float32_t * restrict A, int L );

/* Return the scratch area size, in bytes. */
size_t matinv2x2sf_getScratchSize ( int L );

/* Stream Order, Floating-Point, 3x3, SA=9
   Restrictions: 
     L must be a multiple of 8
*/
void matinv3x3sf ( void * restrict pScr, float32_t * restrict A, int L );

/* Return the scratch area size, in bytes. */
size_t matinv3x3sf_getScratchSize ( int L );

/* Stream Order, Floating-Point, 4x4, SA=16
   Restrictions: 
     L must be a multiple of 8
*/
void matinv4x4sf ( void * restrict pScr, float32_t * restrict A, int L );

/* Return the scratch area size, in bytes. */
size_t matinv4x4sf_getScratchSize ( int L );

/*-------------------------------------------------------------------------
Direct Matrix Inversion For Hermitian Matrices

Description: perform in-place inversion of 2x2 and 4x4 complex Hermitian 
matrices. 2x2 matrices are inverted by Cramer's rule. For 4x4 matrices we
employ the blockwise inversion algorithm encompassed with a suboptimal 
row/column permutation that gains better conditioning of the block structure.
Due to special properties of Hermitian matrices, the complexity of direct
matrix inversion is much lower than for a generic matrix of the same size.

Data format: block ordered 16-bit fixed-point data

NOTES:
1. In general, accuracy of a matrix inversion algorithm implementation is a
   function of input matrix condition number. Thus it is user's responsibility
   to qualify the reliability of numeric results. Refer to NatureDSP Baseband 
   Library Reference for details. 
2. For blockwise inversion of 4x4 fixed-point matrices, it is reasonable to
   limit the dynamic range of input data by 11..13 significant bits. This
   measure reduces the possibility of an overflow at internal computations.

Parameters:
Temporary:
pScr       Scratch memory area. To determine the scratch area size required by
           a function <fun>, use the respective helper function 
           <fun>_getScratchSize()
Input:
L          Number of matrices
qA         Number of fractional bits for input/output data
Input/Output:
A[L][N*N]  Sequence of L NxN complex input/result matrices

Restrictions:
pScr,A     Aligned on 32-byte boundary

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/
/* Block Order, Fixed-Point, 2x2
   Restrictions: 
     L must be even
*/
void chermmatinv2x2n ( void * restrict pScr, complex_fract16 * restrict A, int qA, int L );

/* Return the scratch area size, in bytes. */
size_t chermmatinv2x2n_getScratchSize ( int L );

/* Block Order, Fixed-Point, 4x4
   Restrictions:
     None
*/
void chermmatinv4x4n ( void * restrict pScr, complex_fract16 * restrict A, int qA, int L );

/* Return the scratch area size, in bytes. */
size_t chermmatinv4x4n_getScratchSize ( int L );

/*===========================================================================
Matrix Inversion By Gauss-Jordan algortihm
[c]matinvgjn   Gauss-Jordan algortihm for real/complex block ordered matrices
[c]matinvgjs   Gauss-Jordan algortihm for real/complex stream ordered matrices
===========================================================================*/

/*-------------------------------------------------------------------------
Inversion of Block Ordered Matrices By Gauss-Jordan Algortihm

Description: perform in-place inversion of real/complex matrices by Gauss-
Jordan elimination method, with full pivoting. The algorithm is applied to
a sequence of input matrices stored in block order.

Inversion result is not defined for a close to singular input matrix.

Storage size SA denote the number of data elements required to store a
matrix an NxN matrix A in block order. If matrix size is less than the
SIMD vector size, then the storage_size(matrix_size) equals the matrix_size
rounded up to the next power of two, otherwise it is matrix_size rounded up
to the next multiple of the SIMD vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(N*N)

Data format: IEEE-754 Std single precision floating-point

Temporary:
  pScr      Scratch area. Required size (in bytes) is defined by 
            functions [c]matinvgj<size>nf_getScratchSize(N,L)
Input:
  N         Matrix size
  L         Number of matrices
Input/Output:
  A[L][SA]  Input matrices, inverted matrices on output
Restrictions:
  pScr,A    Must not overlap and must be aligned on 32-byte boundary 
  N         Must be a positive multiple of 4
---------------------------------------------------------------------------*/
/* Real-valued functions */
void matinvgj8x8nf   ( void * pScr, float32_t * restrict A, int L );
void matinvgj16x16nf ( void * pScr, float32_t * restrict A, int L );
void matinvgjnxnnf   ( void * pScr, float32_t * restrict A, int N, int L );

/* Return the scratch area size, in bytes. */
size_t matinvgj8x8nf_getScratchSize   ( int N, int L );
size_t matinvgj16x16nf_getScratchSize ( int N, int L );
size_t matinvgjnxnnf_getScratchSize   ( int N, int L );

/* Complex-valued functions */
void cmatinvgj8x8nf   ( void * pScr, complex_float * restrict A, int L );
void cmatinvgj16x16nf ( void * pScr, complex_float * restrict A, int L );
void cmatinvgjnxnnf   ( void * pScr, complex_float * restrict A, int N, int L );

/* Return the scratch area size, in bytes. */
size_t cmatinvgj8x8nf_getScratchSize   ( int N, int L );
size_t cmatinvgj16x16nf_getScratchSize ( int N, int L );
size_t cmatinvgjnxnnf_getScratchSize   ( int N, int L );

/*-------------------------------------------------------------------------
Inversion of Stream Ordered Matrices By Gauss-Jordan Algortihm

Description: perform in-place inversion of real/complex matrices by Gauss-
Jordan elimination method, with full pivoting. The algorithm is applied to
a sequence of input matrices stored in stream order.

Inversion result is not defined for a close to singular input matrix.

Data format: IEEE-754 Std single precision floating-point

Temporary:
  pScr       Scratch area. Required size (in bytes) is defined by 
             functions [c]matinvgj<size>sf_getScratchSize(N,L)
Input:
  N          Matrix size
  L          Number of matrices
Input/Output:
  A[N*N][L]  Input matrices, inverted matrices on output
Restrictions:
  pScr,A     Must not overlap and must be aligned on 32-byte boundary 
  N          Must be greater than 1
  L          Must be a multiple of 8 for real-valued functions, or a multiple
             of 4 for complex-valued functions
---------------------------------------------------------------------------*/
/* Real-valued functions */
void matinvgj2x2sf ( void * pScr, float32_t * restrict A, int L );
void matinvgj3x3sf ( void * pScr, float32_t * restrict A, int L );
void matinvgj4x4sf ( void * pScr, float32_t * restrict A, int L );
void matinvgjnxnsf ( void * pScr, float32_t * restrict A, int N, int L );

/* Return the scratch area size, in bytes. */
size_t matinvgj2x2sf_getScratchSize ( int N, int L );
size_t matinvgj3x3sf_getScratchSize ( int N, int L );
size_t matinvgj4x4sf_getScratchSize ( int N, int L );
size_t matinvgjnxnsf_getScratchSize ( int N, int L );

/* Complex-valued functions */
void cmatinvgj2x2sf ( void * pScr, complex_float * restrict A, int L );
void cmatinvgj3x3sf ( void * pScr, complex_float * restrict A, int L );
void cmatinvgj4x4sf ( void * pScr, complex_float * restrict A, int L );
void cmatinvgjnxnsf ( void * pScr, complex_float * restrict A, int N, int L );

/* Return the scratch area size, in bytes. */
size_t cmatinvgj2x2sf_getScratchSize ( int N, int L );
size_t cmatinvgj3x3sf_getScratchSize ( int N, int L );
size_t cmatinvgj4x4sf_getScratchSize ( int N, int L );
size_t cmatinvgjnxnsf_getScratchSize ( int N, int L );

/*===========================================================================
LU Decomposition
[c]lun    LU decomposition for real/complex block ordered matrices
[c]lus    LU decomposition for real/complex stream ordered matrices
===========================================================================*/

/*-------------------------------------------------------------------------
LU Decomposition For Block Ordered Matrices

Description: compute LU decomposition of a square matrix using partial pivoting
with row interchanges: P*A = L*U, where P is the permutation matrix, L is the 
lower triangular factor (with diagnoal elements equal to 1), U is the upper 
triangular factor.

Algorithm is applied in-place to a sequence of real/complex matrices stored
in block order. For each input matrix A[k], the resutling factor U[k] replaces
the upper triangle of input matrix A[k], and subdiagonal elements of factor L[k]
replace the lower triangle of matrix A[k]. Ones on the main diagonal of resulting
factor L[k] are discarded.

Row ordering used by LU decomposition algorithm for input matrix A[k] is stored
to the vector of permutation indices P[k][0..N-1]: i-th row of the matrix was 
interchanged with row P[k][i], i=0..N-1.

Decomposition result is not defined for a close to singular input matrix.

Storage sizes SA,SP denote the number of data elements required to store a
matrix or a vector in block order. If matrix size is less than the SIMD vector
size, then the storage_size(matrix_size) equals the matrix_size rounded up to
the next power of two, otherwise it is matrix_size rounded up to the next
multiple of the SIMD vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(N*N)
SP = storage_size(N)

Data format: IEEE-754 Std single precision floating-point

Temporary:
  pScr      Scratch area. Required size (in bytes) is defined by 
            functions [c]lu<size>nf_getScratchSize(N,L)
Input:
  N         Matrix size
  L         Number of matrices
Input/Output:
  A[L][SA]  Input matrices, packed L and U factors on output
Output:
  P[L][SP]  Permutation index vectors
Restrictions:
  pScr,A,P  Must not overlap and must be aligned on 32-byte boundary 
  N         Must be a positive multiple of 4
---------------------------------------------------------------------------*/
/* Real-valued functions */
void lu8x8nf ( 
            void * pScr,
            float32_t * restrict A, 
            int16_t   * restrict P,
            int L );
void lu16x16nf ( 
            void * pScr,
            float32_t * restrict A,
            int16_t   * restrict P,
            int L );
void lunxnnf ( 
            void * pScr,
            float32_t * restrict A, 
            int16_t   * restrict P,
            int N, int L );

/* Return the scratch area size, in bytes. */
size_t lu8x8nf_getScratchSize   ( int N, int L );
size_t lu16x16nf_getScratchSize ( int N, int L );
size_t lunxnnf_getScratchSize   ( int N, int L );

/* Complex-valued functions */
void clu8x8nf ( 
            void * pScr,
            complex_float * restrict A, 
            int16_t       * restrict P,
            int L );
void clu16x16nf (
            void * pScr,
            complex_float * restrict A,
            int16_t       * restrict P,
            int L );
void clunxnnf ( 
            void * pScr,
            complex_float * restrict A, 
            int16_t       * restrict P,
            int N, int L );

/* Return the scratch area size, in bytes. */
size_t clu8x8nf_getScratchSize   ( int N, int L );
size_t clu16x16nf_getScratchSize ( int N, int L );
size_t clunxnnf_getScratchSize   ( int N, int L );

/*-------------------------------------------------------------------------
LU Decomposition For Stream Ordered Matrices

Description: compute LU decomposition of a square matrix using partial pivoting
with row interchanges: P*A = L*U, where P is the permutation matrix, L is the 
lower triangular factor (with diagnoal elements equal to 1), R is the upper 
triangular factor.

Algorithm is applied in-place to a sequence of real/complex matrices stored
in stream order. For each input matrix A[k], the resutling factor U[k] replaces
the upper triangle of input matrix A[k], and subdiagonal elements of factor L[k]
replace the lower triangle of matrix A[k]. Ones on the main diagonal of resulting
factor L[k] are discarded.

Row ordering used by LU decomposition algorithm for input matrix A[k] is stored
to the vector of permutation indices P[0..N-1][k]: i-th row of the matrix was 
interchanged with row P[i][k], i=0..N-1.

Decomposition result is not defined for a close to singular input matrix.

Data format: IEEE-754 Std single precision floating-point

Temporary:
  pScr       Scratch area. Required size (in bytes) is defined by 
             functions [c]lu<size>sf_getScratchSize(N,L)
Input:
  N          Matrix size
  L          Number of matrices
Input/Output:
  A[N*N][L]  Input matrices, packed L and U factors on output
Output:
  P[N][L]    Permutation index vectors
Restrictions:
  pScr,A,P   Must not overlap and must be aligned on 32-byte boundary 
  N          Must be greater than 1
  L          Must be a multiple of 8 for real-valued functions, or a multiple
             of 4 for complex-valued functions
---------------------------------------------------------------------------*/
/* Real-valued functions */
void lu2x2sf ( 
            void *pScr,
            float32_t * restrict A, 
            int16_t   * restrict P,
            int L );
void lu3x3sf ( 
            void *pScr,
            float32_t * restrict A,
            int16_t   * restrict P,
            int L );
void lu4x4sf ( 
            void *pScr,
            float32_t * restrict A,
            int16_t   * restrict P,
            int L );
void lunxnsf ( 
            void *pScr,
            float32_t * restrict A,
            int16_t   * restrict P,
            int N, int L );

void clu2x2sf ( 
            void *pScr,
            complex_float * restrict A, 
            int16_t       * restrict P,
            int L );
void clu3x3sf ( 
            void *pScr,
            complex_float * restrict A,
            int16_t       * restrict P,
            int L );
void clu4x4sf ( 
            void *pScr,
            complex_float * restrict A,
            int16_t       * restrict P,
            int L );
void clunxnsf ( 
            void *pScr,
            complex_float * restrict A,
            int16_t       * restrict P,
            int N, int L );

/* Return the scratch area size, in bytes. */
size_t lu2x2sf_getScratchSize ( int N, int L );
size_t lu3x3sf_getScratchSize ( int N, int L );
size_t lu4x4sf_getScratchSize ( int N, int L );
size_t lunxnsf_getScratchSize ( int N, int L );

/* Return the scratch area size, in bytes. */
size_t clu2x2sf_getScratchSize ( int N, int L );
size_t clu3x3sf_getScratchSize ( int N, int L );
size_t clu4x4sf_getScratchSize ( int N, int L );
size_t clunxnsf_getScratchSize ( int N, int L );

/*===========================================================================
Matrix determinant
[c]detn         Determinant for real/complex block ordered matrices
[c]dets         Determinant for real/complex stream ordered matrices
===========================================================================*/

/*-------------------------------------------------------------------------
Determinant For Block Ordered Matrices

Description: compute determinant of a real/complex matrix from its LU 
decomposition (see lu<size>nf() and clu<size>nf() functions) by multiplying
together diagonal elements of the upper triangular factor U. This operation
is accomplished for a sequence of LU matrices stored in block order.

Storage size SLU denotes the number of data elements required for each NxN LU 
matrix when stored in block order. If matrix size N*N is less than the SIMD
vector size, then the storage size equals N*N rounded up to the next power of
two, otherwise storage size is N*N rounded up to the next multiple of the SIMD
vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

Data format: IEEE-754 Std single precision floating-point

Input:
  N           Matrix size
  L           Number of matrices
Input/Output:
  LU[L][SLU]  Packed L and U factors computed by [c]lu<size>nf()
Output:
  D[L]        Determinant values
Restrictions:
  D,LU        Must not overlap and must be aligned on 32-byte boundary 
  N           Must be a positive multiple of 4
---------------------------------------------------------------------------*/
/* Real-valued functions */
void det8x8nf   ( float32_t * restrict D, const float32_t * restrict LU, int L );
void det16x16nf ( float32_t * restrict D, const float32_t * restrict LU, int L );
void detnxnnf   ( float32_t * restrict D, const float32_t * restrict LU, int N, int L );

/* Complex-valued functions */
void cdet8x8nf   ( complex_float * restrict D, const complex_float * restrict LU, int L );
void cdet16x16nf ( complex_float * restrict D, const complex_float * restrict LU, int L );
void cdetnxnnf   ( complex_float * restrict D, const complex_float * restrict LU, int N, int L );

/*-------------------------------------------------------------------------
Determinant For Stream Ordered Matrices

Description: compute determinant of a real/complex matrix from its LU 
decomposition (see lu<size>sf() and clu<size>sf() functions) by multiplying
together diagonal elements of the upper triangular factor U. This operation
is accomplished for a sequence of LU matrices stored in stream order.

Data format: IEEE-754 Std single precision floating-point

Input:
  N           Matrix size
  L           Number of matrices
Input/Output:
  LU[N*N][L]  Packed L and U factors computed by [c]lu<size>sf()
Output:
  D[L]        Determinant values
Restrictions:
  D,LU        Must not overlap and must be aligned on 32-byte boundary 
  N           Must be greater than 1
  L           Must be a multiple of 8 for real-valued functions, or a mutiple
              of 4 for complex-valued functions.
---------------------------------------------------------------------------*/
/* Real-valued functions */
void det2x2sf ( float32_t * restrict D, const float32_t * restrict LU, int L );
void det3x3sf ( float32_t * restrict D, const float32_t * restrict LU, int L );
void det4x4sf ( float32_t * restrict D, const float32_t * restrict LU, int L );
void detnxnsf ( float32_t * restrict D, const float32_t * restrict LU, int N, int L );

/* Complex-valued functions */
void cdet2x2sf ( complex_float * restrict D, const complex_float * restrict LU, int L );
void cdet3x3sf ( complex_float * restrict D, const complex_float * restrict LU, int L );
void cdet4x4sf ( complex_float * restrict D, const complex_float * restrict LU, int L );
void cdetnxnsf ( complex_float * restrict D, const complex_float * restrict LU, int N, int L );

/*===========================================================================
Cholesky Decomposition
===========================================================================*/
/*===========================================================================
  Cholesky decomposition for a complex-valued fixed-point pseudo-inversion:
  chol2x2s  2x2 complex matrices
  chol4x4s  4x4 complex matrices
  chol8x8s  8x8 complex matrices
  cholmxns  MxN complex matrices
===========================================================================*/

/*-------------------------------------------------------------------------
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex-valued least squares problem: A*X=B, where A is
an MxN coefficient matrix with M >= N; X is an NxP matrix of unknowns; and
B is an MxP right hand matrix.

The decomposition results in an upper triangular complex NxN matrix R with
real and positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A + sigma2*I,
where adj(...) denotes the conjugate transpose of a matrix, and sigma2*I is
the NxN identity matrix multiplied with the regularization term.

The decomposition algorithm is applied to a few matrices per single call,
with input/output matrix sequences being stored in the streaming order.

Fixed-point data type of upper triangular matrices R is the same as the
data type of input matrices A. Fixed point position for the regularization
term sigma2 must match the scale of product adj(A)*A. If, for instance,
matrix A is represented as Q15, then Q30 is expected for sigma2.

In order to find the solution to the above-mentioned least squares problem
A*X=B, one has to follow the Cholesky decomposition with forward and backward
substitution procedures; see cholfwdmxnxps() and cholbkwnxps(), respectively.

Input:
  M, N           Dimensional parameters
  L              Number of matrices
  sigma2[L]      Regularization term; fixed point position is twice the
                 number of fractional bits for matrices A, R
  A[M*N][L]      sequence of L complex matrices A
Output:
  R[N*N][L]      Sequence of L upper triangular complex matrices R
  D[L/4][N][8]   Reciprocal of main diagonal (mantissa, exponent) in the 
                 special format
Restrictions:
  1. A, R, D, sigma2 must not overlap
  2. A, R, D, sigma2 must be aligned on 32-byte boundary
  3. Number of matrices L must be a multiple of 8
  4. Matrix sizes must be greater than 1
  5. Number of columns for input matrices A must not exceed the number
     of rows: N <= M.
---------------------------------------------------------------------------*/
void chol2x2s(
            complex_fract16 * restrict R,
            complex_fract16 * restrict D,
            const complex_fract16 * restrict A, 
            const int32_t * restrict sigma2,
            int L);
void chol4x4s(
            complex_fract16 * restrict R, 
            complex_fract16 * restrict D,
            const complex_fract16 * restrict A, 
            const int32_t * restrict sigma2,
            int L);
void chol8x8s(
            complex_fract16 * restrict R, 
            complex_fract16 * restrict D,
            const complex_fract16 * restrict A, 
            const int32_t * restrict sigma2,
            int L);
void cholmxns(
            complex_fract16 * restrict R, 
            complex_fract16 * restrict D,
            const complex_fract16 * restrict A, 
            const int32_t * restrict sigma2,
            int M,int N,
            int L);

/*===========================================================================
  Cholesky forward recursion for pseudo-inversion API (complex fixed-point data)
  cholfwd2x2x1s   matrices 2x2x1
  cholfwd4x4x1s   matrices 4x4x1
  cholfwd8x8x1s   matrices 8x8x1
  cholfwdmxnxps   matrices MxNxP
===========================================================================*/

/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Input:
M               Matrix dimension (number of rows in matrices A)
N               Matrix dimension (number of columns and rows in matrices 
                R)
P               Number of columns in right-side matrices B
L               Number of matrices
R[N*N][L]       Cholesky upper triangular matrices R
A[M*N][L]       Original left-side matrices A
B[M*P][L]       Original right-side matrices B
D[L/4][N][8]    Reciprocal of main diagonal (mantissa, exponent) in the 
                special format
qA,qB,qY        Fixed point representation of matrices A (or R which is 
                the same), B and y
Output:
y[N*P][L]       Decision matrix y

Restrictions:
1. All matrices must not overlap and be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8
3. Matrix sizes (M,N) must be greater than 1, P must be >=1
4. M >= N
---------------------------------------------------------------------------*/
void cholfwd2x2x1s(
            complex_fract16* restrict y,
            const complex_fract16* restrict R,
            const complex_fract16* restrict D,
            const complex_fract16* restrict A, 
            const complex_fract16* restrict B, 
            int qA,int qB,int qY,
            int L);
void cholfwd4x4x1s(
            complex_fract16* restrict y,
            const complex_fract16* restrict R, 
            const complex_fract16* restrict D,
            const complex_fract16* restrict A, 
            const complex_fract16* restrict B, 
            int qA,int qB,int qY,
            int L);
void cholfwd8x8x1s(
            complex_fract16* restrict y,
            const complex_fract16* restrict R, 
            const complex_fract16* restrict D,
            const complex_fract16* restrict A, 
            const complex_fract16* restrict B, 
            int qA,int qB,int qY,
            int L);
void cholfwdmxnxps(
            complex_fract16* restrict y,
            const complex_fract16* restrict R, 
            const complex_fract16* restrict D,
            const complex_fract16* restrict A, 
            const complex_fract16* restrict B, 
            int qA,int qB,int qY,
            int M,int N, int P,
            int L);


/*===========================================================================
  Cholesky backward recursion for pseudo-inversion API (complex fixed-point data)
  cholbkw2x1s vectors 2x1
  cholbkw4x1s vectors 4x1
  cholbkw8x1s vectors 8x1
  cholbkwnxps vectors NxP
===========================================================================*/

/*-------------------------------------------------------------------------
   These functions make backward recursion stage of pseudo-inversion. They
   use Cholesky decomposition of original matrices and results of forward 
   recursion. 
   NOTE:
   Data layout for matrices is selected as for other matrices written in a 
   streaming order. 

   Input:
   N                Matrix dimension (number of columns and rows in 
                    matrices R)
   P                Number of columns in right-side matrices B
   L                Number of matrices
   R[N*N][L]        Cholesky upper triangular matrices R
   D[L/4][N][8]     reciprocal of main diagonal (mantissa, exponent) in the 
                    special format
   y[N*P][L]        Results of forward recursion stage
   qA,qX,qY         fixed point representation of matrices A(or R which is 
                    the same), x and y
   Output:
   x[N*P][L]        Decision matrix x


   Restrictions:
   1. All matrices must not overlap and be aligned on 32-byte boundary 
   2. Number of matrices L must be a multiple of 8
   3. Matrix sizes (M,N) must be greater than 1, P must be >=1
   4. qX+qA-qY must be <=16 
---------------------------------------------------------------------------*/
void cholbkw2x1s(
            complex_fract16* restrict x, 
            const complex_fract16* restrict R,
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int L);
void cholbkw4x1s(
            complex_fract16* restrict x, 
            const complex_fract16* restrict R,
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int L);
void cholbkw8x1s(
            complex_fract16* restrict x, 
            const complex_fract16* restrict R,
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int L);
void cholbkwnxps(
            complex_fract16* restrict x, 
            const complex_fract16* restrict R,
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int N, int P, int L);

/*===========================================================================
  Cholesky decomposition for block ordered fixed-point matrices:
  cholmxnn   MxN matrices
  chol8x8n   8x8
  chol16x16n 16x16
  chol32x32n 32x32
===========================================================================*/

/*-------------------------------------------------------------------------
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex-valued least squares problem: A*X=B, where A is
an MxN coefficient matrix with M >= N; X is an NxP matrix of unknowns; and
B is an MxP right hand matrix.

The decomposition results in an upper triangular complex NxN matrix R with
real and positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A 
where adj(...) denotes the conjugate transpose of a matrix

The decomposition algorithm is applied to a few matrices per single call,
with input/output matrix sequences being stored in the block order.

Fixed-point data type of upper triangular matrices R is the same as the
data type of input matrices A. 
Matrix R is stored in special format: only upper-diagonal elements are 
stored and they are written column by column. So, total number of elements
in one matrix R is the sum of arithmetic progression 1,2...N == ((N+1)*N)/2

In order to find the solution to the above-mentioned least squares problem
A*X=B, one has to follow the Cholesky decomposition with forward and backward
substitution procedures; see cholfwdmxnxpn() and cholbkwnxpn(), respectively.

Matrix sizes SA,SR,SD are selected as usual for complex block ordered matrix 
sequencies, i.e. total size is rounded up to the closest bigger multiple of 
BBE_SIMD_WIDTH/2==8 elements. 
SA=size(M*N)
SR=size(((N+1)*N)/2)
SD=size(N)
Scratch size in bytes is defined by scratch allocation functions

Input:
 M, N         Dimensional parameters
 L            Number of matrices
 A[L][SA]     Sequence of L complex matrices A
Output:
 R[L][SR]     Sequence of L upper triangular complex matrices R
 D[L][SD]     Reciprocal of main diagonal (mantissa, exponent) 
              in the special format

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. Number of matrices L must be positive
3. M and N must be positive multiples of 4
4. Number of columns for input matrices A must not exceed the number
   of rows: N<=M.
---------------------------------------------------------------------------*/
void cholmxnn(
            void *pScr,
            complex_fract16 * restrict R, 
            complex_fract16 * restrict D,
            const complex_fract16 * restrict A, 
            int M,int N,
            int L);
void chol8x8n(
            void *pScr,
            complex_fract16 * restrict R, 
            complex_fract16 * restrict D,
            const complex_fract16 * restrict A, 
            int L);
void chol16x16n(
            void *pScr,
            complex_fract16 * restrict R, 
            complex_fract16 * restrict D,
            const complex_fract16 * restrict A, 
            int L);
void chol32x32n(
            void *pScr,
            complex_fract16 * restrict R, 
            complex_fract16 * restrict D,
            const complex_fract16 * restrict A, 
            int L);
/* scratch allocation functions */
size_t cholmxnn_getScratchSize  (int M,int N, int L);
size_t chol8x8n_getScratchSize  (int M,int N, int L);
size_t chol16x16n_getScratchSize(int M,int N, int L);
size_t chol32x32n_getScratchSize(int M,int N, int L);

/*===========================================================================
  Cholesky forward recursion for block ordered fixed-point matrices:
  cholfwdmxnxpn     MxNxP
  cholfwd8x8x1n     8x8x1
  cholfwd16x16x1n   16x16x1
  cholfwd32x32x1n   32x32x1
===========================================================================*/

/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. 

Matrix sizes SA,SR,SD,SB,SY are selected as usual for complex block ordered 
matrix sequencies, i.e. total size is rounded up to the closest bigger 
multiple of BBE_SIMD_WIDTH/2==8 elements or, if it is lesser, to the 
closest bigger multiple of degree of 2. 
SA=size(M*N)
SR=size(((N+1)*N)/2)
SD=size(N)
SB=size(M*P)
SY=size(N*P)
Scratch size in bytes is defined by cholfwdxxx_getScratchSize()

Input:
 M            Matrix dimension (number of rows in matrices A)
 N            Matrix dimension (number of columns and rows in 
              matrices R)
 P            Number of columns in right-side matrices B
 L            Number of matrices
 R[L][SR]     Sequence of L upper triangular complex matrices R
 A[L][SA]     Sequence of L complex matrices A
 B[L][SB]     Sequence original right-side matrices B
 D[L][SD]     Reciprocal of main diagonal (mantissa, exponent) 
              in the special format
qB,qY         Fixed point representation of matrices B and y
Output:
y[L][SY]      Sequence of intermediate decision matrices y

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
3. Matrix sizes M,N,P must be positive
4. M and N must be multiples of 4 
5. M>=N
---------------------------------------------------------------------------*/
void cholfwdmxnxpn(
            void *pScr,
            complex_fract16* restrict y,
            const complex_fract16* restrict R, 
            const complex_fract16* restrict D,
            const complex_fract16* restrict A, 
            const complex_fract16* restrict B, 
            int qB,int qY,
            int M,int N, int P,
            int L);
void cholfwd8x8x1n(
            void *pScr,
            complex_fract16* restrict y,
            const complex_fract16* restrict R, 
            const complex_fract16* restrict D,
            const complex_fract16* restrict A, 
            const complex_fract16* restrict B, 
            int qB,int qY,
            int L);
void cholfwd16x16x1n(
            void *pScr,
            complex_fract16* restrict y,
            const complex_fract16* restrict R, 
            const complex_fract16* restrict D,
            const complex_fract16* restrict A, 
            const complex_fract16* restrict B, 
            int qB,int qY,
            int L);
void cholfwd32x32x1n(
            void *pScr,
            complex_fract16* restrict y,
            const complex_fract16* restrict R, 
            const complex_fract16* restrict D,
            const complex_fract16* restrict A, 
            const complex_fract16* restrict B, 
            int qB,int qY,
            int L);
/* scratch allocation functions */
size_t cholfwdmxnxpn_getScratchSize(int M,int N, int P,int L);
size_t cholfwd8x8x1n_getScratchSize(int M,int N, int P,int L);
size_t cholfwd16x16x1n_getScratchSize(int M,int N, int P,int L);
size_t cholfwd32x32x1n_getScratchSize(int M,int N, int P,int L);

/*===========================================================================
  Cholesky backward recursion for block ordered fixed-point matrices:
  cholbkwnxps nxp
===========================================================================*/

/*-------------------------------------------------------------------------
These functions make backward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices and results of forward recursion. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. 

Matrix sizes SR,SD,SY,SX are selected as usual for complex block ordered 
matrix sequencies, i.e. total size is rounded up to the closest bigger 
multiple of BBE_SIMD_WIDTH/2==8 elements or, if it is lesser, to the 
closest bigger multiple of degree of 2. 
SR=size(((N+1)*N)/2)
SD=size(N)
SY=size(N*P)
SX=size(N*P)
Scratch size in bytes is defined by cholbkwxxx_getScratchSize()

Input:
 N            Matrix dimension (number of columns and rows in 
              matrices R)
 P            Number of columns in right-side matrices B
 L            Number of matrices
 R[L][SR]     Sequence of L upper triangular complex matrices R
 D[L][SD]     Reciprocal of main diagonal (mantissa, exponent) 
              in the special format
 y[L][SY]     Sequence of intermediate decision matrices y
 qA,qX,qY     Fixed point representation of matrices A(or R which 
              is the same), x and y
Output:         
x[L][SX]      sequence of decision matrix x

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. Number of matrices L must be positive
3. Matrix sizes M,N,P must be positive
4. M and N must be multiples of 4
5. qX+qA-qY must be <=16 
---------------------------------------------------------------------------*/
void cholbkwnxpn(
            void *pScr,
            complex_fract16* restrict x, 
            const complex_fract16* restrict R,
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int N,int P, int L);
void cholbkw8x1n(
            void *pScr,
            complex_fract16* restrict x, 
            const complex_fract16* restrict R,
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int L);
void cholbkw16x1n(
            void *pScr,
            complex_fract16* restrict x, 
            const complex_fract16* restrict R,
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int L);
void cholbkw32x1n(
            void *pScr,
            complex_fract16* restrict x, 
            const complex_fract16* restrict R,
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int L);
/* scratch allocation functions */
size_t cholbkwnxpn_getScratchSize(int N,int P,int L);
size_t cholbkw8x1n_getScratchSize(int N,int P,int L);
size_t cholbkw16x1n_getScratchSize(int N,int P,int L);
size_t cholbkw32x1n_getScratchSize(int N,int P,int L);

/*===========================================================================
  Cholesky decomposition for real/complex stream ordered floating-point matrices
  [r]chol2x2sf  2x2 real/complex matrices
  [r]chol3x3sf  3x3 real/complex matrices
  [r]chol4x4sf  4x4 real/complex matrices
  [r]cholmxnsf  MxN real/complex matrices
===========================================================================*/

/*-------------------------------------------------------------------------
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex- or real-valued least squares problem: A*X=B, 
where A is an MxN coefficient matrix with M >= N; X is an NxP matrix of
unknowns; and B is an MxP right hand matrix.

The decomposition results in an upper triangular NxN matrix R with real and
positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A + sigma2*I,
where adj(...) denotes the (conjugate) transpose of a matrix, and sigma2*I is
an NxN identity matrix multiplied by the regularization term.

The decomposition algorithm is applied to a few matrices per single call,
with input/output matrix sequences being stored in the stream order.

In order to find the solution to the above-mentioned least squares problem
A*X=B, one has to follow the Cholesky decomposition with forward and backward
substitution procedures; see [r]cholfwdmxnxpsf() and [r]cholbkwnxpsf(), 
respectively.

Data format: IEEE-754 Std. single precision floating-point

Input:
  M, N           Dimensional parameters
  L              Number of matrices
  sigma2[L]      Regularization term
  A[M*N][L]      Sequence of L real/complex matrices A
Output:
  R[N*N][L]      Sequence of L upper triangular real/complex matrices R
  D[L*N]         Reciprocal of main diagonal written in a special format
Restrictions:
  1. A, R, D, sigma2 must not overlap
  2. A, R, D, sigma2 must be aligned on 32-byte boundary
  3. Number of matrices L must be a multiple of 4 for complex-valued 
     functions, or a multiple of 8 for real-valued functions.
  4. Matrix sizes must be greater than 1
  5. Number of columns for input matrices A must not exceed the number
     of rows: N <= M.
---------------------------------------------------------------------------*/
/* Complex-valued functions */
void chol2x2sf (
            complex_float * restrict R,
            complex_float * restrict D,
      const complex_float * restrict A, 
      const float32_t     * restrict sigma2,
      int L );
void chol3x3sf (
            complex_float * restrict R,
            complex_float * restrict D,
      const complex_float * restrict A, 
      const float32_t     * restrict sigma2,
      int L );
void chol4x4sf (
            complex_float * restrict R,
            complex_float * restrict D,
      const complex_float * restrict A, 
      const float32_t     * restrict sigma2,
      int L );
void cholmxnsf (
            complex_float * restrict R,
            complex_float * restrict D,
      const complex_float * restrict A, 
      const float32_t     * restrict sigma2,
      int M, int N, int L );

/* Real-valued functions */
void rchol2x2sf (
            float32_t * restrict R,
            float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict sigma2,
      int L );
void rchol3x3sf (
            float32_t * restrict R,
            float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict sigma2,
      int L );
void rchol4x4sf (
            float32_t * restrict R,
            float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict sigma2,
      int L );
void rcholmxnsf (
            float32_t * restrict R,
            float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict sigma2,
      int M, int N, int L );

/*===========================================================================
  Cholesky forward recursion for real/complex stream ordered floating-point
  matrices
  [r]cholfwd2x2x1sf   2x2x1 real/complex least squares problem
  [r]cholfwd3x3x1sf   3x3x1 real/complex least squares problem
  [r]cholfwd4x4x1sf   4x4x1 real/complex least squares problem
  [r]cholfwdmxnxpsf   MxNxP real/complex least squares problem
===========================================================================*/

/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Data format: IEEE-754 Std. single precision floating-point

Input:
M               Matrix dimension (number of rows in matrices A)
N               Matrix dimension (number of columns and rows in matrices 
                R)
P               Number of columns in right-side matrices B
L               Number of matrices
R[N*N][L]       Cholesky upper triangular matrices R
A[M*N][L]       Original left-side matrices A
B[M*P][L]       Original right-side matrices B
D[L*N]          Reciprocal of main diagonal written in a special format
Output:
y[N*P][L]       Decision matrix y

Restrictions:
1. All matrices must not overlap and be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 4 for complex-valued 
   functions, or a multiple of 8 for real-valued functions.
3. Matrix sizes (M,N) must be greater than 1, P must be >=1
4. M >= N
---------------------------------------------------------------------------*/
/* Complex-valued functions */
void cholfwd2x2x1sf (
            complex_float * restrict y,
      const complex_float * restrict R,
      const complex_float * restrict D,
      const complex_float * restrict A, 
      const complex_float * restrict B, 
      int L );
void cholfwd3x3x1sf (
            complex_float * restrict y,
      const complex_float * restrict R,
      const complex_float * restrict D,
      const complex_float * restrict A, 
      const complex_float * restrict B, 
      int L );
void cholfwd4x4x1sf (
            complex_float * restrict y,
      const complex_float * restrict R,
      const complex_float * restrict D,
      const complex_float * restrict A, 
      const complex_float * restrict B, 
      int L );
void cholfwdmxnxpsf (
            complex_float * restrict y,
      const complex_float * restrict R,
      const complex_float * restrict D,
      const complex_float * restrict A, 
      const complex_float * restrict B, 
      int M, int N, int P, int L );

/* Real-valued functions */
void rcholfwd2x2x1sf (
            float32_t * restrict y,
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict B, 
      int L );
void rcholfwd3x3x1sf (
            float32_t * restrict y,
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict B, 
      int L );
void rcholfwd4x4x1sf (
            float32_t * restrict y,
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict B, 
      int L );
void rcholfwdmxnxpsf (
            float32_t * restrict y,
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict B, 
      int M, int N, int P, int L );

/*===========================================================================
  Cholesky backward recursion for real/complex stream ordered floating-point
  matrices
  [r]cholbkw2x1sf   2x1 real/complex solution vectors
  [r]cholbkw3x1sf   3x1 real/complex solution vectors
  [r]cholbkw4x1sf   4x1 real/complex solution vectors
  [r]cholbkwnxpsf   NxP real/complex solution vectors
===========================================================================*/

/*-------------------------------------------------------------------------
These functions make backward recursion stage of pseudo-inversion. They
use Cholesky decomposition of original matrices and results of forward 
recursion. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Data format: IEEE-754 Std. single precision floating-point

Input:
N                Matrix dimension (number of columns and rows in 
                 matrices R)
P                Number of columns in right-side matrices B
L                Number of matrices
R[N*N][L]        Cholesky upper triangular matrices R
D[L*N]           Reciprocal of main diagonal written in a special format
y[N*P][L]        Results of forward recursion stage
Output:
x[N*P][L]        Decision matrix x

Restrictions:
1. All matrices must not overlap and be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 4 for complex-valued 
   functions, or a multiple of 8 for real-valued functions.
3. N must be greater than 1, P must be >=1
---------------------------------------------------------------------------*/
/* Complex-valued functions */
void cholbkw2x1sf (
            complex_float * restrict x, 
      const complex_float * restrict R,
      const complex_float * restrict D,
      const complex_float * restrict y, 
      int L );
void cholbkw3x1sf (
            complex_float * restrict x, 
      const complex_float * restrict R,
      const complex_float * restrict D,
      const complex_float * restrict y, 
      int L );
void cholbkw4x1sf (
            complex_float * restrict x, 
      const complex_float * restrict R,
      const complex_float * restrict D,
      const complex_float * restrict y, 
      int L );
void cholbkwnxpsf (
            complex_float * restrict x, 
      const complex_float * restrict R,
      const complex_float * restrict D,
      const complex_float * restrict y, 
      int N, int P, int L );

/* Real-valued functions */
void rcholbkw2x1sf (
            float32_t * restrict x, 
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict y, 
      int L );
void rcholbkw3x1sf (
            float32_t * restrict x, 
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict y, 
      int L );
void rcholbkw4x1sf (
            float32_t * restrict x, 
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict y, 
      int L );
void rcholbkwnxpsf (
            float32_t * restrict x, 
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict y, 
      int N, int P, int L );

/*===========================================================================
  Cholesky MMSE solver for real/complex stream ordered floating-point data
  [r]cholmmse2x2x1sf   2x2x1 real/complex MMSE problem
  [r]cholmmse3x3x1sf   3x3x1 real/complex MMSE problem
  [r]cholmmse4x4x1sf   4x4x1 real/complex MMSE problem
  [r]cholmmsemxnxpsf   MxNxP real/complex MMSE problem
===========================================================================*/

/*-------------------------------------------------------------------------
Compute the MMSE solution for a system of linear equations A*x=b, where A
is an MxN real (complex) matrix with M>=N and rank(A)==N, x is an Nx1 vector
of unknowns, and b is an Mx1 right hand side vector. This task is accomplished
in 3 steps:
  1. Cholesky decomposition is applied to the matrix of normal equations
     system, which results in an upper triangular NxN matrix R with real
     and positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A + sigma2*I,
     where adj(...) denotes the (conjugate) transpose of a matrix, and
     sigma2*I is an NxN identity matrix multiplied by the regularization
     term.
  2. Forward recursion step: the system adj(R)*y=adj(A)*b is solved for an 
     Nx1 vector y.
  3. Backward recursion step: the system R*x=y is solved for the Nx1 vector
     of unknows x.

For a single MxN matrix A, these 3 steps may be done simultaneously for P
variants of Mx1 right hand side column vectors b gathered into an MxP input
matrix B. MMSE solution is computed independently for each of P columns,
with resulting column vectors forming the solution matrix X of size NxP.

[r]cholmmse<size>sf() functions process L pairs of MxN matrices A and MxP 
matrices B in a single call, which results in L solution matrices X of
size NxP. L matrices of each kind are stored as stream ordered sequences.

Data format: IEEE-754 Std. single precision floating-point

Temporary:
  pScr       Scratch area. Required size (in bytes) is defined by 
             functions [r]cholmmse<size>sf_getScratchSize(M,N,P,L)
Input:
  M,N,P      Dimensional parameters
  L          Number of matrices
  sigma2[L]  Regularization term
  A[M*N][L]  Sequence of L matrices A
  B[M*P][L]  Sequence of L right hand side matrices B
Output:
  x[N*P][L]  Sequence of L solution matrices X
Restrictions:
  1. pScr,x,A,B,sigma2 must not overlap
  2. pScr,x,A,B,sigma2 must be aligned on 32-byte boundary
  3. Number of matrices L must be a multiple of 4 for complex-valued 
     functions, or a multiple of 8 for real-valued functions.
  4. Matrix sizes M,N must be greater than 1
  5. Number of columns for input matrices A must not exceed the number
     of rows: N <= M.
---------------------------------------------------------------------------*/
/* Complex-valued functions */
void cholmmse2x2x1sf (
            void * pScr,
            complex_float * restrict x,
      const complex_float * restrict A,
      const complex_float * restrict B,
      const float32_t     * restrict sigma2,
      int L );
void cholmmse3x3x1sf (
            void * pScr,
            complex_float * restrict x,
      const complex_float * restrict A,
      const complex_float * restrict B,
      const float32_t     * restrict sigma2,
      int L );
void cholmmse4x4x1sf (
            void * pScr,
            complex_float * restrict x,
      const complex_float * restrict A,
      const complex_float * restrict B,
      const float32_t     * restrict sigma2,
      int L );
void cholmmsemxnxpsf (
            void * pScr,
            complex_float * restrict x,
      const complex_float * restrict A,
      const complex_float * restrict B,
      const float32_t     * restrict sigma2,
      int M, int N, int P, int L );

void rcholmmse2x2x1sf (
            void * pScr,
            float32_t * restrict x,
      const float32_t * restrict A,
      const float32_t * restrict B,
      const float32_t * restrict sigma2,
      int L );
void rcholmmse3x3x1sf (
            void * pScr,
            float32_t * restrict x,
      const float32_t * restrict A,
      const float32_t * restrict B,
      const float32_t * restrict sigma2,
      int L );
void rcholmmse4x4x1sf (
            void * pScr,
            float32_t * restrict x,
      const float32_t * restrict A,
      const float32_t * restrict B,
      const float32_t * restrict sigma2,
      int L );
void rcholmmsemxnxpsf (
            void * pScr,
            float32_t * restrict x,
      const float32_t * restrict A,
      const float32_t * restrict B,
      const float32_t * restrict sigma2,
      int M, int N, int P, int L );

/* Return the scratch area size, in bytes. */
size_t cholmmse2x2x1sf_getScratchSize ( int M, int N, int P, int L );
size_t cholmmse3x3x1sf_getScratchSize ( int M, int N, int P, int L );
size_t cholmmse4x4x1sf_getScratchSize ( int M, int N, int P, int L );
size_t cholmmsemxnxpsf_getScratchSize ( int M, int N, int P, int L );
/* Return the scratch area size, in bytes. */
size_t rcholmmse2x2x1sf_getScratchSize ( int M, int N, int P, int L );
size_t rcholmmse3x3x1sf_getScratchSize ( int M, int N, int P, int L );
size_t rcholmmse4x4x1sf_getScratchSize ( int M, int N, int P, int L );
size_t rcholmmsemxnxpsf_getScratchSize ( int M, int N, int P, int L );

/*===========================================================================
  Cholesky decomposition for real/complex block ordered floating-point matrices
  [r]chol8x8nf    8x8 real/complex matrices
  [r]chol16x16nf  16x16 real/complex matrices
  [r]cholmxnnf    MxN real/complex matrices
===========================================================================*/

/*-------------------------------------------------------------------------
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex- or real-valued least squares problem: A*X=B, 
where A is an MxN coefficient matrix with M >= N; X is an NxP matrix of
unknowns; and B is an MxP right hand matrix.

The decomposition results in an upper triangular NxN matrix R with real and
positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A + sigma2*I,
where adj(...) denotes the (conjugate) transpose of a matrix, and sigma2*I is
an NxN identity matrix multiplied by the regularization term.

The decomposition algorithm is applied to a few matrices per single call,
with input/output matrix sequences being stored in the block order.

Matrix R is stored in special format: only upper-diagonal elements are 
stored and they are written column by column. So, total number of elements
in one matrix R is the sum of arithmetic progression 1,2...N == ((N+1)*N)/2

In order to find the solution to the above-mentioned least squares problem
A*X=B, one has to follow the Cholesky decomposition with forward and backward
substitution procedures; see [r]cholfwdmxnxpsf() and [r]cholbkwnxpsf(), 
respectively.

Storage sizes SA,SR,SD denote the number of data elements required to store a
matrix in block order. If matrix size is less than the SIMD vector size, then
the storage_size(matrix_size) equals the matrix_size rounded up to the next
power of two, otherwise it is matrix_size rounded up to the next multiple of
the SIMD vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(M*N)
SR = storage_size((N+1)*N/2)
SD = storage_size(N)

Scratch size in bytes is defined by scratch allocation functions

Data format: IEEE-754 Std. single precision floating-point

Input:
 M, N      Dimensional parameters
 L         Number of matrices
 A[L][SA]  Sequence of L complex matrices A
 sigma2[L] regularization term
Output:
 R[L][SR]  Sequence of L upper triangular complex matrices R
 D[L][SD]  Reciprocal of main diagonal

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. M and N must be positive multiples of 4
3. Number of columns for input matrices A must not exceed the number
   of rows: N<=M.
---------------------------------------------------------------------------*/

/* Complex-valued functions */
void chol8x8nf (
            void * pScr,
            complex_float * restrict R, 
            complex_float * restrict D,
      const complex_float * restrict A, 
      const float32_t     * restrict sigma2,
      int L );
void chol16x16nf (
            void * pScr,
            complex_float * restrict R, 
            complex_float * restrict D,
      const complex_float * restrict A, 
      const float32_t     * restrict sigma2,
      int L );
void cholmxnnf (
            void * pScr,
            complex_float * restrict R, 
            complex_float * restrict D,
      const complex_float * restrict A, 
      const float32_t     * restrict sigma2,
      int M, int N, int L );

void rchol8x8nf (
            void * pScr,
            float32_t * restrict R, 
            float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict sigma2,
      int L );
void rchol16x16nf (
            void * pScr,
            float32_t * restrict R, 
            float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict sigma2,
      int L );
void rcholmxnnf (
            void * pScr,
            float32_t * restrict R, 
            float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict sigma2,
      int M, int N, int L );

/* Scratch area size, in bytes. */
size_t cholmxnnf_getScratchSize  (int M,int N, int L);
size_t chol8x8nf_getScratchSize  (int M,int N, int L);
size_t chol16x16nf_getScratchSize(int M,int N, int L);
/* Scratch area size, in bytes. */
size_t rcholmxnnf_getScratchSize  (int M,int N, int L);
size_t rchol8x8nf_getScratchSize  (int M,int N, int L);
size_t rchol16x16nf_getScratchSize(int M,int N, int L);

/*===========================================================================
  Cholesky forward recursion for real/complex block ordered floating-point
  matrices
  [r]cholfwd8x8x1nf    8x8x1 real/complex least squares problem
  [r]cholfwd16x16x1nf  16x16x1 real/complex least squares problem
  [r]cholfwdmxnxpnf    MxNxP real/complex least squares problem
===========================================================================*/

/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. 

Storage sizes SA,SR,SD,SB,SY denote the number of data elements required to store
a matrix in block order. If matrix size is less than the SIMD vector size, then
the storage_size(matrix_size) equals the matrix_size rounded up to the next power
of two, otherwise it is matrix_size rounded up to the next multiple of the SIMD 
vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(M*N)
SR = storage_size((N+1)*N/2)
SD = storage_size(N)
SB = storage_size(M*P)
SY = storage_size(N*P)

Scratch size in bytes is defined by [r]cholfwd<...>nf_getScratchSize()

Data format: IEEE-754 Std. single precision floating-point

Input:
 M         Matrix dimension (number of rows in matrices A)
 N         Matrix dimension (number of columns and rows in 
           matrices R)
 P         Number of columns in right-side matrices B
 L         Number of matrices
 R[L][SR]  Sequence of L upper triangular complex matrices R
 A[L][SA]  Sequence of L complex matrices A
 B[L][SB]  Sequence of original right-side matrices B
 D[L][SD]  Reciprocal of main diagonal
Output:
 y[L][SY]   Sequence of intermediate decision matrices y

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. Matrix sizes M,N,P must be positive
3. M and N must be multiples of 4 
4. M>=N
---------------------------------------------------------------------------*/
/* Complex-valued functions */
void cholfwd8x8x1nf (
            void * pScr,
            complex_float * restrict y,
      const complex_float * restrict R, 
      const complex_float * restrict D,
      const complex_float * restrict A, 
      const complex_float * restrict B, 
      int L );
void cholfwd16x16x1nf (
            void * pScr,
            complex_float * restrict y,
      const complex_float * restrict R, 
      const complex_float * restrict D,
      const complex_float * restrict A, 
      const complex_float * restrict B, 
      int L );
void cholfwdmxnxpnf (
            void * pScr,
            complex_float * restrict y,
      const complex_float * restrict R, 
      const complex_float * restrict D,
      const complex_float * restrict A, 
      const complex_float * restrict B, 
      int M, int N, int P, int L );

void rcholfwd8x8x1nf (
            void * pScr,
            float32_t * restrict y,
      const float32_t * restrict R, 
      const float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict B, 
      int L );
void rcholfwd16x16x1nf (
            void * pScr,
            float32_t * restrict y,
      const float32_t * restrict R, 
      const float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict B, 
      int L );
void rcholfwdmxnxpnf (
            void * pScr,
            float32_t * restrict y,
      const float32_t * restrict R, 
      const float32_t * restrict D,
      const float32_t * restrict A, 
      const float32_t * restrict B, 
      int M, int N, int P, int L );

/* Return the scratch area size, in bytes. */
size_t cholfwd8x8x1nf_getScratchSize   ( int M, int N, int P, int L );
size_t cholfwd16x16x1nf_getScratchSize ( int M, int N, int P, int L );
size_t cholfwdmxnxpnf_getScratchSize   ( int M, int N, int P, int L );
/* Return the scratch area size, in bytes. */
size_t rcholfwd8x8x1nf_getScratchSize   ( int M, int N, int P, int L );
size_t rcholfwd16x16x1nf_getScratchSize ( int M, int N, int P, int L );
size_t rcholfwdmxnxpnf_getScratchSize   ( int M, int N, int P, int L );

/*===========================================================================
  Cholesky backward recursion for real/complex block ordered floating-point
  matrices
  [r]cholbkw8x1nf   8x1 real/complex solution vectors
  [r]cholbkw16x1nf  16x1 real/complex solution vectors
  [r]cholbkwnxpnf   NxP real/complex solution vectors
===========================================================================*/

/*-------------------------------------------------------------------------
These functions make backward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices and results of forward recursion. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. 

Storage sizes SR,SD,SY,SX denote the number of data elements required to store a
matrix in block order. If matrix size is less than the SIMD vector size, then the
storage_size(matrix_size) equals the matrix_size rounded up to the next power of
two, otherwise it is matrix_size rounded up to the next multiple of the SIMD
vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SR = storage_size((N+1)*N/2)
SD = storage_size(N)
SY = storage_size(N*P)
SX = storage_size(N*P)

Scratch size in bytes is defined by [r]cholbkw<...>nf_getScratchSize()

Data format: IEEE-754 Std. single precision floating-point

Input:
 N         Matrix dimension (number of columns and rows in matrices R)
 P         Number of columns in right-side matrices B
 L         Number of matrices
 R[L][SR]  Sequence of L upper triangular complex matrices R
 D[L][SD]  Reciprocal of main diagonal
 y[L][SY]  Sequence of intermediate decision matrices y
Output:         
 x[L][SX]  Sequence of decision matrix x

Restrictions:
1. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
2. Matrix sizes M,N,P must be positive
3. M and N must be multiples of 4
---------------------------------------------------------------------------*/
/* Complex-valued functions */
void cholbkw8x1nf (
            void * pScr,
            complex_float * restrict x, 
      const complex_float * restrict R,
      const complex_float * restrict D,
      const complex_float * restrict y, 
      int L );
void cholbkw16x1nf (
            void * pScr,
            complex_float * restrict x, 
      const complex_float * restrict R,
      const complex_float * restrict D,
      const complex_float * restrict y, 
      int L );
void cholbkwnxpnf (
            void * pScr,
            complex_float * restrict x, 
      const complex_float * restrict R,
      const complex_float * restrict D,
      const complex_float * restrict y, 
      int N, int P, int L );

void rcholbkw8x1nf (
            void * pScr,
            float32_t * restrict x, 
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict y, 
      int L );
void rcholbkw16x1nf (
            void * pScr,
            float32_t * restrict x, 
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict y, 
      int L );
void rcholbkwnxpnf (
            void * pScr,
            float32_t * restrict x, 
      const float32_t * restrict R,
      const float32_t * restrict D,
      const float32_t * restrict y, 
      int N, int P, int L );

/* Return the scratch area size, in bytes. */
size_t cholbkw8x1nf_getScratchSize  ( int N, int P, int L );
size_t cholbkw16x1nf_getScratchSize ( int N, int P, int L );
size_t cholbkwnxpnf_getScratchSize  ( int N, int P, int L );
/* Return the scratch area size, in bytes. */
size_t rcholbkw8x1nf_getScratchSize  ( int N, int P, int L );
size_t rcholbkw16x1nf_getScratchSize ( int N, int P, int L );
size_t rcholbkwnxpnf_getScratchSize  ( int N, int P, int L );

/*===========================================================================
  Cholesky MMSE solver for real/complex block ordered floating-point data
  [r]cholmmse8x8x1nf   8x8x1 real/complex MMSE problem
  [r]cholmmse16x16x1nf 16x16x1 real/complex MMSE problem
  [r]cholmmsemxnxpnf   MxNxP real/complex MMSE problem
===========================================================================*/

/*-------------------------------------------------------------------------
Compute the MMSE solution for a system of linear equations A*x=b, where A
is an MxN real (complex) matrix with M>=N and rank(A)==N, x is an Nx1 vector
of unknowns, and b is an Mx1 right hand side vector. This task is accomplished
in 3 steps:
  1. Cholesky decomposition is applied to the matrix of normal equations
     system, which results in an upper triangular NxN matrix R with real
     and positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A + sigma2*I,
     where adj(...) denotes the (conjugate) transpose of a matrix, and
     sigma2*I is an NxN identity matrix multiplied by the regularization
     term.
  2. Forward recursion step: the system adj(R)*y=adj(A)*b is solved for an 
     Nx1 vector y.
  3. Backward recursion step: the system R*x=y is solved for the Nx1 vector
     of unknows x.

For a single MxN matrix A, these 3 steps may be done simultaneously for P
variants of Mx1 right hand side column vectors b gathered into an MxP input
matrix B. MMSE solution is computed independently for each of P columns,
with resulting column vectors forming the solution matrix X of size NxP.

Storage sizes SA,SB,SX denote the number of data elements required to store a
matrix in block order. If matrix size is less than the SIMD vector size, then the
storage_size(matrix_size) equals the matrix_size rounded up to the next power of
two, otherwise it is matrix_size rounded up to the next multiple of the SIMD
vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(M*N)
SB = storage_size(M*P)
SX = storage_size(N*P)

[r]cholmmse<size>nf() functions process L pairs of MxN matrices A and MxP 
matrices B in a single call, which results in L solution matrices X of
size NxP. L matrices of each kind are stored as block ordered sequences.

Data format: IEEE-754 Std. single precision floating-point

Temporary:
  pScr       Scratch area. Required size (in bytes) is defined by 
             functions [r]cholmmse<size>sf_getScratchSize(M,N,P,L)
Input:
  M,N,P      Dimensional parameters
  L          Number of matrices
  sigma2[L]  Regularization term
  A[L][SA]  Sequence of L matrices A
  B[L][SB]  Sequence of L right hand side matrices B
Output:
  x[L][SZ]  Sequence of L solution matrices X
Restrictions:
  1. pScr,x,A,B,sigma2 must not overlap
  2. pScr,x,A,B,sigma2 must be aligned on 32-byte boundary
  3. M and N should be a multiple of 4
  4. Matrix sizes M,N must be greater than 1
  5. Number of columns for input matrices A must not exceed the number
     of rows: N <= M.
---------------------------------------------------------------------------*/
/* Complex-valued functions */
void cholmmse8x8x1nf (
            void * pScr,
            complex_float * restrict x,
      const complex_float * restrict A,
      const complex_float * restrict B,
      const float32_t     * restrict sigma2,
      int L );
void cholmmse16x16x1nf (
            void * pScr,
            complex_float * restrict x,
      const complex_float * restrict A,
      const complex_float * restrict B,
      const float32_t     * restrict sigma2,
      int L );
void cholmmsemxnxpnf (
            void * pScr,
            complex_float * restrict x,
      const complex_float * restrict A,
      const complex_float * restrict B,
      const float32_t     * restrict sigma2,
      int M, int N, int P, int L );

void rcholmmse8x8x1nf (
            void * pScr,
            float32_t * restrict x,
      const float32_t * restrict A,
      const float32_t * restrict B,
      const float32_t * restrict sigma2,
      int L );
void rcholmmse16x16x1nf (
            void * pScr,
            float32_t * restrict x,
      const float32_t * restrict A,
      const float32_t * restrict B,
      const float32_t * restrict sigma2,
      int L );
void rcholmmsemxnxpnf (
            void * pScr,
            float32_t * restrict x,
      const float32_t * restrict A,
      const float32_t * restrict B,
      const float32_t * restrict sigma2,
      int M, int N, int P, int L );

/* Return the scratch area size, in bytes. */
size_t cholmmse8x8x1nf_getScratchSize   ( int M, int N, int P, int L );
size_t cholmmse16x16x1nf_getScratchSize ( int M, int N, int P, int L );
size_t cholmmsemxnxpnf_getScratchSize   ( int M, int N, int P, int L );
/* Return the scratch area size, in bytes. */
size_t rcholmmse8x8x1nf_getScratchSize   ( int M, int N, int P, int L );
size_t rcholmmse16x16x1nf_getScratchSize ( int M, int N, int P, int L );
size_t rcholmmsemxnxpnf_getScratchSize   ( int M, int N, int P, int L );

/*===========================================================================
  Banded Cholesky decomposition for a complex-valued fixed-point pseudo-inversion:
  bchol4x16n  19x16 complex matrices (band width 4)
  bchol8x16n  23x16 complex matrices (band width 8)
  bchol4x32n  35x32 complex matrices (band width 4)
  bchol8x32n  39x32 complex matrices (band width 8)
  bchol16x32n 47x32 complex matrices (band width 8)
  bcholwxnn  (W+N-1)xN complex matrices, (band width W)
===========================================================================*/

/*-------------------------------------------------------------------------
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex-valued least squares problem: A*X=B, where A is
an MxN coefficient matrix with M >= N; X is an NxP matrix of unknowns; and
B is an MxP right hand matrix.

The decomposition results in an upper triangular complex NxN matrix R with
real and positive numbers on the main diagonal, such that 
                     adj(R)*R = adj(A)*A + sigma2*I,
where adj(...) denotes the conjugate transpose of a matrix, and sigma2*I is
the NxN identity matrix multiplied with the regularization term.

The decomposition algorithm is applied to a few matrices per single call,
with input/output matrix sequences being stored in the streaming order.

Fixed-point data type of upper triangular matrices R is the same as the
data type of input matrices A. Fixed point position for the regularization
term sigma2 must match the scale of product adj(A)*A. If, for instance,
matrix A is represented as Q15, then Q30 is expected for sigma2.

In order to find the solution to the above-mentioned least squares problem
A*X=B, one has to follow the Cholesky decomposition with forward and backward
substitution procedures; see cholfwdmxnxps() and cholbkwnxps(), respectively.

The code for banded matrices is intended for cases where matrix A contains 
W non-zero elements on the main diagonal and below. So, size M is N+W-1. 
Matrix A may be stored in the compact form of size WxN. Cholesky matrix R 
also has WxN non-zero elements

NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. Specifically, matrix sizes SA,SR,SD are selected as usual for 
complex block ordered matrix sequencies, i.e. total size is rounded up to 
the closest bigger multiple of BBE_SIMD_WIDTH/2==8 elements. 
SA=size(W*N)
SR=size(W*N)
SD=size(N)


Input:
  W             Band width
  N             Dimensional parameters
  L             Number of matrices
  sigma2[L]     Regularization term; fixed point position is twice the
                number of fractional bits for matrices A, R
  At[L][SA]     Sequence of L complex matrices A represented in the 
                compact form (only band)
Output:
  Rt[L][SR]     Sequence of L upper triangular complex matrices R 
                represented in the compact form (saved only elements on 
                the main diagonal and above in such a way that diagoanal
                elements are in the last raw)
  D[L][SD]      Sequence of L reciprocals of main diagonal A represented 
                in the  block floating point (mantissa and exponent).

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. For banded matrices, W must be 4, 8, 12 or 16, N must be a positive multiple of 4
Note:
  Function may speculatively read up to (W-2)*W complex elements
  beyond the upper bound of At and Rt.
---------------------------------------------------------------------------*/
#define BCHOLN_MAXW 16  // maximum band width
void bcholwxnn  (complex_fract16 * restrict Rt, complex_fract16 * restrict D, const complex_fract16 * restrict At, const int32_t * restrict sigma2,int W, int N, int L);
void bchol4x16n (complex_fract16 * restrict Rt, complex_fract16 * restrict D, const complex_fract16 * restrict At, const int32_t * restrict sigma2,int L);
void bchol8x16n (complex_fract16 * restrict Rt, complex_fract16 * restrict D, const complex_fract16 * restrict At, const int32_t * restrict sigma2,int L);
void bchol4x32n (complex_fract16 * restrict Rt, complex_fract16 * restrict D, const complex_fract16 * restrict At, const int32_t * restrict sigma2,int L);
void bchol8x32n (complex_fract16 * restrict Rt, complex_fract16 * restrict D, const complex_fract16 * restrict At, const int32_t * restrict sigma2,int L);
void bchol16x32n(complex_fract16 * restrict Rt, complex_fract16 * restrict D, const complex_fract16 * restrict At, const int32_t * restrict sigma2,int L);

/*===========================================================================
  Banded Cholesky forward recursion for pseudo-inversion API (complex fixed-point data)
  bcholfwd4x16x1s   matrices 19x16x1 (band width 4)
  bcholfwd8x16x1n   matrices 23x16x1 (band width 8)
  bcholfwd4x32x1s   matrices 35x32x1 (band width 4)
  bcholfwd8x32x1s   matrices 39x32x1 (band width 8)
  bcholfwd16x32x1n  matrices 47x32xP (band width 16)
  bcholfwdwxnxpn    matrices (W+N-1)xNxP (band width W)
===========================================================================*/

/*-------------------------------------------------------------------------
These functions make forward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. Specifically, matrix sizes SA,SR,SD,SB,SY are selected as 
usual for complex block ordered matrix sequencies, i.e. total size is 
rounded up to the closest bigger multiple of BBE_SIMD_WIDTH/2==8 elements. 
SA=size(W*N)
SR=size(W*N)
SD=size(N)
SB=size((W+N-1)*P)
SY=size(N*P)

Input:
W             Band width
N             Matrix dimension (number of columns in matrices R)
P             Number of columns in right-side matrices B
L             Number of matrices
Rt[L][SR]     Cholesky upper triangular matrices R represented in the compact
              form (saved only elements on the main diagonal and above in 
              such a way that diagoanal elements are in the last raw)
D[L][SD]      Sequence of L reciprocals of main diagonal A represented in the  
              block floating point (mantissa and exponent). N' is computed as 
              for complex block ordered matrices of size N
At[L][SA]     Original left-side matrices A represented in the compact 
              form (only band)
Bt[L][SB]     Original right-side matrices B. SB is computed as for complex 
              block ordered matrices of size (W+N-1)*P
qA,qB,qY      Fixed point representation of matrices A (or R which is the 
              same),B and y

Output:
Yt[L][SY]     Decision matrix y. SY is computed as for complex 
              block ordered matrices of size N*P

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. For banded matrices, W must be 4, 8, 12 or 16, N must be a positive multiple of 4
3. P>=1
---------------------------------------------------------------------------*/
void bcholfwdwxnxpn   (complex_fract16* restrict Yt,const complex_fract16* restrict Rt, const complex_fract16* restrict D, const complex_fract16* restrict At, const complex_fract16* restrict Bt, int qB,int qY,int W, int N, int P, int L);
void bcholfwd4x16x1n  (complex_fract16* restrict Yt,const complex_fract16* restrict Rt, const complex_fract16* restrict D, const complex_fract16* restrict At, const complex_fract16* restrict Bt, int qB,int qY,int L);
void bcholfwd8x16x1n  (complex_fract16* restrict Yt,const complex_fract16* restrict Rt, const complex_fract16* restrict D, const complex_fract16* restrict At, const complex_fract16* restrict Bt, int qB,int qY,int L);
void bcholfwd4x32x1n  (complex_fract16* restrict Yt,const complex_fract16* restrict Rt, const complex_fract16* restrict D, const complex_fract16* restrict At, const complex_fract16* restrict Bt, int qB,int qY,int L);
void bcholfwd8x32x1n  (complex_fract16* restrict Yt,const complex_fract16* restrict Rt, const complex_fract16* restrict D, const complex_fract16* restrict At, const complex_fract16* restrict Bt, int qB,int qY,int L);
void bcholfwd16x32x1n (complex_fract16* restrict Yt,const complex_fract16* restrict Rt, const complex_fract16* restrict D, const complex_fract16* restrict At, const complex_fract16* restrict Bt, int qB,int qY,int L);
#define BCHOLFWDNSCRATCH(W,N,P,L) ( 0 )

/*===========================================================================
  Banded Cholesky backward recursion for pseudo-inversion API (complex fixed-point data)

  bcholbkw4x16x1n   matrices 19x16x1 (band width 4)
  bcholbkw8x16x1n   matrices 23x16x1 (band width 8)
  bcholbkw4x32x1n   matrices 35x32x1 (band width 4)
  bcholbkw8x32x1n   matrices 39x32x1 (band width 8)
  bcholbkw16x32x1n  matrices 47x32x1 (band width 16)
  bcholbkwwxnxpn    matrices (W+N-1)xNxP (band width W)
===========================================================================*/

/*-------------------------------------------------------------------------
These functions make backward recursion stage of pseudo-inversion. They use 
Cholesky decomposition of original matrices and results of forward recursion. 
NOTE:
Data layout for matrices is selected as for other matrices written in a 
block order. Specifically, matrix sizes SR,SD,SY,SX are selected as usual for 
complex block ordered matrix sequencies, i.e. total size is rounded up to 
the closest bigger multiple of BBE_SIMD_WIDTH/2==8 elements. 
SR=size(W*N)
SD=size(N)
SY=size(N*P)
SX=size(N*P)

Input:
W             Band width
N             Matrix dimension (number of columns in matrices R)
P             Number of columns in right-side matrices B
L             Number of matrices
Rt[L][SR]     Cholesky upper triangular matrices R represented in the compact
              form (saved only elements on the main diagonal and above in 
              such a way that diagoanal elements are in the last raw)
D[L][SD]      Sequence of L reciprocals of main diagonal R represented in the  
              block floating point (mantissa and exponent). N' is computed as 
              for complex block ordered matrices of size N
Yt[L][SY]     Results of forward recursion stage. SY is computed as for complex 
              block ordered matrices of size N*P
qA,qX,qY      Fixed point representation of matrices A(or R which is the same), 
              x and y
Output:
Xt[L][SX]     Decision matrix x. SX is computed as for complex block ordered 
              matrices of size N*P

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. For banded matrices, W must be 4, 8, 12 or 16, N must be a positive multiple of 4
3. P>=1
---------------------------------------------------------------------------*/
void bcholbkwwxnxpn   (complex_fract16* restrict Xt, const complex_fract16* restrict Rt, const complex_fract16* restrict D, const complex_fract16* restrict Yt, int qA, int qY, int qX,int W,int N, int P, int L);
void bcholbkw4x16x1n  (complex_fract16* restrict Xt, const complex_fract16* restrict Rt, const complex_fract16* restrict D, const complex_fract16* restrict Yt, int qA, int qY, int qX,int L);
void bcholbkw8x16x1n  (complex_fract16* restrict Xt, const complex_fract16* restrict Rt, const complex_fract16* restrict D, const complex_fract16* restrict Yt, int qA, int qY, int qX,int L);
void bcholbkw4x32x1n  (complex_fract16* restrict Xt, const complex_fract16* restrict Rt, const complex_fract16* restrict D, const complex_fract16* restrict Yt, int qA, int qY, int qX,int L);
void bcholbkw8x32x1n  (complex_fract16* restrict Xt, const complex_fract16* restrict Rt, const complex_fract16* restrict D, const complex_fract16* restrict Yt, int qA, int qY, int qX,int L);
void bcholbkw16x32x1n (complex_fract16* restrict Xt, const complex_fract16* restrict Rt, const complex_fract16* restrict D, const complex_fract16* restrict Yt, int qA, int qY, int qX,int L);
#define BCHOLBKWNSCRATCH(W,N,P,L) ( 0 )

/*===========================================================================
QR-based matrix decomposition and inversion for stream ordered fixed-point 
data
[c]qrMxNs             QR decomposition of MxN complex/real matrices by 
                      Householder reflections with explicit computation of
                      Q and R factors
[c]qr_build_rMxNs     Compute Householder reflection vectors for MxN complex/
                      real matrices with no explicit computation of Q and R
                      factors
[c]qr_calc_qbMxNxPs   Apply Householder reflections to the right hand side of
                      a real- or complex-valued system of linear equations
[c]cqr_bkwNxPs        Perform back substitution for the right hand side of
                      a real- or complex-valued system of linear equations
                      transformed by Householder reflections
===========================================================================*/

/*-------------------------------------------------------------------------
These functions apply QR decomposition procedure to the sequence of complex 
matrices written in a streaming order. The transformation is done in-place 
so the result replaces the original input.
Rotation matrix Q is calculated in Q15 fixed point representation. Fixed 
point representation of upper-diagonal matrix R is the same as of input. 

Functions return nonzero if overflow is detected 

NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Input:
R[M*N][L][C]   input matrices
C              1 for real, 2 for complex data
Output:
Q[M*M][L][C]   output rotation matrices (L matrices of size MxM)
R[M*N][L][C]   output upper triangular matrices (L matrices of size MxN)

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Matrix sizes M,N,L must be greater than 1
4. Scratch memory must be aligned on 32-byte boundary. Its size  (in bytes)
   is defined by cqr2x2s_getScratchSize, qr2x2s_getScratchSize
---------------------------------------------------------------------------*/
// complex data
void cqr2x2s (void* pScr, complex_fract16 * restrict Q, complex_fract16 * restrict R, int L);
// real data
void qr2x2s (void* pScr, int16_t * restrict Q, int16_t * restrict R, int L);

// sizeof for Q matrices
#define CQR2x2_SIZE_Q(L) ( 4*2*(L)*sizeof(int16_t) )
#define  QR2x2_SIZE_Q(L) ( 4*  (L)*sizeof(int16_t) )

// scratch memory needed for QR decomposition
size_t cqr2x2s_getScratchSize(int M, int N,int L);
size_t  qr2x2s_getScratchSize(int M, int N,int L);

/*-----------------------------------------------------------------------
[c]qr_build_rMxNs

QR decomposition of MxN complex matrices.
Instead of direct computation of Q factors, these functions produce a
set of N Householder vectors V for each of input matrices A. This
approach allow us to save CPU cycles and memory when solving a system
of linear equations: it is cheaper to perform N elementary reflections
for a right hand side vector if compared to explicit multiplication of
that vector by matrix Q.

Fixed point representation of output matrices R is the same as for
input matrices, but Householder vectors V are always Q14.

Data transform is performed in-place: upper triangular matrices R replace
input matrices A.

NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Input:
R[M*N][L]                  matrices A (L matrices of size MxN)
Output:
V[((M*N+((N-1)*N)/2+M)*L]  L sets of Householder vectors
R[M*N][L]                  upper triangular matrices (L matrices
                           of size MxN)

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,N,L)
4. M must greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
-------------------------------------------------------------------------*/
// complex data
void cqr_build_r3x3s   ( void* pScr, complex_fract16 * restrict V, complex_fract16 * restrict R, int L); 
void cqr_build_r4x2s   ( void* pScr, complex_fract16 * restrict V, complex_fract16 * restrict R, int L); 
void cqr_build_r4x4s   ( void* pScr, complex_fract16 * restrict V, complex_fract16 * restrict R, int L); 
void cqr_build_r5x5s   ( void* pScr, complex_fract16 * restrict V, complex_fract16 * restrict R, int L); 
void cqr_build_r6x6s   ( void* pScr, complex_fract16 * restrict V, complex_fract16 * restrict R, int L); 
void cqr_build_r7x7s   ( void* pScr, complex_fract16 * restrict V, complex_fract16 * restrict R, int L); 
void cqr_build_r8x4s   ( void* pScr, complex_fract16 * restrict V, complex_fract16 * restrict R, int L); 
void cqr_build_r8x8s   ( void* pScr, complex_fract16 * restrict V, complex_fract16 * restrict R, int L); 
void cqr_build_r16x8s  ( void* pScr, complex_fract16 * restrict V, complex_fract16 * restrict R, int L); 
void cqr_build_r16x16s ( void* pScr, complex_fract16 * restrict V, complex_fract16 * restrict R, int L); 
void cqr_build_rmxns   ( void* pScr, complex_fract16 * restrict V, complex_fract16 * restrict R, int M, int N,int L);
// real data
void qr_build_r4x4s   ( void* pScr, int16_t * restrict V, int16_t * restrict R, int L); 
void qr_build_r8x8s   ( void* pScr, int16_t * restrict V, int16_t * restrict R, int L); 
void qr_build_rmxns   ( void* pScr, int16_t * restrict V, int16_t * restrict R, int M, int N,int L);

/* scratch memory needed for buildr functions */
size_t cqr_build_r3x3s_getScratchSize  (int M, int N,int L);
size_t cqr_build_r4x2s_getScratchSize  (int M, int N,int L);
size_t cqr_build_r4x4s_getScratchSize  (int M, int N,int L);
size_t cqr_build_r5x5s_getScratchSize  (int M, int N,int L);
size_t cqr_build_r6x6s_getScratchSize  (int M, int N,int L);
size_t cqr_build_r7x7s_getScratchSize  (int M, int N,int L);
size_t cqr_build_r8x4s_getScratchSize  (int M, int N,int L);
size_t cqr_build_r8x8s_getScratchSize  (int M, int N,int L);
size_t cqr_build_r16x8s_getScratchSize (int M, int N,int L);
size_t cqr_build_r16x16s_getScratchSize(int M, int N,int L);
size_t cqr_build_rmxns_getScratchSize  (int M, int N,int L);

size_t qr_build_r4x4s_getScratchSize   (int M, int N,int L);
size_t qr_build_r8x8s_getScratchSize   (int M, int N,int L);
size_t qr_build_rmxns_getScratchSize   (int M, int N,int L);

/*-------------------------------------------------------------------------
cqr_calc_qbMxNxPs/qr_calc_qbMxNxPs

These functions apply Householder reflections to L MxP matrices B in the
course of solving a set of complex-valued linear problems A*X=B through
the QR decomposition of matrices A: A*X=B, A=Q*R => Q*R*X=B => R*X=Q'*B.
Instead of direct multiplication of each matrix B by conjugate transpose
of the corresponding matrix Q, we use a set of vectors V to perform a
sequence of Householder  reflections (see QR decomposition routines
cqr_build_rMxN/qr_build_rMxN).

Fixed point representation of output matrices is the same as for input
matrices.

Data transform is performed in-place.

NOTE:
Data layout for matrices is selected as for other matrices written in a stream 
order. 

Input:
B[M*P]L]                   Matrices B (L matrices of size MxP)
V[((M*N+((N-1)*N)/2+M)*L]  L sets of Householder vectors
Output:
B[M*P][L]                  Matrices Q'*B (L matrices of size MxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,P,L)
4. M must be greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
---------------------------------------------------------------------------*/
// complex data
void cqr_calc_qb3x3x1s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb3x3x3s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb4x2x1s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb4x4x1s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb4x4x4s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb5x5x1s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb5x5x5s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb6x6x1s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb6x6x6s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb7x7x1s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb7x7x7s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb8x4x1s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb8x8x1s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb8x8x8s   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb16x8x1s  (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb16x16x1s (void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qb16x16x16s(void *pScr, complex_fract16 *B, const complex_fract16 *V , int L);
void cqr_calc_qbmxnxps   (void *pScr, complex_fract16 *B, const complex_fract16 *V , int M, int N, int P, int L);
// real data
void  qr_calc_qb4x4x1s   (void *pScr, int16_t *B, const int16_t *V , int L);
void  qr_calc_qb4x4x4s   (void *pScr, int16_t *B, const int16_t *V , int L);
void  qr_calc_qb8x8x1s   (void *pScr, int16_t *B, const int16_t *V , int L);
void  qr_calc_qb8x8x8s   (void *pScr, int16_t *B, const int16_t *V , int L);
void  qr_calc_qbmxnxps   (void *pScr, int16_t *B, const int16_t *V , int M, int N, int P, int L);

/* scratch memory needed for calc_qb functions */
size_t cqr_calc_qb3x3x1s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb3x3x3s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb4x2x1s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb4x4x1s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb4x4x4s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb5x5x1s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb5x5x5s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb6x6x1s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb6x6x6s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb7x7x1s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb7x7x7s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb8x4x1s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb8x8x1s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb8x8x8s_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb16x8x1s_getScratchSize  (int M, int P, int L);
size_t cqr_calc_qb16x16x1s_getScratchSize (int M, int P, int L);
size_t cqr_calc_qb16x16x16s_getScratchSize(int M, int P, int L);
size_t cqr_calc_qbmxnxps_getScratchSize   (int M, int P, int L);

size_t  qr_calc_qb4x4x1s_getScratchSize   (int M, int P, int L);
size_t  qr_calc_qb4x4x4s_getScratchSize   (int M, int P, int L);
size_t  qr_calc_qb8x8x1s_getScratchSize   (int M, int P, int L);
size_t  qr_calc_qb8x8x8s_getScratchSize   (int M, int P, int L);
size_t  qr_calc_qbmxnxps_getScratchSize   (int M, int P, int L);

/*-------------------------------------------------------------------------
cqr_bkwNxPs/qr_bkwNxPs

Last stage of solving a set of L complex-valued linear problems A*X=B
through the QR decomposition by Householder reflections: back substitution
process for L systems of complex-valued linear equations R*X=QB, where R is
an MxM upper triangular matrix, X is an MxP matrix of unknowns, QB is an MxP
matrix resulting from Householder reflections being applied to the right
hand matrix B of the original linear problem: QB=Q'*B.

Fixed-point representation for output data is a function of fixed-point
format of input data: FPP(X) = FPP(QB)-FPP(R)+10, where FPP(x) stands for
the Fixed-Point Position of data item x.

Data transform is performed in-place.

NOTE:
1. Data layout for matrices is selected as for other matrices written 
   in a stream order. So, shorter dimension of output matrix B (NxP 
   instead of MxP as on input) does not require special management - 
   remaining (M-N)*P*L elements are kept unchanged

Input
B[M*P][L]  Matrices QB=Q'*B (L matrices of size MxP)
R[M*N][L]  upper triangular matrices R (L matrices of size MxN)
Output:
B[N*P][L]  Matrices X (L matrices of size NxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(N,P,L)
4. Matrix sizes N,L must be greater than 1
---------------------------------------------------------------------------*/
// complex data: for cqr_build_r+cqr_calc_qb 
void cqr_bkw2x1s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw3x1s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw3x3s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw4x1s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw4x4s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw5x1s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw5x5s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw6x1s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw6x6s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw7x1s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw7x7s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw8x1s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw8x8s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw16x1s (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkw16x16s(void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L); 
void cqr_bkwnxps  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int N, int P, int L);
// real data: for cqr_build_r+qr_calc_qb 
void  qr_bkw4x1s  (void* pScr, int16_t* restrict B, const int16_t* restrict R, int L); 
void  qr_bkw4x4s  (void* pScr, int16_t* restrict B, const int16_t* restrict R, int L); 
void  qr_bkw8x1s  (void* pScr, int16_t* restrict B, const int16_t* restrict R, int L); 
void  qr_bkw8x8s  (void* pScr, int16_t* restrict B, const int16_t* restrict R, int L); 
void  qr_bkwnxps  (void* pScr, int16_t* restrict B, const int16_t* restrict R, int N, int P, int L);

// complex data: special variant for cqr2x2s
void cqr_bkw2x2x1s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict Q, const complex_fract16* restrict  R, int L); 
void cqr_bkw2x2x2s  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict Q, const complex_fract16* restrict  R, int L); 
// real data: special variant for qr2x2s
void  qr_bkw2x2x1s  (void* pScr, int16_t* restrict B, const int16_t* restrict Q, const int16_t* restrict  R, int L); 
void  qr_bkw2x2x2s  (void* pScr, int16_t* restrict B, const int16_t* restrict Q, const int16_t* restrict  R, int L); 

// get scratch size required for bkw functions
size_t cqr_bkw2x1s_getScratchSize  (int N, int P, int L);
size_t cqr_bkw3x1s_getScratchSize  (int N, int P, int L);
size_t cqr_bkw3x3s_getScratchSize  (int N, int P, int L);
size_t cqr_bkw4x1s_getScratchSize  (int N, int P, int L);
size_t cqr_bkw4x4s_getScratchSize  (int N, int P, int L);
size_t cqr_bkw5x1s_getScratchSize  (int N, int P, int L);
size_t cqr_bkw5x5s_getScratchSize  (int N, int P, int L);
size_t cqr_bkw6x1s_getScratchSize  (int N, int P, int L);
size_t cqr_bkw6x6s_getScratchSize  (int N, int P, int L);
size_t cqr_bkw7x1s_getScratchSize  (int N, int P, int L);
size_t cqr_bkw7x7s_getScratchSize  (int N, int P, int L);
size_t cqr_bkw8x1s_getScratchSize  (int N, int P, int L);
size_t cqr_bkw8x8s_getScratchSize  (int N, int P, int L);
size_t cqr_bkw16x1s_getScratchSize (int N, int P, int L);
size_t cqr_bkw16x16s_getScratchSize(int N, int P, int L);
size_t cqr_bkwnxps_getScratchSize  (int N, int P, int L);
size_t  qr_bkw4x1s_getScratchSize  (int N, int P, int L);
size_t  qr_bkw4x4s_getScratchSize  (int N, int P, int L);
size_t  qr_bkw8x1s_getScratchSize  (int N, int P, int L);
size_t  qr_bkw8x8s_getScratchSize  (int N, int P, int L);
size_t  qr_bkwnxps_getScratchSize  (int N, int P, int L);

size_t cqr_bkw2x2x1s_getScratchSize(int N, int P, int L);
size_t cqr_bkw2x2x2s_getScratchSize(int N, int P, int L);
size_t  qr_bkw2x2x1s_getScratchSize(int N, int P, int L);
size_t  qr_bkw2x2x2s_getScratchSize(int N, int P, int L);

/*===========================================================================
   Apply the QR decomposition to the matrix of normal equations system
   associated with a fixed-point complex-valued least squares problem: A*X=B,
   where A is an MxN coefficient matrix with M >= N; X is an NxP matrix of
   unknowns; and B is an MxP right hand matrix.

   Decomposition decomposes matrix A to the upper triangle matrix R and 
   square rotation (unitary) matrix Q 
   A=QR
   where matrix Q is given implicitly via sequence of Housholder vectors V 
   and diagonal rotation vector Fi.
   So, least squares solution of original problem is given in 3 steps:
   - find decomposition : A->R,V,Fi
   - update right side of equation: V,Fi,B->Z=Q'*B
   - find MMSE solution of upper triangle system via backward recursion:
     R,Z->X=R\Z

   The decomposition algorithm is applied to a few matrices per single call,
   with input/output matrix sequences being stored in the block order.
   QR decomposition makes in-place trasformations for R matrix replaces A in 
   the cqr_buildrmxnn functions, Z matrix replaces B matrix in cqr_calcqbmxnxpn 
   and X matrix replaces Z in cqr_bkwmxnxpn

   User selects fixed point representation for A,B and X matrices. It is assumed 
   matrix Z has the same fixed point representation as matrix B. Matrix R has the 
   same representation as matrix A. Housholder vectors V are scaled by sqrt(2) 
   and stored in Q14 format. Rotation vector Fi goes in Q15 representation.
   Housholder and rotation vectors are stored into a single array in special format:
   L Housholder vectors for 1-st column (of length M), L Housholder vectors for 
   2-nd column (of length M-1)..., L Housholder vectors for N-th column (of length 
   M-N+1), L 1-st elements of rotation vector... L N-th elements of rotation 
   vector. Totally, they require ((2*M-N+1)*N/2+N)*L complex elements

   Matrix sizes SA,SB,SV are selected as usual for complex block ordered matrix 
   sequencies, i.e. total size is rounded up to the closest bigger multiple of 
   BBE_SIMD_WIDTH/2==8 elements or, if it is lesser, to the closest bigger 
   multiple of degree of 2. 
   SA=size(M*N)
   SB=size(M*P)
   SV=size(((2*M-N+1)*N/2+N)*L)
===========================================================================*/
/*-------------------------------------------------------------------------
Make QR decomposition for block ordered matrices.
Matrix sizes SA,SV are selected as usual for complex block ordered matrix 
sequencies, i.e. total size is rounded up to the closest bigger multiple of 
BBE_SIMD_WIDTH/2==8 elements or, if it is lesser, to the closest bigger 
multiple of degree of 2. 
SA=size(M*N)
SV=size(((2*M-N+1)*N/2+N)*L)
SD=size(N)
Scratch size in bytes is defined by functions cqr_build_rmxnn_getScratchSize()

Input:
 M, N         Dimensional parameters
 L            Number of matrices
Input/output:
 A[L][SA]     On input it is the sequence of L complex matrices A. 
              At the end of the process, matrices R replace input
              matrices A. In a case of non-square matrices (N!=M),
              only N*N elements of each output matrix are valid.
Output:
 V[SV]        Sequence of L Housholder rotation vectors 
 D[L][SD]     Reciprocal of main diagonal (mantissa, exponent) 
              in a special format

Restrictions:
1. A, V, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. Number of columns for input matrices A must not exceed the number
   of rows: N <= M.
---------------------------------------------------------------------------*/
void cqr_build_rmxnn(void *pScr,
                    complex_fract16* A,
                    complex_fract16* V,
                    complex_fract16* D,
                    int M, int N,
                    int L);
void cqr_build_r8x8n(void *pScr,
                    complex_fract16* A,
                    complex_fract16* V,
                    complex_fract16* D,
                    int L);
void cqr_build_r16x16n(void *pScr,
                    complex_fract16* A,
                    complex_fract16* V,
                    complex_fract16* D,
                    int L);
void cqr_build_r32x32n(void *pScr,
                    complex_fract16* A,
                    complex_fract16* V,
                    complex_fract16* D,
                    int L);

/* scratch memory needed for build_r functions */
size_t cqr_build_rmxnn_getScratchSize(int M, int N,int L);
size_t cqr_build_r8x8n_getScratchSize(int M, int N,int L);
size_t cqr_build_r16x16n_getScratchSize(int M, int N,int L);
size_t cqr_build_r32x32n_getScratchSize(int M, int N,int L);

/*-------------------------------------------------------------------------
Update right side of equations for QR process for block ordered matrices.
Matrix sizes SB,SV are selected as usual for complex block ordered matrix 
sequencies, i.e. total size is rounded up to the closest bigger multiple of 
BBE_SIMD_WIDTH/2==8 elements or, if it is lesser, to the closest bigger 
multiple of degree of 2. 
SB=size(M*P)
SV=size(((2*M-N+1)*N/2+N)*L)
Scratch size in bytes is defined by cqr_calc_qbmxnn_getScratchSize(M,N,P,L)
functions

Input:
 M, N, P      dimensional parameters
 L            Number of matrices
Input/output:
 B[L][SB]     On input it is the sequence of L complex matrices B. 
              At the end of the process, matrices Z replace input
              matrices A. In a case of non-square matrices (N!=M), 
              only N*P elements of each output matrix will be valid.
Input:
 V[SV]        Sequence of L Housholder rotation vectors 

Restrictions:
1. B, V, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. N <= M
---------------------------------------------------------------------------*/
void cqr_calc_qbmxnxpn(void *pScr,
                    complex_fract16* B,
                    const complex_fract16* V,
                    int M, int N, int P,
                    int L);
void cqr_calc_qb8x8x1n(void *pScr,
                    complex_fract16* B,
                    const complex_fract16* V,
                    int L);
void cqr_calc_qb8x8x8n(void *pScr,
                    complex_fract16* B,
                    const complex_fract16* V,
                    int L);
void cqr_calc_qb16x16x1n(void *pScr,
                    complex_fract16* B,
                    const complex_fract16* V,
                    int L);
void cqr_calc_qb16x16x16n(void *pScr,
                    complex_fract16* B,
                    const complex_fract16* V,
                    int L);
void cqr_calc_qb32x32x1n(void *pScr,
                    complex_fract16* B,
                    const complex_fract16* V,
                    int L);
void cqr_calc_qb32x32x32n(void *pScr,
                    complex_fract16* B,
                    const complex_fract16* V,
                    int L);
/* scratch memory needed for calc_qb functions */
size_t cqr_calc_qbmxnxpn_getScratchSize(int M, int N,int P,int L);
size_t cqr_calc_qb8x8x1n_getScratchSize(int M, int N,int P,int L);
size_t cqr_calc_qb8x8x8n_getScratchSize(int M, int N,int P,int L);
size_t cqr_calc_qb16x16x1n_getScratchSize(int M, int N,int P,int L);
size_t cqr_calc_qb16x16x16n_getScratchSize(int M, int N,int P,int L);
size_t cqr_calc_qb32x32x1n_getScratchSize(int M, int N,int P,int L);
size_t cqr_calc_qb32x32x32n_getScratchSize(int M, int N,int P,int L);

/*-------------------------------------------------------------------------
Apply backward recursion process for QR decomposition for block ordered 
matrices.
Matrix sizes SA,SB are selected as usual for complex block ordered matrix 
sequencies, i.e. total size is rounded up to the closest bigger multiple of 
BBE_SIMD_WIDTH/2==8 elements or, if it is lesser, to the closest bigger 
multiple of degree of 2. 
SA=size(M*N)
SB=size(M*P)
SD=size(N)
Scratch size in bytes is defined by cqr_bkwmxnxpn_getScratchSize(M,N,P,L)
functions

Input:
 M, N, P      Dimensional parameters
 L            Number of matrices
 qABX         qA-qB+qX where qA,qB,qX - fixed point representations of 
              matrices A,B,X
Input/output:
 X[L][SB]     On input it is the sequence of L updated right parts Z=Q'B.
              They will be replaced with MMSE solution vectors X (only N*P 
              elements are used)
Input:
 R[L][SA]     Upper triangular matrices R (only N*N 
              elements of each matrix are used)
 D[L][SD]     Reciprocal of main diagonal (mantissa, exponent) 
              in the special format

Restrictions:
1. X, R, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. N <= M
---------------------------------------------------------------------------*/
void cqr_bkwmxnxpn(void *pScr,
                    complex_fract16* X,
                    const complex_fract16* R,
                    const complex_fract16* D,
                    int qABX,
                    int M, int N, int P,
                    int L);
void cqr_bkw8x8x1n(void *pScr,
                    complex_fract16* X,
                    const complex_fract16* R,
                    const complex_fract16* D,
                    int qABX,
                    int L);
void cqr_bkw8x8x8n(void *pScr,
                    complex_fract16* X,
                    const complex_fract16* R,
                    const complex_fract16* D,
                    int qABX,
                    int L);
void cqr_bkw16x16x1n(void *pScr,
                    complex_fract16* X,
                    const complex_fract16* R,
                    const complex_fract16* D,
                    int qABX,
                    int L);
void cqr_bkw16x16x16n(void *pScr,
                    complex_fract16* X,
                    const complex_fract16* R,
                    const complex_fract16* D,
                    int qABX,
                    int L);
void cqr_bkw32x32x1n(void *pScr,
                    complex_fract16* X,
                    const complex_fract16* R,
                    const complex_fract16* D,
                    int qABX,
                    int L);
void cqr_bkw32x32x32n(void *pScr,
                    complex_fract16* X,
                    const complex_fract16* R,
                    const complex_fract16* D,
                    int qABX,
                    int L);

/* scratch memory needed for bkw functions */
size_t cqr_bkwmxnxpn_getScratchSize    (int M, int N,int P,int L);
size_t cqr_bkw8x8x1n_getScratchSize    (int M, int N,int P,int L);
size_t cqr_bkw8x8x8n_getScratchSize    (int M, int N,int P,int L);
size_t cqr_bkw16x16x1n_getScratchSize  (int M, int N,int P,int L);
size_t cqr_bkw16x16x16n_getScratchSize (int M, int N,int P,int L);
size_t cqr_bkw32x32x1n_getScratchSize  (int M, int N,int P,int L);
size_t cqr_bkw32x32x32n_getScratchSize (int M, int N,int P,int L);

/*===========================================================================
  QR decomposition/solution for real/complex stream ordered floating-point 
  matrices
===========================================================================*/

/*-----------------------------------------------------------------------
[c]qr_build_rMxNsf

QR decomposition of MxN real/complex matrices.
Instead of direct computation of Q factors, these functions produce a
set of N Householder vectors V for each of input matrices A. This
approach allow us to save CPU cycles and memory when solving a system
of linear equations: it is cheaper to perform N elementary reflections
for a right hand side vector if compared to explicit multiplication of
that vector by matrix Q.

Data transform is performed in-place: upper triangular matrices R replace
input matrices A.

Data format: IEEE-754 Std single precision floating-point

NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Input:
R[M*N][L]               matrices A (L matrices of size MxN)
Output:
V[((2*M-N+1)*N/2+N)*L]  L sets of Householder vectors
R[M*N][L]               upper triangular matrices (L matrices
                        of size MxN)
D[N*L]                  reciprocals of main diagonal

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 4 for complex data and 
   8 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,N,L)
4. M must greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
-------------------------------------------------------------------------*/
// complex data
void cqr_build_r2x2sf   ( void* pScr, complex_float* restrict V, complex_float * restrict R, complex_float * restrict D, int L); 
void cqr_build_r3x3sf   ( void* pScr, complex_float* restrict V, complex_float * restrict R, complex_float * restrict D, int L); 
void cqr_build_r4x4sf   ( void* pScr, complex_float* restrict V, complex_float * restrict R, complex_float * restrict D, int L); 
void cqr_build_rmxnsf   ( void* pScr, complex_float* restrict V, complex_float * restrict R, complex_float * restrict D, int M, int N,int L);
// real data
void qr_build_r2x2sf   ( void* pScr, float32_t* restrict V, float32_t * restrict R, float32_t * restrict D, int L); 
void qr_build_r3x3sf   ( void* pScr, float32_t* restrict V, float32_t * restrict R, float32_t * restrict D, int L); 
void qr_build_r4x4sf   ( void* pScr, float32_t* restrict V, float32_t * restrict R, float32_t * restrict D, int L); 
void qr_build_rmxnsf   ( void* pScr, float32_t* restrict V, float32_t * restrict R, float32_t * restrict D, int M, int N,int L);

/* scratch memory needed for buildr functions */
size_t cqr_build_r2x2sf_getScratchSize  (int M, int N,int L);
size_t cqr_build_r3x3sf_getScratchSize  (int M, int N,int L);
size_t cqr_build_r4x4sf_getScratchSize  (int M, int N,int L);
size_t cqr_build_rmxnsf_getScratchSize  (int M, int N,int L);
size_t qr_build_r2x2sf_getScratchSize  (int M, int N,int L);
size_t qr_build_r3x3sf_getScratchSize  (int M, int N,int L);
size_t qr_build_r4x4sf_getScratchSize  (int M, int N,int L);
size_t qr_build_rmxnsf_getScratchSize  (int M, int N,int L);

/*-------------------------------------------------------------------------
[c]qr_calc_qbMxNxPsf

These functions apply Householder reflections to L MxP matrices B in the
course of solving a set of complex-valued linear problems A*X=B through
the QR decomposition of matrices A: A*X=B, A=Q*R => Q*R*X=B => R*X=Q'*B.
Instead of direct multiplication of each matrix B by conjugate transpose
of the corresponding matrix Q, we use a set of vectors V to perform a
sequence of Householder  reflections (see QR decomposition routines
[c]qr_build_rMxNsf.

Data format: IEEE-754 Std single precision floating-point

Data transform is performed in-place.

NOTE:
Data layout for matrices is selected as for other matrices written in a stream 
order. 

Input:
B[M*P]L]                Matrices B (L matrices of size MxP)
V[((2*M-N+1)*N/2+N)*L]  L sets of Householder vectors
Output:
B[M*P][L]               Matrices Q'*B (L matrices of size MxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 4 for complex data and 
   8 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,P,L)
4. M must be greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
---------------------------------------------------------------------------*/
// complex data
void cqr_calc_qb2x2x1sf  (void *pScr, complex_float *B, const complex_float *V , int L);
void cqr_calc_qb3x3x1sf  (void *pScr, complex_float *B, const complex_float *V , int L);
void cqr_calc_qb4x4x1sf  (void *pScr, complex_float *B, const complex_float *V , int L);
void cqr_calc_qbmxnxpsf  (void *pScr, complex_float *B, const complex_float *V , int M, int N, int P, int L);
// real data
void qr_calc_qb2x2x1sf  (void *pScr, float32_t *B, const float32_t *V , int L);
void qr_calc_qb3x3x1sf  (void *pScr, float32_t *B, const float32_t *V , int L);
void qr_calc_qb4x4x1sf  (void *pScr, float32_t *B, const float32_t *V , int L);
void qr_calc_qbmxnxpsf  (void *pScr, float32_t *B, const float32_t *V , int M, int N, int P, int L);

/* scratch memory needed for calc_qb functions */
size_t cqr_calc_qb2x2x1sf_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb3x3x1sf_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qb4x4x1sf_getScratchSize   (int M, int P, int L);
size_t cqr_calc_qbmxnxpsf_getScratchSize   (int M, int P, int L);

size_t qr_calc_qb2x2x1sf_getScratchSize   (int M, int P, int L);
size_t qr_calc_qb3x3x1sf_getScratchSize   (int M, int P, int L);
size_t qr_calc_qb4x4x1sf_getScratchSize   (int M, int P, int L);
size_t qr_calc_qbmxnxpsf_getScratchSize   (int M, int P, int L);

/*-------------------------------------------------------------------------
[c]qr_bkwNxPsf

Last stage of solving a set of L complex-valued linear problems A*X=B
through the QR decomposition by Householder reflections: back substitution
process for L systems of complex-valued linear equations R*X=QB, where R is
an MxM upper triangular matrix, X is an MxP matrix of unknowns, QB is an MxP
matrix resulting from Householder reflections being applied to the right
hand matrix B of the original linear problem: QB=Q'*B.

Data transform is performed in-place.

NOTE:
1. Data layout for matrices is selected as for other matrices written 
   in a stream order. So, shorter dimension of output matrix B (NxP 
   instead of MxP as on input) does not require special management - 
   remaining (M-N)*P*L elements are kept unchanged

Input
B[M*P][L]  Matrices QB=Q'*B (L matrices of size MxP)
R[M*N][L]  upper triangular matrices R (L matrices of size MxN)
D[N*L]     reciprocals of main diagonal written in a special format
Output:
B[N*P][L]  Matrices X (L matrices of size NxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 4 for complex data and 
   8 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(N,P,L)
4. Matrix sizes N,L must be greater than 1
---------------------------------------------------------------------------*/
// complex data: for cqr_build_r+cqr_calc_qb 
void cqr_bkw2x1sf  (void* pScr, complex_float* restrict B, const complex_float* restrict R, const complex_float* restrict D,int L); 
void cqr_bkw3x1sf  (void* pScr, complex_float* restrict B, const complex_float* restrict R, const complex_float* restrict D,int L); 
void cqr_bkw4x1sf  (void* pScr, complex_float* restrict B, const complex_float* restrict R, const complex_float* restrict D,int L); 
void cqr_bkwnxpsf  (void* pScr, complex_float* restrict B, const complex_float* restrict R, const complex_float* restrict D,int N, int P, int L);
// real data: for cqr_build_r+qr_calc_qb 
void qr_bkw2x1sf  (void* pScr, float32_t* restrict B, const float32_t* restrict R, const float32_t* restrict D, int L); 
void qr_bkw3x1sf  (void* pScr, float32_t* restrict B, const float32_t* restrict R, const float32_t* restrict D, int L); 
void qr_bkw4x1sf  (void* pScr, float32_t* restrict B, const float32_t* restrict R, const float32_t* restrict D, int L); 
void qr_bkwnxpsf  (void* pScr, float32_t* restrict B, const float32_t* restrict R, const float32_t* restrict D, int N, int P, int L);

/* get scratch size required for bkw functions */
size_t cqr_bkw2x1sf_getScratchSize  (int N, int P, int L);
size_t cqr_bkw3x1sf_getScratchSize  (int N, int P, int L);
size_t cqr_bkw4x1sf_getScratchSize  (int N, int P, int L);
size_t cqr_bkwnxpsf_getScratchSize  (int N, int P, int L);
size_t  qr_bkw2x1sf_getScratchSize  (int N, int P, int L);
size_t  qr_bkw3x1sf_getScratchSize  (int N, int P, int L);
size_t  qr_bkw4x1sf_getScratchSize  (int N, int P, int L);
size_t  qr_bkwnxpsf_getScratchSize  (int N, int P, int L);
/*===========================================================================
  QR decomposition/solution for real/complex block ordered floating-point 
  matrices
===========================================================================*/

/*-------------------------------------------------------------------------
Make QR decomposition for block ordered matrices.
Matrix sizes SA,SV are selected as usual for block ordered matrix 
sequencies of corresponding type, i.e. total size is rounded up to the 
closest bigger multiple of 
- BBE_SIMD_WIDTH/2==8 elements for float32_t
- BBE_SIMD_WIDTH/4==4 elements for complex_float
or, if it is less, to the closest bigger 
multiple of degree of 2. 
SA=size(M*N)
SV=size(((2*M-N+1)*N/2+N)*L)
SD=size(N)
Scratch size in bytes is defined by functions xxxx_getScratchSize()

Input:
 M, N         Dimensional parameters
 L            Number of matrices
Input/output:
 A[L][SA]     On input it is the sequence of L matrices A. 
              At the end of the process, matrices R replace input
              matrices A. In a case of non-square matrices (N!=M),
              only N*N elements of each output matrix are valid.
Output:
 V[SV]        Sequence of L Housholder rotation vectors 
 D[L*SD]      Reciprocals of main diagonal in a special format

Restrictions:
1. A, V, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. Number of columns for input matrices A must not exceed the number
   of rows: N <= M.
---------------------------------------------------------------------------*/
void cqr_build_rmxnnf  (void *pScr,complex_float* A,complex_float* V,complex_float* D,int M, int N,int L);
void cqr_build_r8x8nf  (void *pScr,complex_float* A,complex_float* V,complex_float* D,int L);
void cqr_build_r16x16nf(void *pScr,complex_float* A,complex_float* V,complex_float* D,int L);

void  qr_build_rmxnnf  (void *pScr,float32_t* A,float32_t* V,float32_t* D,int M, int N,int L);
void  qr_build_r8x8nf  (void *pScr,float32_t* A,float32_t* V,float32_t* D,int L);
void  qr_build_r16x16nf(void *pScr,float32_t* A,float32_t* V,float32_t* D,int L);

/* scratch allocation functions */
size_t cqr_build_rmxnnf_getScratchSize  (int M, int N,int L);
size_t cqr_build_r8x8nf_getScratchSize  (int M, int N,int L);
size_t cqr_build_r16x16nf_getScratchSize(int M, int N,int L);
size_t  qr_build_rmxnnf_getScratchSize  (int M, int N,int L);
size_t  qr_build_r8x8nf_getScratchSize  (int M, int N,int L);
size_t  qr_build_r16x16nf_getScratchSize(int M, int N,int L);

/*-------------------------------------------------------------------------
Update right side of equations for QR process for block ordered matrices.
Matrix sizes SB,SV are selected as usual for block ordered matrix 
sequencies of corresponding type, i.e. total size is rounded up to the 
closest bigger multiple of 
- BBE_SIMD_WIDTH/2==8 elements for float32_t
- BBE_SIMD_WIDTH/4==4 elements for complex_float
or, if it is less, to the closest bigger 
multiple of degree of 2.  
SB=size(M*P)
SV=size(((2*M-N+1)*N/2+N)*L)
Scratch size in bytes is defined by cqr_calc_qbmxnn_getScratchSize(M,N,P,L)
functions

Input:
 M, N, P      dimensional parameters
 L            Number of matrices
Input/output:
 B[L][SB]     On input it is the sequence of L matrices B. 
              At the end of the process, matrices Z replace input
              matrices A. In a case of non-square matrices (N!=M), 
              only N*P elements of each output matrix will be valid.
Input:
 V[SV]        Sequence of L Housholder rotation vectors 

Restrictions:
1. B, V, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. N <= M
---------------------------------------------------------------------------*/
void cqr_calc_qbmxnxpnf(void *pScr,
                    complex_float* B,const complex_float* V,int M, int N, int P,int L);
void cqr_calc_qb8x8x1nf(void *pScr,
                    complex_float* B,const complex_float* V,int L);
void cqr_calc_qb16x16x1nf(void *pScr,
                    complex_float* B,const complex_float* V,int L);
void qr_calc_qbmxnxpnf(void *pScr,
                    float32_t* B,const float32_t* V,int M, int N, int P,int L);
void qr_calc_qb8x8x1nf(void *pScr,
                    float32_t* B,const float32_t* V,int L);
void qr_calc_qb16x16x1nf(void *pScr,
                    float32_t* B,const float32_t* V,int L);
/* scratch memory needed for calc_qb functions */
size_t cqr_calc_qbmxnxpnf_getScratchSize  (int M, int N,int P,int L);
size_t cqr_calc_qb8x8x1nf_getScratchSize  (int M, int N,int P,int L);
size_t cqr_calc_qb16x16x1nf_getScratchSize(int M, int N,int P,int L);
size_t qr_calc_qbmxnxpnf_getScratchSize   (int M, int N,int P,int L);
size_t qr_calc_qb8x8x1nf_getScratchSize   (int M, int N,int P,int L);
size_t qr_calc_qb16x16x1nf_getScratchSize (int M, int N,int P,int L);
/*-------------------------------------------------------------------------
Apply backward recursion process for QR decomposition for block ordered 
matrices.
Matrix sizes SA,SB are selected as usual for block ordered matrix 
sequencies of corresponding type, i.e. total size is rounded up to the 
closest bigger multiple of 
- BBE_SIMD_WIDTH/2==8 elements for float32_t
- BBE_SIMD_WIDTH/4==4 elements for complex_float
or, if it is less, to the closest bigger 
multiple of degree of 2. 
SA=size(M*N)
SB=size(M*P)
SD=size(N)
Scratch size in bytes is defined by cqr_bkwmxnxpn_getScratchSize(M,N,P,L)
functions

Input:
 M, N, P      Dimensional parameters
 L            Number of matrices
Input/output:
 X[L][SB]     On input it is the sequence of L updated right parts Z=Q'B.
              They will be replaced with MMSE solution vectors X (only N*P 
              elements are used)
Input:
 R[L][SA]     Upper triangular matrices R (only N*N 
              elements of each matrix are used)
 D[L*SD]      Reciprocals of main diagonal in a special format

Restrictions:
1. X, R, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. N <= M
---------------------------------------------------------------------------*/
void cqr_bkwmxnxpnf(void *pScr,
                          complex_float* X,
                    const complex_float* R,
                    const complex_float* D,
                    int M, int N, int P,
                    int L);
void cqr_bkw8x8x1nf(void *pScr,
                          complex_float* X,
                    const complex_float* R,
                    const complex_float* D,
                    int L);
void cqr_bkw16x16x1nf(void *pScr,
                          complex_float* X,
                    const complex_float* R,
                    const complex_float* D,
                    int L);

void  qr_bkwmxnxpnf(void *pScr,
                          float32_t* X,
                    const float32_t* R,
                    const float32_t* D,
                    int M, int N, int P,
                    int L);
void  qr_bkw8x8x1nf(void *pScr,
                          float32_t* X,
                    const float32_t* R,
                    const float32_t* D,
                    int L);
void  qr_bkw16x16x1nf(void *pScr,
                          float32_t* X,
                    const float32_t* R,
                    const float32_t* D,
                    int L);
/* scratch memory needed for bkw functions */
size_t cqr_bkwmxnxpnf_getScratchSize    (int M, int N,int P,int L);
size_t cqr_bkw8x8x1nf_getScratchSize    (int M, int N,int P,int L);
size_t cqr_bkw16x16x1nf_getScratchSize  (int M, int N,int P,int L);
size_t  qr_bkwmxnxpnf_getScratchSize    (int M, int N,int P,int L);
size_t  qr_bkw8x8x1nf_getScratchSize    (int M, int N,int P,int L);
size_t  qr_bkw16x16x1nf_getScratchSize  (int M, int N,int P,int L);

/*===========================================================================
Singular Value Decomposition (SVD)
[r]svdn         Thin SVD for real/complex block ordered matrices
[r]svds         Thin SVD for real/complex stream ordered matrices
===========================================================================*/

/*---------------------------------------------------------------------------
For a real (complex) matrix A of size MxN with M>=N, SVD functions compute:
  - the ordered sequence of N singular values s[N], s[0]>=s[1]>=...>=s[N-1]
  - N orthogonal (orthonormal) left-singular vectors comprising the MxN 
    column space matrix U
  - N orthogonal (orthonormal) right-singular vectors comprising the NxN 
    row space matrix V
Left- and right-singular vectors form othonormal (orthogonal) bases of the
column (U) and row (V) spaces of matrix A. These bases, together with the
singular values sequence, constitute the Thin Singular Value Decomposition
of A: A = U*S*V'. Here S is an NxN diagonal matrix with s[0..N-1] on its main
diagonal, and V' is the (conjugate) transpose of V.
-----------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
Thin SVD For Real/Complex Block Ordered Matrices

Description: compute the Thin Singular Value Decomposition of L complex 
(real) MxN matrices, with the number of rows greater than or equal to the
number of columns: M>=N. Input and output matrices are stored in block order.

Data format: IEEE-754 Std single precision floating-point

Storage sizes SA,Ss,SV denote the number of data elements required to store
a matrix or a vector in block order. If matrix size is less than the SIMD vector
size, then the storage_size(matrix_size) equals the matrix_size rounded up to
the next power of two, otherwise it is matrix_size rounded up to the next
multiple of the SIMD vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(M*N)
Ss = storage_size(N)
SV = storage_size(N*N)

Notes:
1. SVD implementation may perform in-place transformations of input matrices,
   so INPUT DATA MAY APPEAR DAMAGED after the call.
2. Once U or V matrix is not required, set the corresponding output pointer
   parameter to zero to allow for a lower complexity implementation of the
   SVD algorithm.
3. Floating-point functions assume that input data are reasonably scaled. That
   is, the base-2 exponent e of the maximum absolute value over an input matrix
   belongs to the range -E<e<E, where E = 63-log2(N)/2.

Temporary:
  pScr          Scratch area. Required size (in bytes) is defined by 
                functions [r]svd<size>nf_getScratchSize(N,L)
Input:
  M,N           Matrix dimensions
  L             Number of matrices
  A[L][SA]      MxN input matrices
Output:
  U[L][SA]      MxN matrices comprised of M left-singular column 
                vectors (optional)
  s[L][Ss]      Nx1 vectors of singular values in descending order. In an 
                exceptional case when the iterative algorithm fails to 
                converge for a particular matrix, all elements of the
                respective vector are set to NaN.
  V[L][SV]      NxN matrices comprised of N right-singular column
                vectors (optional)
Restrictions:
  pScr,U,s,V,A  Must not overlap and must be aligned on 32-byte boundary 
  M,N           Must be positive multiples of 4 such that M>=N
---------------------------------------------------------------------------*/
/* Complex-valued functions */
void svd8x8nf (
            void * pScr,
            complex_float * restrict U,
            float32_t     * restrict s,
            complex_float * restrict V,
            complex_float * restrict A,
            int L );
void svd16x16nf (
            void * pScr,
            complex_float * restrict U,
            float32_t     * restrict s,
            complex_float * restrict V,
            complex_float * restrict A,
            int L );
void svdmxnnf (
            void * pScr,
            complex_float * restrict U,
            float32_t     * restrict s,
            complex_float * restrict V,
            complex_float * restrict A,
            int M, int N, int L );


/* Real-valued functions */
void rsvd8x8nf (
            void * pScr,
            float32_t * restrict U,
            float32_t * restrict s,
            float32_t * restrict V,
            float32_t * restrict A,
            int L );
void rsvd16x16nf (
            void * pScr,
            float32_t * restrict U,
            float32_t * restrict s,
            float32_t * restrict V,
            float32_t * restrict A,
            int L );
void rsvdmxnnf (
            void * pScr,
            float32_t * restrict U,
            float32_t * restrict s,
            float32_t * restrict V,
            float32_t * restrict A,
            int M, int N, int L );

/* Return the scratch area size, in bytes. */
size_t svd8x8nf_getScratchSize   ( int M, int N, int L );
size_t svd16x16nf_getScratchSize ( int M, int N, int L );
size_t svdmxnnf_getScratchSize   ( int M, int N, int L );
/* Return the scratch area size, in bytes. */
size_t rsvd8x8nf_getScratchSize   ( int M, int N, int L );
size_t rsvd16x16nf_getScratchSize ( int M, int N, int L );
size_t rsvdmxnnf_getScratchSize   ( int M, int N, int L );

/*-------------------------------------------------------------------------
Thin SVD For Real/Complex Stream Ordered Matrices

Description: compute the Thin Singular Value Decomposition of L complex 
(real) MxN matrices, with the number of rows greater than or equal to the
number of columns: M>=N. Input and output matrices are stored in stream order.

Data format: IEEE-754 Std single precision floating-point

Notes:
1. SVD implementation may perform in-place transformations of input matrices,
   so INPUT DATA MAY APPEAR DAMAGED after the call.
2. Once U or V matrix is not required, set the corresponding output pointer
   parameter to zero to allow for a lower complexity implementation of the
   SVD algorithm.
3. Floating-point functions assume that input data are reasonably scaled. That
   is, the base-2 exponent e of the maximum absolute value over an input matrix
   belongs to the range -E<e<E, where E = 63-log2(N)/2.

Temporary:
  pScr          Scratch area. Required size (in bytes) is defined by 
                functions [r]svd<size>nf_getScratchSize(N,L)
Input:
  M,N           Matrix dimensions
  L             Number of matrices
  A[M*N][L]     MxN input matrices
Output:
  U[M*N][L]     MxN matrices comprised of M left-singular column 
                vectors If NULL, this agrument is assumed to be optional 
                and will not be computed
  s[N][L]       Nx1 vectors of singular values in descending order. In an 
                exceptional case when the iterative algorithm fails to 
                converge for a particular matrix, all elements of the
                respective vector are set to NaN.
  V[N*N][L]     NxN matrices comprised of N right-singular column
                vectors. If NULL, this agrument is assumed to be optional 
                and will not be computed
Restrictions:
  pScr,U,s,V,A  Must not overlap and must be aligned on 32-byte boundary 
  M,N           M>=N>1
  L             Must be a multiple of 8
---------------------------------------------------------------------------*/
/* Complex-valued functions */
void svd2x2sf (
            void * pScr,
            complex_float * restrict U,
            float32_t     * restrict s,
            complex_float * restrict V,
            complex_float * restrict A,
            int L );
void svd3x3sf (
            void * pScr,
            complex_float * restrict U,
            float32_t     * restrict s,
            complex_float * restrict V,
            complex_float * restrict A,
            int L );
void svd4x4sf (
            void * pScr,
            complex_float * restrict U,
            float32_t     * restrict s,
            complex_float * restrict V,
            complex_float * restrict A,
            int L );
void svdmxnsf (
            void * pScr,
            complex_float * restrict U,
            float32_t     * restrict s,
            complex_float * restrict V,
            complex_float * restrict A,
            int M, int N, int L );

/* Real-valued functions */
void rsvd2x2sf (
            void * pScr,
            float32_t * restrict U,
            float32_t * restrict s,
            float32_t * restrict V,
            float32_t * restrict A,
            int L );
void rsvd3x3sf (
            void * pScr,
            float32_t * restrict U,
            float32_t * restrict s,
            float32_t * restrict V,
            float32_t * restrict A,
            int L );
void rsvd4x4sf (
            void * pScr,
            float32_t * restrict U,
            float32_t * restrict s,
            float32_t * restrict V,
            float32_t * restrict A,
            int L );
void rsvdmxnsf (
            void * pScr,
            float32_t * restrict U,
            float32_t * restrict s,
            float32_t * restrict V,
            float32_t * restrict A,
            int M, int N, int L );

/* Return the scratch area size, in bytes. */
size_t svd2x2sf_getScratchSize ( int M, int N, int L );
size_t svd3x3sf_getScratchSize ( int M, int N, int L );
size_t svd4x4sf_getScratchSize ( int M, int N, int L );
size_t svdmxnsf_getScratchSize ( int M, int N, int L );

/* Return the scratch area size, in bytes. */
size_t rsvd2x2sf_getScratchSize ( int M, int N, int L );
size_t rsvd3x3sf_getScratchSize ( int M, int N, int L );
size_t rsvd4x4sf_getScratchSize ( int M, int N, int L );
size_t rsvdmxnsf_getScratchSize ( int M, int N, int L );

/*===========================================================================
Eigenvalues and eigenvectors
[c]eigenn   Eigenvalues and eigenvectors of real/complex block ordered matrices
[c]eigens   Eigenvalues and eigenvectors of real/complex stream ordered matrices
===========================================================================*/

/*---------------------------------------------------------------------------
For a complex (real) matrix A of size NxN, eigen (reigen) functions compute:
  - N eigenvalues e[N], with possible repetitions
  - N right eigenvectors V[N], normalized such that the Euclidian norm of an
    eigenvector is 1
Right eigenvectors satisfy the equation A*V[n] = s[n]*V[n], n=0..N-1.
For complex matrices, eigenvalues and eigenvectors are also complex. For real
matrices, eigenvalues and corresponding eigenvectors are either real, or come
in conjugate pairs.
-----------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
Eigenvalues And Eigenvectors Of Real/Complex Block Ordered Matrices

Description: for each complex/real input matrix A of size NxN, compute N
(possibly repeated) eigenvalues s[N], and (optonally) N right eigenvectors
of size Nx1 V[N]. Input and output data are stored in block order.

Data format: IEEE-754 Std single precision floating-point

Storage sizes SA and Se denote the number of data elements required to store
a matrix or a vector in block order. If matrix size is less than the SIMD vector
size, then the storage_size(matrix_size) equals the matrix_size rounded up to
the next power of two, otherwise it is matrix_size rounded up to the next
multiple of the SIMD vector size.

SIMD vector size:
  - for complex floating-point data 2*BBE_SIMD_WIDTH/sizeof(complex_float) == 4
  - for real floating-point data 2*BBE_SIMD_WIDTH/sizeof(float32_t) == 8

SA = storage_size(N*N)
Ss = storage_size(N)

Notes:
1. Functions may perform in-place transformations of input matrices, so that
   INPUT DATA MAY APPEAR DAMAGED after the call.
2. Once the eigenvectors are not required, set the corresponding output pointer
   V to zero, so that a lower complexity algorithm will be used.
3. Floating-point functions assume that input data are reasonably scaled. That
   is, the base-2 exponent e of the maximum absolute value over an input matrix
   belongs to the range -E<e<E, where E = 63-log2(N)/2.
4. In order to reduce the computational complexity, a preprocessing step known
   as "matrix balancing" is omitted from the implementation.

Temporary:
  pScr        Scratch area. Required size (in bytes) is defined by 
              functions [r]eigen<size>nf_getScratchSize(N,L)
Input:
  N           Matrix size
  L           Number of matrices
  A[L][SA]    NxN input matrices
Output:
  e[L][Se]    Nx1 vectors of eigenvalues. In an exceptional case when the
              iterative algorithm fails to converge for a particular matrix,
              all elements of the respective vector are set to NaN.
  V[L][SA]    NxN matrices comprised of N column eigenvectors (optional)
Restrictions:
  pScr,e,V,A  Must not overlap and must be aligned on 32-byte boundary 
  N           Must be a positive multipe of 4
---------------------------------------------------------------------------*/
/* Complex-valued functions */
void eigen8x8nf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            complex_float * restrict A,
            int L );
void eigen16x16nf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            complex_float * restrict A,
            int L );
void eigennxnnf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            complex_float * restrict A,
            int N, int L );

/* Real-valued functions */
void reigen8x8nf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            float32_t     * restrict A,
            int L );
void reigen16x16nf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            float32_t     * restrict A,
            int L );
void reigennxnnf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            float32_t     * restrict A,
            int N, int L );

/* Return the scratch area size, in bytes. */
size_t eigen8x8nf_getScratchSize   ( int N, int L );
size_t eigen16x16nf_getScratchSize ( int N, int L );
size_t eigennxnnf_getScratchSize   ( int N, int L );
/* Return the scratch area size, in bytes. */
size_t reigen8x8nf_getScratchSize   ( int N, int L );
size_t reigen16x16nf_getScratchSize ( int N, int L );
size_t reigennxnnf_getScratchSize   ( int N, int L );

/*-------------------------------------------------------------------------
Eigenvalues And Eigenvectors Of Real/Complex Stream Ordered Matrices

Description: for each complex/real input matrix A of size NxN, compute N
(possibly repeated) eigenvalues s[N], and (optonally) N right eigenvectors
of size Nx1 V[N]. Input and output data are stored in stream order.

Data format: IEEE-754 Std single precision floating-point

Notes:
1. Functions may perform in-place transformations of input matrices, so that
   INPUT DATA MAY APPEAR DAMAGED after the call.
2. Once the eigenvectors are not required, set the corresponding output pointer
   V to zero, so that a lower complexity algorithm will be used.
3. Floating-point functions assume that input data are reasonably scaled. That
   is, the base-2 exponent e of the maximum absolute value over an input matrix
   belongs to the range -E<e<E, where E = 63-log2(N)/2.
4. In order to reduce the computational complexity, a preprocessing step known
   as "matrix balancing" is omitted from the implementation.

Temporary:
  pScr        Scratch area. Required size (in bytes) is defined by 
              functions [r]eigen<size>sf_getScratchSize(N,L)
Input:
  N           Matrix size
  L           Number of matrices
  A[N*N][L]   NxN input matrices
Output:
  e[N][L]     Nx1 vectors of eigenvalues. In an exceptional case when the
              iterative algorithm fails to converge for a particular matrix,
              all elements of the respective vector are set to NaN.
  V[N*N][L]   NxN matrices comprised of N column eigenvectors (optional)
Restrictions:
  pScr,e,V,A  Must not overlap and must be aligned on 32-byte boundary 
  N           N>1
  L           Must be a multiple of 8 for real data, or a multiple of 4 for
              complex data
---------------------------------------------------------------------------*/
/* Complex-valued functions */
void eigen2x2sf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            complex_float * restrict A,
            int L );
void eigen3x3sf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            complex_float * restrict A,
            int L );
void eigen4x4sf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            complex_float * restrict A,
            int L );
void eigennxnsf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            complex_float * restrict A,
            int N, int L );

/* Real-valued functions */
void reigen2x2sf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            float32_t     * restrict A,
            int L );
void reigen3x3sf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            float32_t     * restrict A,
            int L );
void reigen4x4sf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            float32_t     * restrict A,
            int L );
void reigennxnsf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            float32_t     * restrict A,
            int N, int L );

/* Return the scratch area size, in bytes. */
size_t eigen2x2sf_getScratchSize ( int N, int L );
size_t eigen3x3sf_getScratchSize ( int N, int L );
size_t eigen4x4sf_getScratchSize ( int N, int L );
size_t eigennxnsf_getScratchSize ( int N, int L );

/* Return the scratch area size, in bytes. */
size_t reigen2x2sf_getScratchSize ( int N, int L );
size_t reigen3x3sf_getScratchSize ( int N, int L );
size_t reigen4x4sf_getScratchSize ( int N, int L );
size_t reigennxnsf_getScratchSize ( int N, int L );

#ifdef __cplusplus
};
#endif

#endif /* __NATUREDSP_BASEBAND_MATINV_H */
