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
/*          Copyright (C) 2009-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
  NatureDSP_Baseband library. Cholesky decomposition for a complex-valued pseudo-inversion:
    Apply the Cholesky decomposition to the matrix of normal equations system
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"

#if !(HAVE_VSAMATH && HAVE_NSAENX40 && 1)
DISCARD_FUN(void,chol2x2s,(
                  complex_fract16 * restrict R,
                  complex_fract16 * restrict D,
            const complex_fract16 * restrict A, 
            const int32_t * restrict sigma2,
            int L))
#else

#define VARIANT 1

/*
    Input:
    A - vector Nx40
    _6 - vector of all sixes
    _17 - vector of all 17-s
    Output:
    mag - magnitudes
    sh6,sh17 - exponents
    pd  - stored magnitudes and exponents
*/
#if VARIANT==0
#define INV_SQRT1(mag,sh6,sh17,pd,pdstep,_8,A)                  \
{                                                               \
    xb_vecNx16   b_vec, xn_vec,cn_vec;                          \
    vsaN c_vec;                                                 \
    A=BBE_ADDNX40(A,A);                                         \
    c_vec=BBE_NSAENX40(A);                                      \
    A=BBE_SLLNX40(A,c_vec);                                     \
    BBE_RSQRTLUNX40_0(A,b_vec, cn_vec, A);                 \
    BBE_MULUUSNX16( A, cn_vec,  b_vec);                        \
    /* note: not possible to use packv because */               \
    /* output is uint16_t                      */               \
    A=BBE_SRAINX40(A,23);                                       \
    xn_vec=BBE_PACKLNX40(A);                                    \
    magu=BBE_SELNX16I(xn_vec,xn_vec,BBE_SELI_INTERLEAVE_1_EVEN);\
    BBE_SVNX16_IP(magu,pd,(2*BBE_SIMD_WIDTH));                  \
    exp=BBE_MOVVVS(c_vec);                                      \
    exp=BBE_SRAINX16(exp,1);                                    \
    exp=BBE_SELNX16I(exp,exp,BBE_SELI_INTERLEAVE_1_EVEN);       \
    exp=BBE_SUBNX16(exp,_8);                                    \
    BBE_SVNX16_IP(exp,pd,pdstep);                               \
    sh6=BBE_MOVVSV(exp,0);                                      \
    sh17=BBE_SUBSAVSN(11,sh6);                               \
}
// another variant
#define INV_SQRT2(mag,sh6,sh17,pd,pdstep,_17,A)                 \
{                                                               \
    xb_vecNx16   b_vec, xn_vec,cn_vec;                          \
    vsaN c_vec;                                                 \
    A=BBE_ADDNX40(A,A);                                         \
    c_vec=BBE_NSAENX40(A);                                      \
    A=BBE_SLLNX40(A,c_vec);                                     \
    BBE_RSQRTLUNX40_0(A,b_vec, cn_vec, A);                 \
    BBE_MULUUSNX16( A, cn_vec,  b_vec);                         \
    A=BBE_SRAINX40(A,23);                                       \
    xn_vec=BBE_PACKLNX40(A);                                    \
    magu=BBE_SELNX16I(xn_vec,xn_vec,BBE_SELI_INTERLEAVE_1_EVEN);\
    BBE_SVNX16_IP(magu,pd,(2*BBE_SIMD_WIDTH));                  \
    exp=BBE_MOVVVS(c_vec);                                      \
    exp=BBE_SELNX16I(exp,exp,BBE_SELI_INTERLEAVE_1_EVEN);       \
    exp= BBE_SUBSR1RNX16(exp,_17);                              \
    BBE_SVNX16_IP(exp,pd,pdstep);                               \
    sh6=BBE_MOVVSV(exp,0);                                      \
    sh17=BBE_SUBSAVSN(11,sh6);                                  \
}
#elif VARIANT==1
#define INV_SQRT3(mag,sh6,sh17,pd,pdstep,A)                     \
{                                                               \
    xb_vecNx16   b_vec, xn_vec,cn_vec;                          \
    vsaN c_vec;                                                 \
    A=BBE_ADDNX40(A,A);                                         \
    c_vec=BBE_NSAENX40(A);                                      \
    A=BBE_SLLNX40(A,c_vec);                                     \
    BBE_RSQRTLUNX40_0(A,b_vec, cn_vec, A);                      \
    BBE_MULUUSNX16( A, cn_vec,  b_vec);                         \
    /* note: not possible to use packv because */               \
    /* output is uint16_t                      */               \
    A=BBE_SRAINX40(A,23);                                       \
    xn_vec=BBE_PACKLNX40(A);                                    \
    magu=BBE_SELNX16I(xn_vec,xn_vec,BBE_SELI_INTERLEAVE_1_EVEN);\
    BBE_SVNX16_IP(magu,pd,(2*BBE_SIMD_WIDTH));                  \
    c_vec=BBE_SHFLVSNI(c_vec,BBE_VSA_SHFLI_DUPLICATE_1_EVEN);   \
    sh17=BBE_SUBSR1SAVSN(19,c_vec);                             \
    sh6=BBE_SUBSAVSN(11,sh17);                                  \
    exp=BBE_MOVVVS(sh6);                                        \
    BBE_SVNX16_IP(exp,pd,pdstep);                               \
}

#define INV_SQRT4(mag,sh6,sh17,pd,pdstep,A)                     \
{                                                               \
    xb_vecNx16   b_vec, xn_vec,cn_vec;                          \
    vsaN c_vec;                                                 \
    A=BBE_ADDNX40(A,A);                                         \
    c_vec=BBE_NSAENX40(A);                                      \
    A=BBE_SLLNX40(A,c_vec);                                     \
    BBE_RSQRTLUNX40_0(A,b_vec, cn_vec, A);                      \
    BBE_MULUUSNX16( A, cn_vec,  b_vec);                         \
    A=BBE_SRAINX40(A,23);                                       \
    xn_vec=BBE_PACKLNX40(A);                                    \
    magu=BBE_SELNX16I(xn_vec,xn_vec,BBE_SELI_INTERLEAVE_1_EVEN);\
    BBE_SVNX16_IP(magu,pd,(2*BBE_SIMD_WIDTH));                  \
    c_vec=BBE_SHFLVSNI(c_vec,BBE_VSA_SHFLI_DUPLICATE_1_EVEN);   \
    sh17=BBE_SUBSR1SAVSN(19,c_vec);                             \
    sh6=BBE_SUBSAVSN(11,sh17);                                  \
    exp=BBE_MOVVVS(sh6);                                        \
    BBE_SVNX16_IP(exp,pd,pdstep);                               \
}
#endif

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

void chol2x2s (
                  complex_fract16 * restrict _R, 
                  complex_fract16 * restrict _D,
            const complex_fract16 * restrict _A, 
            const int32_t * restrict sigma2,
            int L)
{
    int16_t *       restrict R=(int16_t *      )_R;
    int16_t *       restrict D=(int16_t *      )_D;
    const int16_t * restrict A=(const int16_t *)_A;
  int l,_4L;
  const xb_vecNx16    * restrict pa0;
  const xb_vecNx16    * restrict ps;
  const xb_vecNx16    * restrict prr; // read only pointer for R
        xb_vecNx16    * restrict prw; // write only pointer for R
        xb_vecNx16    * restrict pd ;
  xb_vecNx40 A10,C10;
  xb_vecNx16 t0,a0,a1,y0,am0,am1;
#if VARIANT==0
  xb_vecNx16 _8,_17;
#endif
  xb_vecNx16 zero;
  vsaN sh16,sh17,sh6;
  xb_vecNx16 res,exp,a10_hi,a10_lo,magu;

  NASSERT_ALIGN(R     ,(2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(A     ,(2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(D     ,(2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(sigma2,(2*BBE_SIMD_WIDTH));
  NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
  NASSERT(L>0);

#if VARIANT==0
  _8=BBE_MOVVINT16(8);
  _17=BBE_MOVVINT16(17);
#endif
  sh16=BBE_MOVVSA32(16);

  zero = 0;
  _4L=L<<2;
  ps  = (const xb_vecNx16   *) sigma2 ;
  prw = (xb_vecNx16   *) (R);
  pa0 = (const xb_vecNx16   *) (A);
  pd  = (xb_vecNx16   *) D ;
  if (L<=BBE_SIMD_WIDTH*2)
  {
    //--------------------------------------------
    // regular algorithm: all computations in one 
    // loop: preferrable on small L
    //--------------------------------------------
    __Pragma("loop_count min=1");
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        // m=0
        BBE_LVNX16_IP(t0,ps,0);    // load sigma and compute diagonal element
        A10 = BBE_MOVSWV(t0,t0);
        BBE_LVNX16_XP(a0,pa0,4*L);
        BBE_LVNX16_XP(am0,pa0,4*L);
        BBE_LVNX16_XP(a1,pa0,4*L);
        BBE_LVNX16_XP(am1,pa0,-3*4*L+2*BBE_SIMD_WIDTH);

        A10 = BBE_MOVSWV(t0,t0);
        BBE_MAGIANX16C(A10,a0,a0);
        BBE_MAGIANX16C(A10,a1,a1);

#if VARIANT==0
        INV_SQRT2(magu,sh6,sh17,pd,(2*BBE_SIMD_WIDTH), _17,A10)
#elif VARIANT==1
        INV_SQRT4(magu,sh6,sh17,pd,(2*BBE_SIMD_WIDTH), A10)
#endif
        A10=BBE_MULNX16J(a0,a0);
        BBE_MULANX16J(A10,a1,a1);
        a10_lo=BBE_PACKVNX40(A10,sh17);
        A10=BBE_MULUUNX16(a10_lo,magu);
        res=BBE_PACKQNX40(A10);
        BBE_SVNX16_XP(res , prw, 4*L);   // r(0,0)

        A10=BBE_MULNX16J (am0,a0);
        BBE_MULANX16J(A10,am1,a1);
        A10=BBE_SLSNX40(A10,sh6);
        a10_lo=BBE_PACKLNX40(A10);
        a10_hi=BBE_PACKVNX40(A10,sh16);
        A10=BBE_MULUUNX16(a10_lo,magu);
        A10=BBE_SRAINX40(A10,16);
        BBE_MULUSANX16(A10,magu,a10_hi);
        y0=BBE_PACKPNX40(A10);
        #ifdef COMPILER_XTENSA
        #pragma no_reorder
        #endif
        BBE_SVNX16_XP(y0  , prw, 4*L);   // r(0,1)
        BBE_SVNX16_XP(zero, prw, 4*L);   // r(1,0);

        // m=1;
        BBE_LVNX16_IP(t0,ps,(2*BBE_SIMD_WIDTH));    // reload sigma and compute diagonal element
        A10 = BBE_MOVSWV(t0,t0);
        BBE_MAGIANX16C(A10,am0,am0);
        BBE_MAGIANX16C(A10,am1,am1);
        C10=BBE_MAGINX16C(y0,y0);
        A10 = BBE_SUBNX40(A10,C10);
#if VARIANT==0
        INV_SQRT2(magu,sh6,sh17,pd,(2*BBE_SIMD_WIDTH), _17,A10)
#elif VARIANT==1
        INV_SQRT4(magu,sh6,sh17,pd,(2*BBE_SIMD_WIDTH), A10);
#endif

        // compute element r11
        A10=BBE_MULNX16J(am0,am0);
        BBE_MULANX16J(A10,am1,am1);
        BBE_MULSNX16J(A10 ,y0,y0);
        a10_lo=BBE_PACKVNX40(A10,sh17);
        A10=BBE_MULUUNX16(a10_lo,magu);
        res=BBE_PACKQNX40(A10);

        BBE_SVNX16_XP(res, prw, -3*4*L+2*BBE_SIMD_WIDTH); // r(1,1)
    }
  }
  else
  {
    //--------------------------------------------
    // Alternative algo: splitted by 3 loops
    // beneficial on bigger L
    //--------------------------------------------
    __Pragma("loop_count min=1");
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        // m=0
        BBE_LVNX16_IP(t0,ps,2*BBE_SIMD_WIDTH);    // load sigma and compute diagonal element
        a1=BBE_LVNX16_X (pa0, 2*_4L);
        BBE_LVNX16_IP(a0,pa0, 2*BBE_SIMD_WIDTH);
    #if 1
        // longer but better
        A10 = BBE_MOVSWV(t0,t0);
        BBE_MAGIANX16C(A10,a1,a1);
        BBE_MAGIANX16C(A10,a0,a0);
#if VARIANT==0
        INV_SQRT2(magu,sh6,sh17,pd ,3*2*BBE_SIMD_WIDTH, _17,A10)
#elif VARIANT==1
        INV_SQRT4(magu,sh6,sh17,pd ,3*2*BBE_SIMD_WIDTH,A10);
#endif
        A10=BBE_MULNX16J (a1,a1);
        BBE_MULANX16J(A10,a0,a0);
    #else
        C10 = BBE_MOVSWV(t0,t0);
        A10=BBE_MULNX16J(a0,a0);
        BBE_MULANX16J(A10,a1,a1);
        C10=BBE_ADDNX40(A10,C10);
        INV_SQRT(magu,sh6,sh17,pd,3*2*BBE_SIMD_WIDTH, _6,_17,C10)
    #endif
        a10_lo=BBE_PACKVNX40(A10,sh17);
        A10=BBE_MULUUNX16(a10_lo,magu);
        res=BBE_PACKQNX40(A10);
        BBE_SVNX16_IP(res, prw, 2*BBE_SIMD_WIDTH);   // r(0,0)
    }
    // compute // r(0,1),r(1,0)
    pd  = (      xb_vecNx16 *)(D);
    pa0 = (const xb_vecNx16 *)(A);
    prw = (      xb_vecNx16 *)(R+1*2*L);
    
    __Pragma("loop_count min=1");
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        xb_vecNx16 exp;
        BBE_LVNX16_XP(a0 ,pa0, _4L);
        BBE_LVNX16_XP(am0,pa0, _4L);
        BBE_LVNX16_XP(a1 ,pa0, _4L);
        BBE_LVNX16_XP(am1,pa0,-3*_4L+2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(magu,pd,  2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(exp ,pd,3*2*BBE_SIMD_WIDTH);
        sh6=BBE_MOVVSV(exp,0);
        A10=BBE_MULNX16J (am0,a0);
        BBE_MULANX16J(A10,am1,a1);
        A10=BBE_SLSNX40(A10,sh6);
        a10_lo=BBE_PACKLNX40(A10);
        a10_hi=BBE_PACKVNX40(A10,sh16);
        A10=BBE_MULUUNX16(a10_lo,magu);
        A10=BBE_SRAINX40(A10,16);
        BBE_MULUSANX16(A10,magu,a10_hi);
        y0=BBE_PACKPNX40(A10);
        BBE_SVNX16_X (zero, prw, _4L);                // r(1,0)
        BBE_SVNX16_IP(y0  , prw, 2*BBE_SIMD_WIDTH);   // r(0,1)
    }

    // compute r(1,1)
    pa0 = (const xb_vecNx16   *) (A+2*L);
    ps  = (const xb_vecNx16   *) sigma2 ;
    pd  = (xb_vecNx16   *) (D+2*BBE_SIMD_WIDTH) ;
    prw = (xb_vecNx16   *) (R+3*2*L);
    prr = (xb_vecNx16   *) (R+1*2*L);

    __Pragma("loop_count min=1");
    for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
    {
        BBE_LVNX16_IP(y0, prr, 2*BBE_SIMD_WIDTH);   // r(0,1)
        am1=BBE_LVNX16_X (pa0, 2*_4L);
        BBE_LVNX16_IP(am0,pa0, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(t0 ,ps , 2*BBE_SIMD_WIDTH);    // reload sigma and compute diagonal element
    #if 1   
        // longer but better
        A10 = BBE_MOVSWV(t0,t0);
        BBE_MAGIANX16C(A10,am1,am1);
        BBE_MAGIANX16C(A10,am0,am0);
        C10=BBE_MAGINX16C(y0,y0);
        A10 = BBE_SUBNX40(A10,C10);
#if VARIANT==0
        INV_SQRT1(magu,sh6,sh17,pd,3*2*BBE_SIMD_WIDTH, _8,A10)
#else
        INV_SQRT3(magu,sh6,sh17,pd,3*2*BBE_SIMD_WIDTH, A10)
#endif
        // compute element r11
        A10=BBE_MULNX16J (am1,am1);
        BBE_MULANX16J(A10,am0,am0);
        BBE_MULSNX16J(A10,y0 ,y0 );
    #else
        A10=BBE_MULNX16J(am0,am0);
        BBE_MULANX16J(A10,am1,am1);
        BBE_MULSNX16J(A10 ,y0,y0);
        C10 = BBE_MOVSWV(t0,t0);
        C10 = BBE_ADDNX40(A10,C10);
        INV_SQRT(magu,sh6,sh17,pd,3*2*BBE_SIMD_WIDTH, _6,_17,C10)
    #endif
        // compute element r11
        a10_lo=BBE_PACKVNX40(A10,sh17);
        A10=BBE_MULUUNX16(a10_lo,magu);
        res=BBE_PACKQNX40(A10);
        BBE_SVNX16_IP(res, prw, 2*BBE_SIMD_WIDTH); 
    }
  }
} /* chol2x2s() */
#endif
