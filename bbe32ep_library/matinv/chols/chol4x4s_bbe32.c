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
DISCARD_FUN(void,chol4x4s,(
                  complex_fract16 * restrict R,
                  complex_fract16 * restrict D,
            const complex_fract16 * restrict A, 
    const int32_t * restrict sigma2,
    int L))
#else

/*
Reference Matlab code:
function R=chol2(A,sigma2)
sz=size(A); M=sz(1); N=sz(2);
R=zeros(N,N);
for m=1:N
Rm=R(:,m);  % take m-th column of original and decomposing matrix
Am=A(:,m);
Amm=Am'*Am+sigma2;
Rmm=Rm'*Rm;
Rmm=sqrt(real(Amm-Rmm));
x(1,1:m)=[zeros(1,m-1) Rmm];
for k=m+1:N
Akm=A(:,k)'*Am;
Rkm=R(:,k)'*Rm;
x(1,k)=(Akm-Rkm)/Rmm;
end
R(m,:)=conj(x);
end
*/

#define VARIANT 4   // 2 or 4

#define INN0(_pa)                   \
  BBE_LVNX16_XP(am0,_pa, 2*8*L);    \
  A10=BBE_MULNX16J(am0,a0);         \
  BBE_LVNX16_XP(am1,_pa, 2*8*L );   \
  BBE_MULANX16J(A10,am1,a1);        \
  BBE_LVNX16_XP(am2,_pa, 2*8*L );   \
  BBE_MULANX16J(A10,am2,a2);        \
  BBE_LVNX16_XP(am3, _pa, -2*22*L); \
  BBE_MULANX16J(A10,am3,a3)


#define SAV0(_y,_A10,_sh6,_sh16,_magu)\
{                                    \
    xb_vecNx16 a10_lo,a10_hi;        \
    _A10=BBE_SLSNX40(_A10,_sh6);     \
    a10_lo=BBE_PACKLNX40(_A10);      \
    a10_hi=BBE_PACKVNX40(_A10,_sh16);\
    _A10=BBE_MULUUNX16(a10_lo,_magu);\
    _A10=BBE_SRAINX40(_A10,16);      \
    BBE_MULUSANX16(_A10,magu,a10_hi);\
    _y=BBE_PACKPNX40(_A10);          \
}

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
    sh17=BBE_SUBSAVSN(11,sh6);                               \
}

// 2 different variants of calculation the main diagonal element
#define DIAG1(/*result:*/y0,                               \
              /*internals:*/sh6,sh17,pd,_17,t0,            \
              /*a0...a3:*/a,/*r0...r2:*/r,                 \
              /* iteration number (0...3):*/iteration)     \
{                                                          \
    xb_vecNx40 A,C;                                        \
    xb_vecNx16 a10_lo,magu;                                \
    A=BBE_MULNX16J (a##0,a##0);                            \
    BBE_MULANX16J(A,a##1,a##1);                            \
    BBE_MULANX16J(A,a##2,a##2);                            \
    BBE_MULANX16J(A,a##3,a##3);                            \
    if(iteration>0)    BBE_MULSNX16J(A,r##0,r##0);         \
    if(iteration>1)    BBE_MULSNX16J(A,r##1,r##1);         \
    if(iteration>2)    BBE_MULSNX16J(A,r##2,r##2);         \
    C = BBE_MOVSWV(t0,t0);                                 \
    C = BBE_ADDNX40(C,A);                                  \
    INV_SQRT2(magu,sh6,sh17,pd,7*2*BBE_SIMD_WIDTH,_17,C);  \
    a10_lo=BBE_PACKVNX40(A,sh17);                          \
    A=BBE_MULUUNX16(a10_lo,magu);                          \
    y0=BBE_PACKQNX40(A);                                   \
}

#define DIAG2(/*result:*/y0,                               \
              /*internals:*/sh6,sh17,pd,_17,t0,            \
              /*a0...a3:*/a,/*r0...r2:*/r,                 \
              /* iteration number (0...3):*/iteration)     \
{                                                          \
    xb_vecNx40 A,C;                                        \
    xb_vecNx16 a10_lo,magu;                                \
    C = BBE_MOVSWV(t0,t0);                                 \
    BBE_MAGIANX16C(C,a##0,a##0);                          \
    BBE_MAGIANX16C(C,a##1,a##1);                          \
    BBE_MAGIANX16C(C,a##2,a##2);                          \
    BBE_MAGIANX16C(C,a##3,a##3);                          \
    if(iteration>0)    A=BBE_MAGINX16C (r##0,r##0);       \
    if(iteration>1)    BBE_MAGIANX16C(A,r##1,r##1);       \
    if(iteration>2)    BBE_MAGIANX16C(A,r##2,r##2);       \
    if(iteration>0)    C = BBE_SUBNX40(C,A);               \
    INV_SQRT2(magu,sh6,sh17,pd,7*2*BBE_SIMD_WIDTH,_17,C);  \
    A=BBE_MULNX16J (a##0,a##0);                            \
    BBE_MULANX16J(A,a##1,a##1);                            \
    BBE_MULANX16J(A,a##2,a##2);                            \
    BBE_MULANX16J(A,a##3,a##3);                            \
    if(iteration>0)    BBE_MULSNX16J(A,r##0,r##0);         \
    if(iteration>1)    BBE_MULSNX16J(A,r##1,r##1);         \
    if(iteration>2)    BBE_MULSNX16J(A,r##2,r##2);         \
    a10_lo=BBE_PACKVNX40(A,sh17);                          \
    A=BBE_MULUUNX16(a10_lo,magu);                          \
    y0=BBE_PACKQNX40(A);                                   \
}

// another variants
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

#define DIAG3(/*result:*/y0,                               \
              /*internals:*/sh6,sh17,pd,t0,                \
              /*a0...a3:*/a,/*r0...r2:*/r,                 \
              /* iteration number (0...3):*/iteration)     \
{                                                          \
    xb_vecNx40 A,C;                                        \
    xb_vecNx16 a10_lo,magu;                                \
    A=BBE_MULNX16J (a##0,a##0);                            \
    BBE_MULANX16J(A,a##1,a##1);                            \
    BBE_MULANX16J(A,a##2,a##2);                            \
    BBE_MULANX16J(A,a##3,a##3);                            \
    if(iteration>0)    BBE_MULSNX16J(A,r##0,r##0);         \
    if(iteration>1)    BBE_MULSNX16J(A,r##1,r##1);         \
    if(iteration>2)    BBE_MULSNX16J(A,r##2,r##2);         \
    C = BBE_MOVSWV(t0,t0);                                 \
    C = BBE_ADDNX40(C,A);                                  \
    INV_SQRT4(magu,sh6,sh17,pd,7*2*BBE_SIMD_WIDTH,C);      \
    a10_lo=BBE_PACKVNX40(A,sh17);                          \
    A=BBE_MULUUNX16(a10_lo,magu);                          \
    y0=BBE_PACKQNX40(A);                                   \
}

#define DIAG4(/*result:*/y0,                               \
              /*internals:*/sh6,sh17,pd,t0,                \
              /*a0...a3:*/a,/*r0...r2:*/r,                 \
              /* iteration number (0...3):*/iteration)     \
{                                                          \
    xb_vecNx40 A,C;                                        \
    xb_vecNx16 a10_lo,magu;                                \
    C = BBE_MOVSWV(t0,t0);                                 \
    BBE_MAGIANX16C(C,a##0,a##0);                           \
    BBE_MAGIANX16C(C,a##1,a##1);                           \
    BBE_MAGIANX16C(C,a##2,a##2);                           \
    BBE_MAGIANX16C(C,a##3,a##3);                           \
    if(iteration>0)    A=BBE_MAGINX16C (r##0,r##0);        \
    if(iteration>1)    BBE_MAGIANX16C(A,r##1,r##1);        \
    if(iteration>2)    BBE_MAGIANX16C(A,r##2,r##2);        \
    if(iteration>0)    C = BBE_SUBNX40(C,A);               \
    INV_SQRT4(magu,sh6,sh17,pd,7*2*BBE_SIMD_WIDTH,C);      \
    A=BBE_MULNX16J (a##0,a##0);                            \
    BBE_MULANX16J(A,a##1,a##1);                            \
    BBE_MULANX16J(A,a##2,a##2);                            \
    BBE_MULANX16J(A,a##3,a##3);                            \
    if(iteration>0)    BBE_MULSNX16J(A,r##0,r##0);         \
    if(iteration>1)    BBE_MULSNX16J(A,r##1,r##1);         \
    if(iteration>2)    BBE_MULSNX16J(A,r##2,r##2);         \
    a10_lo=BBE_PACKVNX40(A,sh17);                          \
    A=BBE_MULUUNX16(a10_lo,magu);                          \
    y0=BBE_PACKQNX40(A);                                   \
}

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

void chol4x4s (
                  complex_fract16 * restrict _R, 
                  complex_fract16 * restrict _D,
            const complex_fract16 * restrict _A, 
            const int32_t * restrict sigma2,
            int L)
#if 0
{
    int offset[1];
    int incw;
    int idx;
    int _16L;   // 16*L
    int _44L;   // -44*L
    int _4L;    // 4*L
    xb_vecNx40 A10,C10;
    xb_vecNx16 t0,a0,a1,a2,a3,r0,r1,r2,y0,am0,am1,am2,am3,_6,_17,y1,y2,magu;
    xb_vecNx16 yx;
    vsaN sh16,sh17,sh6;
    xb_vecNx16 exp,a10_lo;
    const xb_vecNx16    * restrict pa = (xb_vecNx16   *) A ;
    const xb_vecNx16    * restrict ps = (xb_vecNx16   *) sigma2 ;
    xb_vecNx16    * restrict pd = (xb_vecNx16   *) D ;
    const xb_vecNx16    * restrict pa1 = (xb_vecNx16   *) A ;
    const xb_vecNx16    * restrict prr; // read only pointer for R
    xb_vecNx16    * restrict prw; // write only pointer for R

    int l;
    NASSERT_ALIGN(R     ,*2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(A     ,*2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D     ,*2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(sigma2,*2*BBE_SIMD_WIDTH));
    NASSERT((L%(BBE_SIMD_WIDTH/2)==0);
    NASSERT(L>0);

    _6=BBE_MOVVINT16(6);
    _17=BBE_MOVVINT16(17);
    sh16=BBE_MOVVSA32(16);
    offset[0]=L<<3;
    // make parallel processing by 8 matrices
    prr=(const xb_vecNx16 *)R;
    prw=(xb_vecNx16 *)R;
    pa1=(const xb_vecNx16 *)A;
    _16L =L<<4;
    _44L=-3*2*2*4*L+2*2*L;
    _4L =4*L;
    l=L>>4;
    do
    {
        incw = _4L;
        BBE_LVNX16_IP(t0,ps,0);    // load sigma and compute diagonal element
        // m=0
        BBE_LVNX16_XP(a0,pa1,_16L);
        BBE_LVNX16_XP(a1,pa1,_16L);
        BBE_LVNX16_XP(a2,pa1,_16L);
        BBE_LVNX16_XP(a3,pa1,_44L);
        pa = pa1;
        A10=BBE_MULNX16J (a0,a0);
        BBE_MULANX16J(A10,a1,a1);
        BBE_MULANX16J(A10,a2,a2);
        BBE_MULANX16J(A10,a3,a3);
        C10 = BBE_MOVSWV(t0,t0);
        C10 = BBE_ADDNX40(C10,A10);
        //INV_SQRT(magu,sh6,sh17,pd, _6,_17,C10);
        INV_SQRT2(magu,sh6,sh17,pd,2*BBE_SIMD_WIDTH,_17,C10);
        a10_lo=BBE_PACKVNX40(A10,sh17);
        A10=BBE_MULUUNX16(a10_lo,magu);
        y0=BBE_PACKQNX40(A10);
        BBE_SVNX16_XP(y0, prw, _16L);  // diagonal element r(0,0)
        // write zeros in r(1,0),r(2,0),r(3,0)
        {
            xb_vecNx16 zero=0;
            BBE_SVNX16_XP(zero, prw, _16L); 
            BBE_SVNX16_XP(zero, prw, _16L); 
            BBE_SVNX16_XP(zero, prw, _44L); 
        }
        // compute elements r(0,1)...r(0,3)
        BBE_LVNX16_XP(am0 ,pa , _16L);
        BBE_LVNX16_XP(am1 ,pa , _16L);
        BBE_LVNX16_XP(am2 ,pa , _16L);
        BBE_LVNX16_XP(am3 ,pa , _44L);
        A10=BBE_MULNX16J (am0,a0);
        BBE_MULANX16J(A10,am1,a1);
        BBE_MULANX16J(A10,am2,a2);
        BBE_MULANX16J(A10,am3,a3);
        SAV0(yx,A10,sh6,sh16,magu);
        BBE_SVNX16_XP(yx, prw, _4L); // r(0,1)
        BBE_LVNX16_XP(am0 ,pa , _16L);
        BBE_LVNX16_XP(am1 ,pa , _16L);
        BBE_LVNX16_XP(am2 ,pa , _16L);
        BBE_LVNX16_XP(am3 ,pa , _44L);
        A10=BBE_MULNX16J (am0,a0);
        BBE_MULANX16J(A10,am1,a1);
        BBE_MULANX16J(A10,am2,a2);
        BBE_MULANX16J(A10,am3,a3);
        SAV0(y1,A10,sh6,sh16,magu);
        BBE_LVNX16_XP(am0 ,pa , _16L);
        BBE_LVNX16_XP(am1 ,pa , _16L);
        BBE_LVNX16_XP(am2 ,pa , _16L);
        BBE_LVNX16_XP(am3 ,pa , _44L);
        A10=BBE_MULNX16J (am0,a0);
        BBE_MULANX16J(A10,am1,a1);
        BBE_MULANX16J(A10,am2,a2);
        BBE_MULANX16J(A10,am3,a3);
        SAV0(y0,A10,sh6,sh16,magu);
        // m=1
        BBE_LVNX16_XP(a0,pa1, _16L);
        BBE_LVNX16_XP(a1,pa1, _16L);
        BBE_LVNX16_XP(a2,pa1, _16L);
        BBE_LVNX16_XP(a3,pa1, _44L);
        pa = pa1;
        A10=BBE_MULNX16J (a0,a0);
        BBE_MULANX16J(A10,a1,a1);
        BBE_MULANX16J(A10,a2,a2);
        BBE_MULANX16J(A10,a3,a3);
        r0=yx; // r(0,1)
        BBE_LVNX16_IP(t0,ps,0);    // reload sigma and compute diagonal element
        BBE_MULSNX16J(A10,r0,r0);
        C10 = BBE_MOVSWV(t0,t0);
        C10 = BBE_ADDNX40(C10,A10);
        //INV_SQRT(magu,sh6,sh17,pd, _6,_17,C10);
        INV_SQRT2(magu,sh6,sh17,pd,2*BBE_SIMD_WIDTH,_17,C10);
        a10_lo=BBE_PACKVNX40(A10,sh17);
        A10=BBE_MULUUNX16(a10_lo,magu);
        yx=BBE_PACKQNX40(A10);
        // compute elements r(1,2),r(1,3)
        BBE_LVNX16_XP(am0 ,pa , _16L);
        BBE_LVNX16_XP(am1 ,pa , _16L);
        BBE_LVNX16_XP(am2 ,pa , _16L);
        BBE_LVNX16_XP(am3 ,pa , _44L);
        A10=BBE_MULNX16J (am0,a0);
        BBE_MULANX16J(A10,am1,a1);
        BBE_MULANX16J(A10,am2,a2);
        BBE_MULANX16J(A10,am3,a3);
        am0=y1; // r(0,2)
        BBE_MULSNX16J(A10,am0,r0);
        BBE_SVNX16_XP(y1, prw, _4L); // r(0,2)
        incw+=_4L;
        BBE_SVNX16_XP(y0, prw, incw); // r(0,3)
        // write zeros in r(2,1),r(3,1)
        {
            xb_vecNx16 zero=0;
            BBE_SVNX16_X (zero, prw, _16L); 
            BBE_SVNX16_X (zero, prw, (2*16)*L); 
        }
        BBE_SVNX16_XP(yx, prw, _4L);  // diagonal element r(1,1)
        SAV0(yx,A10,sh6,sh16,magu);
        BBE_LVNX16_XP(am0 ,pa , _16L);
        BBE_LVNX16_XP(am1 ,pa , _16L);
        BBE_LVNX16_XP(am2 ,pa , _16L);
        BBE_LVNX16_XP(am3 ,pa , _44L);
        A10=BBE_MULNX16J (am0,a0);
        BBE_MULANX16J(A10,am1,a1);
        BBE_MULANX16J(A10,am2,a2);
        BBE_MULANX16J(A10,am3,a3);
        am0=y0; // r(0,3)
        y1=y0;
        BBE_MULSNX16J(A10,am0,r0);
        SAV0(y0,A10,sh6,sh16,magu);
        // m=2
        BBE_LVNX16_XP(a0,pa1, _16L);
        BBE_LVNX16_XP(a1,pa1, _16L);
        BBE_LVNX16_XP(a2,pa1, _16L);
        BBE_LVNX16_XP(a3,pa1, _44L);
        pa = pa1;
        A10=BBE_MULNX16J (a0,a0);
        BBE_MULANX16J(A10,a1,a1);
        BBE_MULANX16J(A10,a2,a2);
        BBE_MULANX16J(A10,a3,a3);
        idx=XT_L32I_N(offset,0);
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
        r0=BBE_LVNX16_X (prr,idx); // r(0,2)
        r1=yx; // r(1,2)
        BBE_SVNX16_XP(yx, prw, _4L);  // r(1,2)
        BBE_MULSNX16J(A10,r0,r0);
        BBE_MULSNX16J(A10,r1,r1);
        BBE_LVNX16_IP(t0,ps,(2*BBE_SIMD_WIDTH));    // reload sigma and compute diagonal element
        C10 = BBE_MOVSWV(t0,t0);
        C10 = BBE_ADDNX40(C10,A10);
        //INV_SQRT(magu,sh6,sh17,pd, _6,_17,C10);
        INV_SQRT2(magu,sh6,sh17,pd,2*BBE_SIMD_WIDTH,_17,C10);
        a10_lo=BBE_PACKVNX40(A10,sh17);
        A10=BBE_MULUUNX16(a10_lo,magu);
        yx=BBE_PACKQNX40(A10);
        // compute element r(2,3)
        BBE_LVNX16_XP(am0 ,pa , _16L);
        BBE_LVNX16_XP(am1 ,pa , _16L);
        BBE_LVNX16_XP(am2 ,pa , _16L);
        BBE_LVNX16_XP(am3 ,pa , _44L);
        A10=BBE_MULNX16J (am0,a0);
        BBE_MULANX16J(A10,am1,a1);
        BBE_MULANX16J(A10,am2,a2);
        BBE_MULANX16J(A10,am3,a3);
        am0=y1; // r(0,3)
        am1=y0; // r(1,3)
        y2=y0;
        incw+=_4L;
        BBE_SVNX16_XP(y0, prw, incw);   // r(1,3)
        BBE_MULSNX16J(A10,am0,r0);
        BBE_MULSNX16J(A10,am1,r1);
        SAV0(y0,A10,sh6,sh16,magu);
        //        m==3
        BBE_LVNX16_XP(a0,pa1, _16L);
        BBE_LVNX16_XP(a1,pa1, _16L);
        BBE_LVNX16_XP(a2,pa1, _16L);
        BBE_LVNX16_XP(a3,pa1,  _44L);
        pa=pa1;
        A10=BBE_MULNX16J (a0,a0);
        BBE_MULANX16J(A10,a1,a1);
        BBE_MULANX16J(A10,a2,a2);
        BBE_MULANX16J(A10,a3,a3);
        r0=y1; // r(0,3)
        r1=y2; // r(1,3)
        r2=y0; // r(2,3)
        BBE_MULSNX16J(A10,r0,r0);
        BBE_MULSNX16J(A10,r1,r1);
        BBE_MULSNX16J(A10,r2,r2);
        C10 = BBE_MOVSWV(t0,t0);
        C10 = BBE_ADDNX40(C10,A10);
        //INV_SQRT(magu,sh6,sh17,pd, _6,_17,C10);
        INV_SQRT2(magu,sh6,sh17,pd,2*BBE_SIMD_WIDTH,_17,C10);
        a10_lo=BBE_PACKVNX40(A10,sh17);
        A10=BBE_MULUUNX16(a10_lo,magu);
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
        // write zeros in r(3,2)
        {
            xb_vecNx16 zero=0;
            BBE_SVNX16_X(zero, prw, _16L); 
        }
        BBE_SVNX16_XP(yx, prw, _4L);  // diagonal element r(2,2)
        incw+=_4L;
        BBE_SVNX16_XP(y0, prw, incw); // r(2,3)
        y0=BBE_PACKQNX40(A10);
        BBE_SVNX16_XP(y0, prw, _4L);  // diagonal element r(3,3)
        // go to next matrices
        pa1 = (xb_vecNx16*)(((uintptr_t)pa1)-_16L);
        BBE_LVNX16_IP(t0,pa1,(2*BBE_SIMD_WIDTH));
        BBE_LVNX16_IP(t0,prr,(2*BBE_SIMD_WIDTH));
        prw=(xb_vecNx16*)prr;
    }
    while(--l);

    return 0;
} /* chol4x4s() */
#else
{
    int16_t *       restrict R=(int16_t *      )_R;
    int16_t *       restrict D=(int16_t *      )_D;
    const int16_t * restrict A=(const int16_t *)_A;
    xb_vecNx16 zero = 0;
    int incw;
    int _16L;   // 16*L
    int _44L;   // -44*L
    int _4L;    // 4*L
    xb_vecNx40 A10, C10;
    xb_vecNx16 t0, a0, a1, a2, a3, r0, r1, r2, y0, r3, am0, am1, am2, am3, y1, magu;
#if VARIANT==2
    xb_vecNx16 _17;
#elif VARIANT==4
#endif
    xb_vecNx16 yx;
    vsaN sh16, sh17, sh6;
    xb_vecNx16 exp, a10_lo;
    const xb_vecNx16    * restrict pa;
    const xb_vecNx16    * restrict pa1;
    const xb_vecNx16    * restrict ps;
    xb_vecNx16    * restrict pd;
    const xb_vecNx16    * restrict prr; // read only pointer for R
    xb_vecNx16    * restrict prw; // write only pointer for R

    int l;
    NASSERT_ALIGN(R, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(A, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(sigma2, (2 * BBE_SIMD_WIDTH));
    NASSERT(L % (BBE_SIMD_WIDTH / 2) == 0);
    NASSERT(L > 0);

#if VARIANT==2
    _17=BBE_MOVVINT16(17);
#endif
    sh16 = BBE_MOVVSA32(16);
    _16L = L << 4;
    _4L = L << 2;
    _44L = -11 * _4L;
    if (L <= BBE_SIMD_WIDTH)
    {
        //====================================================================================
        // algorithm for smaller L
        //====================================================================================
        int idx, offset[1];
        xb_vecNx16 y2;
        offset[0] = L << 3;
        // make parallel processing by 8 matrices
        prr = (const xb_vecNx16 *)R;
        prw = (xb_vecNx16 *)R;
        pa = pa1 = (const xb_vecNx16 *)A;
        ps = (xb_vecNx16   *)sigma2;
        pd = (xb_vecNx16   *)D;
        _16L = L << 4;
        _44L = -3 * 2 * 2 * 4 * L + 2 * 2 * L;
        _4L = 4 * L;
        l = L >> (LOG2_BBE_SIMD_WIDTH - 1);
        do
        {
            incw = _4L;
            BBE_LVNX16_IP(t0, ps, 0);    // load sigma and compute diagonal element
            // m=0
            BBE_LVNX16_XP(a0, pa1, _16L);
            BBE_LVNX16_XP(a1, pa1, _16L);
            BBE_LVNX16_XP(a2, pa1, _16L);
            BBE_LVNX16_XP(a3, pa1, _44L);
            pa = pa1;
            A10 = BBE_MULNX16J(a0, a0);
            BBE_MULANX16J(A10, a1, a1);
            BBE_MULANX16J(A10, a2, a2);
            BBE_MULANX16J(A10, a3, a3);
            C10 = BBE_MOVSWV(t0, t0);
            C10 = BBE_ADDNX40(C10, A10);
#if VARIANT==2
            INV_SQRT2(magu,sh6,sh17,pd,2*BBE_SIMD_WIDTH,_17,C10);
#elif VARIANT==4
            INV_SQRT4(magu, sh6, sh17, pd, 2 * BBE_SIMD_WIDTH, C10);
#endif
            a10_lo = BBE_PACKVNX40(A10, sh17);
            A10 = BBE_MULUUNX16(a10_lo, magu);
            y0 = BBE_PACKQNX40(A10);
            BBE_SVNX16_XP(y0, prw, _16L);  // diagonal element r(0,0)
            // write zeros in r(1,0),r(2,0),r(3,0)
            {
                xb_vecNx16 zero = 0;
                BBE_SVNX16_XP(zero, prw, _16L);
                BBE_SVNX16_XP(zero, prw, _16L);
                BBE_SVNX16_XP(zero, prw, _44L);
            }
            // compute elements r(0,1)...r(0,3)
            BBE_LVNX16_XP(am0, pa, _16L);
            BBE_LVNX16_XP(am1, pa, _16L);
            BBE_LVNX16_XP(am2, pa, _16L);
            BBE_LVNX16_XP(am3, pa, _44L);
            A10 = BBE_MULNX16J(am0, a0);
            BBE_MULANX16J(A10, am1, a1);
            BBE_MULANX16J(A10, am2, a2);
            BBE_MULANX16J(A10, am3, a3);
            SAV0(yx, A10, sh6, sh16, magu);
            BBE_SVNX16_XP(yx, prw, _4L); // r(0,1)
            BBE_LVNX16_XP(am0, pa, _16L);
            BBE_LVNX16_XP(am1, pa, _16L);
            BBE_LVNX16_XP(am2, pa, _16L);
            BBE_LVNX16_XP(am3, pa, _44L);
            A10 = BBE_MULNX16J(am0, a0);
            BBE_MULANX16J(A10, am1, a1);
            BBE_MULANX16J(A10, am2, a2);
            BBE_MULANX16J(A10, am3, a3);
            SAV0(y1, A10, sh6, sh16, magu);
            BBE_LVNX16_XP(am0, pa, _16L);
            BBE_LVNX16_XP(am1, pa, _16L);
            BBE_LVNX16_XP(am2, pa, _16L);
            BBE_LVNX16_XP(am3, pa, _44L);
            A10 = BBE_MULNX16J(am0, a0);
            BBE_MULANX16J(A10, am1, a1);
            BBE_MULANX16J(A10, am2, a2);
            BBE_MULANX16J(A10, am3, a3);
            SAV0(y0, A10, sh6, sh16, magu);
            // m=1
            BBE_LVNX16_XP(a0, pa1, _16L);
            BBE_LVNX16_XP(a1, pa1, _16L);
            BBE_LVNX16_XP(a2, pa1, _16L);
            BBE_LVNX16_XP(a3, pa1, _44L);
            pa = pa1;
            A10 = BBE_MULNX16J(a0, a0);
            BBE_MULANX16J(A10, a1, a1);
            BBE_MULANX16J(A10, a2, a2);
            BBE_MULANX16J(A10, a3, a3);
            r0 = yx; // r(0,1)
            BBE_LVNX16_IP(t0, ps, 0);    // reload sigma and compute diagonal element
            BBE_MULSNX16J(A10, r0, r0);
            C10 = BBE_MOVSWV(t0, t0);
            C10 = BBE_ADDNX40(C10, A10);
            //INV_SQRT(magu,sh6,sh17,pd, _6,_17,C10);
#if VARIANT==2
            INV_SQRT2(magu,sh6,sh17,pd,2*BBE_SIMD_WIDTH,_17,C10);
#elif VARIANT==4
            INV_SQRT4(magu, sh6, sh17, pd, 2 * BBE_SIMD_WIDTH, C10);
#endif
            a10_lo = BBE_PACKVNX40(A10, sh17);
            A10 = BBE_MULUUNX16(a10_lo, magu);
            yx = BBE_PACKQNX40(A10);
            // compute elements r(1,2),r(1,3)
            BBE_LVNX16_XP(am0, pa, _16L);
            BBE_LVNX16_XP(am1, pa, _16L);
            BBE_LVNX16_XP(am2, pa, _16L);
            BBE_LVNX16_XP(am3, pa, _44L);
            A10 = BBE_MULNX16J(am0, a0);
            BBE_MULANX16J(A10, am1, a1);
            BBE_MULANX16J(A10, am2, a2);
            BBE_MULANX16J(A10, am3, a3);
            am0 = y1; // r(0,2)
            BBE_MULSNX16J(A10, am0, r0);
            BBE_SVNX16_XP(y1, prw, _4L); // r(0,2)
            incw += _4L;
            BBE_SVNX16_XP(y0, prw, incw); // r(0,3)
            // write zeros in r(2,1),r(3,1)
            {
                xb_vecNx16 zero = 0;
                BBE_SVNX16_X(zero, prw, _16L);
                BBE_SVNX16_X(zero, prw, (2 * 16)*L);
            }
            BBE_SVNX16_XP(yx, prw, _4L);  // diagonal element r(1,1)
            SAV0(yx, A10, sh6, sh16, magu);
            BBE_LVNX16_XP(am0, pa, _16L);
            BBE_LVNX16_XP(am1, pa, _16L);
            BBE_LVNX16_XP(am2, pa, _16L);
            BBE_LVNX16_XP(am3, pa, _44L);
            A10 = BBE_MULNX16J(am0, a0);
            BBE_MULANX16J(A10, am1, a1);
            BBE_MULANX16J(A10, am2, a2);
            BBE_MULANX16J(A10, am3, a3);
            am0 = y0; // r(0,3)
            y1 = y0;
            BBE_MULSNX16J(A10, am0, r0);
            SAV0(y0, A10, sh6, sh16, magu);
            // m=2
            BBE_LVNX16_XP(a0, pa1, _16L);
            BBE_LVNX16_XP(a1, pa1, _16L);
            BBE_LVNX16_XP(a2, pa1, _16L);
            BBE_LVNX16_XP(a3, pa1, _44L);
            pa = pa1;
            A10 = BBE_MULNX16J(a0, a0);
            BBE_MULANX16J(A10, a1, a1);
            BBE_MULANX16J(A10, a2, a2);
            BBE_MULANX16J(A10, a3, a3);
            idx = XT_L32I_N(offset, 0);
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
            r0 = BBE_LVNX16_X(prr, idx); // r(0,2)
            r1 = yx; // r(1,2)
            BBE_SVNX16_XP(yx, prw, _4L);  // r(1,2)
            BBE_MULSNX16J(A10, r0, r0);
            BBE_MULSNX16J(A10, r1, r1);
            BBE_LVNX16_IP(t0, ps, (2 * BBE_SIMD_WIDTH));    // reload sigma and compute diagonal element
            C10 = BBE_MOVSWV(t0, t0);
            C10 = BBE_ADDNX40(C10, A10);
#if VARIANT==2
            INV_SQRT2(magu,sh6,sh17,pd,2*BBE_SIMD_WIDTH,_17,C10);
#elif VARIANT==4
            INV_SQRT4(magu, sh6, sh17, pd, 2 * BBE_SIMD_WIDTH, C10);
#endif
            a10_lo = BBE_PACKVNX40(A10, sh17);
            A10 = BBE_MULUUNX16(a10_lo, magu);
            yx = BBE_PACKQNX40(A10);
            // compute element r(2,3)
            BBE_LVNX16_XP(am0, pa, _16L);
            BBE_LVNX16_XP(am1, pa, _16L);
            BBE_LVNX16_XP(am2, pa, _16L);
            BBE_LVNX16_XP(am3, pa, _44L);
            A10 = BBE_MULNX16J(am0, a0);
            BBE_MULANX16J(A10, am1, a1);
            BBE_MULANX16J(A10, am2, a2);
            BBE_MULANX16J(A10, am3, a3);
            am0 = y1; // r(0,3)
            am1 = y0; // r(1,3)
            y2 = y0;
            incw += _4L;
            BBE_SVNX16_XP(y0, prw, incw);   // r(1,3)
            BBE_MULSNX16J(A10, am0, r0);
            BBE_MULSNX16J(A10, am1, r1);
            SAV0(y0, A10, sh6, sh16, magu);
            //        m==3
            BBE_LVNX16_XP(a0, pa1, _16L);
            BBE_LVNX16_XP(a1, pa1, _16L);
            BBE_LVNX16_XP(a2, pa1, _16L);
            BBE_LVNX16_XP(a3, pa1, _44L);
            pa = pa1;
            A10 = BBE_MULNX16J(a0, a0);
            BBE_MULANX16J(A10, a1, a1);
            BBE_MULANX16J(A10, a2, a2);
            BBE_MULANX16J(A10, a3, a3);
            r0 = y1; // r(0,3)
            r1 = y2; // r(1,3)
            r2 = y0; // r(2,3)
            BBE_MULSNX16J(A10, r0, r0);
            BBE_MULSNX16J(A10, r1, r1);
            BBE_MULSNX16J(A10, r2, r2);
            C10 = BBE_MOVSWV(t0, t0);
            C10 = BBE_ADDNX40(C10, A10);
#if VARIANT==2
            INV_SQRT2(magu,sh6,sh17,pd,2*BBE_SIMD_WIDTH,_17,C10);
#elif VARIANT==4
            INV_SQRT4(magu, sh6, sh17, pd, 2 * BBE_SIMD_WIDTH, C10);
#endif
            a10_lo = BBE_PACKVNX40(A10, sh17);
            A10 = BBE_MULUUNX16(a10_lo, magu);
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
            // write zeros in r(3,2)
            {
                xb_vecNx16 zero = 0;
                BBE_SVNX16_X(zero, prw, _16L);
            }
            BBE_SVNX16_XP(yx, prw, _4L);  // diagonal element r(2,2)
            incw += _4L;
            BBE_SVNX16_XP(y0, prw, incw); // r(2,3)
            y0 = BBE_PACKQNX40(A10);
            BBE_SVNX16_XP(y0, prw, _4L);  // diagonal element r(3,3)
            // go to next matrices
            pa1 = (xb_vecNx16*)(((uintptr_t)pa1) - _16L);
            BBE_LVNX16_IP(t0, pa1, (2 * BBE_SIMD_WIDTH));
            BBE_LVNX16_IP(t0, prr, (2 * BBE_SIMD_WIDTH));
            prw = (xb_vecNx16*)prr;
        } while (--l);
    }
    else
    {
        //====================================================================================
        // algorithm for bigger L
        //====================================================================================
        incw = -15 * _4L + 2 * BBE_SIMD_WIDTH;
        //---------------------------------------
        // compute diagonal element r(0,0)
        //---------------------------------------
        prw = (xb_vecNx16 *)R;
        pa = (const xb_vecNx16 *)A;
        ps = (xb_vecNx16   *)sigma2;
        pd = (xb_vecNx16   *)D;
#ifdef COMPILER_XTENSA
#pragma loop_count min=1
#endif
        for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
        {
            BBE_LVNX16_IP(t0, ps, 2 * BBE_SIMD_WIDTH);    // load sigma
            BBE_LVNX16_XP(a0, pa, _16L);
            BBE_LVNX16_XP(a1, pa, _16L);
            BBE_LVNX16_XP(a2, pa, _16L);
            BBE_LVNX16_XP(a3, pa, -3 * _16L + 2 * BBE_SIMD_WIDTH);
#if VARIANT==2
            DIAG2(y0,sh6,sh17,pd,_17,t0,a,r,0);
#elif VARIANT==4
            DIAG4(y0, sh6, sh17, pd, t0, a, r, 0);
#endif
            // write zeros in r(1,0),r(2,0),r(3,0)
            BBE_SVNX16_X(zero, prw, 1 * _16L);
            BBE_SVNX16_X(zero, prw, 2 * _16L);
            BBE_SVNX16_X(zero, prw, 3 * _16L);
            BBE_SVNX16_IP(y0, prw, 2 * BBE_SIMD_WIDTH);  // diagonal element r(0,0)
        }

        //---------------------------------------
        // compute r(0,1),r(0,2),r(0,3)
        //---------------------------------------
        prw = (xb_vecNx16 *)(R + (0 * 4 + 1) * 2 * L);
        pa = (const xb_vecNx16 *)(A);    // A(0,0)
        pd = (xb_vecNx16   *)D;
#ifdef COMPILER_XTENSA
#pragma loop_count min=1
#endif
        for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
        {
            BBE_LVNX16_IP(magu, pd, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(exp, pd, 7 * 2 * BBE_SIMD_WIDTH);
            sh6 = BBE_MOVVSV(exp, 0);
            BBE_LVNX16_XP(a0, pa, _16L);
            BBE_LVNX16_XP(a1, pa, _16L);
            BBE_LVNX16_XP(a2, pa, _16L);
            BBE_LVNX16_XP(a3, pa, _44L);
            BBE_LVNX16_XP(am0, pa, _16L);
            BBE_LVNX16_XP(am1, pa, _16L);
            BBE_LVNX16_XP(am2, pa, _16L);
            BBE_LVNX16_XP(am3, pa, _44L);
            A10 = BBE_MULNX16J(am0, a0);
            BBE_MULANX16J(A10, am1, a1);
            BBE_MULANX16J(A10, am2, a2);
            BBE_MULANX16J(A10, am3, a3);
            SAV0(yx, A10, sh6, sh16, magu);
            BBE_LVNX16_XP(am0, pa, _16L);
            BBE_LVNX16_XP(am1, pa, _16L);
            BBE_LVNX16_XP(am2, pa, _16L);
            BBE_LVNX16_XP(am3, pa, _44L);
            A10 = BBE_MULNX16J(am0, a0);
            BBE_MULANX16J(A10, am1, a1);
            BBE_MULANX16J(A10, am2, a2);
            BBE_MULANX16J(A10, am3, a3);
            SAV0(y1, A10, sh6, sh16, magu);
            BBE_LVNX16_XP(am0, pa, _16L);
            BBE_LVNX16_XP(am1, pa, _16L);
            BBE_LVNX16_XP(am2, pa, _16L);
            BBE_LVNX16_XP(am3, pa, incw);
            A10 = BBE_MULNX16J(am0, a0);
            BBE_MULANX16J(A10, am1, a1);
            BBE_MULANX16J(A10, am2, a2);
            BBE_MULANX16J(A10, am3, a3);
            SAV0(y0, A10, sh6, sh16, magu);
            BBE_SVNX16_X(y1, prw, 1 * _4L); // r(0,2)
            BBE_SVNX16_X(y0, prw, 2 * _4L); // r(0,3)
            BBE_SVNX16_IP(yx, prw, 2 * BBE_SIMD_WIDTH);   //r(0,1)
        }
        incw += _4L;

        //---------------------------------------
        // compute r(1,1), r(2,1),r(3,1)
        //---------------------------------------
        prr = (const xb_vecNx16 *)(R + (0 * 4 + 1) * 2 * L); // r(0,1)
        prw = (xb_vecNx16 *)(R + (1 * 4 + 1) * 2 * L); // r(1,1)
        pa = (const xb_vecNx16 *)(A + (0 * 4 + 1) * 2 * L); // A(0,1)
        ps = (xb_vecNx16   *)sigma2;
        pd = (xb_vecNx16   *)(D + 1 * 2 * BBE_SIMD_WIDTH);
#ifdef COMPILER_XTENSA
#pragma loop_count min=1
#endif
        for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
        {
            BBE_LVNX16_IP(t0, ps, 2 * BBE_SIMD_WIDTH);    // load sigma 
            BBE_LVNX16_XP(a0, pa, _16L);
            BBE_LVNX16_XP(a1, pa, _16L);
            BBE_LVNX16_XP(a2, pa, _16L);
            BBE_LVNX16_XP(a3, pa, -3 * _16L + 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(r0, prr, 2 * BBE_SIMD_WIDTH); // r(0,1)
#if VARIANT==2
            DIAG2(y0,sh6,sh17,pd,_17,t0,a,r,1);
#elif VARIANT==4
            DIAG4(y0, sh6, sh17, pd, t0, a, r, 1);
#endif
            BBE_SVNX16_X(zero, prw, 1 * _16L);
            BBE_SVNX16_X(zero, prw, 2 * _16L);
            BBE_SVNX16_IP(y0, prw, 2 * BBE_SIMD_WIDTH);  // diagonal element r(1,1)
        }
        //---------------------------------------
        // compute elements r(1,2),r(1,3)
        //---------------------------------------
        prr = (const xb_vecNx16 *)(R + (0 * 4 + 1) * 2 * L); // r(0,1)
        prw = (xb_vecNx16 *)(R + (1 * 4 + 2) * 2 * L);       // r(1,2)
        pa = (const xb_vecNx16 *)(A + (0 * 4 + 1) * 2 * L); // A(0,1)
        pd = (xb_vecNx16   *)(D + 1 * 2 * BBE_SIMD_WIDTH);
#ifdef COMPILER_XTENSA
#pragma loop_count min=1
#endif
        for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
        {
            BBE_LVNX16_IP(magu, pd, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(exp, pd, 7 * 2 * BBE_SIMD_WIDTH);
            sh6 = BBE_MOVVSV(exp, 0);
            r2 = BBE_LVNX16_X(prr, 2 * _4L);            // r(0,3)
            r1 = BBE_LVNX16_X(prr, 1 * _4L);            // r(0,2)
            BBE_LVNX16_IP(r0, prr, 2 * BBE_SIMD_WIDTH); // r(0,1)
            BBE_LVNX16_XP(a0, pa, _16L);
            BBE_LVNX16_XP(a1, pa, _16L);
            BBE_LVNX16_XP(a2, pa, _16L);
            BBE_LVNX16_XP(a3, pa, _44L);
            BBE_LVNX16_XP(am0, pa, _16L);
            BBE_LVNX16_XP(am1, pa, _16L);
            BBE_LVNX16_XP(am2, pa, _16L);
            BBE_LVNX16_XP(am3, pa, _44L);
            A10 = BBE_MULNX16J(am0, a0);
            BBE_MULANX16J(A10, am1, a1);
            BBE_MULANX16J(A10, am2, a2);
            BBE_MULANX16J(A10, am3, a3);
            BBE_MULSNX16J(A10, r1, r0);
            SAV0(y0, A10, sh6, sh16, magu);
            BBE_LVNX16_XP(am0, pa, _16L);
            BBE_LVNX16_XP(am1, pa, _16L);
            BBE_LVNX16_XP(am2, pa, _16L);
            BBE_LVNX16_XP(am3, pa, incw);
            A10 = BBE_MULNX16J(am0, a0);
            BBE_MULANX16J(A10, am1, a1);
            BBE_MULANX16J(A10, am2, a2);
            BBE_MULANX16J(A10, am3, a3);
            BBE_MULSNX16J(A10, r2, r0);
            SAV0(y1, A10, sh6, sh16, magu);
            BBE_SVNX16_X(y1, prw, 1 * _4L); // r(1,3)
            BBE_SVNX16_IP(y0, prw, 2 * BBE_SIMD_WIDTH);   //r(1,2)
        }
        incw += _4L;
        //---------------------------------------
        // compute r(2,2), r(3,2)
        //---------------------------------------
        prr = (const xb_vecNx16 *)(R + (0 * 4 + 2) * 2 * L); // r(0,2)
        prw = (xb_vecNx16 *)(R + (2 * 4 + 2) * 2 * L); // r(2,2)
        pa = (const xb_vecNx16 *)(A + (0 * 4 + 2) * 2 * L); // A(0,2)
        ps = (xb_vecNx16   *)sigma2;
        pd = (xb_vecNx16   *)(D + 2 * 2 * BBE_SIMD_WIDTH);
#ifdef COMPILER_XTENSA
#pragma loop_count min=1
#endif
        for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
        {
            BBE_LVNX16_IP(t0, ps, 2 * BBE_SIMD_WIDTH);    // load sigma 
            BBE_LVNX16_XP(a0, pa, _16L);
            BBE_LVNX16_XP(a1, pa, _16L);
            BBE_LVNX16_XP(a2, pa, _16L);
            BBE_LVNX16_XP(a3, pa, -3 * _16L + 2 * BBE_SIMD_WIDTH);
            r1 = BBE_LVNX16_X(prr, _16L);           // r(1,2)
            BBE_LVNX16_IP(r0, prr, 2 * BBE_SIMD_WIDTH); // r(0,2)
#if VARIANT==2
            DIAG1(y0,sh6,sh17,pd,_17,t0,a,r,2);
#elif VARIANT==4
            //            DIAG3(y0,sh6,sh17,pd,t0,a,r,2);
            {
                xb_vecNx40 A, C;
                xb_vecNx16 a10_lo, magu, cn_vec;
                A = BBE_MULNX16J(a0, a0);
                BBE_MULANX16J(A, a1, a1);
                BBE_MULANX16J(A, a2, a2);
                BBE_MULANX16J(A, a3, a3);
                BBE_MULSNX16J(A, r0, r0);
                BBE_MULSNX16J(A, r1, r1);
                C = BBE_MOVSWV(t0, t0);
                C = BBE_ADDNX40(C, C);
                A = BBE_ADDNX40(A, A);
                C = BBE_ADDNX40(C, A);
                sh6 = BBE_NSAENX40(C);
                C = BBE_SLLNX40(C, sh6);
                BBE_RSQRTLUNX40_0(C, magu, cn_vec, C);
                BBE_MULUUSNX16(C, cn_vec, magu);
                C = BBE_SRAINX40(C, 23);
                magu = BBE_PACKLNX40(C);
                magu = BBE_SELNX16I(magu, magu, BBE_SELI_INTERLEAVE_1_EVEN);
                BBE_SVNX16_IP(magu, pd, (2 * BBE_SIMD_WIDTH));
                sh6 = BBE_SHFLVSNI(sh6, BBE_VSA_SHFLI_DUPLICATE_1_EVEN);
                sh17 = BBE_SUBSR1SAVSN(20, sh6);
                sh6 = BBE_SUBSAVSN(12, sh17);
                exp = BBE_MOVVVS(sh6);
                BBE_SVNX16_IP(exp, pd, 7 * 2 * BBE_SIMD_WIDTH);
                a10_lo = BBE_PACKVNX40(A, sh17);
                A = BBE_MULUUNX16(a10_lo, magu);
                y0 = BBE_PACKQNX40(A);
            }
#endif

            BBE_SVNX16_X(zero, prw, 1 * _16L);
            BBE_SVNX16_IP(y0, prw, 2 * BBE_SIMD_WIDTH);  // diagonal element r(2,2)
        }
        //---------------------------------------
        // compute element r(2,3)
        //---------------------------------------
        prr = (const xb_vecNx16 *)(R + (0 * 4 + 2) * 2 * L); // r(0,2)
        prw = (xb_vecNx16 *)(R + (2 * 4 + 3) * 2 * L); // r(2,3)
        pa = (const xb_vecNx16 *)(A + (0 * 4 + 2) * 2 * L); // A(0,2)
        pd = (xb_vecNx16   *)(D + 2 * 2 * BBE_SIMD_WIDTH);
#ifdef COMPILER_XTENSA
#pragma loop_count min=1
#endif
        for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
        {
            BBE_LVNX16_IP(magu, pd, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(exp, pd, 7 * 2 * BBE_SIMD_WIDTH);
            sh6 = BBE_MOVVSV(exp, 0);
            r3 = BBE_LVNX16_X(prr, _16L + _4L);         // r(1,3)
            r2 = BBE_LVNX16_X(prr, _16L);             // r(1,2)
            r1 = BBE_LVNX16_X(prr, 1 * _4L);            // r(0,3)
            BBE_LVNX16_IP(r0, prr, 2 * BBE_SIMD_WIDTH); // r(0,2)
            BBE_LVNX16_XP(a0, pa, _16L);
            BBE_LVNX16_XP(a1, pa, _16L);
            BBE_LVNX16_XP(a2, pa, _16L);
            BBE_LVNX16_XP(a3, pa, _44L);
            BBE_LVNX16_XP(am0, pa, _16L);
            BBE_LVNX16_XP(am1, pa, _16L);
            BBE_LVNX16_XP(am2, pa, _16L);
            BBE_LVNX16_XP(am3, pa, incw);
            A10 = BBE_MULNX16J(am0, a0);
            BBE_MULANX16J(A10, am1, a1);
            BBE_MULANX16J(A10, am2, a2);
            BBE_MULANX16J(A10, am3, a3);
            BBE_MULSNX16J(A10, r1, r0);
            BBE_MULSNX16J(A10, r3, r2);
            SAV0(y0, A10, sh6, sh16, magu);
            BBE_SVNX16_IP(y0, prw, 2 * BBE_SIMD_WIDTH);   //r(2,3)
        }
        //---------------------------------------
        // compute r(2,2), r(3,2)
        //---------------------------------------
        prr = (const xb_vecNx16 *)(R + (0 * 4 + 3) * 2 * L); // r(0,3)
        prw = (xb_vecNx16 *)(R + (3 * 4 + 3) * 2 * L); // r(3,3)
        pa = (const xb_vecNx16 *)(A + (0 * 4 + 3) * 2 * L); // A(0,3)
        ps = (xb_vecNx16   *)sigma2;
        pd = (xb_vecNx16   *)(D + 3 * 2 * BBE_SIMD_WIDTH);
#ifdef COMPILER_XTENSA
#pragma loop_count min=1
#endif
        for (l = 0; l < (L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
        {
            BBE_LVNX16_IP(t0, ps, 2 * BBE_SIMD_WIDTH);    // load sigma 
            BBE_LVNX16_XP(a0, pa, _16L);
            BBE_LVNX16_XP(a1, pa, _16L);
            BBE_LVNX16_XP(a2, pa, _16L);
            BBE_LVNX16_XP(a3, pa, -3 * _16L + 2 * BBE_SIMD_WIDTH);
            r2 = BBE_LVNX16_X(prr, 2 * _16L);           // r(2,3)
            r1 = BBE_LVNX16_X(prr, 1 * _16L);           // r(1,3)
            BBE_LVNX16_IP(r0, prr, 2 * BBE_SIMD_WIDTH); // r(0,3)
#if VARIANT==2
            DIAG1(y0,sh6,sh17,pd,_17,t0,a,r,3);
#elif VARIANT==4
            DIAG3(y0, sh6, sh17, pd, t0, a, r, 3);
#endif
            BBE_SVNX16_IP(y0, prw, 2 * BBE_SIMD_WIDTH);  // diagonal element r(3,3)
        }
    }
} /* chol4x4s() */
#endif

#endif
