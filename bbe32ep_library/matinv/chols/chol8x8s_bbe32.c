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
DISCARD_FUN(void,chol8x8s,(
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

#define VARIANT 4// 2 or 4
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
#define INV_SQRT2(magu,sh6,sh17,pd,pdstep,_17,A)                \
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
    exp=BBE_MOVVVS(c_vec);                                      \
    exp=BBE_SELNX16I(exp,exp,BBE_SELI_INTERLEAVE_1_EVEN);       \
    exp= BBE_SUBSR1RNX16(exp,_17);                              \
    BBE_SVNX16_XP(exp,pd,pdstep);                               \
    sh6=BBE_MOVVSV(exp,0);                                      \
    sh17=BBE_SUBSAVSN(11,sh6);                                  \
}

// another variant
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
    BBE_SVNX16_XP(exp,pd,pdstep);                               \
}

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

void chol8x8s (
                  complex_fract16 * restrict _R, 
                  complex_fract16 * restrict _D,
            const complex_fract16 * restrict _A, 
            const int32_t * restrict sigma2,
            int L)
{
    int16_t *       restrict R=(int16_t *      )_R;
    int16_t *       restrict D=(int16_t *      )_D;
    const int16_t * restrict A=(const int16_t *)_A;
    int _32L,_4L,stepa,stepr;  
    xb_vecNx40 A10,C10,D10;
    xb_vecNx16 t0,r0,rm0,y0,am0,_16;
#if VARIANT==2
    xb_vecNx16 _17;
#endif
    vsaN sh16,sh17,sh6;
    xb_vecNx16 exp,magu,a10_lo;

    const xb_vecNx16    * restrict pa = (const xb_vecNx16   *) A ;
    const xb_vecNx16    * restrict ps = (const xb_vecNx16   *) sigma2 ;
          xb_vecNx16    * restrict pr = (xb_vecNx16   *) R ;
          xb_vecNx16    * restrict pd = (xb_vecNx16   *) D ;
          xb_vecNx16    * restrict prm;
    const xb_vecNx16    * restrict pam;

    int k,l;
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    NASSERT(L>0);

    _32L = L<<5;
    _4L = L<<2;
    _16=BBE_MOVVINT16(16);
#if VARIANT==2
    _17=BBE_MOVVINT16(17);
#endif
    sh16=BBE_MOVVSV(_16,0);

    {
        int k;
        xb_vecNx16 z=BBE_MOVVINT16(0);
        xb_vecNx16 * pR=(xb_vecNx16 *)R;
        for (k=0; k<(8*8*(L>>(LOG2_BBE_SIMD_WIDTH-1))); k++) BBE_SVNX16_IP(z,pR,(2*BBE_SIMD_WIDTH));
    }
    //------------------------------------
    // compute r00,r01...r07
    //------------------------------------
    pa = (xb_vecNx16 *)(A);
    pr = prm= (xb_vecNx16 *)(R);
    stepa=2*BBE_SIMD_WIDTH-7*_32L;
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;

        BBE_LVNX16_IP(t0,ps,2*BBE_SIMD_WIDTH);    // load sigma and compute diagonal element
        BBE_LVNX16_XP(a0,pa, _32L); A10=BBE_MULNX16J( a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa,stepa); BBE_MULANX16J(A10,a0,a0); 

        C10 = BBE_MOVSWV(t0,t0);
        C10 = BBE_ADDNX40(C10,A10);
#if VARIANT==2
        INV_SQRT2(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,_17,C10);
#elif VARIANT==4
        INV_SQRT4(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,C10);
#endif
        A10=BBE_SRANX40(A10,sh17);
        a10_lo=BBE_PACKLNX40(A10);
        A10=BBE_MULUUNX16(a10_lo,magu);
        y0=BBE_PACKQNX40(A10);
        BBE_SVNX16_IP(y0, pr, 2*BBE_SIMD_WIDTH);  // r00
    }
    //------------------------------------
    // compute r01...r07
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    pd = (xb_vecNx16 *)(((uintptr_t)pd)-64*L); // next 32 D
    pa = (const xb_vecNx16 *)(A);
    prm= (xb_vecNx16 *)(R);
    stepa-=2*BBE_SIMD_WIDTH;
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0,a1,a2,a3,a4,a5,a6,a7;
        BBE_LVNX16_IP(magu,pd,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(exp,pd,15*2*BBE_SIMD_WIDTH);
        sh6=BBE_MOVVSV(exp,0);   
        pr = (xb_vecNx16 *)(((uintptr_t)prm)+_4L);
        pam = (xb_vecNx16 *)(((uintptr_t)pa)+_4L);
        BBE_LVNX16_XP(a0,pa, _32L);
        BBE_LVNX16_XP(a1,pa, _32L);
        BBE_LVNX16_XP(a2,pa, _32L);
        BBE_LVNX16_XP(a3,pa, _32L);
        BBE_LVNX16_XP(a4,pa, _32L);
        BBE_LVNX16_XP(a5,pa, _32L);
        BBE_LVNX16_XP(a6,pa, _32L);
        BBE_LVNX16_XP(a7,pa,stepa);
        BBE_LVNX16_XP(am0 ,pam, _32L); A10=BBE_MULNX16J (am0,a0);
        BBE_LVNX16_XP(am0 ,pam, _32L); BBE_MULANX16J(A10,am0,a1);
        BBE_LVNX16_XP(am0 ,pam, _32L); BBE_MULANX16J(A10,am0,a2);
        BBE_LVNX16_XP(am0 ,pam, _32L); BBE_MULANX16J(A10,am0,a3);
        BBE_LVNX16_XP(am0 ,pam, _32L); BBE_MULANX16J(A10,am0,a4);
        BBE_LVNX16_XP(am0 ,pam, _32L); BBE_MULANX16J(A10,am0,a5);
        BBE_LVNX16_XP(am0 ,pam, _32L); BBE_MULANX16J(A10,am0,a6);
        BBE_LVNX16_XP(am0 ,pam,stepa); BBE_MULANX16J(A10,am0,a7);
        SAV0(y0,A10,sh6,sh16,magu);
        BBE_SVNX16_XP(y0, pr,_4L);
        pam = (xb_vecNx16 *)(((uintptr_t)pam)+_4L);
        for (k=2; k<8; k+=2)
        {
            am0=BBE_LVNX16_X(pam ,_4L); C10=BBE_MULNX16J (am0,a0);BBE_LVNX16_XP(am0 ,pam, _32L); A10=BBE_MULNX16J (am0,a0);
            am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a1);BBE_LVNX16_XP(am0 ,pam, _32L); BBE_MULANX16J(A10,am0,a1);
            am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a2);BBE_LVNX16_XP(am0 ,pam, _32L); BBE_MULANX16J(A10,am0,a2);
            am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a3);BBE_LVNX16_XP(am0 ,pam, _32L); BBE_MULANX16J(A10,am0,a3);
            am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a4);BBE_LVNX16_XP(am0 ,pam, _32L); BBE_MULANX16J(A10,am0,a4);
            am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a5);BBE_LVNX16_XP(am0 ,pam, _32L); BBE_MULANX16J(A10,am0,a5);
            am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a6);BBE_LVNX16_XP(am0 ,pam, _32L); BBE_MULANX16J(A10,am0,a6);
            am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a7);BBE_LVNX16_XP(am0 ,pam,stepa); BBE_MULANX16J(A10,am0,a7);
            SAV0(y0,C10,sh6,sh16,magu);
            BBE_SVNX16_X (y0, pr,_4L);
            SAV0(y0,A10,sh6,sh16,magu);
            BBE_SVNX16_XP(y0, pr, 8*L);
            pam = (xb_vecNx16 *)(((uintptr_t)pam)+8*L);
        }
        BBE_LVNX16_IP(a0,pa,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(a0,prm,2*BBE_SIMD_WIDTH);
    }
    stepa+=2*BBE_SIMD_WIDTH;

    //------------------------------------
    // compute r11
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    ps = (const xb_vecNx16   *)(((uintptr_t)ps)-_4L);   // back to sigma2
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L+4*BBE_SIMD_WIDTH); // next 32 D
    A=(const int16_t *)(((uintptr_t)A)+_4L);  // go to the next column
    R=(      int16_t *)(((uintptr_t)R)+_4L);
    pa = (xb_vecNx16 *)((int16_t *)A);
    prm= (xb_vecNx16 *)(R);
    pr = (xb_vecNx16 *)(((uintptr_t)R)+1*_32L);
    stepr=2*BBE_SIMD_WIDTH;
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;
        BBE_LVNX16_IP(t0,ps,2*BBE_SIMD_WIDTH);    // load sigma and compute diagonal element
        BBE_LVNX16_XP(a0,pa, _32L); A10=BBE_MULNX16J( a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, stepr);BBE_MULSNX16J(A10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa,stepa); BBE_MULANX16J(A10,a0,a0); 

        C10 = BBE_MOVSWV(t0,t0);
        C10 = BBE_ADDNX40(C10,A10);
#if VARIANT==2
        INV_SQRT2(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,_17,C10);
#elif VARIANT==4
        INV_SQRT4(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,C10);
#endif
        A10=BBE_SRANX40(A10,sh17);
        a10_lo=BBE_PACKLNX40(A10);
        A10=BBE_MULUUNX16(a10_lo,magu);
        y0=BBE_PACKQNX40(A10);
        BBE_SVNX16_IP(y0, pr, 2*BBE_SIMD_WIDTH);  // r11
    }    
    //------------------------------------
    // compute r12,r13,r14
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L); // next 32 D
    pa = (xb_vecNx16 *)((int16_t *)A);
    pr = (xb_vecNx16 *)((int16_t *)R);
    pam=(const xb_vecNx16*)(((uintptr_t)pa)+_4L);  // go to the next column
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;
        BBE_LVNX16_IP(magu,pd,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(exp,pd,15*2*BBE_SIMD_WIDTH);
        sh6=BBE_MOVVSV(exp,0);   
        prm=(      xb_vecNx16*)(((uintptr_t)pr)+_4L);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); D10=BBE_MULNX16J (am0,a0);am0=BBE_LVNX16_X(pam ,_4L); C10=BBE_MULNX16J (am0,a0);BBE_LVNX16_XP(am0,pam , _32L); A10=BBE_MULNX16J (am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,stepa);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam ,stepa); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(rm0,pr,stepr); r0=BBE_LVNX16_X (prm , 8*L); BBE_MULSNX16J(D10,r0,rm0);r0=BBE_LVNX16_X (prm ,_4L); BBE_MULSNX16J(C10,r0,rm0);BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        SAV0(y0,D10,sh6,sh16,magu);
        BBE_SVNX16_X(y0, prm, 8*L); // r14
        SAV0(y0,C10,sh6,sh16,magu);
        BBE_SVNX16_X(y0, prm,_4L); // r13
        SAV0(y0,A10,sh6,sh16,magu);
        BBE_SVNX16_I(y0, prm, 0);   // r12
    }
    //------------------------------------
    // compute r15,r16,r17
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L); // next 32 D
    pa = (xb_vecNx16 *)((int16_t *)A);
    pr = (xb_vecNx16 *)((int16_t *)R);
    pam=(const xb_vecNx16*)XT_ADDX4(_4L,(uintptr_t)pa);  // go to the next column
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;
        BBE_LVNX16_IP(magu,pd,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(exp,pd,15*2*BBE_SIMD_WIDTH);
        sh6=BBE_MOVVSV(exp,0);   
        prm=(      xb_vecNx16*)XT_ADDX4(_4L,(uintptr_t)pr);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); D10=BBE_MULNX16J (am0,a0);am0=BBE_LVNX16_X(pam ,_4L); C10=BBE_MULNX16J (am0,a0);BBE_LVNX16_XP(am0,pam , _32L); A10=BBE_MULNX16J (am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,stepa);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam ,stepa); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(rm0,pr,stepr); r0=BBE_LVNX16_X (prm , 8*L); BBE_MULSNX16J(D10,r0,rm0);r0=BBE_LVNX16_X (prm ,_4L); BBE_MULSNX16J(C10,r0,rm0);BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        SAV0(y0,D10,sh6,sh16,magu);
        BBE_SVNX16_X(y0, prm, 8*L); // r17
        SAV0(y0,C10,sh6,sh16,magu);
        BBE_SVNX16_X(y0, prm,_4L); // r16
        SAV0(y0,A10,sh6,sh16,magu);
        BBE_SVNX16_I(y0, prm, 0);   // r15
    }
    //------------------------------------
    // compute r22
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    ps = (const xb_vecNx16   *)(((uintptr_t)ps)-_4L);   // back to sigma2
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L+4*BBE_SIMD_WIDTH); // next 32 D
    A=(const int16_t *)(((uintptr_t)A)+_4L);  // go to the next column
    R=(      int16_t *)(pr);
    pa = (xb_vecNx16 *)((int16_t *)A);
    prm= (xb_vecNx16 *)(R);
    pr = (xb_vecNx16 *)XT_ADDX2(_32L,(uintptr_t)R);
    stepr-=_32L;
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;
        BBE_LVNX16_IP(t0,ps,2*BBE_SIMD_WIDTH);    // load sigma and compute diagonal element
        BBE_LVNX16_XP(a0,pa, _32L); A10=BBE_MULNX16J( a0,a0); BBE_LVNX16_XP(r0,prm, _32L);C10=BBE_MULNX16J( r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm,stepr);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa,stepa); BBE_MULANX16J(A10,a0,a0); 

        A10=BBE_SUBNX40(A10,C10);
        C10 = BBE_MOVSWV(t0,t0);
        C10 = BBE_ADDNX40(C10,A10);
#if VARIANT==2
        INV_SQRT2(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,_17,C10);
#elif VARIANT==4
        INV_SQRT4(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,C10);
#endif
        A10=BBE_SRANX40(A10,sh17);
        a10_lo=BBE_PACKLNX40(A10);
        A10=BBE_MULUUNX16(a10_lo,magu);
        y0=BBE_PACKQNX40(A10);
        BBE_SVNX16_IP(y0, pr, 2*BBE_SIMD_WIDTH);  // diagonal element r(m,m)
    }    
    //------------------------------------
    // compute r23,r24,r25
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L); // next 32 D
    pa = (xb_vecNx16 *)((int16_t *)A);
    pr = (xb_vecNx16 *)((int16_t *)R);
    pam=(const xb_vecNx16*)(((uintptr_t)pa)+_4L);  // go to the next column
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;
        BBE_LVNX16_IP(magu,pd,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(exp,pd,15*2*BBE_SIMD_WIDTH);
        sh6=BBE_MOVVSV(exp,0);   
        prm=(      xb_vecNx16*)(((uintptr_t)pr)+_4L);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); D10=BBE_MULNX16J (am0,a0);am0=BBE_LVNX16_X(pam ,_4L); C10=BBE_MULNX16J (am0,a0);BBE_LVNX16_XP(am0,pam , _32L); A10=BBE_MULNX16J (am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,stepa);  am0=BBE_LVNX16_X(pam , 8*L); BBE_MULANX16J(D10,am0,a0);am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam ,stepa); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(rm0,pr, _32L); r0=BBE_LVNX16_X (prm , 8*L); BBE_MULSNX16J(D10,r0,rm0);r0=BBE_LVNX16_X (prm ,_4L); BBE_MULSNX16J(C10,r0,rm0);BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        BBE_LVNX16_XP(rm0,pr,stepr); r0=BBE_LVNX16_X (prm , 8*L); BBE_MULSNX16J(D10,r0,rm0);r0=BBE_LVNX16_X (prm ,_4L); BBE_MULSNX16J(C10,r0,rm0);BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        SAV0(y0,D10,sh6,sh16,magu);
        BBE_SVNX16_X(y0, prm, 8*L); // r25
        SAV0(y0,C10,sh6,sh16,magu);
        BBE_SVNX16_X(y0, prm,_4L); // r24
        SAV0(y0,A10,sh6,sh16,magu);
        BBE_SVNX16_I(y0, prm, 0);   // r23
    }
    //------------------------------------
    // compute r26,r27
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L); // next 32 D
    pa = (xb_vecNx16 *)((int16_t *)A);
    pr = (xb_vecNx16 *)((int16_t *)R);
    pam=(const xb_vecNx16*)XT_ADDX4(_4L,(uintptr_t)pa);  // go to the next column
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;
        BBE_LVNX16_IP(magu,pd,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(exp,pd,15*2*BBE_SIMD_WIDTH);
        sh6=BBE_MOVVSV(exp,0);   
        prm=(      xb_vecNx16*)XT_ADDX4(_4L,(uintptr_t)pr);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); C10=BBE_MULNX16J (am0,a0);BBE_LVNX16_XP(am0,pam , _32L); A10=BBE_MULNX16J (am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,stepa);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam ,stepa); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(rm0,pr, _32L); r0=BBE_LVNX16_X (prm ,_4L); BBE_MULSNX16J(C10,r0,rm0);BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        BBE_LVNX16_XP(rm0,pr,stepr); r0=BBE_LVNX16_X (prm ,_4L); BBE_MULSNX16J(C10,r0,rm0);BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        SAV0(y0,C10,sh6,sh16,magu);
        BBE_SVNX16_X(y0, prm,_4L); // r27
        SAV0(y0,A10,sh6,sh16,magu);
        BBE_SVNX16_I(y0, prm, 0);   // r26
    }
    //------------------------------------
    // compute r33
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    ps = (const xb_vecNx16   *)(((uintptr_t)ps)-_4L);   // back to sigma2
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L+4*BBE_SIMD_WIDTH); // next 32 D
    A=(const int16_t *)(((uintptr_t)A)+_4L);  // go to the next column
    R=(      int16_t *)(pr);
    pa = (xb_vecNx16 *)((int16_t *)A);
    prm= (xb_vecNx16 *)(R);
    pr = (xb_vecNx16 *)(((uintptr_t)R)+3*_32L);
    stepr-=_32L;
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;
        BBE_LVNX16_IP(t0,ps,2*BBE_SIMD_WIDTH);    // load sigma and compute diagonal element
        BBE_LVNX16_XP(a0,pa, _32L); A10=BBE_MULNX16J( a0,a0); BBE_LVNX16_XP(r0,prm, _32L);C10=BBE_MULNX16J( r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm,stepr);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa,stepa); BBE_MULANX16J(A10,a0,a0); 

        A10=BBE_SUBNX40(A10,C10);
        C10 = BBE_MOVSWV(t0,t0);
        C10 = BBE_ADDNX40(C10,A10);
#if VARIANT==2
        INV_SQRT2(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,_17,C10);
#elif VARIANT==4
        INV_SQRT4(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,C10);
#endif
        A10=BBE_SRANX40(A10,sh17);
        a10_lo=BBE_PACKLNX40(A10);
        A10=BBE_MULUUNX16(a10_lo,magu);
        y0=BBE_PACKQNX40(A10);
        BBE_SVNX16_IP(y0, pr, 2*BBE_SIMD_WIDTH);  // diagonal element r(m,m)
    }    
    //------------------------------------
    // compute r34,r35
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L); // next 32 D
    pa = (xb_vecNx16 *)((int16_t *)A);
    pr = (xb_vecNx16 *)((int16_t *)R);
    pam=(const xb_vecNx16*)(((uintptr_t)pa)+_4L);  // go to the next column
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;
        BBE_LVNX16_IP(magu,pd,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(exp,pd,15*2*BBE_SIMD_WIDTH);
        sh6=BBE_MOVVSV(exp,0);   
        prm=(      xb_vecNx16*)(((uintptr_t)pr)+_4L);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); C10=BBE_MULNX16J (am0,a0);BBE_LVNX16_XP(am0,pam , _32L); A10=BBE_MULNX16J (am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,stepa);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam ,stepa); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(rm0,pr, _32L); r0=BBE_LVNX16_X (prm ,_4L); BBE_MULSNX16J(C10,r0,rm0);BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        BBE_LVNX16_XP(rm0,pr, _32L); r0=BBE_LVNX16_X (prm ,_4L); BBE_MULSNX16J(C10,r0,rm0);BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        BBE_LVNX16_XP(rm0,pr,stepr); r0=BBE_LVNX16_X (prm ,_4L); BBE_MULSNX16J(C10,r0,rm0);BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        SAV0(y0,C10,sh6,sh16,magu);
        BBE_SVNX16_X(y0, prm,_4L); // r35
        SAV0(y0,A10,sh6,sh16,magu);
        BBE_SVNX16_I(y0, prm, 0);   // r34
    }
    //------------------------------------
    // compute r36,r37
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L); // next 32 D
    pa = (xb_vecNx16 *)((int16_t *)A);
    pr = (xb_vecNx16 *)((int16_t *)R);
    pam=(const xb_vecNx16*)(((uintptr_t)pa)+12*L);  // go to the next column
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;
        BBE_LVNX16_IP(magu,pd,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(exp,pd,15*2*BBE_SIMD_WIDTH);
        sh6=BBE_MOVVSV(exp,0);   
        prm=(      xb_vecNx16*)(((uintptr_t)pr)+12*L);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); C10=BBE_MULNX16J (am0,a0);BBE_LVNX16_XP(am0,pam , _32L); A10=BBE_MULNX16J (am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa, _32L);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam , _32L); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,stepa);  am0=BBE_LVNX16_X(pam ,_4L); BBE_MULANX16J(C10,am0,a0);BBE_LVNX16_XP(am0,pam ,stepa); BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(rm0,pr, _32L); r0=BBE_LVNX16_X (prm ,_4L); BBE_MULSNX16J(C10,r0,rm0);BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        BBE_LVNX16_XP(rm0,pr, _32L); r0=BBE_LVNX16_X (prm ,_4L); BBE_MULSNX16J(C10,r0,rm0);BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        BBE_LVNX16_XP(rm0,pr,stepr); r0=BBE_LVNX16_X (prm ,_4L); BBE_MULSNX16J(C10,r0,rm0);BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        SAV0(y0,C10,sh6,sh16,magu);
        BBE_SVNX16_X(y0, prm,_4L); //r37
        SAV0(y0,A10,sh6,sh16,magu); 
        BBE_SVNX16_I (y0, prm, 0);  //r36
    }
    //------------------------------------
    // compute r44
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    ps = (const xb_vecNx16   *)(((uintptr_t)ps)-_4L);   // back to sigma2
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L+4*BBE_SIMD_WIDTH); // next 32 D
    A=(const int16_t *)(((uintptr_t)A)+_4L);  // go to the next column
    R=(      int16_t *)(pr);
    pa = (xb_vecNx16 *)((int16_t *)A);
    prm= (xb_vecNx16 *)(R);
    pr = (xb_vecNx16 *)XT_ADDX4(_32L,(uintptr_t)R);
    stepr-=_32L;
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;
        BBE_LVNX16_IP(t0,ps,2*BBE_SIMD_WIDTH);    // load sigma and compute diagonal element
        BBE_LVNX16_XP(a0,pa, _32L); A10=BBE_MULNX16J( a0,a0); BBE_LVNX16_XP(r0,prm, _32L);C10=BBE_MULNX16J( r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm,stepr);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULSNX16J(C10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa,stepa); BBE_MULSNX16J(C10,a0,a0); 

        A10=BBE_SUBNX40(A10,C10);
        C10 = BBE_MOVSWV(t0,t0);
        C10 = BBE_ADDNX40(C10,A10);
#if VARIANT==2
        INV_SQRT2(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,_17,C10);
#elif VARIANT==4
        INV_SQRT4(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,C10);
#endif
        A10=BBE_SRANX40(A10,sh17);
        a10_lo=BBE_PACKLNX40(A10);
        A10=BBE_MULUUNX16(a10_lo,magu);
        y0=BBE_PACKQNX40(A10);
        BBE_SVNX16_IP(y0, pr, 2*BBE_SIMD_WIDTH);  // diagonal element r(m,m)
    }    
    //------------------------------------
    // compute r45,r46,r47
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L); // next 32 D
    pa = (const xb_vecNx16 *)((int16_t *)A);
    pr = (xb_vecNx16 *)((int16_t *)R);
    pam = (const xb_vecNx16   *)(((uintptr_t)pa)+_4L);   
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0,exp;
        BBE_LVNX16_IP(magu,pd,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(exp,pd,15*2*BBE_SIMD_WIDTH);
        sh6=BBE_MOVVSV(exp,0);   
        prm = (      xb_vecNx16   *)(((uintptr_t)pr)+_4L);   
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam, 8*L); D10=BBE_MULNX16J (am0,a0); am0=BBE_LVNX16_X(pam,_4L); C10=BBE_MULNX16J (am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  A10=BBE_MULNX16J (am0,a0);
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam, 8*L); BBE_MULANX16J(D10,am0,a0); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam, 8*L); BBE_MULANX16J(D10,am0,a0); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam, 8*L); BBE_MULANX16J(D10,am0,a0); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam, 8*L); BBE_MULANX16J(D10,am0,a0); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam, 8*L); BBE_MULANX16J(D10,am0,a0); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam, 8*L); BBE_MULANX16J(D10,am0,a0); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,  stepa); am0=BBE_LVNX16_X(pam, 8*L); BBE_MULANX16J(D10,am0,a0); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam,stepa);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(rm0,pr, _32L);  r0 =BBE_LVNX16_X(prm, 8*L); BBE_MULSNX16J(D10,r0,rm0); r0 =BBE_LVNX16_X(prm,_4L); BBE_MULSNX16J(C10,r0,rm0); BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        BBE_LVNX16_XP(rm0,pr, _32L);  r0 =BBE_LVNX16_X(prm, 8*L); BBE_MULSNX16J(D10,r0,rm0); r0 =BBE_LVNX16_X(prm,_4L); BBE_MULSNX16J(C10,r0,rm0); BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        BBE_LVNX16_XP(rm0,pr, _32L);  r0 =BBE_LVNX16_X(prm, 8*L); BBE_MULSNX16J(D10,r0,rm0); r0 =BBE_LVNX16_X(prm,_4L); BBE_MULSNX16J(C10,r0,rm0); BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        BBE_LVNX16_XP(rm0,pr, stepr); r0 =BBE_LVNX16_X(prm, 8*L); BBE_MULSNX16J(D10,r0,rm0); r0 =BBE_LVNX16_X(prm,_4L); BBE_MULSNX16J(C10,r0,rm0); BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        SAV0(y0,D10,sh6,sh16,magu);
        BBE_SVNX16_X (y0, prm, 8*L); //r47
        SAV0(y0,C10,sh6,sh16,magu);
        BBE_SVNX16_X (y0, prm,_4L); //r46
        SAV0(y0,A10,sh6,sh16,magu);
        BBE_SVNX16_I (y0, prm,   0); //r45
    }
    //------------------------------------
    // compute r55
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    ps = (const xb_vecNx16   *)(((uintptr_t)ps)-_4L);   // back to sigma2
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L+4*BBE_SIMD_WIDTH); // next 32 D
    A=(const int16_t *)(((uintptr_t)A)+_4L);  // go to the next column
    R=(      int16_t *)(pr);
    pa = (xb_vecNx16 *)((int16_t *)A);
    prm= (xb_vecNx16 *)(R);
    pr = (xb_vecNx16 *)(((uintptr_t)R)+5*_32L);
    stepr-=_32L;
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;
        BBE_LVNX16_IP(t0,ps,2*BBE_SIMD_WIDTH);    // load sigma and compute diagonal element
        BBE_LVNX16_XP(a0,pa, _32L); A10=BBE_MULNX16J( a0,a0); BBE_LVNX16_XP(r0,prm, _32L);C10=BBE_MULNX16J( r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm,stepr);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa,stepa); BBE_MULANX16J(A10,a0,a0); 

        A10=BBE_SUBNX40(A10,C10);
        C10 = BBE_MOVSWV(t0,t0);
        C10 = BBE_ADDNX40(C10,A10);
#if VARIANT==2
        INV_SQRT2(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,_17,C10);
#elif VARIANT==4
        INV_SQRT4(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,C10);
#endif
        A10=BBE_SRANX40(A10,sh17);
        a10_lo=BBE_PACKLNX40(A10);
        A10=BBE_MULUUNX16(a10_lo,magu);
        y0=BBE_PACKQNX40(A10);
        BBE_SVNX16_IP(y0, pr, 2*BBE_SIMD_WIDTH);  // diagonal element r(m,m)
    }
    //------------------------------------
    // compute r56,r57
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L); // next 32 D
    pr = (xb_vecNx16 *)((int16_t *)R);
    pa = (const xb_vecNx16   *)(((uintptr_t)pam)-_4L);   
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0,exp;
        BBE_LVNX16_IP(magu,pd,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(exp,pd,15*2*BBE_SIMD_WIDTH);
        sh6=BBE_MOVVSV(exp,0);   
        prm = (      xb_vecNx16   *)(((uintptr_t)pr)+_4L);   
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam,_4L); C10=BBE_MULNX16J (am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  A10=BBE_MULNX16J (am0,a0);
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,   _32L); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam, _32L);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(a0,pa,  stepa); am0=BBE_LVNX16_X(pam,_4L); BBE_MULANX16J(C10,am0,a0); BBE_LVNX16_XP(am0,pam,stepa);  BBE_MULANX16J(A10,am0,a0);
        BBE_LVNX16_XP(rm0,pr, _32L);  r0 =BBE_LVNX16_X(prm,_4L);  BBE_MULSNX16J(C10,r0,rm0); BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        BBE_LVNX16_XP(rm0,pr, _32L);  r0 =BBE_LVNX16_X(prm,_4L);  BBE_MULSNX16J(C10,r0,rm0); BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        BBE_LVNX16_XP(rm0,pr, _32L);  r0 =BBE_LVNX16_X(prm,_4L);  BBE_MULSNX16J(C10,r0,rm0); BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        BBE_LVNX16_XP(rm0,pr, _32L);  r0 =BBE_LVNX16_X(prm,_4L);  BBE_MULSNX16J(C10,r0,rm0); BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        BBE_LVNX16_XP(rm0,pr, stepr); r0 =BBE_LVNX16_X(prm,_4L);  BBE_MULSNX16J(C10,r0,rm0); BBE_LVNX16_XP(r0 ,prm , _32L); BBE_MULSNX16J(A10,r0,rm0);
        SAV0(y0,C10,sh6,sh16,magu);
        BBE_SVNX16_X (y0, prm,_4L);    //r57
        SAV0(y0,A10,sh6,sh16,magu);
        BBE_SVNX16_I (y0, prm,   0); //r56
    }
    //------------------------------------
    // compute r66
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    ps = (const xb_vecNx16   *)(((uintptr_t)ps)-_4L);   // back to sigma2
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L+4*BBE_SIMD_WIDTH); // next 32 D
    A=(const int16_t *)(((uintptr_t)A)+_4L);  // go to the next column
    R=(      int16_t *)(pr);
    pa = (xb_vecNx16 *)((int16_t *)A);
    prm= (xb_vecNx16 *)(R);
    pr = (xb_vecNx16 *)(((uintptr_t)R)+6*_32L);
    stepr-=_32L;
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;
        BBE_LVNX16_IP(t0,ps,2*BBE_SIMD_WIDTH);    // load sigma and compute diagonal element
        BBE_LVNX16_XP(a0,pa, _32L); A10=BBE_MULNX16J( a0,a0); BBE_LVNX16_XP(r0,prm, _32L);C10=BBE_MULNX16J( r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm,stepr);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); 
        BBE_LVNX16_XP(a0,pa,stepa); BBE_MULANX16J(A10,a0,a0); 
        A10=BBE_SUBNX40(A10,C10);
        C10 = BBE_MOVSWV(t0,t0);
        C10 = BBE_ADDNX40(C10,A10);
#if VARIANT==2
        INV_SQRT2(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,_17,C10);
#elif VARIANT==4
        INV_SQRT4(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,C10);
#endif
        A10=BBE_SRANX40(A10,sh17);
        a10_lo=BBE_PACKLNX40(A10);
        A10=BBE_MULUUNX16(a10_lo,magu);
        y0=BBE_PACKQNX40(A10);
        BBE_SVNX16_IP(y0, pr, 2*BBE_SIMD_WIDTH);  // diagonal element r(m,m)
    }
    //------------------------------------
    // compute r67
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L); // next 32 D
    pa  = (xb_vecNx16 *)((int16_t *)A);
    prm = (xb_vecNx16 *)((int16_t *)R);
    pr  = (xb_vecNx16 *)(((uintptr_t)R)+6*_32L+_4L);
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0,a1,a2,a3,a4,a5,a6,a7,exp;
        BBE_LVNX16_IP(magu,pd,2*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(exp,pd,15*2*BBE_SIMD_WIDTH);
        sh6=BBE_MOVVSV(exp,0);   
        BBE_LVNX16_XP(a0,pa, _4L);BBE_LVNX16_XP(am0,pa , _32L-_4L); A10=BBE_MULNX16J (am0,a0);
        BBE_LVNX16_XP(a1,pa, _4L);BBE_LVNX16_XP(am0,pa , _32L-_4L); BBE_MULANX16J(A10,am0,a1);
        BBE_LVNX16_XP(a2,pa, _4L);BBE_LVNX16_XP(am0,pa , _32L-_4L); BBE_MULANX16J(A10,am0,a2);
        BBE_LVNX16_XP(a3,pa, _4L);BBE_LVNX16_XP(am0,pa , _32L-_4L); BBE_MULANX16J(A10,am0,a3);
        BBE_LVNX16_XP(a4,pa, _4L);BBE_LVNX16_XP(am0,pa , _32L-_4L); BBE_MULANX16J(A10,am0,a4);
        BBE_LVNX16_XP(a5,pa, _4L);BBE_LVNX16_XP(am0,pa , _32L-_4L); BBE_MULANX16J(A10,am0,a5);
        BBE_LVNX16_XP(a6,pa, _4L);BBE_LVNX16_XP(am0,pa , _32L-_4L); BBE_MULANX16J(A10,am0,a6);
        BBE_LVNX16_XP(a7,pa, _4L);BBE_LVNX16_XP(am0,pa , _32L-_4L-_32L*8+2*BBE_SIMD_WIDTH); BBE_MULANX16J(A10,am0,a7);
        // compute elements r(m,m+1)...r(m,N)
        BBE_LVNX16_XP(rm0,prm,_4L); BBE_LVNX16_XP(r0 ,prm, _32L-_4L);  C10=BBE_MULNX16J( r0,rm0);
        BBE_LVNX16_XP(rm0,prm,_4L); BBE_LVNX16_XP(r0 ,prm, _32L-_4L);  BBE_MULANX16J(C10,r0,rm0);
        BBE_LVNX16_XP(rm0,prm,_4L); BBE_LVNX16_XP(r0 ,prm, _32L-_4L);  BBE_MULANX16J(C10,r0,rm0);
        BBE_LVNX16_XP(rm0,prm,_4L); BBE_LVNX16_XP(r0 ,prm, _32L-_4L);  BBE_MULANX16J(C10,r0,rm0);
        BBE_LVNX16_XP(rm0,prm,_4L); BBE_LVNX16_XP(r0 ,prm, _32L-_4L);  BBE_MULANX16J(C10,r0,rm0);
        BBE_LVNX16_XP(rm0,prm,_4L); BBE_LVNX16_XP(r0 ,prm, _32L-_4L-_32L*6+2*BBE_SIMD_WIDTH);  BBE_MULANX16J(C10,r0,rm0);
        A10=BBE_SUBNX40(A10,C10);
        SAV0(y0,A10,sh6,sh16,magu);
        BBE_SVNX16_IP(y0, pr, 2*BBE_SIMD_WIDTH);
    }
    //------------------------------------
    // compute r77
    //------------------------------------
#ifdef COMPILER_XTENSA
#pragma no_reorder
#endif
    ps = (const xb_vecNx16   *)(((uintptr_t)ps)-_4L);   // back to sigma2
    pd = ( xb_vecNx16   *)(((uintptr_t)pd)-64*L+4*BBE_SIMD_WIDTH); // next 32 D
    A=(const int16_t *)(((uintptr_t)A)+_4L);  // go to the next column
    R=(      int16_t *)(((uintptr_t)R)+_4L);
    pa = (xb_vecNx16 *)((int16_t *)A);
    prm= (xb_vecNx16 *)(R);
    pr = (xb_vecNx16 *)(((uintptr_t)R)+7*_32L);
    stepr-=_32L;
#ifdef COMPILER_XTENSA
  #pragma loop_count min=1
#endif
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx16 a0;
        BBE_LVNX16_IP(t0,ps,2*BBE_SIMD_WIDTH);    // load sigma and compute diagonal element
        BBE_LVNX16_XP(a0,pa, _32L); A10=BBE_MULNX16J( a0,a0); BBE_LVNX16_XP(r0,prm, _32L);C10=BBE_MULNX16J( r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm, _32L);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa, _32L); BBE_MULANX16J(A10,a0,a0); BBE_LVNX16_XP(r0,prm,stepr);BBE_MULANX16J(C10,r0,r0);
        BBE_LVNX16_XP(a0,pa,stepa); BBE_MULANX16J(A10,a0,a0); 
        A10=BBE_SUBNX40(A10,C10);
        C10 = BBE_MOVSWV(t0,t0);
        C10 = BBE_ADDNX40(C10,A10);
#if VARIANT==2
        INV_SQRT2(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,_17,C10);
#elif VARIANT==4
        INV_SQRT4(magu,sh6,sh17,pd,15*2*BBE_SIMD_WIDTH,C10);
#endif
        A10=BBE_SRANX40(A10,sh17);
        a10_lo=BBE_PACKLNX40(A10);
        A10=BBE_MULUUNX16(a10_lo,magu);
        y0=BBE_PACKQNX40(A10);
        BBE_SVNX16_IP(y0, pr, 2*BBE_SIMD_WIDTH);  // diagonal element r(m,m)
    }
} /* chol8x8s() */
#endif
