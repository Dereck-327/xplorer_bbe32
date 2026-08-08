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
  NatureDSP_Baseband library. Cholesky decomposition for block ordered matrices:
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
/* Common utility declarations. */
#include "choln_common.h"

#if !HAVE_CHOLN
DISCARD_FUN(void,chol16x16n,( void    * pScr,
                              complex_fract16 * restrict R, 
                              complex_fract16 * restrict D,
                        const complex_fract16 * restrict A, 
                              int L ) )
#else

/*---------------------------------------------------
   compute n-th column of A'*A for all L matrices

   Input:
   A[L][SA]     L matrices MxN
   n            number of column
   Output:
   Z[L][n+1][2] results
---------------------------------------------------*/
// n=0...7
static void conv0_7(int32_t* Z,const int16_t* A,int n,int M,int N,int L,int SA)
#if 0
{
    int l,k,m;
    int64_t B_re,B_im;
    int16_t amk_re,amk_im,amn_re,amn_im;
    NASSERT(n>=0 && n<=3);
    for (l=0; l<L; l++)
    {
        for (k=0; k<n+1; k++)
        {
            B_re=B_im=0;
            for (m=0; m<M; m++)
            {
                amk_re=A[2*(m*N+k)+0];amk_im=A[2*(m*N+k)+1];
                amn_re=A[2*(m*N+n)+0];amn_im=A[2*(m*N+n)+1];
                B_re+=amk_re*amn_re+amk_im*amn_im;
                B_im+=amk_re*amn_im-amk_im*amn_re;
            }
            Z[0]=(int32_t)B_re;Z[1]=(int32_t)B_im;
            Z+=2;
        }
        A+=SA;
    }
}
#else
{
    int l;
    const xb_vecNx16 *restrict pA = (const xb_vecNx16*)A;
    xb_vecNx16 *restrict pZ = (xb_vecNx16*)Z;

    xb_vecNx16 A0;
    xb_vecNx16 repA;
    xb_vecNx16 Zl, Zh;
    xb_vecNx40 acc0;
    valign align;
    const vselN sel= BBE_MOVVSELNX16(BBE_MOVVA16C(((2*n+1) << 16) + 2*n),0);

    NASSERT(n>=0 && n<=7);
    NASSERT(N==16 && M==16 && SA==2*(16*16));
    (void)M;
    (void)N;
    (void)SA;

    align = BBE_ZALIGN();
    __Pragma("loop_count min=1");
    for (l=0; l<L; l++)
    {
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); acc0= BBE_MULNX16J(repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        BBE_LVNX16_IP(A0, pA, 4*BBE_SIMD_WIDTH); repA = BBE_SHFLNX16(A0,sel); BBE_MULANX16J(acc0,repA,A0);
        // store res
        Zl=BBE_MOVSVWL(acc0);
        Zh=BBE_MOVSVWH(acc0);
        BBE_SAVNX16_XP(Zl,align,pZ,(n+1)*2*4);
        BBE_SAVNX16_XP(Zh,align,pZ,(n+1)*2*4 - 2*BBE_SIMD_WIDTH);
    }
    BBE_SAPOS_FP(align,pZ);
}
#endif

/* --------------------------------------------------
   make forward recursion to update n column elements of 1 matrix
   Input:
   Z[L][SZ]  convolutions in N-th column
   D[L][SD]  reciprocals of main diagonal
   Output:
   y[L][SY]  result of recursion (N elements filled)
--------------------------------------------------*/
// N>0 && N<=8, L==1
static void cholnFwdrec8L1(int16_t* y,const int16_t* R,const int16_t* D,const int32_t* Z,int N,int L,int SR,int SD,int SY,int SZ)
{
    int n;
    const int32_t   * restrict pZrd;
    const xb_vecNx16*          pYrd;
    const xb_vecNx16*          pRrd;
    const xb_vecNx16* restrict pDrd;
          int16_t   *          pYwr;
    xb_vecNx16 dd,d0,rr0,yy0,tl,th,b_res,b;
    xb_vecNx40 B;
    xb_c40 r_summ;
    valign align,ay;
    vsaN q;
    const vsaN sh16= BBE_MOVVSA32(16);
    NASSERT(N>0 && N<=8);
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Z,(2*BBE_SIMD_WIDTH));

    pYwr= (int16_t*)y;
    pDrd= (const xb_vecNx16*)D;
    pZrd= (const int32_t*)Z;
    __Pragma("loop_count min=1");
    for(n=0; n<N; n++)
    {
        // calculate A(:,n)'*B-Rn'*Y, 1xP
        // load 32 bit complex value (32 bit re and 32 bit im)
        // load B
        b = BBE_LV4X16_I(pZrd,0);
        B=BBE_MOVSWVL(b);
        // pZrd+=2
        pZrd=(const int32_t*)XT_ADDX4(2,(uintptr_t)pZrd);
        // load rr (load only n elements. others = 0)
        pRrd = (const xb_vecNx16*)(R+(n*(n+1)));
        align = BBE_LA_PP(pRrd);
        BBE_LAVNX16_XP(rr0,align,pRrd,n*4);
        // load y
        pYrd = (const xb_vecNx16*)y;
        ay=BBE_LA_PP(pYrd);
        BBE_LANX16_IP(yy0,ay,pYrd);
        //B-= rr*yy
        BBE_MULSNX16J(B,yy0,rr0);
        // reduced add
        r_summ = BBE_RADDNX40C(B);
        // type conversion
        B=BBE_MOVNX40_FROMC40(r_summ);
        // get low 16 bits of B
        tl = BBE_PACKLNX40(B);
        // get high 16 bits of B
        th = BBE_PACKVNX40(B,sh16);
        // load D
        BBE_LPNX16_IP(dd,pDrd,2*sizeof(int16_t));
        // replicate dd[0]
        d0 = BBE_REPNX16(dd,0);
        // mul low (unsigned*unsigned)
        B = BBE_MULUUNX16(d0,tl);
        // shift low left
        B = BBE_SRAINX40(B,16);
        // Bres+= mul high (unsigned*signed)
        BBE_MULUSANX16(B,d0,th);
        // make q
        d0= BBE_REPNX16(dd,1);
        q = BBE_MOVVSV(d0,0);
        // result rounding
        B = BBE_RNDADJNX40(B,q);
        // store res
        b_res = BBE_PACKVNX40(B,q);
        BBE_SPNX16_IP(b_res,pYwr,2*sizeof(int16_t));
    }
}

/* --------------------------------------------------
   make forward recursion to update n column elements of 1 matrix
   Input:
   Z[L][SZ]  convolutions in N-th column
   D[L][SD]  reciprocals of main diagonal
   Output:
   y[L][SY]  result of recursion (N elements filled)
--------------------------------------------------*/
// N>8 && N<=16, L==1
void cholnFwdrec16L1(int16_t* y,const int16_t* R,const int16_t* D,const int32_t* Z,int N,int L,int SR,int SD,int SY,int SZ)
{
    int n;
    const int32_t   * restrict pZrd;
    const xb_vecNx16*          pYrd;
    const xb_vecNx16*          pRrd;
    const xb_vecNx16* restrict pDrd;
          int16_t   *          pYwr;
    xb_vecNx16 dd,d0,rr0,rr1,yy0,yy1,tl,th,b_res,b;
    xb_vecNx40 B;
    xb_c40 r_summ;
    valign align,ay;
    vsaN q;
    const vsaN sh16= BBE_MOVVSA32(16);
    NASSERT(N>8 && N<=16);
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Z,(2*BBE_SIMD_WIDTH));

    cholnFwdrec8L1(y,R,D,Z,8,L,SR,SD,SY,SZ);
    pZrd= (const int32_t*)(Z+2*8);
    pYwr= (int16_t*)(y+2*8);
    pDrd= (const xb_vecNx16*)(D+2*8);
    __Pragma("loop_count min=1");
    for(n=8; n<N; n++)
    {
        // calculate A(:,n)'*B-Rn'*Y, 1xP
        // load 32 bit complex value (32 bit re and 32 bit im)
        // load B
        b = BBE_LV4X16_I(pZrd,0);
        B=BBE_MOVSWVL(b);
        // pZrd+=2
        pZrd=(const int32_t*)XT_ADDX4(2,(uintptr_t)pZrd);
        // load rr (load only n elements. others = 0)
        pRrd = (const xb_vecNx16*)(R+(n*(n+1)));
        align = BBE_LA_PP(pRrd);
        BBE_LANX16_IP (rr0,align,pRrd);
        BBE_LAVNX16_XP(rr1,align,pRrd,n*4-2*BBE_SIMD_WIDTH);
        // load y
        pYrd = (const xb_vecNx16*)y;
        ay=BBE_LA_PP(pYrd);
        BBE_LANX16_IP(yy0,ay,pYrd);
        BBE_LANX16_IP(yy1,ay,pYrd);
        //B-= rr*yy
        BBE_MULSNX16J(B,yy0,rr0);
        BBE_MULSNX16J(B,yy1,rr1);
        // reduced add
        r_summ = BBE_RADDNX40C(B);
        // type conversion
        B=BBE_MOVNX40_FROMC40(r_summ);
        // get low 16 bits of B
        tl = BBE_PACKLNX40(B);
        // get high 16 bits of B
        th = BBE_PACKVNX40(B,sh16);
        // load D
        BBE_LPNX16_IP(dd,pDrd,2*sizeof(int16_t));
        // replicate dd[0]
        d0 = BBE_REPNX16(dd,0);
        // mul low (unsigned*unsigned)
        B = BBE_MULUUNX16(d0,tl);
        // shift low left
        B = BBE_SRAINX40(B,16);
        // Bres+= mul high (unsigned*signed)
        BBE_MULUSANX16(B,d0,th);
        // make q
        d0= BBE_REPNX16(dd,1);
        q = BBE_MOVVSV(d0,0);
        // result rounding
        B = BBE_RNDADJNX40(B,q);
        // store res
        b_res = BBE_PACKVNX40(B,q);
        BBE_SPNX16_IP(b_res,pYwr,2*sizeof(int16_t));
    }
}

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

void chol16x16n ( void    * pScr,
                  complex_fract16 * restrict _R, 
                  complex_fract16 * restrict _D,
            const complex_fract16 * restrict _A, 
                  int L )
{
    int16_t *         restrict R=(int16_t *      )_R;
    int16_t *         restrict D=(int16_t *      )_D;
    const int16_t * restrict   A=(const int16_t *)_A;
    fnConv    conv0, conv1;
    fnFwdrec  fwdrec0, fwdrec1;
    fnDiagUpd diagUpd0, diagUpd1;
    int32_t* Z=(int32_t*)pScr; // L columns of A'*A
    int16_t* y; // pointer to the new column in R
    int n;
    int SA=2*(16*16);
    int SR=2*(8*17);
    int SD=2*(16);

    NASSERT(L>0);
    if (L<=0) return;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);

    y=R;
    conv0 = conv0_7;
    conv1 = cholnConv8_15;
    diagUpd0 = cholnDiagUpd8;
    diagUpd1 = cholnDiagUpd16;
    if(1==L)
    {
        fwdrec0 = cholnFwdrec8L1;
        fwdrec1 = cholnFwdrec16L1;
    }
    else
    {
        fwdrec0 = cholnFwdrec8;
        fwdrec1 = cholnFwdrec16;
    }
    
    // 1-st stage
    // compute n-th column of A'*A for all L matrices
    conv0(Z,A,0,16,16,L,SA);
    // update n-th diagonal element
    diagUpd0(y,D,Z,0,L,SR,SD);
    for (n=1; n<8; n++)
    {
        y+=2*n; // go to the next column
        // compute n-th column of A'*A for all L matrices
        conv0(Z,A,n,16,16,L,SA);
        // make forward recursion to update n new column elements
        fwdrec0(y,R,D,Z,n,L,SR,SD,SR,2*(n+1));
        // update n-th diagonal element
        diagUpd0(y,D+2*n,Z,n,L,SR,SD);
    }

    // 2-nd stage
    y+=2*8; // go to the next (8-th) column
    conv1(Z,A,8,16,16,L,SA);
    fwdrec0(y,R,D,Z,8,L,SR,SD,SR,2*(8+1));
    diagUpd0(y,D+2*8,Z,8,L,SR,SD);
    for (n=9; n<16; n++)
    {
        y+=2*n; // go to the next column
        // compute n-th column of A'*A for all L matrices
        conv1(Z,A,n,16,16,L,SA);
        // make forward recursion to update n new column elements
        fwdrec1(y,R,D,Z,n,L,SR,SD,SR,2*(n+1));
        // update n-th diagonal element
        diagUpd1(y,D+2*n,Z,n,L,SR,SD);
    }
} /* chol16x16n() */

#endif
#define CHOLN_SCRATCH(M,N,L) (2*(N)*(L)*sizeof(int32_t))
size_t chol16x16n_getScratchSize  (int M,int N, int L)
{
    NASSERT(N==16 && N==16);
    if  (L<0) L=0;
    return CHOLN_SCRATCH(M,N,L);
}
