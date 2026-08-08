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
  NatureDSP_Baseband library. Cholesky forward recursion for block ordered matrices:
    These functions make forward recursion stage of pseudo-inversion. They use
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

#if HAVE_CHOLN

/*
    compute L of matrix product Z[L][NxP]=A[L][MxN]'*B[L][MxP]
    Input:
    A[L][SA]    L complex matrices MxN
    B[L][SB]    L complex matrices MxP
    Output:
    Z[L][N*P]   L complex matrices NxP
*/
static void computeAB32x32x1(int32_t* Z,const int16_t* A,const int16_t* B, int qYB, int L)
#if 0
{
    int64_t B_re,B_im;
    int n,m,l;
    int SA=2*32*32;
    int SB=2*32;

    for(n=0; n<32; n++)
    for(l=0; l<L; l++)
    {
        B_re=B_im=0;
        for (m=0; m<32; m++) 
        {
            int16_t a_re,a_im,b_re,b_im;
            a_re=A[l*SA+2*n+m*32*2+0];a_im=A[l*SA+2*n+m*32*2+1];
            b_re=B[l*SB+m*2+0];       b_im=B[l*SB+m*2+1];
            B_re+=L_mul_ss(a_re,b_re)+L_mul_ss(a_im,b_im);  // representation qA+qB
            B_im+=L_mul_ss(a_re,b_im)-L_mul_ss(a_im,b_re);
        }
        B_re=(qYB>=0) ? B_re<<qYB : B_re>>(-qYB);// representation qA+qB->qA+qY-16
        B_im=(qYB>=0) ? B_im<<qYB : B_im>>(-qYB);
        Z[2*l*32+2*n+0]=(int32_t)B_re;
        Z[2*l*32+2*n+1]=(int32_t)B_im;
    }
}
#else
{
    int l,j;
    xb_vecNx16 a, b, b1, zh, zl;
    xb_vecNx40 acc;
    vsaN vqYB=BBE_MOVVSA32(qYB);
    const xb_vecNx16* restrict pa0=(const xb_vecNx16*)A;
    const xb_vecNx16* restrict pb0=(const xb_vecNx16*)B;
          xb_vecNx16* restrict pz0=(      xb_vecNx16*)Z;
    j=0;
    /* set up circular pointers to the boundaries of B */
    WUR_CBEGIN((uintptr_t)(B));
    WUR_CEND  ((uintptr_t)(B+2*32*L));
    __Pragma("loop_count min=4");
    for(l=0; l<L*4; l++)
    {
        int a_stride;
        int z_stride;
        j=BBE_ADDMOD16U(j,(L<<16)|1);
        a_stride=4*2*BBE_SIMD_WIDTH;
        z_stride=2*BBE_SIMD_WIDTH + 6*2*BBE_SIMD_WIDTH;
        XT_MOVEQZ(a_stride,5*2*BBE_SIMD_WIDTH-32*4*2*BBE_SIMD_WIDTH*L,j);
        XT_MOVEQZ(z_stride,9*2*BBE_SIMD_WIDTH-8*2*BBE_SIMD_WIDTH*L,j);
        //computing quarter of Z[l]------------------------------------------------
        //--------part0 of MULACC--------------------------------------------------
        // load B[l] row
        BBE_LVNX16_IP(b, pb0, 2*BBE_SIMD_WIDTH);

        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,0);acc=BBE_MULNX16J( b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,1);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,2);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,3);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,4);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,5);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,6);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,7);BBE_MULANX16J(acc,b1,a);
        //--------part1 of MULACC--------------------------------------------------
        // load B[l] row
        BBE_LVNX16_IP(b, pb0, 2*BBE_SIMD_WIDTH);

        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,0);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,1);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,2);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,3);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,4);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,5);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,6);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,7);BBE_MULANX16J(acc,b1,a);
        //--------part2 of MULACC--------------------------------------------------
        // load B[l] row
        BBE_LVNX16_IP(b, pb0, 2*BBE_SIMD_WIDTH);

        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,0);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,1);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,2);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,3);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,4);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,5);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,6);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,7);BBE_MULANX16J(acc,b1,a);
        //--------part3 of MULACC--------------------------------------------------
        // load B[l] row
        BBE_LVNX16_IC(b, pb0);  /* use circular addressing for automatic returning to the beginning of B */

        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,0);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,1);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,2);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,3);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,4);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,5);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 4*2*BBE_SIMD_WIDTH);b1 = BBE_REPNX16C(b,6);BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_XP(a, pa0, a_stride          );b1 = BBE_REPNX16C(b,7);BBE_MULANX16J(acc,b1,a);
        
        //prepare res
        acc=BBE_SLSNX40(acc,vqYB);
        // moves 8 40-bit elements from the lower half of acc. Each 40-bit element of acc is truncated to 32 bits before the move.
        zl=BBE_MOVSVWL(acc);
        // moves 8 40-bit elements from the upper half of acc. Each 40-bit element of acc is truncated to 32 bits before the move.
        zh=BBE_MOVSVWH(acc);
        // store quarter of Z[l]
        BBE_SVNX16_IP(zl, pz0, 2*BBE_SIMD_WIDTH);
        // skip 8*6 32 bit elements
        BBE_SVNX16_XP(zh, pz0, z_stride);
    }
}
#endif

/*
    make forward recursion (P==1)
   Input:
    R[L][SR][2]     sequence of L upper triangular complex matrices R
    D[L][SD][2]     reciprocal of main diagonal (mantissa, exponent) 
                    in the special format
    Z[L][32*1]      L complex matrices 32x1
    L               Number of matrices
   Output:
   y[L][SY][2]      sequence of intermediate decision matrices y
*/
static void fwd32x1(
            int16_t* restrict y,
            const int16_t* restrict R, 
            const int16_t* restrict D,
            const int32_t* restrict Z, 
            int L)
{
    #define N 32
    #define SR 1056
    
    int l;
    int n;
    const int32_t *pZrd;
    const xb_vecNx16 *pYrd;
    const xb_vecNx16* restrict pDrd;
    int16_t* restrict pYwr;
    const xb_vecNx16 *pRrd;
    xb_vecNx16 rr0,rr1,rr2,rr3;
    xb_vecNx16 yy0,yy1,yy2,yy3;
    xb_vecNx16 dd,d0,tl,th,b_res,b;
    xb_vecNx40 B;
    xb_c40 r_summ;
    valign align;
    vsaN q;

    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Z,(2*BBE_SIMD_WIDTH));

    const vsaN sh16= BBE_MOVVSA32(16);
    __Pragma("loop_count min=4");
    for(n=0; n<N; n++){
        pYrd = (const xb_vecNx16*)y;
        //pYwr = &y[n+n];
        pYwr= (int16_t*)XT_ADDX4(n,(uintptr_t)y);
        //pZrd = &Z[n+n];
        pZrd= (const int32_t*)XT_ADDX8(n,(uintptr_t)Z);
        pRrd = (const xb_vecNx16*)(R+(n*(n+1)));
        //pDrd = (const xb_vecNx16*)&D[n+n];
        pDrd= (const xb_vecNx16*)XT_ADDX4(n,(uintptr_t)D); 
        __Pragma("loop_count min=1");
        for(l=0; l<L; l++){
            // calculate A(:,n)'*B-Rn'*Y, 1xP
            // load 32 bit complex value (32 bit re and 32 bit im)
            // load B
            b = BBE_LV4X16_I(pZrd,0);
            B=BBE_MOVSWVL(b);
            // pZrd+=64 (2*N = 64)
            pZrd=(const int32_t*)XT_ADDX4(64,(uintptr_t)pZrd);
            // load rr
            align = BBE_LA_PP(pRrd);
            BBE_LAVNX16_XP(rr0,align,pRrd,n*4);
            BBE_LAVNX16_XP(rr1,align,pRrd,(n*4)-2*BBE_SIMD_WIDTH);
            BBE_LAVNX16_XP(rr2,align,pRrd,(n*4)-4*BBE_SIMD_WIDTH);
            BBE_LAVNX16_XP(rr3,align,pRrd,(n*4)-6*BBE_SIMD_WIDTH);
            // pRrd+=SR-n*2
            pRrd=(const xb_vecNx16*)XT_ADDX2(SR-n*2,(uintptr_t)pRrd);
            // load y
            BBE_LVNX16_IP(yy0,pYrd,2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(yy1,pYrd,2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(yy2,pYrd,2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(yy3,pYrd,2*BBE_SIMD_WIDTH);
            //B-= rr*yy
            BBE_MULSNX16J(B,yy0,rr0);
            BBE_MULSNX16J(B,yy1,rr1);
            BBE_MULSNX16J(B,yy2,rr2);
            BBE_MULSNX16J(B,yy3,rr3);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            // get low 16 bits of B
            tl = BBE_PACKLNX40(B);
            // get high 16 bits of B
            th = BBE_PACKVNX40(B,sh16);
            // load D
            BBE_LPNX16_XP(dd,pDrd,2*2*2*BBE_SIMD_WIDTH);
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
            BBE_SPNX16_XP(b_res,pYwr,2*2*2*BBE_SIMD_WIDTH);
        }
        __Pragma("no_reorder")
    }
    #undef N
    #undef SR
}

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

void cholfwd32x32x1n ( void    * pScr,
                     complex_fract16 * restrict _y,
               const complex_fract16 * restrict _R, 
               const complex_fract16 * restrict _D,
               const complex_fract16 * restrict _A, 
               const complex_fract16 * restrict _B, 
                       int qB,int qY,int L )
{
          int16_t * restrict y=(      int16_t *)_y;
    const int16_t * restrict R=(const int16_t *)_R;
    const int16_t * restrict D=(const int16_t *)_D;
    const int16_t * restrict A=(const int16_t *)_A;
    const int16_t * restrict B=(const int16_t *)_B;
    int32_t* Z=(int32_t* )pScr;
    int qYB=qY-qB;

    NASSERT(L>0);
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);

    // compute A'*B
    computeAB32x32x1(Z,A,B,qYB,L);
    //cholnFwdrec32(y,R,D,Z,32,L,2*16*33,2*32,2*32,2*32);
    fwd32x1(y,R,D,Z,L);
} /* cholfwd32x32x1n() */

size_t cholfwd32x32x1n_getScratchSize (int M,int N, int P,int L)
{
    size_t Z_size;
    NASSERT(M==32 && N==32 && P==1 && L>0);
    Z_size= (2*L*N*sizeof(int32_t));
    return (Z_size);
} /* cholfwd32x32x1n_getScratchSize() */

#else
DISCARD_FUN(void,cholfwd32x32x1n,( void    * pScr,
                                 complex_fract16 * restrict _y,
                           const complex_fract16 * restrict _R, 
                           const complex_fract16 * restrict _D,
                           const complex_fract16 * restrict _A, 
                           const complex_fract16 * restrict _B, 
                                   int qB,int qY,int L))
size_t cholfwd32x32x1n_getScratchSize  (int M,int N, int P,int L) { (void)M,(void)N,(void)P,(void)L; return 0; }
#endif
