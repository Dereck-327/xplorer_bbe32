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

#define GET_ALIGNED_BBE_SIMD_WIDTH_SIZE(size) ((size)+2*BBE_SIMD_WIDTH-1)&(~(2*BBE_SIMD_WIDTH-1))

/*
    compute L of matrix product Z[L][NxP]=A[L][MxN]'*B[L][MxP]
    Input:
    A[L][SA]    L complex matrices MxN
    B[L][SB]    L complex matrices MxP
    Output:
    Z[L][N*P]   L complex matrices NxP
*/
static void computeAB8x8x1(int32_t* Z,const int16_t* A,const int16_t* B, int qYB, int L)
#if 0
{
    int64_t B_re,B_im;
    int n,m,l;
    int SA=2*8*8;
    int SB=2*8;

    for(n=0; n<8; n++)
    for(l=0; l<L; l++)
    {
        B_re=B_im=0;
        for (m=0; m<8; m++) 
        {
            int16_t a_re,a_im,b_re,b_im;
            a_re=A[l*SA+2*n+m*8*2+0];a_im=A[l*SA+2*n+m*8*2+1];
            b_re=B[l*SB+m*2+0];      b_im=B[l*SB+m*2+1];
            B_re+=L_mul_ss(a_re,b_re)+L_mul_ss(a_im,b_im);  // representation qA+qB
            B_im+=L_mul_ss(a_re,b_im)-L_mul_ss(a_im,b_re);
        }
        B_re=(qYB>=0) ? B_re<<qYB : B_re>>(-qYB);// representation qA+qB->qA+qY-16
        B_im=(qYB>=0) ? B_im<<qYB : B_im>>(-qYB);
        Z[2*l*8+2*n+0]=(int32_t)B_re;
        Z[2*l*8+2*n+1]=(int32_t)B_im;
    }
}
#else
{
    int l;
    xb_vecNx16 a, b, b1,zh,zl;
    xb_vecNx40 acc;
    vsaN vqYB=BBE_MOVVSA32(qYB);
    const xb_vecNx16* restrict pa0=(const xb_vecNx16*)A;
    const xb_vecNx16* restrict pb0=(const xb_vecNx16*)B;
          xb_vecNx16* restrict pz0=(      xb_vecNx16*)Z;
    __Pragma("loop_count min=1");
    for(l=0; l<L; l++)
    {
        //computing Z[l]---------------------------------------------------
        // load B[l] row
        BBE_LVNX16_IP(b, pb0, 2*BBE_SIMD_WIDTH);
        // load A[l] matrix and compute Z[l]
        BBE_LVNX16_IP(a, pa0, 2*BBE_SIMD_WIDTH); b1 = BBE_REPNX16C(b,0); acc= BBE_MULNX16J(b1,a);
        BBE_LVNX16_IP(a, pa0, 2*BBE_SIMD_WIDTH); b1 = BBE_REPNX16C(b,1); BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 2*BBE_SIMD_WIDTH); b1 = BBE_REPNX16C(b,2); BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 2*BBE_SIMD_WIDTH); b1 = BBE_REPNX16C(b,3); BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 2*BBE_SIMD_WIDTH); b1 = BBE_REPNX16C(b,4); BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 2*BBE_SIMD_WIDTH); b1 = BBE_REPNX16C(b,5); BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 2*BBE_SIMD_WIDTH); b1 = BBE_REPNX16C(b,6); BBE_MULANX16J(acc,b1,a);
        BBE_LVNX16_IP(a, pa0, 2*BBE_SIMD_WIDTH); b1 = BBE_REPNX16C(b,7); BBE_MULANX16J(acc,b1,a);
        //-------------------------------------------------------------------------
        // shift acc
        acc=BBE_SLSNX40(acc,vqYB);
        // moves 8 40-bit elements from the lower half of acc. Each 40-bit element of acc is truncated to 32 bits before the move.
        zl=BBE_MOVSVWL(acc);
        // moves 8 40-bit elements from the upper half of acc. Each 40-bit element of acc is truncated to 32 bits before the move.
        zh=BBE_MOVSVWH(acc);
        // save
        BBE_SVNX16_IP(zl, pz0, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(zh, pz0, 2*BBE_SIMD_WIDTH);
    }
}
#endif

/*
    make forward recursion (P==1)
   Input:
    R[L][SR][2]     sequence of L upper triangular complex matrices R
    D[L][SD][2]     reciprocal of main diagonal (mantissa, exponent) 
                    in the special format
    Z[L][8*1]       L complex matrices 8x1
    L               Number of matrices
   Output:
   y[L][SY][2]      sequence of intermediate decision matrices y
   Restrictions:
   L is multiple of 8
*/
static void fwd8x1_8(
            int16_t* restrict y,
            const int16_t* restrict R, 
            const int16_t* restrict Dtr,
            const int32_t* restrict Z,
            int L,
            void *pScr)
{
    #define N 8
    #define SR 80
    #define SY 16
    #define SD 16

    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Dtr,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Z,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr,(2*BBE_SIMD_WIDTH));
    NASSERT((L%8)==0);

    int l;
    int n;
    const int32_t* restrict pZrd;
    const xb_vecNx16* restrict pYrd;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd= (const xb_vecNx16*)Dtr;
    int16_t* restrict pYwr;
    xb_vecNx16* restrict pYrdwr;

    xb_vecNx16 d0,rr,yy,tl,th,b,vTmp0,vTmp1,vTmp2,vTmp3;

    xb_vecNx40 B;
    xb_c40 r_summ;
    valign align;
    valign Yrdwr_align;
    vsaN q;

    //const vsaN sh16= BBE_MOVVSA32(16);
    Yrdwr_align= BBE_ZALIGN();

    {
        //pYwr = &y[n+n];
        pYwr= y;
        //pZrd = &Z[n+n];
        pZrd= Z;
        pYrdwr= (xb_vecNx16*)pScr;

        __Pragma("loop_count min=2 factor=2");
        for(l=0; l<(L>>2); l++)
        {
            //0----------------------------------------------------------
            vTmp0 = BBE_LV4X16_I(pZrd,0);
            //1----------------------------------------------------------
            vTmp1 = BBE_LV4X16_I(pZrd,1*4*N*2);
            //1----------------------------------------------------------
            vTmp2 = BBE_LV4X16_I(pZrd,2*4*N*2);
            //1----------------------------------------------------------
            vTmp3 = BBE_LV4X16_I(pZrd,3*4*N*2);

            vTmp0= BBE_SELNX16I(vTmp1,vTmp0,BBE_SELI_PACK_4);
            vTmp2= BBE_SELNX16I(vTmp3,vTmp2,BBE_SELI_PACK_4);
            vTmp0= BBE_SELNX16I(vTmp2,vTmp0,BBE_SELI_PACK_8);
            BBE_SVNX16_IP(vTmp0,pYrdwr,2*BBE_SIMD_WIDTH);
            // pZrd+=N*2
            pZrd=(const int32_t*)XT_ADDX8(2*N*2,(uintptr_t)pZrd);
        }
        pYrdwr= (xb_vecNx16*)pScr;

        __Pragma("loop_count min=1");
        for(l=0; l<(L>>3); l++)
        {
            BBE_LVNX16_IP(vTmp0,pYrdwr,2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(vTmp1,pYrdwr,2*BBE_SIMD_WIDTH);
            // get low 16 and high 16 bits of B
            BBE_DSELNX16I(th, tl, vTmp1, vTmp0, BBE_DSELI_DEINTERLEAVE_1);
            // load D
            BBE_LVNX16_IP(d0,pDrd,2*BBE_SIMD_WIDTH);
            vTmp0= BBE_SHFLNX16I(d0,BBE_SHFLI_DUPLICATE_1_ODD);
            d0= BBE_SHFLNX16I(d0,BBE_SHFLI_DUPLICATE_1_EVEN);
            q = BBE_MOVVSV(vTmp0,0);
            // mul low (unsigned*unsigned)
            B = BBE_MULUUNX16(d0,tl);
            // shift low left
            B = BBE_SRAINX40(B,16);
            // Bres+= mul high (unsigned*signed)
            BBE_MULUSANX16(B,d0,th);
            // result rounding
            B = BBE_RNDADJNX40(B,q);
            // store res
            vTmp0 = BBE_PACKVNX40(B,q);
            BBE_SPNX16_IP(vTmp0,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,1);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,2);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,3);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,4);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,5);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,6);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,7);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
        }
    }
    for(n=1; n<N; n++)
    {
        pYrd = (const xb_vecNx16*)y;
        //pYwr = &y[n+n];
        pYwr= (int16_t*)XT_ADDX4(n,(uintptr_t)y);
        //pZrd = &Z[n+n];
        pZrd= (const int32_t*)XT_ADDX8(n,(uintptr_t)Z);
        pRrd = (const xb_vecNx16*)(R+(n*(n+1)));
        pYrdwr= (xb_vecNx16*)pScr;

        __Pragma("loop_count min=1 factor=4");
        for(l=0; l<(L>>1); l++)
        {
            //0----------------------------------------------------------
            b = BBE_LV4X16_I(pZrd,0);
            B=BBE_MOVSWVL(b);
            // load rr (load only n elements. others = 0)
            align = BBE_LA_PP(pRrd);
            BBE_LAVNX16_XP(rr,align,pRrd,n*4);
            // pRrd+=SR-n*2
            pRrd=(const xb_vecNx16*)XT_ADDX2(SR-n*2,(uintptr_t)pRrd);
            // load y, pYrd+= SY*2
            BBE_LVNX16_IP(yy,pYrd,2*SY);
            //B-= rr*yy
            BBE_MULSNX16J(B,yy,rr);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            vTmp0= BBE_MOVVWL(B);
            //1----------------------------------------------------------
            b = BBE_LV4X16_I(pZrd,4*N*2);
            B=BBE_MOVSWVL(b);
            // load rr (load only n elements. others = 0)
            align = BBE_LA_PP(pRrd);
            BBE_LAVNX16_XP(rr,align,pRrd,n*4);
            // pRrd+=SR-n*2
            pRrd=(const xb_vecNx16*)XT_ADDX2(SR-n*2,(uintptr_t)pRrd);
            // load y, pYrd+= SY*2
            BBE_LVNX16_IP(yy,pYrd,2*SY);
            //B-= rr*yy
            BBE_MULSNX16J(B,yy,rr);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            vTmp1= BBE_MOVVWL(B);

            vTmp0= BBE_SELNX16I(vTmp1,vTmp0,BBE_SELI_PACK_4);
            BBE_SAVNX16_XP(vTmp0,Yrdwr_align,pYrdwr,2*8);
            // pZrd+=N*2
            pZrd=(const int32_t*)XT_ADDX4(2*N*2,(uintptr_t)pZrd);
        }
        BBE_SAPOS_FP(Yrdwr_align,pYrdwr);
        pYrdwr= (xb_vecNx16*)pScr;

        __Pragma("loop_count min=1");
        for(l=0; l<(L>>3); l++)
        {
            BBE_LVNX16_IP(vTmp0,pYrdwr,2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(vTmp1,pYrdwr,2*BBE_SIMD_WIDTH);
            // get low 16 and high 16 bits of B
            BBE_DSELNX16I(th, tl, vTmp1, vTmp0, BBE_DSELI_DEINTERLEAVE_1);
            // load D
            BBE_LVNX16_IP(d0,pDrd,2*BBE_SIMD_WIDTH);
            vTmp0= BBE_SHFLNX16I(d0,BBE_SHFLI_DUPLICATE_1_ODD);
            d0= BBE_SHFLNX16I(d0,BBE_SHFLI_DUPLICATE_1_EVEN);
            q = BBE_MOVVSV(vTmp0,0);
            // mul low (unsigned*unsigned)
            B = BBE_MULUUNX16(d0,tl);
            // shift low left
            B = BBE_SRAINX40(B,16);
            // Bres+= mul high (unsigned*signed)
            BBE_MULUSANX16(B,d0,th);
            // result rounding
            B = BBE_RNDADJNX40(B,q);
            // store res
            vTmp0 = BBE_PACKVNX40(B,q);
            BBE_SPNX16_IP(vTmp0,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,1);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,2);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,3);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,4);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,5);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,6);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
            vTmp1= BBE_REPNX16C(vTmp0,7);
            BBE_SPNX16_IP(vTmp1,pYwr,2*SY);
        }
        //__Pragma("no_reorder")
    }
    #undef N
    #undef SR
    #undef SY
    #undef SD
}

/*
    Convert Diagonal D from block to streaming format
    
    Input:
    D_in            D[L][8][2] diagonal elements in block format
    D_out           D[8][L][2] diagonal element in streaming format
    L               number of matrices
    restrictions:
    L is multiple of 8
*/
static
void transformD8( int16_t* restrict D_out,
            const int16_t* restrict D_in,
            const int L)
{
    #define N 8
    #define SD 16

    NASSERT((L%8)==0);

    int l;
    const int stride= (L >> 3)*2*BBE_SIMD_WIDTH;
    const xb_vecNx16* restrict pD_in_rd= (const xb_vecNx16*)D_in;
    xb_vecNx16* restrict pD_out_wr= (xb_vecNx16*)D_out;
    xb_vecNx16 vX0, vX1, vX2, vX3, vX4, vX5, vX6, vX7;

    __Pragma("loop_count min=1");
    for (l=0; l<(L>>3); l++)
    {
        // read 8 D*8*2 elements
        BBE_LVNX16_IP(vX0, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX1, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX2, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX3, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX4, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX5, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX6, pD_in_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vX7, pD_in_rd, 2*BBE_SIMD_WIDTH);

        // complex transpose
        BBE_DSELNX16I(vX1,vX0,vX1,vX0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX3,vX2,vX3,vX2,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX5,vX4,vX5,vX4,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX7,vX6,vX7,vX6,BBE_DSELI_DEINTERLEAVE_2);

        BBE_DSELNX16I(vX2,vX0,vX2,vX0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX3,vX1,vX3,vX1,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX6,vX4,vX6,vX4,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX7,vX5,vX7,vX5,BBE_DSELI_DEINTERLEAVE_2);

        BBE_DSELNX16I(vX4,vX0,vX4,vX0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX5,vX1,vX5,vX1,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX6,vX2,vX6,vX2,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vX7,vX3,vX7,vX3,BBE_DSELI_DEINTERLEAVE_2);
        // now transposed matrix in vX0, vX1, ..., vX7
        
        // store res
        BBE_SVNX16_XP(vX0,pD_out_wr,stride);
        BBE_SVNX16_XP(vX1,pD_out_wr,stride);
        BBE_SVNX16_XP(vX2,pD_out_wr,stride);
        BBE_SVNX16_XP(vX3,pD_out_wr,stride);
        BBE_SVNX16_XP(vX4,pD_out_wr,stride);
        BBE_SVNX16_XP(vX5,pD_out_wr,stride);
        BBE_SVNX16_XP(vX6,pD_out_wr,stride);
        BBE_SVNX16_XP(vX7,pD_out_wr,-7*stride+2*BBE_SIMD_WIDTH);
    }
    #undef N
    #undef SD
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

void cholfwd8x8x1n ( void    * pScr,
                     complex_fract16 * restrict _y,
               const complex_fract16 * restrict _R, 
               const complex_fract16 * restrict _D,
               const complex_fract16 * restrict _A, 
               const complex_fract16 * restrict _B, 
                     int qB,
                     int qY,
                     int L )
{
          int16_t * restrict y=(      int16_t *)_y;
    const int16_t * restrict R=(const int16_t *)_R;
    const int16_t * restrict D=(const int16_t *)_D;
    const int16_t * restrict A=(const int16_t *)_A;
    const int16_t * restrict B=(const int16_t *)_B;
#define N 8
#define SD 16
#define SY 16
#define SR 80
#define SZ 16
    int L8, Lmod;
    int32_t *Z=(int32_t *)pScr;
    int qYB=qY-qB;

    NASSERT(L>0);
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);

    // compute A'*B
    computeAB8x8x1(Z, A, B, qYB, L);
    
    Lmod = L&(BBE_SIMD_WIDTH/2-1);
    L8 = L - Lmod;
    // calculate values that less or equal than L and multiple of SIMD_WIDTH/2 (L8)
    if (0!=L8)
    {
        size_t Z_size    = (2*L*N*sizeof(int32_t));
        size_t D_tr_size = GET_ALIGNED_BBE_SIMD_WIDTH_SIZE(2*L*N*sizeof(int16_t));
        int16_t *D_tr = (int16_t*)((char*)pScr + Z_size);
        pScr= ((char*)pScr + Z_size + D_tr_size);
        transformD8(D_tr,D,L8);
        fwd8x1_8(y,R,D_tr,Z,L8,pScr);
    }
    if (0!=Lmod)
    {
              int16_t* y2 = y + L8*SY;
        const int16_t* R2 = R + L8*SR;
        const int16_t* D2 = D + L8*SD;
        const int32_t* Z2 = Z + L8*SZ;
        cholnFwdrec8(y2,R2,D2,Z2,N,Lmod,SR,SD,SY,SZ);
    }
#undef SZ
#undef SR
#undef SY
#undef SD
#undef N
} /* cholfwd8x8x1n() */

size_t cholfwd8x8x1n_getScratchSize (int M, int N, int P,int L)
{
    size_t Z_size;
    size_t D_tr_size = 0;
    size_t Intrmd_size = 0;
    NASSERT(M==8 && N==8 && P==1 && L>0);
    Z_size= (2*L*N*sizeof(int32_t));
    if (L>7)
    {
        D_tr_size = 2*L*N*sizeof(int16_t);
        Intrmd_size = GET_ALIGNED_BBE_SIMD_WIDTH_SIZE(L*2*sizeof(int32_t));
    }
    return (Z_size + D_tr_size + Intrmd_size);
} /* cholfwd8x8x1n_getScratchSize() */

#else
DISCARD_FUN(void,cholfwd8x8x1n,(
            void    * pScr,
                     complex_fract16 * restrict _y,
               const complex_fract16 * restrict _R, 
               const complex_fract16 * restrict _D,
               const complex_fract16 * restrict _A, 
               const complex_fract16 * restrict _B, 
                     int qB,
                     int qY,
                     int L ))
size_t cholfwd8x8x1n_getScratchSize  (int M,int N, int P,int L) { (void)M,(void)N,(void)P,(void)L; return 0; }
#endif
