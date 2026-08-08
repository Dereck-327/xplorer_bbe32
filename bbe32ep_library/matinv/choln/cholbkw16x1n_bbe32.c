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
  NatureDSP_Baseband library. Cholesky backward recursion for block ordered matrices:
    These functions make backward recursion stage of pseudo-inversion. They use
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
    backward recursion: P==1
    Input:
    Rt[L][16*16][2]    L transformed R matrices
    D[L][SD][2]        reciprocal of main diagonal (mantissa, exponent) 
                       in the special format
    y[L][SY][2]        sequence of intermediate decision matrices y
    qA,qX,qY           fixed point representation of matrices A(or R which 
                       is the same), x and y
    L                  number of matrices
    Output:
    x[L][SX][2]        sequence of decision matrix x
*/
static void bkw16x1(
            int16_t* restrict x, 
            const int16_t* restrict Rt,
            const int16_t* restrict D,
            const int16_t* restrict y, 
            int qXYA,
            int L)

{
    #define N 16
    #define SX 32
    #define SD 32

    int k;
    int l;
    
    xb_vecNx16 tl,th,dd,d0,b_res,b;

    xb_vecNx40 B;
    xb_c40 r_summ;
    vsaN q;

    const int16_t* restrict pYrd;
    const xb_vecNx16* restrict pXrd;
    xb_vecNx16* restrict pXwr;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd;
    
    xb_vecNx16 coef;
    xb_vecNx16 xx0, xx1;
    xb_vecNx16 rr0, rr1;

    NASSERT_ALIGN(x,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Rt,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));

    // load shift amount
    coef = 1;
    coef = BBE_SLSNX16(coef,qXYA);

    const vsaN sh16=BBE_MOVVSA32(16);

    for (k=N-1; k>=0; k--)
    {
        pXrd= (const xb_vecNx16*)x;
        //pXwr= (xb_vecNx16*)&x[k+k];
        pXwr= (xb_vecNx16*)XT_ADDX4(k,(uintptr_t)x);
        //pYrd= &y[k+k];
        pYrd= (const int16_t*)XT_ADDX4(k,(uintptr_t)y);
        //pDrd= (const xb_vecNx16*)&D[k+k];
        pDrd= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);
        pRrd= (const xb_vecNx16*)&Rt[2*N*k];
        __Pragma("loop_count min=1");
        for(l=0; l<L; l++)
        {
            // calculate y(m,:)-R(m,:)*X, 1xP

            // load 16 bit complex value (16 bit re and 16 bit im)
            // load B and shl with saturation
            BBE_LPNX16_XP(b,pYrd,SX*2);
            B=BBE_MULNX16(b,coef);
            // load xx
            BBE_LVNX16_IP(xx0, pXrd, 2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(xx1, pXrd, 2*BBE_SIMD_WIDTH);
            // load rr
            BBE_LVNX16_IP(rr0, pRrd, 2*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(rr1, pRrd, 2*2*N*N - 2*BBE_SIMD_WIDTH);
            //B-= xx*rr
            BBE_MULSNX16C(B,xx0,rr0);
            BBE_MULSNX16C(B,xx1,rr1);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            // get low 16 bits of B
            tl = BBE_PACKLNX40(B);
            // get high 16 bits of B
            th = BBE_PACKVNX40(B,sh16);
            // load D
            BBE_LPNX16_XP(dd,pDrd,SD*2);
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
            BBE_SPNX16_XP(b_res,pXwr,SX*2);
        }
    }
    #undef SX
    #undef N
    #undef SD
}

/*----------------------------------------------------------------------------------
   reversing R matrices for easier readings by rows (diagonal elements are omitted):
   original R    transformed R
   0 1 3 6 a     d x x x
     2 4 7 b     8 c x x
       5 8 c     4 7 b x
         9 d     1 3 6 a
           e

   Input:
   R[L][SR]        L input matrices
   Rt[L*N*(N-1)]   stream of L trasposed matrices
----------------------------------------------------------------------------------*/
static void transformR16x16(int16_t* Rt, void *pScr, const int16_t* R,int L)
{
    int l,k,j;
    int wstride,zstride;

    const xb_vecNx16 * restrict pR  = (const xb_vecNx16*)R;
          xb_vecNx16 * restrict pUR = (      xb_vecNx16*)pScr;
          xb_vecNx16 * restrict pRt = (      xb_vecNx16*)Rt;

    xb_vecNx16 r0, r1;
    xb_vecNx16 x0,x1,x2,x3,x4,x5,x6,x7;

    valign align_rd;

    NASSERT_ALIGN(R,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Rt,(2*BBE_SIMD_WIDTH));

    {
        /* Clear part of scratch that should be filled with zeros */
        r1 = BBE_ZERONX16();
        BBE_SVNX16_IP(r1,pUR,2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1,pUR,4*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1,pUR,4*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1,pUR,4*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1,pUR,4*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1,pUR,4*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1,pUR,4*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1,pUR,4*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(r1,pUR,4*BBE_SIMD_WIDTH);
    }

    __Pragma("loop_count min=1");
    for(l=0; l<L; l++)
    {
        pUR=(xb_vecNx16*)((uintptr_t)pScr+2*2*BBE_SIMD_WIDTH);

        // skip R[0]
        pR=(const xb_vecNx16*)XT_ADDX2(2,(uintptr_t)pR);
        // unpack R
        align_rd = BBE_LA_PP(pR);
        for (k = 2; k <= 8; k++)
        {
            BBE_LAVNX16_XP(r0,align_rd,pR,k*4);

            BBE_SVNX16_IP(r0,pUR,4*BBE_SIMD_WIDTH);
        }
        for (k = 9; k <= 16; k++)
        {
            BBE_LANX16_IP(r0,align_rd,pR);
            BBE_LAVNX16_XP(r1,align_rd,pR,(k*4)-2*BBE_SIMD_WIDTH);

            BBE_SVNX16_IP(r0,pUR,2*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(r1,pUR,2*BBE_SIMD_WIDTH);
        }
        // restore pUR to UR
        pUR=(xb_vecNx16*)pScr;
        // transpose R
        for(k=0; k<4; k++)
        {
            wstride = 2*2*BBE_SIMD_WIDTH;
            zstride = 2*BBE_SIMD_WIDTH;
            j=XT_AND(k,1);
            XT_MOVNEZ(wstride,-15*2*2*BBE_SIMD_WIDTH+2*BBE_SIMD_WIDTH,j);
            XT_MOVEQZ(zstride,-7*2*2*BBE_SIMD_WIDTH+2*BBE_SIMD_WIDTH,j);

            BBE_LVNX16_IP( x0, pUR, 2*2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP( x1, pUR, 2*2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP( x2, pUR, 2*2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP( x3, pUR, 2*2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP( x4, pUR, 2*2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP( x5, pUR, 2*2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP( x6, pUR, 2*2*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP( x7, pUR, wstride);

            BBE_DSELNX16I(x1,x0,x1,x0,BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(x3,x2,x3,x2,BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(x5,x4,x5,x4,BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(x7,x6,x7,x6,BBE_DSELI_DEINTERLEAVE_2);

            BBE_DSELNX16I(x2,x0,x2,x0,BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(x3,x1,x3,x1,BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(x6,x4,x6,x4,BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(x7,x5,x7,x5,BBE_DSELI_DEINTERLEAVE_2);

            BBE_DSELNX16I(x4,x0,x4,x0,BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(x5,x1,x5,x1,BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(x6,x2,x6,x2,BBE_DSELI_DEINTERLEAVE_2);
            BBE_DSELNX16I(x7,x3,x7,x3,BBE_DSELI_DEINTERLEAVE_2);

            BBE_SVNX16_IP( x0, pRt, 2*2*BBE_SIMD_WIDTH); 
            BBE_SVNX16_IP( x1, pRt, 2*2*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP( x2, pRt, 2*2*BBE_SIMD_WIDTH); 
            BBE_SVNX16_IP( x3, pRt, 2*2*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP( x4, pRt, 2*2*BBE_SIMD_WIDTH); 
            BBE_SVNX16_IP( x5, pRt, 2*2*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP( x6, pRt, 2*2*BBE_SIMD_WIDTH); 
            BBE_SVNX16_XP( x7, pRt, zstride);
        }
    }
}

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

void cholbkw16x1n ( void    * pScr,
              complex_fract16* restrict _x, 
        const complex_fract16* restrict _R,
        const complex_fract16* restrict _D,
        const complex_fract16* restrict _y, 
                    int qA, int qY, int qX,
                    int L)
{
          int16_t * restrict x=(      int16_t *)_x;
    const int16_t * restrict R=(const int16_t *)_R;
    const int16_t * restrict D=(const int16_t *)_D;
    const int16_t * restrict y=(const int16_t *)_y;
    #define N 16
    int16_t* Rt=(int16_t*)pScr; // transformed R
    int16_t* pUR = &((int16_t*)pScr)[2*N*N*L]; // unpacked R

    xb_vecNx16 ZeroVec;
    xb_vecNx16 *restrict px = (xb_vecNx16 *)x;
    int l;
    int qXYA=qX-qY+qA;

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(x,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);

    // zeroize x
    ZeroVec=  BBE_MOVVINT16(0);
    __Pragma("loop_count min=1");
    for (l = 0; l < L; l++)
    {
        BBE_SVNX16_IP(ZeroVec,px,2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(ZeroVec,px,2*BBE_SIMD_WIDTH);
    }
    transformR16x16(Rt,pUR,R,L);
    bkw16x1(x,Rt,D,y,qXYA,L);
    #undef N
} /* cholbkw16x1n() */

size_t cholbkw16x1n_getScratchSize (int N,int P,int L)
{
    NASSERT(N==16 && P==1 && L>0);
    return (2*N*N*L + 2*N*N)*sizeof(int16_t);
} /* cholbkw16x1n_getScratchSize() */

#else
DISCARD_FUN(void,cholbkw16x1n,(
            void *pScr,
                  complex_fract16* restrict x, 
            const complex_fract16* restrict R,
            const complex_fract16* restrict D,
            const complex_fract16* restrict y, 
            int qA, int qY, int qX,
            int L))
size_t cholbkw16x1n_getScratchSize(int N,int P,int L) { (void)N,(void)P,(void)L; return 0; }

#endif
