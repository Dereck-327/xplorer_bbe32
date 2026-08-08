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
    NatureDSP_Baseband library. Direct Matrix Inversion
    Direct inversion of 2x2 complex floating point matrices 
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
#include "matinv2x2sf_inv.h"

#if HAVE_VFPU

#if 0
#include <complex.h>

static complex_float subc(complex_float x,complex_float y)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    z.s.re=crealf(x)-crealf(y);
    z.s.im=cimagf(x)-cimagf(y);
    return z.u;
}

static complex_float mulc(complex_float x,complex_float y)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    z.s.re=crealf(x)*crealf(y) - cimagf(x)*cimagf(y);
    z.s.im=crealf(x)*cimagf(y) + cimagf(x)*crealf(y);
    return z.u;
}

static complex_float negc(complex_float x)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    z.s.re=-crealf(x);
    z.s.im=-cimagf(x);
    return z.u;
}

static complex_float recipc(complex_float x)
{
    union {complex_float u; struct {float32_t re,im; } s;} z;
    float32_t d;
    d=crealf(x)*crealf(x) + cimagf(x)*cimagf(x);
    d=1.0f/d;
    z.s.re= crealf(x)*d;
    z.s.im=-cimagf(x)*d;
    return z.u;
}
#endif

static void cmatinv2x2f_cbs2x2_inplace(complex_float *X,int L,eLayout layout)
#if 0
{
    int l,n,p;
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    (void)layout;
    NASSERT(layout==e2x2_block)
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/4)
    {
        complex_float x[4][BBE_SIMD_WIDTH/4];
        for (n=0; n<4; n++)
        for (p=0; p<BBE_SIMD_WIDTH/4; p++)
        {
            x[n][p]=X[n+p*4];
        }
        for (p=0; p<BBE_SIMD_WIDTH/4; p++)
        {
            X[p+0*BBE_SIMD_WIDTH/4]=x[0][p];
            X[p+1*BBE_SIMD_WIDTH/4]=x[1][p];
            X[p+2*BBE_SIMD_WIDTH/4]=x[2][p];
            X[p+3*BBE_SIMD_WIDTH/4]=x[3][p];
        }
        X+=4*BBE_SIMD_WIDTH/4;
    }
}
#else
{
    const xb_vecNx16* restrict pX = (const xb_vecNx16*)X;
          xb_vecNx16* restrict pY = (xb_vecNx16*)X;
    xb_vecNx16 x0,x1,x2,x3;
    int l;
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/4)==0);
    (void)layout;
    NASSERT(layout==e2x2_block)
    __Pragma("loop_count min=1")
    for (l = 0; l<(L >> (LOG2_BBE_SIMD_WIDTH - 2)); l++)
    {
        BBE_LVNX16_IP(x0 , pX, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1 , pX, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(x2 , pX, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(x3 , pX, 2*BBE_SIMD_WIDTH); 

        BBE_DSELNX16I(x3, x1, x3, x1 ,BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(x2, x0, x2, x0 ,BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(x3, x2, x3, x2 ,BBE_DSELI_INTERLEAVE_4);
        BBE_DSELNX16I(x1, x0, x1, x0 ,BBE_DSELI_INTERLEAVE_4);

        BBE_SVNX16_IP(x0 , pY, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x1 , pY, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x2 , pY, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x3 , pY, 2*BBE_SIMD_WIDTH);
    }
}
#endif

inline_ void  cmatinv2x2nf_inv_small(complex_float *A,int L)
#if 0
{
    NASSERT_ALIGN(A, (2 * BBE_SIMD_WIDTH));
    NASSERT(L<BBE_SIMD_WIDTH/4);

    int l;
    for (l=0; l<L; l++,A+=4)
    {
        complex_float a,b,c,d,det,rdet;
        a=A[0]; b=A[1]; c=A[2]; d=A[3];
        det=subc(mulc(a,d),mulc(b,c));
        rdet=recipc(det);
        A[0]=mulc( d,rdet);      A[1]=negc(mulc(b,rdet)); 
        A[2]=negc(mulc(c,rdet)); A[3]= mulc(a,rdet); 
    }
}
#else
{
    const xb_vecN_4xcf32 *pArd=(const xb_vecN_4xcf32 *)A;
          xb_vecN_4xcf32 *pAwr=(      xb_vecN_4xcf32 *)A;
    xb_vecN_4xcf32 a,b,c,d,det,t;
    int nbytesrd,nbyteswr;
    valign ar,aw;
    NASSERT(L>0 && L<BBE_SIMD_WIDTH/2);
    nbytesrd=nbyteswr=L*4*sizeof(complex_float);
    ar=BBE_LAN_4XCF32_PP(pArd);
    BBE_LAVN_4XCF32_XP(a,ar,pArd,nbytesrd); nbytesrd-=2*BBE_SIMD_WIDTH;
    BBE_LAVN_4XCF32_XP(b,ar,pArd,nbytesrd); nbytesrd-=2*BBE_SIMD_WIDTH;
    BBE_LAVN_4XCF32_XP(c,ar,pArd,nbytesrd); nbytesrd-=2*BBE_SIMD_WIDTH;
    BBE_LAVN_4XCF32_XP(d,ar,pArd,nbytesrd); 

    BBE_DSELN_4XCF32I(c, a, c, a ,BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(d, b, d, b ,BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(d, c, d, c ,BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(b, a, b, a ,BBE_DSELI_INTERLEAVE_4);
    det=BBE_SUBN_4XCF32(BBE_MULN_4XCF32(a,d),BBE_MULN_4XCF32(b,c));
    det=BBE_RECIPN_4XCF32(det);
    a=BBE_MULN_4XCF32 (a,det);
    t=BBE_MULMN_4XCF32(b,det,3,4); BBE_MULMASN_4XCF32(t,b,det,2,11); b=t;
    t=BBE_MULMN_4XCF32(c,det,3,4); BBE_MULMASN_4XCF32(t,c,det,2,11); c=t;
    d=BBE_MULN_4XCF32 (d,det);
    t=a;a=d;d=t;
    BBE_DSELN_4XCF32I(c, a, c, a ,BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(d, b, d, b ,BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(d, c, d, c ,BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_4XCF32I(b, a, b, a ,BBE_DSELI_INTERLEAVE_4);
    aw=BBE_ZALIGN();
    BBE_SAVN_4XCF32_XP(a,aw,pAwr,nbyteswr); nbyteswr-=2*BBE_SIMD_WIDTH;
    BBE_SAVN_4XCF32_XP(b,aw,pAwr,nbyteswr); nbyteswr-=2*BBE_SIMD_WIDTH;
    BBE_SAVN_4XCF32_XP(c,aw,pAwr,nbyteswr); nbyteswr-=2*BBE_SIMD_WIDTH;
    BBE_SAVN_4XCF32_XP(d,aw,pAwr,nbyteswr); 
    BBE_SAN_4XCF32POS_FP(aw,pAwr);
}
#endif

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
void cmatinv2x2nf(void * restrict pScr, complex_float * restrict A, int L)
{
    int L0;
    NASSERT_ALIGN(A, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr, (2 * BBE_SIMD_WIDTH));
    if (L<=0) return;
    L0=(L)&~(BBE_SIMD_WIDTH/4-1);
    if (L0)
    {
        cmatinv2x2f_cbs2x2_inplace(A,L0,e2x2_block);
        cmatinv2x2f_inv(pScr,A,L0,e2x2_block);
        cmatinv2x2f_cbs2x2_inplace(A,L0,e2x2_block);    // complex 2x2 block to stream and stream to block for bbe32 are equivalent!
    }
    if (L0!=L)
    {
        A+=4*L0;
        L-=L0;
        cmatinv2x2nf_inv_small(A,L);
    }
}

/* Return the scratch area size, in bytes. */
size_t cmatinv2x2nf_getScratchSize ( int L )
{
    size_t sz=0;
    int L0;
    if (L<=0) return 0;
    L0=L&~(BBE_SIMD_WIDTH/4-1);
    sz+=cmatinv2x2f_inv_getScratchSize(L0); /* for matinv2x2f_inv */
    return sz;
}


#else
DISCARD_FUN(void, cmatinv2x2nf,(void * pScr, complex_float * A, int L))
size_t cmatinv2x2nf_getScratchSize ( int L )
{
  return (0);
}
#endif
