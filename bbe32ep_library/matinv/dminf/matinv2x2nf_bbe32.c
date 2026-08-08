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
    Direct inversion of 2x2 floating point matrices 
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

static void matinv2x2f_rbs2x2_inplace(float32_t *X,int L,eLayout layout)
#if 0
{
    int l,n,p;
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    (void)layout;
    NASSERT(layout==e2x2_block)
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        float32_t x[4][BBE_SIMD_WIDTH/2];
        for (n=0; n<4; n++)
        for (p=0; p<BBE_SIMD_WIDTH/2; p++)
        {
            x[n][p]=X[n+p*4];
        }
        for (p=0; p<BBE_SIMD_WIDTH/2; p++)
        {
            X[p+0*BBE_SIMD_WIDTH/2]=x[0][p];
            X[p+1*BBE_SIMD_WIDTH/2]=x[1][p];
            X[p+2*BBE_SIMD_WIDTH/2]=x[2][p];
            X[p+3*BBE_SIMD_WIDTH/2]=x[3][p];
        }
        X+=4*BBE_SIMD_WIDTH/2;
    }
}
#else
{
    const xb_vecNx16* restrict pX = (const xb_vecNx16*)X;
          xb_vecNx16* restrict pY = (xb_vecNx16*)X;
    xb_vecNx16 x0,x1,x2,x3;
    int l;
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    (void)layout;
    NASSERT(layout==e2x2_block)
    __Pragma("loop_count min=1")
    for (l = 0; l<(L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        BBE_LVNX16_IP(x0 , pX, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1 , pX, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(x2 , pX, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(x3 , pX, 2*BBE_SIMD_WIDTH); 

        BBE_DSELNX16I(x1, x0, x1, x0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(x3, x2, x3, x2,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(x2, x0, x2, x0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(x3, x1, x3, x1,BBE_DSELI_DEINTERLEAVE_2);

        BBE_SVNX16_IP(x0 , pY, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x1 , pY, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x2 , pY, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x3 , pY, 2*BBE_SIMD_WIDTH);
    }
}
#endif

static void matinv2x2f_rsb2x2_inplace(float32_t *X,int L,eLayout layout)
#if 0
{
    int l,n,p;
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    (void)layout;
    NASSERT(layout==e2x2_block)
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        float32_t x[4][BBE_SIMD_WIDTH/2];
        for (p=0; p<BBE_SIMD_WIDTH/2; p++)
        {
            x[0][p]=X[p+0*BBE_SIMD_WIDTH/2];
            x[1][p]=X[p+1*BBE_SIMD_WIDTH/2];
            x[2][p]=X[p+2*BBE_SIMD_WIDTH/2];
            x[3][p]=X[p+3*BBE_SIMD_WIDTH/2];
        }
        for (n=0; n<4; n++)
        for (p=0; p<BBE_SIMD_WIDTH/2; p++)
        {
            X[n+p*4]=x[n][p];
        }
        X+=4*BBE_SIMD_WIDTH/2;
    }
}
#else
{
    const xb_vecNx16* restrict pX = (const xb_vecNx16*)X;
          xb_vecNx16* restrict pY = (xb_vecNx16*)X;
    xb_vecNx16 x0,x1,x2,x3;
    int l;
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    (void)layout;
    NASSERT(layout==e2x2_block)
    __Pragma("loop_count min=1")
    for (l = 0; l<(L >> (LOG2_BBE_SIMD_WIDTH - 1)); l++)
    {
        BBE_LVNX16_IP(x0 , pX, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(x1 , pX, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(x2 , pX, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(x3 , pX, 2*BBE_SIMD_WIDTH); 

        BBE_DSELNX16I(x3, x1, x3, x1 ,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(x2, x0, x2, x0 ,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(x3, x2, x3, x2 ,BBE_DSELI_INTERLEAVE_2);
        BBE_DSELNX16I(x1, x0, x1, x0 ,BBE_DSELI_INTERLEAVE_2);

        BBE_SVNX16_IP(x0 , pY, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x1 , pY, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x2 , pY, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(x3 , pY, 2*BBE_SIMD_WIDTH);
    }
}
#endif

/*
    invert small number of matrices
*/
inline_ void matinv2x2nf_inv_small(float32_t *A,int L)
#if 0
{
    int l;
    NASSERT(L>0 && L<BBE_SIMD_WIDTH/2);
    for (l=0; l<L; l++,A+=4)
    {
        float32_t a,b,c,d,det,rdet;
        a=A[0]; b=A[1]; c=A[2]; d=A[3];
        det=a*d-b*c;
        rdet=1.0f/det;
        A[0]= d*rdet; A[1]=-b*rdet; 
        A[2]=-c*rdet; A[3]= a*rdet; 
    }
}
#else
{
    const xb_vecN_2xf32 *pArd=(const xb_vecN_2xf32 *)A;
          xb_vecN_2xf32 *pAwr=(      xb_vecN_2xf32 *)A;
    xb_vecN_2xf32 a,b,c,d,det,t;
    int nbytesrd,nbyteswr;
    valign ar,aw;
    NASSERT(L>0 && L<BBE_SIMD_WIDTH/2);
    nbytesrd=nbyteswr=L*4*sizeof(float32_t);
    ar=BBE_LAN_2XF32_PP(pArd);
    BBE_LAVN_2XF32_XP(a,ar,pArd,nbytesrd); nbytesrd-=2*BBE_SIMD_WIDTH;
    BBE_LAVN_2XF32_XP(b,ar,pArd,nbytesrd); nbytesrd-=2*BBE_SIMD_WIDTH;
    BBE_LAVN_2XF32_XP(c,ar,pArd,nbytesrd); nbytesrd-=2*BBE_SIMD_WIDTH;
    BBE_LAVN_2XF32_XP(d,ar,pArd,nbytesrd); 

    BBE_DSELN_2XF32I(b, a, b, a,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(d, c, d, c,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(c, a, c, a,BBE_DSELI_DEINTERLEAVE_2);
    BBE_DSELN_2XF32I(d, b, d, b,BBE_DSELI_DEINTERLEAVE_2);

    det=BBE_SUBN_2XF32(BBE_MULN_2XF32(a,d),BBE_MULN_2XF32(b,c));
    det=BBE_RECIPN_2XF32(det);

    a=BBE_MULN_2XF32 (a,det);
    b=BBE_MULMN_2XF32(b,det,3,12);
    c=BBE_MULMN_2XF32(c,det,3,12);
    d=BBE_MULN_2XF32 (d,det);
    t=a;a=d;d=t;

    BBE_DSELN_2XF32I(c, a, c, a ,BBE_DSELI_INTERLEAVE_2);
    BBE_DSELN_2XF32I(d, b, d, b ,BBE_DSELI_INTERLEAVE_2);
    BBE_DSELN_2XF32I(b, a, b, a ,BBE_DSELI_INTERLEAVE_2);
    BBE_DSELN_2XF32I(d, c, d, c ,BBE_DSELI_INTERLEAVE_2);

    aw=BBE_ZALIGN();
    BBE_SAVN_2XF32_XP(a,aw,pAwr,nbyteswr); nbyteswr-=2*BBE_SIMD_WIDTH;
    BBE_SAVN_2XF32_XP(b,aw,pAwr,nbyteswr); nbyteswr-=2*BBE_SIMD_WIDTH;
    BBE_SAVN_2XF32_XP(c,aw,pAwr,nbyteswr); nbyteswr-=2*BBE_SIMD_WIDTH;
    BBE_SAVN_2XF32_XP(d,aw,pAwr,nbyteswr); 
    BBE_SAN_2XF32POS_FP(aw,pAwr);
}
#endif

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
void matinv2x2nf(void * restrict pScr, float32_t* restrict A, int L)
{
    int L0;
    NASSERT_ALIGN(A, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(pScr, (2 * BBE_SIMD_WIDTH));
    NASSERT(L%2==0);
    if (L<=0) return;

    L0=(L)&~(BBE_SIMD_WIDTH/2-1);
    if (L0)
    {
        matinv2x2f_rbs2x2_inplace(A,L0,e2x2_block);
        matinv2x2f_inv(pScr,A,L0,e2x2_block);
        matinv2x2f_rsb2x2_inplace(A,L0,e2x2_block);
    }
    if (L0!=L)
    {
        A+=4*L0;
        L-=L0;
        matinv2x2nf_inv_small(A,L);
    }
}

/* Return the scratch area size, in bytes. */
size_t matinv2x2nf_getScratchSize ( int L )
{
    size_t sz=0;
    int L0;
    if (L<=0) return 0;
    L0=L&~(BBE_SIMD_WIDTH/2-1);
    sz+=matinv2x2f_inv_getScratchSize(L0); /* for matinv2x2f_inv */
    return sz;
}
#else
DISCARD_FUN(void, matinv2x2nf,(void * restrict pScr, float32_t* restrict A, int L))
size_t matinv2x2nf_getScratchSize ( int L )
{
  return (0);
}
#endif
