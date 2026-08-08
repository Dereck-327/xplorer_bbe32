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
  NatureDSP_Baseband library. Apply the QR decomposition to the matrix of normal equations system
    Apply backward recursion process for QR decomposition for block ordered
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
#include "cqrn_common.h"

#if !HAVE_CQRN
DISCARD_FUN(void,cqr_bkw8x8x1n,(void *pScr,
                          complex_fract16* X,
                    const complex_fract16* R,
                    const complex_fract16* D,
                    int qABX,
                    int L))
#else

/*
    Convert Diagonal D from block to streaming format
    
    Input:
    D_in            D[L][8][2] diagonal elements in block format
    D_out           D[8][L][2] diagonal element in streaming format
    L               number of matrices
    restrictions:
    L is multiple of 8
*/
static void transformD8(int16_t* restrict D_out,
                  const int16_t* restrict D_in,
                  const int L)
{
    #define N 8
    #define SD 16

    int l;
    const int stride= (L >> 3)*2*BBE_SIMD_WIDTH;
    const xb_vecNx16* restrict pD_rd= (const xb_vecNx16*)D_in;
          xb_vecNx16* restrict pD_wr= (xb_vecNx16*)D_out;
    xb_vecNx16 vD0, vD1, vD2, vD3, vD4, vD5, vD6, vD7;

    NASSERT_ALIGN(D_in ,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(D_out,(2*BBE_SIMD_WIDTH));
    NASSERT((L%8)==0);

    __Pragma("loop_count min=1");
    for (l=0; l<(L>>3); l++)
    {
        // read 8 D*8*2 elements
        BBE_LVNX16_IP(vD0, pD_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vD1, pD_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vD2, pD_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vD3, pD_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vD4, pD_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vD5, pD_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vD6, pD_rd, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vD7, pD_rd, 2*BBE_SIMD_WIDTH);

        // complex transpose
        BBE_DSELNX16I(vD1,vD0,vD1,vD0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vD3,vD2,vD3,vD2,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vD5,vD4,vD5,vD4,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vD7,vD6,vD7,vD6,BBE_DSELI_DEINTERLEAVE_2);

        BBE_DSELNX16I(vD2,vD0,vD2,vD0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vD3,vD1,vD3,vD1,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vD6,vD4,vD6,vD4,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vD7,vD5,vD7,vD5,BBE_DSELI_DEINTERLEAVE_2);

        BBE_DSELNX16I(vD4,vD0,vD4,vD0,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vD5,vD1,vD5,vD1,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vD6,vD2,vD6,vD2,BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(vD7,vD3,vD7,vD3,BBE_DSELI_DEINTERLEAVE_2);
        // now transposed matrix in vD0, vD1, ..., vD7
        
        // store res
        BBE_SVNX16_XP(vD0,pD_wr,stride);
        BBE_SVNX16_XP(vD1,pD_wr,stride);
        BBE_SVNX16_XP(vD2,pD_wr,stride);
        BBE_SVNX16_XP(vD3,pD_wr,stride);
        BBE_SVNX16_XP(vD4,pD_wr,stride);
        BBE_SVNX16_XP(vD5,pD_wr,stride);
        BBE_SVNX16_XP(vD6,pD_wr,stride);
        BBE_SVNX16_XP(vD7,pD_wr,-7*stride+2*BBE_SIMD_WIDTH);
    }
    #undef N
    #undef SD
}

/*
    Apply backward recursion process
    Restrictions:
    L is multiple of 8
*/
static void cqr_bkw8x8x1n_8(int16_t* restrict X,
                      const int16_t* restrict R,
                      const int16_t* restrict Dtr,
                            int qABX,
                            int L,
                            void* restrict pScr)
{
    #define N 8	
    #define SR 128
    #define SX 16
    #define SD 16

    int k;
    int l;
    
    xb_vecNx16 d0,tl,th,b,vTmp0,vTmp1;

    xb_vecNx40 B;
    xb_c40 r_summ;
    valign R_align;
    valign X_align;
    valign Yrdwr_align;
    vsaN q;

    const xb_vecNx16* restrict pXrd_first;
    const xb_vecNx16* restrict pXrd;
          xb_vecNx16* restrict pXwr;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd= (const xb_vecNx16*)&Dtr[7*L*2];
          xb_vecNx16* restrict pYrdwr= (xb_vecNx16*)pScr;

    xb_vecNx16 xx;
    xb_vecNx16 rr;

    const vsaN vqABX= BBE_MOVVSA32(qABX);
    
    NASSERT_ALIGN(X  ,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(R  ,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(Dtr,(2*BBE_SIMD_WIDTH));
    NASSERT(L>0);
    NASSERT((L%8)==0);

    Yrdwr_align= BBE_ZALIGN();

    for (k=N-1; k>=0; k--)
    {
        //pXrd_first= (const xb_vecNx16*)&X[k+k];
        pXrd_first= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)X);
        pXrd      = (const xb_vecNx16*)XT_ADDI((uintptr_t)pXrd_first,2*2);
        pXwr      = (      xb_vecNx16*)pXrd_first;
        pRrd      = (const xb_vecNx16*)(R + (2*k*N+(k+k+2)));
        X_align = BBE_LA_PP(pXrd);
        pYrdwr= (xb_vecNx16*)pScr;

        __Pragma("loop_count min=2 factor=2");
        for(l=0; l<L; l+=2)
        {
            // calculate y(m,:)-R(m,:)*X, 1xP
            // load 16 bit complex value (16 bit re and 16 bit im)
            BBE_LPNX16_IP(b,pXrd_first,SX*2);
            B= BBE_UNPKSNX16(b);
            B=BBE_SLSNX40(B,vqABX);
            // load xx (load only N-k-1 elements. others = 0)
            X_align = BBE_LA_PP(pXrd);
            BBE_LAVNX16_XP(xx,X_align,pXrd,(N-1-k)*2*2);
            //pXrd++
            pXrd=(const xb_vecNx16*)XT_ADDX2(SX-((N-1-k)*2),(uintptr_t)pXrd);
            
            // load rr (load only N-k-1 elements. others = 0)
            R_align = BBE_LA_PP(pRrd);
            BBE_LAVNX16_XP(rr,R_align,pRrd,(N-1-k)*2*2);
            //pRrd++
            pRrd=(const xb_vecNx16*)XT_ADDX2(SR-((N-1-k)*2),(uintptr_t)pRrd);

            //B-= rr*xx
            BBE_MULSNX16C(B,xx,rr);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            vTmp0= BBE_MOVVWL(B);
            //-----------------------------------------------------------------
            BBE_LPNX16_IP(b,pXrd_first,SX*2);
            B= BBE_UNPKSNX16(b);
            B=BBE_SLSNX40(B,vqABX);
            // load xx (load only N-k-1 elements. others = 0)
            X_align = BBE_LA_PP(pXrd);
            BBE_LAVNX16_XP(xx,X_align,pXrd,(N-1-k)*2*2);
            //pXrd++
            pXrd=(const xb_vecNx16*)XT_ADDX2(SX-((N-1-k)*2),(uintptr_t)pXrd);
            
            // load rr (load only N-k-1 elements. others = 0)
            R_align = BBE_LA_PP(pRrd);
            BBE_LAVNX16_XP(rr,R_align,pRrd,(N-1-k)*2*2);
            //pRrd++
            pRrd=(const xb_vecNx16*)XT_ADDX2(SR-((N-1-k)*2),(uintptr_t)pRrd);

            //B-= rr*xx
            BBE_MULSNX16C(B,xx,rr);
            // reduced add
            r_summ = BBE_RADDNX40C(B);
            // type conversion
            B=BBE_MOVNX40_FROMC40(r_summ);
            vTmp1= BBE_MOVVWL(B);
            //-----------------------------------------------------------------
            vTmp0= BBE_SELNX16I(vTmp1,vTmp0,BBE_SELI_PACK_4);
            BBE_SAVNX16_XP(vTmp0,Yrdwr_align,pYrdwr,2*8);
        }
        BBE_SAPOS_FP(Yrdwr_align,pYrdwr);
        pYrdwr= (xb_vecNx16*)pScr;
        __Pragma("loop_count min=1");
        for(l=0; l<(L>>3); l++)
        {
            BBE_LVNX16_IP(vTmp0,pYrdwr,2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(vTmp1,pYrdwr,2*BBE_SIMD_WIDTH);
            // get high 16 and low 16 bits of B
            BBE_DSELNX16I(th,tl,vTmp1,vTmp0,BBE_DSELI_DEINTERLEAVE_1);
            // load D
            BBE_LVNX16_IP(d0,pDrd,2*BBE_SIMD_WIDTH);
            q = BBE_MOVVSV(d0,0);
            d0= BBE_SHFLNX16I(d0,BBE_SHFLI_DUPLICATE_1_EVEN);
            q = BBE_SHFLVSNI(q, BBE_VSA_SHFLI_DUPLICATE_1_ODD);
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
            BBE_SPNX16_IP(vTmp0,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,1);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,2);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,3);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,4);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,5);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,6);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
            vTmp1= BBE_REPNX16C(vTmp0,7);
            BBE_SPNX16_IP(vTmp1,pXwr,2*SX);
        }
        // set pDtr to prev element
        pDrd= (const xb_vecNx16*)XT_ADDX2(-2*L*2,(uintptr_t)pDrd);
        //__Pragma("no_reorder")
    }
    #undef N
    #undef SR
    #undef SX
    #undef SD
}

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

void cqr_bkw8x8x1n (void *pScr,
                          complex_fract16* _X,
                    const complex_fract16* _R,
                    const complex_fract16* _D,
                    int qABX,
                    int L)
{
          int16_t* X=(      int16_t*)_X;
    const int16_t* R=(const int16_t*)_R;
    const int16_t* D=(const int16_t*)_D;
    #define N 8
    #define M 8
    #define SR 128
    #define SX 16
    #define SD 16

    int k;
    int l, L8, Ltail;
    
    xb_vecNx16 tl,th,dd,d0,b_res,b;
    xb_vecNx40 B;
    xb_vecNx16 xx0, rr0;
    xb_c40 r_summ;
    vsaN q;
    xb_vecNx16 cmp0;
    vboolN mask0;

    const xb_vecNx16* restrict pXrd_first;
    const xb_vecNx16* restrict pXrd;
          xb_vecNx16* restrict pXwr;
    const xb_vecNx16* restrict pRrd;
    const xb_vecNx16* restrict pDrd;

    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);

    const vsaN sh16= BBE_MOVVSA32(16);
    const vsaN vqABX= BBE_MOVVSA32(qABX);

    Ltail = L&7;
    L8 = L - Ltail;
    if (0!=L8)
    {
        size_t D_tr_size= 2*L8*N*sizeof(int16_t);
        int16_t *D_tr= (int16_t*)pScr;
        pScr= ((char*)pScr + D_tr_size);
        transformD8(D_tr,D,L8);
        cqr_bkw8x8x1n_8(X,R,D_tr,qABX,L8,pScr);
    }

    if (0!=Ltail)
    {
        X += SX*L8;
        R += SR*L8;
        D += SD*L8;
        for (k=N-1; k>=0; k--)
        {
            pXrd_first= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)X);
            pXrd= (const xb_vecNx16*)(X);
            pXwr= (      xb_vecNx16*)pXrd_first;
            pRrd= (const xb_vecNx16*)(R+(k*N)*2);
            pDrd= (const xb_vecNx16*)XT_ADDX4(k,(uintptr_t)D);
            // set mask for correct loading of data
            cmp0 = BBE_SEQNX16();
            cmp0 = BBE_SHFLNX16I(cmp0, BBE_SHFLI_DOUBLE_1_LO);
            mask0 = BBE_LTNX16(BBE_MOVVA16(k), cmp0);

            __Pragma("loop_count min=1");
            for(l=0; l<Ltail; l++)
            {
                // calculate y(m,:)-R(m,:)*X, 1xP
                // load 16 bit complex value (16 bit re and 16 bit im)
                BBE_LPNX16_XP(b,pXrd_first,2*SX);
                B= BBE_UNPKSNX16(b);
                B=BBE_SLSNX40(B,vqABX);

                BBE_LVNX16T_XP(rr0,pRrd,2*SR,mask0);
                BBE_LVNX16T_XP(xx0,pXrd,2*SX,mask0);
                //B-= rr*xx
                BBE_MULSNX16C(B,xx0,rr0);

                // reduced add
                r_summ = BBE_RADDNX40C(B);
                // type conversion
                B=BBE_MOVNX40_FROMC40(r_summ);
                // get low 16 bits of B
                tl = BBE_PACKLNX40(B);
                // get high 16 bits of B
                th = BBE_PACKVNX40(B,sh16);
                // load D
                BBE_LPNX16_XP(dd,pDrd,2*SD);
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
                BBE_SPNX16_XP(b_res,pXwr,2*SX);
            }
        }
    }
    #undef N
    #undef M
    #undef SR
    #undef SX
    #undef SD
} /* cqr_bkw8x8x1n() */
#endif

size_t cqr_bkw8x8x1n_getScratchSize (int M,int N,int P,int L)
{
    int L8;
    size_t D_tr_size= 0;
    size_t Intrmd_size= 0;
    NASSERT(M==8 && N==8 && P==1 && L>0);
    L8 = L - (L&7);
    if (0!=L8)
    {
        D_tr_size= 2*L8*N*sizeof(int16_t);
        Intrmd_size= L8*2*sizeof(int32_t);
    }
    return (D_tr_size + Intrmd_size);
    return 0;
} /* cqr_bkw8x8x1n_getScratchSize() */
