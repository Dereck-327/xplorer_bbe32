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
    Real Matrix Gauss-Jordan inversion, floating point real data, block 
    format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"
#include <math.h>
#include <float.h>

/*
    Reference Matlab code:
    % Gauss-Jordan inversion with full pivoting
    function A=gj_inversion1(A)
    sz=size(A);
    N=sz(1);
    swap_row=[];
    swap_col=[];
    for k=1:N
        % pivoting
         T=abs(A);
        T(swap_col,:)=0;
        T(:,swap_col)=0;
        [amax,imax]=max(reshape(T,1,N*N));
        pivot_col=floor((imax-1)/N)+1;
        pivot_row=rem(imax-1,N)+1;
        % swap rows
        t=A(pivot_row,:); A(pivot_row,:)=A(pivot_col,:); A(pivot_col,:)=t;
        swap_row(k)=pivot_row;
        swap_col(k)=pivot_col;
        % process pivot row
        norm=1./A(pivot_col,pivot_col);
        A(pivot_col,pivot_col)=1;
        A(pivot_col,:)=A(pivot_col,:)*norm;
        % elimination
        for i=1:N
            if(i==pivot_col) continue; end
            t=A(i,pivot_col);
            A(i,pivot_col)=0;
            A(i,:)=A(i,:)-A(pivot_col,:)*t;
        end
    end
    % final permutation of columns
    for k = N:-1:1
        t = A(:, swap_row(k));
        A(:, swap_row(k))=A(:, swap_col(k));
        A(:, swap_col(k))=t;
    end

*/

#if HAVE_VFPU

/* swap k-th row with rows given in srow[L] */
static void swapRows(float32_t * z, const int32_t* srow, int N,int k,int L)
#if 0
{
    int l,n;
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(srow,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 && N>1);
    for (l=0; l<L; l++)
    {
        int prow=srow[l];
        /* swap rows */
        for(n=0; n<N; n++) 
        {
            float32_t t;
            t=z[l+L*(prow*N+n)];
            z[l+L*(prow*N+n)]=z[l+L*(k*N+n)];
            z[l+L*(k*N+n   )]=t;
        }
    }
}
#else
{
    int n,l;
    xb_vecN_2xf32* restrict pAw;
    xtfloat* restrict pA0;
    xtfloat* restrict pA1;
    xtfloat* restrict pA2;
    xtfloat* restrict pA3;
    xtfloat* restrict pA4;
    xtfloat* restrict pA5;
    xtfloat* restrict pA6;
    xtfloat* restrict pA7;
    const xb_vecNx16 * restrict pRow;
    pRow=(const xb_vecNx16 *)srow;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx40 w;
        xb_vecNx16 xrow,addr,al,ah; 
        BBE_LVNX16_IP(xrow,pRow,2*BBE_SIMD_WIDTH);
        xrow=BBE_SELNX16I(xrow,xrow,BBE_SELI_EXTRACT_1_OF_2_OFF_0);
        w=BBE_SEQNX40();
        BBE_MULANX16(w,xrow,L*N);
        w=BBE_SLLINX40(w,2);
        w=BBE_ADDNX40(w,(int32_t)(z+l));
        al=BBE_PACKLNX40(w);
        ah=BBE_PACKVNX40(w,16);
        addr=BBE_SELNX16I(ah,al,BBE_SELI_INTERLEAVE_1_LO);
        pA0=(xtfloat *)BBE_EXTRNX16C(addr,0);//(xtfloat *)&z[l+0+L*N*srow[l+0]]; 
        pA1=(xtfloat *)BBE_EXTRNX16C(addr,1);//(xtfloat *)&z[l+1+L*N*srow[l+1]]; 
        pA2=(xtfloat *)BBE_EXTRNX16C(addr,2);//(xtfloat *)&z[l+2+L*N*srow[l+2]]; 
        pA3=(xtfloat *)BBE_EXTRNX16C(addr,3);//(xtfloat *)&z[l+3+L*N*srow[l+3]]; 
        pA4=(xtfloat *)BBE_EXTRNX16C(addr,4);//(xtfloat *)&z[l+4+L*N*srow[l+4]]; 
        pA5=(xtfloat *)BBE_EXTRNX16C(addr,5);//(xtfloat *)&z[l+5+L*N*srow[l+5]]; 
        pA6=(xtfloat *)BBE_EXTRNX16C(addr,6);//(xtfloat *)&z[l+6+L*N*srow[l+6]]; 
        pA7=(xtfloat *)BBE_EXTRNX16C(addr,7);//(xtfloat *)&z[l+7+L*N*srow[l+7]]; 

        pAw=(xb_vecN_2xf32*)&z[l+L*k*N];
        for (n = 0; n<N; n++)
        {
            xb_vecN_2xf32 t0,t1,t2,t3,t4,t5,t6,t7,T0,T1;
            T1=BBE_LVN_2XF32_I(pAw,0);
            t0=BBE_LSN_2XF32_I(pA0,0);
            t1=BBE_LSN_2XF32_I(pA1,0);
            t2=BBE_LSN_2XF32_I(pA2,0);
            t3=BBE_LSN_2XF32_I(pA3,0);
            t4=BBE_LSN_2XF32_I(pA4,0);
            t5=BBE_LSN_2XF32_I(pA5,0);
            t6=BBE_LSN_2XF32_I(pA6,0);
            t7=BBE_LSN_2XF32_I(pA7,0);
            t0=BBE_SELN_2XF32I(t1,t0,BBE_SELI_PACK_2);
            t2=BBE_SELN_2XF32I(t3,t2,BBE_SELI_PACK_2);
            t4=BBE_SELN_2XF32I(t5,t4,BBE_SELI_PACK_2);
            t6=BBE_SELN_2XF32I(t7,t6,BBE_SELI_PACK_2);
            t0=BBE_SELN_2XF32I(t2,t0,BBE_SELI_PACK_4);
            t4=BBE_SELN_2XF32I(t6,t4,BBE_SELI_PACK_4);
            T0=BBE_SELN_2XF32I(t4,t0,BBE_SELI_PACK_8);
            BBE_SVN_2XF32_XP(T0,pAw,L*sizeof(float32_t));
            BBE_SSN_2XF32_XP(T1,pA0,L*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_2 ),pA1,L*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_4 ),pA2,L*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_6 ),pA3,L*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_8 ),pA4,L*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_10),pA5,L*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_12),pA6,L*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_14),pA7,L*sizeof(float32_t));
        }
    }
}
#endif

/* swap k-th column with columns given in scol[L] */
static void swapCols(float32_t * z, const int32_t* scol, int N,int k,int L)
#if 0
{
    int l,n;
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(scol,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 && N>1);
    for (l=0; l<L; l++)
    {
        int pcol;
        pcol=scol[l];
        for(n=0; n<N; n++) 
        {
            float32_t t;
            t=z[l+L*(n*N+pcol)];
            z[l+L*(n*N+pcol)]=z[l+L*(n*N+k)];
            z[l+L*(n*N+k   )]=t;
        }
    }
}
#else
{
    int n,l;
    xb_vecN_2xf32* restrict pAw;
    xtfloat* restrict pA0;
    xtfloat* restrict pA1;
    xtfloat* restrict pA2;
    xtfloat* restrict pA3;
    xtfloat* restrict pA4;
    xtfloat* restrict pA5;
    xtfloat* restrict pA6;
    xtfloat* restrict pA7;
    const xb_vecNx16 * restrict pCol;
    pCol=(const xb_vecNx16 *)scol;
    for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
    {
        xb_vecNx40 w;
        xb_vecNx16 xcol,addr,al,ah; 
        BBE_LVNX16_IP(xcol,pCol,2*BBE_SIMD_WIDTH);
        xcol=BBE_SELNX16I(xcol,xcol,BBE_SELI_EXTRACT_1_OF_2_OFF_0);
        w=BBE_SEQNX40();
        BBE_MULANX16(w,xcol,L);
        w=BBE_SLLINX40(w,2);
        w=BBE_ADDNX40(w,(int32_t)(z+l));
        al=BBE_PACKLNX40(w);
        ah=BBE_PACKVNX40(w,16);
        addr=BBE_SELNX16I(ah,al,BBE_SELI_INTERLEAVE_1_LO);
        pA0=(xtfloat *)BBE_EXTRNX16C(addr,0);//(xtfloat *)&z[l+0+L*scol[l+0]]; 
        pA1=(xtfloat *)BBE_EXTRNX16C(addr,1);//(xtfloat *)&z[l+1+L*scol[l+1]]; 
        pA2=(xtfloat *)BBE_EXTRNX16C(addr,2);//(xtfloat *)&z[l+2+L*scol[l+2]]; 
        pA3=(xtfloat *)BBE_EXTRNX16C(addr,3);//(xtfloat *)&z[l+3+L*scol[l+3]]; 
        pA4=(xtfloat *)BBE_EXTRNX16C(addr,4);//(xtfloat *)&z[l+4+L*scol[l+4]]; 
        pA5=(xtfloat *)BBE_EXTRNX16C(addr,5);//(xtfloat *)&z[l+5+L*scol[l+5]]; 
        pA6=(xtfloat *)BBE_EXTRNX16C(addr,6);//(xtfloat *)&z[l+6+L*scol[l+6]]; 
        pA7=(xtfloat *)BBE_EXTRNX16C(addr,7);//(xtfloat *)&z[l+7+L*scol[l+7]]; 

        pAw=(xb_vecN_2xf32*)&z[l+L*k];
        for (n = 0; n<N; n++)
        {
            xb_vecN_2xf32 t0,t1,t2,t3,t4,t5,t6,t7,T0,T1;
            T1=BBE_LVN_2XF32_I(pAw,0);
            t0=BBE_LSN_2XF32_I(pA0,0);
            t1=BBE_LSN_2XF32_I(pA1,0);
            t2=BBE_LSN_2XF32_I(pA2,0);
            t3=BBE_LSN_2XF32_I(pA3,0);
            t4=BBE_LSN_2XF32_I(pA4,0);
            t5=BBE_LSN_2XF32_I(pA5,0);
            t6=BBE_LSN_2XF32_I(pA6,0);
            t7=BBE_LSN_2XF32_I(pA7,0);
            t0=BBE_SELN_2XF32I(t1,t0,BBE_SELI_PACK_2);
            t2=BBE_SELN_2XF32I(t3,t2,BBE_SELI_PACK_2);
            t4=BBE_SELN_2XF32I(t5,t4,BBE_SELI_PACK_2);
            t6=BBE_SELN_2XF32I(t7,t6,BBE_SELI_PACK_2);
            t0=BBE_SELN_2XF32I(t2,t0,BBE_SELI_PACK_4);
            t4=BBE_SELN_2XF32I(t6,t4,BBE_SELI_PACK_4);
            T0=BBE_SELN_2XF32I(t4,t0,BBE_SELI_PACK_8);
            BBE_SVN_2XF32_XP(T0,pAw,L*N*sizeof(float32_t));
            BBE_SSN_2XF32_XP(T1,pA0,L*N*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_2 ),pA1,L*N*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_4 ),pA2,L*N*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_6 ),pA3,L*N*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_8 ),pA4,L*N*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_10),pA5,L*N*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_12),pA6,L*N*sizeof(float32_t));
            BBE_SSN_2XF32_XP(BBE_SELN_2XF32I(T1,T1,BBE_SELI_ROTATE_RIGHT_14),pA7,L*N*sizeof(float32_t));
        }
    }
}
#endif

/*-------------------------------------------------------------------------
Inversion of Stream Ordered Matrices By Gauss-Jordan Algortihm

Description: perform in-place inversion of real/complex matrices by Gauss-
Jordan elimination method, with full pivoting. The algorithm is applied to
a sequence of input matrices stored in stream order.

Inversion result is not defined for a close to singular input matrix.

Data format: IEEE-754 Std single precision floating-point

Temporary:
  pScr       Scratch area. Required size (in bytes) is defined by 
             functions [c]matinvgj<size>sf_getScratchSize(N,L)
Input:
  N          Matrix size
  L          Number of matrices
Input/Output:
  A[N*N][L]  Input matrices, inverted matrices on output
Restrictions:
  pScr,A     Must not overlap and must be aligned on 32-byte boundary 
  N          Must be greater than 1
  L          Must be a multiple of 8 for real-valued functions, or a multiple
             of 4 for complex-valued functions
---------------------------------------------------------------------------*/
void matinvgjnxnsf ( void* pScr, float32_t * restrict z, int N, int L  )
#if 0
{
    float32_t t;
    float32_t *maxVal;
    int32_t *srow,*scol;
    int l,k,m,n;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 && N>1);

    if(L<=0 || N<=1) return;

    maxVal=(float32_t*)pScr;
    srow=(int32_t*)(maxVal+L);
    scol=srow+N*L;

    for (k=0; k<N; k++)
    {
        /* pivoting: search the absolute maximum and its position in the matrix 
            at all positions excluding previously used */
        for (l=0; l<L; l++)
        {
            float32_t maxval=FLT_MIN;
            int prow,pcol;
            prow=k; pcol=k;
            for (m=k; m<N; m++)
            {
                for(n=k; n<N; n++) 
                {
                    t=fabsf(z[l+L*(m*N+n)]);  
                    if(t>maxval) { maxval=t; prow=m; pcol=n; }
                }
            }
            srow[k*L+l]=prow; scol[k*L+l]=pcol;
        }
        swapRows(z, srow+k*L, N,k,L);
        swapCols(z, scol+k*L, N,k,L);
        /* reciprocal of main diagonal */
        for (l=0; l<L; l++) maxVal[l]=1.0f/z[l+L*(k*N+k)];

        /* process pivot row */
        for(l=0; l<L; l++)
        {
            float32_t t;
            t=maxVal[l];
            for (n=0; n<N; n++) 
            {
                z[l+L*(k*N+n)]*=t;
            }
        }
        for(l=0; l<L; l++) z[l+L*(k*N+k)]=maxVal[l];

            /* GS elimination : for all rows excluding the pivot one */
        for (m=0; m<N; m++)
        {
            if(m==k) continue; 
            for(l=0; l<L; l++)
            {
                float32_t t;
                t=z[l+L*(m*N+k)];
                maxVal[l]=t;
                for (n=0; n<N; n++)
                {
                     z[l+L*(m*N+n)] -= z[l+L*(k*N+n)] * t; 
                }
            }
            for(l=0; l<L; l++) z[l+L*(m*N+k)] = -z[l+L*(k*N+k)] * maxVal[l];
        }
    }
    
    /* final reverse permulation of columns  */
    for (k=N-2; k>=0; k--)
    {
        swapCols(z, srow+k*L, N,k,L);
        swapRows(z, scol+k*L, N,k,L);
    }
}
#else
{
    const xb_vecN_2xf32* restrict pZrd;
    const xb_vecN_2xf32* restrict pZrw;
          xb_vecN_2xf32* restrict pZwr;
    xb_vecN_2xc16 * restrict pRow;
    xb_vecN_2xc16 * restrict pCol;
          xb_vecN_2xf32* restrict pMaxVal;
    float32_t *maxVal;
    int32_t *srow,*scol;
    int l,k,m,n;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(z,2*BBE_SIMD_WIDTH);
    NASSERT(z);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 && N>1);

    if(L<=0 || N<=1) return;

    maxVal=(float32_t*)pScr;
    srow=(int32_t*)(maxVal+L);
    scol=srow+N*L;

    for (k=0; k<N; k++)
    {
        if (k!=N-1)
        {
            /* pivoting: search the absolute maximum and its position in the matrix 
                at all positions excluding previously used */
            pRow=(xb_vecN_2xc16 *)(srow+k*L);
            pCol=(xb_vecN_2xc16 *)(scol+k*L);
            pZrd=(const xb_vecN_2xf32*)(z+L*(k*N+k));
            for (l=0; l<L; l+=BBE_SIMD_WIDTH/2)
            {
                xb_vecN_2xf32 maxval=FLT_MIN,t;
                xb_vecN_2xc16 prow,pcol,vm,vn;
                prow=k; pcol=k;
                vm=k;
                __Pragma("loop_count min=2")
                for (m=k; m<N; m++)
                {
                    vn=k;
                    __Pragma("loop_count min=2")
                    for(n=k; n<N; n++) 
                    {
                        vboolN_2 b;
                        BBE_LVN_2XF32_XP(t,pZrd,L*sizeof(float32_t));
                        t=BBE_ABSN_2XF32(t);
                        b=BBE_OGTN_2XF32(t,maxval);
                        maxval=BBE_MAXNUMN_2XF32(t,maxval);
                        prow=BBE_MOVN_2XC16T(vm,prow,b);
                        pcol=BBE_MOVN_2XC16T(vn,pcol,b);
                        vn=BBE_ADDN_2XC16(vn,1);
                    }
                    vm=BBE_ADDN_2XC16(vm,1);
                    pZrd+= (k*L*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
                }
                BBE_SVN_2XC16_IP(prow,pRow,2*BBE_SIMD_WIDTH);
                BBE_SVN_2XC16_IP(pcol,pCol,2*BBE_SIMD_WIDTH);
                pZrd+= ((BBE_SIMD_WIDTH/2-N*(N-k)*L)*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
            }

            swapRows(z, srow+k*L, N,k,L);
            swapCols(z, scol+k*L, N,k,L);
        }
        /* reciprocal of main diagonal */
        pMaxVal=(xb_vecN_2xf32*)(maxVal);
        pZrd   =(const xb_vecN_2xf32*)&z[L*(k*N+k)];
        __Pragma("loop_count min=1")
        for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
        {
            xb_vecN_2xf32 t;
            BBE_LVN_2XF32_IP(t,pZrd,2*BBE_SIMD_WIDTH);
            t=BBE_RECIPN_2XF32(t);
            BBE_SVN_2XF32_IP(t,pMaxVal,2*BBE_SIMD_WIDTH);
        }
        /* process pivot row */
        __Pragma("no_reorder")
        pMaxVal=(xb_vecN_2xf32*)(maxVal);
        pZrd   =(const xb_vecN_2xf32*)&z[L*N*k];
        pZwr   =(      xb_vecN_2xf32*)&z[L*N*k];
        __Pragma("loop_count min=1")
        for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
        {
            xb_vecN_2xf32 x,t;
            BBE_LVN_2XF32_IP(t,pMaxVal,2*BBE_SIMD_WIDTH);
            __Pragma("loop_count min=2")
            for (n=0; n<N; n++) 
            {
                BBE_LVN_2XF32_XP(x,pZrd,L*sizeof(float32_t));
                x=BBE_MULN_2XF32(x,t);
                BBE_SVN_2XF32_XP(x,pZwr,L*sizeof(float32_t));
            }
            pZrd+= ((BBE_SIMD_WIDTH/2-L*N)*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
            pZwr+= ((BBE_SIMD_WIDTH/2-L*N)*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
        }
        __Pragma("no_reorder")
        pMaxVal=(xb_vecN_2xf32*)(maxVal);
        pZwr   =(xb_vecN_2xf32*)&z[L*(k*N+k)];
        __Pragma("loop_count min=1")
        for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
        {
            xb_vecN_2xf32 t;
            BBE_LVN_2XF32_IP(t,pMaxVal,2*BBE_SIMD_WIDTH);
            BBE_SVN_2XF32_IP(t,pZwr   ,2*BBE_SIMD_WIDTH);
        }

            /* GS elimination : for all rows excluding the pivot one */
        __Pragma("no_reorder")
        for (m=0; m<N; m++)
        {
            xb_vecN_2xf32 x,y,t;
            if(m==k) continue; 
            pZrd=(const xb_vecN_2xf32*)&z[L*N*k];
            pZwr=(      xb_vecN_2xf32*)&z[L*N*m];
            pZrw=(const xb_vecN_2xf32*)pZwr;
            __Pragma("loop_count min=1")
            for (l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH-1)); l++)
            {
                t = BBE_LVN_2XF32_X  (pZrw,k*L*sizeof(float32_t));
                __Pragma("loop_count min=2")
                for (n=0; n<N; n++)
                {
                    BBE_LVN_2XF32_XP (x,pZrd,L*sizeof(float32_t));
                    BBE_LVN_2XF32_XP (y,pZrw,L*sizeof(float32_t));
                    BBE_MULSN_2XF32(y,x,t);
                    BBE_SVN_2XF32_XP (y,pZwr,L*sizeof(float32_t));
                }
                x = BBE_LVN_2XF32_X  (pZrd,(k-N)*L*sizeof(float32_t));
                y=BBE_MULMN_2XF32(x,t,3,12);
                BBE_SVN_2XF32_X  (y,pZwr,(k-N)*L*sizeof(float32_t));
                pZrd+= ((BBE_SIMD_WIDTH/2-N*L)*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
                pZwr+= ((BBE_SIMD_WIDTH/2-N*L)*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
                pZrw+= ((BBE_SIMD_WIDTH/2-N*L)*sizeof(float32_t))/sizeof(xb_vecN_2xf32);
            }
        }
    }
    
    /* final reverse permulation of columns  */
    for (k=N-2; k>=0; k--)
    {
        swapCols(z, srow+k*L, N,k,L);
        swapRows(z, scol+k*L, N,k,L);
    }
}
#endif

size_t matinvgjnxnsf_getScratchSize (int N, int L  )
{
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 && N>1);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return sizeof(int32_t)*N*2*L+sizeof(float32_t)*L;
}
#else
DISCARD_FUN(void, matinvgjnxnsf,( void* pScr, float32_t * restrict z, int N, int L  ))

size_t matinvgjnxnsf_getScratchSize (int N, int L  )
{
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0 && N>1);
    if(L<=0 || N<=0) return 0;
    (void)L; (void)N;
    return 0;
}
#endif
