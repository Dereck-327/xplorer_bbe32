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
    NatureDSP_Baseband library. Fitting and Interpolation Routines
    Polynomial Fitting and Interpolation for Real Data
    C code optimized for BBE32EP core with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Fixed point arithmetic. */
#include "NatureDSP_Math.h"
/* Common utility declarations. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fit.h"
#include "pfitf_common.h"

#if HAVE_VFPU

#include <math.h>
#include <float.h>

/*-------------------------------------------------------------------------
Polynomial Fitting and Interpolation

Description: the pfit functions fit (in least squares sense) a degree N 
polynomial to input data sampled at a points grid of length M, and use
that polynomial to interpolate data at arbitrary query points. Namely,
the pfitN_grid() functions compute the Vandermonde matrix for the sample
points grid and perform the Cholesky decomposition of that matrix, the
pfitN_process() functions calculate the least squares solution for the
polynomial coefficients. Finally, the pfitN_eval() functions evaluate
the polynomial at query points.

Please refer to the NatureDSP Baseband Library Reference for full details
on these functions.

Representation:
pfit_gridN,      16-bit fixed-point data. Parameter specifications denote
pfit_processN,   fixed-point format for various data items
pfit_evalN       
pfitf_gridN,     IEEE-754 Std single precision floating-point data
pfitf_processN,  
pfitf_evalN     

Note:
Number of fractional bits specidied for various input/output arguments below apply
for the fixed-point variant

Parameters:
Input:
N                     Degree of polynomial, 1..6
M                     Number of sample points
P                     Number of query points
maxIter               Number of least squares solution enhancement iterations. Right 
                      choice depends on required accuracy, the ad-hoc value is (N+1)/2
x[M]                  Sample points grid, Q15 or floating point
y[M]                  Sampled data values, Q15 or floating point
xi[P]                 Query points, Q15 or floating point
M'=(M+7)&(~7), N'=8   for floating point API
M'=(M+15)&(~15),N'=16 for fixed-point API

Intermediate:
V[M'*8]               Vandermonde matrix, Q15 or floating point
R[N'*8]               Upper triangular Cholesky factor of matrix V, Q11 or floating point
Output:
yi[P]                 Data values interpolated at query points, Q15 or floating point
p[N+1]                Polynomial coefficients, Q8.23 or floating point
Temporary:
pScr                  Scratch memory area. To determine the scratch area size required by
                      a function pfitN_<fun>, use the respective helper function 
                      pfit_<fun>_getScratchSize(M,N)

Restrictions:
M>N                   The number of sample points must exceed the degree of polynomial
x,y,xi,yi,V,R,p,pScr  Must not overlap
V,R,pScr              Aligned on 32-byte boundary
---------------------------------------------------------------------------*/
#define RSTRIDE 8


/*-------------------------------------------------------------------------
   These functions apply Cholesky decomposition procedure to the sequence of 
   real matrices written in a streaming order. 
   R=pfit_chol(A'*A+sigma2*I)
   Fixed point representation of upper-diagonal matrix R is the same as of input. 

   Functions return nonzero if overflow is detected 

   NOTE:
   Data layout for matrices is selected as for other matrices written in a 
   streaming order. 

   Input:
   A[M*N]  input real matrix, Q15
   N       Matrix dimension (number of columns for MxN)
   M       Matrix dimension (number of rows for MxN)
   sigma2  noise estimate, Q31

   Output:
   R[N*N]  output real upper-triangle matrices , Q11

   Scratch:
   t[M+N]

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
---------------------------------------------------------------------------*/

/*-----------------------------
    Cholesky 
    Input:
    A[M'*N]
    M'=(M+7)&~7
    Input/Output:
    R[N*8]  

    Scratch:
    Am[M]
    Rm[N]
-----------------------------*/
static void choliter(
                const float32_t* restrict A,
                float32_t* restrict R,
                float32_t sigma2,
                int m,int M,int N)
#if 0
{
  float32_t D[1]; 
  int n,k;
  const int Astride = ((M + (BBE_SIMD_WIDTH/2 - 1)) & ~(BBE_SIMD_WIDTH/2 - 1));

  NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
  
  { /* take colunms of A and R and calculate diagonal elements */
    float32_t Acc;
    Acc=sigma2;
    /* TBD: optimize */
    for (k=0; k<M; k++)
    {
      float32_t a = A[m*Astride+k];
      Acc   += a*a;
    }
    for (k=0; k<m; k++) /* NOTE: this loop may be from 0 to N: for (k=0; k<N; k++) */
    {
      float32_t r = R[m*RSTRIDE+k];
      Acc   -= r*r;
    }
    if (Acc<=FLT_MIN) Acc=FLT_MIN;
    D[0]=1.f/sqrtf(Acc);
  }
  // find elements in row from m to N
  for (k=m; k<N; k++)
  {
    float32_t Acc;
    Acc=0.f;
    /* TBD: optimize */
    for (n=0; n<M; n++)
    {
      float32_t ak_re;
      float32_t am_re;
      ak_re = A[k*Astride+n];
      am_re = A[m*Astride+n];
      Acc += ak_re*am_re;
    }
    for (n=0; n<m; n++) 
    {
      float32_t rk_re,rm_re;
      rk_re = R[k*RSTRIDE+n];
      rm_re = R[m*RSTRIDE+n];
      Acc -= rk_re*rm_re;
    }
    R[k*RSTRIDE+m] = Acc*D[0];
  }
}
#else
{
    const xb_vecN_2xf32 *restrict pAk;
    const xb_vecN_2xf32 *restrict pAm;
          xb_vecN_2xf32 *restrict pR;
    valign aAk,aAm;
    xtfloat D; 
    xb_vecN_2xf32 AA,ak,am,rk,rm;
    int n,k;
    int K,nbytes;
    const int Astride = ((M + (BBE_SIMD_WIDTH/2 - 1)) & ~(BBE_SIMD_WIDTH/2 - 1));

    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
  
    pAk=(const xb_vecN_2xf32 *)&A[m*Astride];
    pR =(      xb_vecN_2xf32 *)&R[m*RSTRIDE];
    AA=BBE_ZERON_2XF32();
    K=(M-1)>>(LOG2_BBE_SIMD_WIDTH-1);
    nbytes=(M-(K<<(LOG2_BBE_SIMD_WIDTH-1)))*sizeof(float32_t);
    for (n=0; n<K; n++)
    {
        BBE_LVN_2XF32_XP(ak,pAk,sizeof(xb_vecN_2xf32));
        BBE_MULAN_2XF32(AA,ak,ak);
    }
    aAk=BBE_LAN_2XF32_PP(pAk);
    BBE_LAVN_2XF32_XP(ak,aAk,pAk,nbytes);
    BBE_MULAN_2XF32(AA,ak,ak);
    rm=BBE_LVN_2XF32_I(pR,0);
    BBE_MULSN_2XF32(AA,rm,rm);
    D=XT_RSQRT_S(XT_MAX_S(XT_ADD_S(BBE_RADDN_2XF32(AA),sigma2),FLT_MIN));

    for (k=m; k<N; k++)
    {
        xtfloat Acc;
        pAk=(const xb_vecN_2xf32*)(A+k*Astride);
        pAm=(const xb_vecN_2xf32*)(A+m*Astride);
        AA=BBE_ZERON_2XF32();
        for (n=0; n<K; n++)
        {
            BBE_LVN_2XF32_XP(ak,pAk,sizeof(xb_vecN_2xf32));
            BBE_LVN_2XF32_XP(am,pAm,sizeof(xb_vecN_2xf32));
            BBE_MULAN_2XF32(AA,ak,am);
        }
        aAk=BBE_LAN_2XF32_PP(pAk);
        aAm=BBE_LAN_2XF32_PP(pAm);
        BBE_LAVN_2XF32_XP(ak,aAk,pAk,nbytes);
        BBE_LAVN_2XF32_XP(am,aAm,pAm,nbytes);
        BBE_MULAN_2XF32(AA,ak,am);

        rk=BBE_LVN_2XF32_I(pR,0);
        BBE_MULSN_2XF32(AA,rk,rm);

        Acc=XT_MUL_S(BBE_RADDN_2XF32(AA),D);
        BBE_SSN_2XF32_X(BBE_MOVN_2XF32_FROMF32(Acc),(xtfloat*)pR,m*sizeof(float32_t));
        pR++;
    }
}
#endif

void pfitf_chol(float32_t* t,
            float32_t * restrict R, 
            const float32_t* restrict A, 
            float32_t sigma2,
            int M, int N)
#if 0
{
  int m;
  NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(t,2*BBE_SIMD_WIDTH);
  NASSERT(M>=N);
  NASSERT(M>0 && N>0);
  
  /* auto-vectorizable */
  for (m=0; m<RSTRIDE*N; m++) R[m]=0.f;
  for (m=0; m<N; m++)
  {
    choliter(A,R,sigma2,m,M,N);
  }
}
#else
{
    xb_vecN_2xf32 * pR=(xb_vecN_2xf32 *)R;
    int m;
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(t,2*BBE_SIMD_WIDTH);
    NASSERT(M>=N);
    NASSERT(M>0 && N>0);
  
    /* auto-vectorizable */
    for (m=0; m<8; m++) BBE_SVN_2XF32_IP(BBE_ZERON_2XF32(),pR,sizeof(xb_vecN_2xf32));
    for (m=0; m<N; m++)
    {
        choliter(A,R,sigma2,m,M,N);
    }
}
#endif


/*-----------------------------
    Cholesky forward recursion:
    y=R'\(A'*B)
    Input:
    A[M'*N]  
    B[M]     
    R[8*8]  
    M'       (M+7)&~7
    Output:
    Y[N]     
-----------------------------*/
void pfitf_cholfwd (
            float32_t* restrict y, 
            const float32_t* restrict R, 
            const float32_t* restrict D, 
            const float32_t* restrict A, 
            const float32_t* restrict B, 
            int M, int N)
#if 0
{
  int n;
  int m;
  const int Astride = (M+(BBE_SIMD_WIDTH/2-1))&~(BBE_SIMD_WIDTH/2-1);

  NASSERT(M>1);
  NASSERT(N>=2 && N<=7);
  NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);

  for(n=0; n<N; n++)    y[n]=0.f;

  for (n=0; n<N; n++)
  {
    // calculate A(:,n)'*B-Rn'*Y, 
    float32_t B_=0.f, C_=0.f;
    for (m=0; m<M; m++) /* TBD: optimize */
    {
      B_ += (A[n*Astride+m]*B[m]);
    }
    for (m=0; m<n+1; m++)   // this loop may be to N: for (m=0; m<N; m++) 
    {
      C_ += y[m]*R[n*RSTRIDE+m];
    }
    B_ -= C_;
    y[n]=B_*D[n];
  }
}
#else
{
    xb_vecN_2xf32 YY,RR,BB;
    valign aB;
    int nbytes;
    const xb_vecN_2xf32 *restrict pB;
    const xb_vecN_2xf32 *restrict pA;
    const xb_vecN_2xf32 *restrict pR;
    int n;
    int m;
    const int Astride = (M+(BBE_SIMD_WIDTH/2-1))&~(BBE_SIMD_WIDTH/2-1);

    NASSERT(M>1);
    NASSERT(N>=2 && N<=7);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);

    YY=BBE_ZERON_2XF32();
    pR=(const xb_vecN_2xf32 *)R;
    for (n=0; n<N; n++)
    {
        float32_t B_;
        vboolN_2 bN,bN0,bN1;
        bN0=BBE_LTRN_2(n);
        bN1=BBE_LTRN_2(n+1);
        bN=bN0^bN1; /* this boolean mask is rotating by 1 bit each time */
        // calculate A(:,n)'*B-Rn'*Y, 
        pB=(const xb_vecN_2xf32 *)B;
        pA=(const xb_vecN_2xf32 *)(&A[n*Astride]);
        BB=BBE_ZERON_2XF32();
        aB=BBE_LAN_2XF32_PP(pB);
        nbytes=M*sizeof(float32_t);
        BBE_LVN_2XF32_IP(RR,pR,sizeof(xb_vecN_2xf32));
        for (m=0; m<Astride>>(LOG2_BBE_SIMD_WIDTH-1); m++)
        {
            xb_vecN_2xf32 a,b;
            BBE_LAVN_2XF32_XP(b,aB,pB,nbytes);
            nbytes-=sizeof(xb_vecN_2xf32 );
            BBE_LVN_2XF32_IP(a,pA,sizeof(xb_vecN_2xf32 ));
            BBE_MULAN_2XF32(BB,a,b);
        }
        BBE_MULSN_2XF32(BB,YY,RR);
        B_=(float32_t)BBE_RADDN_2XF32(BB);
        B_*=D[n];
        BB=BBE_REPN_2XF32(B_,0);
        {
            xb_vecN_2xc16 y,b;
            y=BBE_MOVN_2XC16_FROMNX16(BBE_MOVNX16_FROMN_2XF32(YY));
            b=BBE_MOVN_2XC16_FROMNX16(BBE_MOVNX16_FROMN_2XF32(BB));
            y=BBE_MOVN_2XC16T(b,y,bN);
            YY=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_2XC16(y));
        }
    }
    BBE_SVN_2XF32_I(YY,(xb_vecN_2xf32 *)y,0);
}
#endif

/*-----------------------------
    Cholesky backward recursion:
    x=R\y
    Input:
    Rt[8*8]  transposed R
    y[N]    
    Output:
    x[N]    
-----------------------------*/
void  pfitf_cholbkw (
            float32_t* restrict x, 
            const float32_t* restrict Rt, 
            const float32_t* restrict D, 
            const float32_t* restrict y, 
            int N)
#if 0
{
  int n,m;
  NASSERT(N>=2 && N<=7);
  NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);

  /* auto-vectorizable */
  for (m=0; m<N; m++) x[m]=0.f;
  for (m=N-1; m>=0; m--)
  {
    // calculate y(m,:)-R(m,:)*X, 1xP
    float32_t A;
    A =y[m];
    for (n=m+1; n<N; n++)   // NOTE: loop may begin from 0: for (n=0; n<N; n++) 
    {
      A -= (x[n]*Rt[m*RSTRIDE+n]);   
    }
    A *= D[m];
    x[m]=A;
  }
}
#else
{
          xb_vecN_2xf32* restrict pX=(xb_vecN_2xf32* )x;
    const xtfloat*       restrict pY=(const xtfloat      *)(y+N-1);
    const xb_vecN_2xf32* restrict pR=(const xb_vecN_2xf32*)(&Rt[(N-1)*RSTRIDE]);
    valign aX;
    xb_vecN_2xf32 X;
    int m;
    NASSERT(N>=2 && N<=7);
    NASSERT_ALIGN(Rt,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(D,2*BBE_SIMD_WIDTH);

    X=BBE_ZERON_2XF32();
    for (m=N-1; m>=0; m--)
    {
        // calculate y(m,:)-R(m,:)*X, 1xP
        float32_t A;
        xb_vecN_2xf32 Y,R,T;
        vboolN_2 bN,bN0,bN1;
        bN0=BBE_LTRN_2(m);
        bN1=BBE_LTRN_2(m+1);
        bN=bN0^bN1; /* this boolean mask is rotating by 1 bit each time */
        BBE_LSN_2XF32_XP(Y,pY,-(int)sizeof(float32_t));
        T=X;
        BBE_LVN_2XF32_IP(R,pR,-(int)sizeof(xb_vecN_2xf32));
        BBE_MULSN_2XF32(Y,T,R);
        A=BBE_RADDN_2XF32(Y);
        A *= D[m];
        Y=BBE_REPN_2XF32(A,0);
        {
            xb_vecN_2xc16 y,x;
            y=BBE_MOVN_2XC16_FROMNX16(BBE_MOVNX16_FROMN_2XF32(Y));
            x=BBE_MOVN_2XC16_FROMNX16(BBE_MOVNX16_FROMN_2XF32(X));
            x=BBE_MOVN_2XC16T(y,x,bN);
            X=BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_2XC16(x));
        }
    }
    aX=BBE_ZALIGN();
    BBE_SAVN_2XF32_XP(X,aX,pX,N*sizeof(float32_t));
    BBE_SAN_2XF32POS_FP(aX,pX);
}
#endif

#endif
