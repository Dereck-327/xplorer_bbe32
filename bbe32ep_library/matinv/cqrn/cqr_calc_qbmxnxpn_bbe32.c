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
    Update right side of equations for QR process for block ordered matrices.
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

#if HAVE_CQRN

#define sz_i16 sizeof(int16_t)

// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    // compute multiple of next degree of 2 (max multiple is 32 for real, 16 for complex, 16 for real+dbl, 8 for complex+dbl )
    m=30-XT_NSA(S);
    if (m>(LOG2_BBE_SIMD_WIDTH-1)) m=LOG2_BBE_SIMD_WIDTH-1;
    // round up to the  next multiple of 32 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}

/*
    update matrix B[L][SB] by housholder vectors V[L][SV]

    Input:
    M,N,P,L     dimensions
    V[L][SV]    Housholder vectors
    Input/output:
    B[L][SV]    B matrices MxP

    Temporary:
    Z[2*P*L]
   
*/
static void cqrnUpdateB(int16_t *Z, int16_t* B,const int16_t* V,int M,int N,int P,int SB,int L)
{
    const xb_vecNx16 * restrict pV;
    const xb_vecNx16 * restrict pBld;
          xb_vecNx16 * restrict pBst;
    const xb_vecNx16 * restrict pZld;
          xb_vecNx16 * restrict pZst;
    xb_vecNx16 V0, B0, Z0, _1Q14;
    xb_vecNx40 ACC0;
    valign al_V, al_Z, al_Bld, al_Bst;
    int l, m, p;
    vsaN rnd14, shift;
    
    NASSERT_ALIGN(Z,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    
    rnd14 = BBE_MOVVSA32(14);
    _1Q14 = BBE_MOVVINT16(4);
    _1Q14 = BBE_SLLINX16(_1Q14, 12);
    al_Bst = BBE_ZALIGN();
    al_Z   = BBE_ZALIGN();
    
    /* Special case for P=1 */
    m=((M+(BBE_SIMD_WIDTH/2-1))/(BBE_SIMD_WIDTH/2));
    if (P==1 )
    {
        int nbytes=M*4;
        int step_m=0;
        int count=0,modulo=(m*0x010000)+1;

        pBld=(const xb_vecNx16 *)(B);
        pBst=(      xb_vecNx16 *)(B);
        pV  =(const xb_vecNx16 *)(V);
        pZst=(      xb_vecNx16 *)(Z);

        al_Bld=BBE_LA_PP(pBld);
        al_V=BBE_LA_PP(pV);
        xb_vecNx40 Z1=0;
        for (l=0; l<L*m; l++)
        {   
            count=BBE_ADDMOD16U(count,modulo); //count=(count+1)&m;

            BBE_LAVNX16_XP(V0, al_V, pV, nbytes);            
            al_Bld=BBE_LA_PP(pBld);
            BBE_LAVNX16_XP(B0, al_Bld, pBld, nbytes);

            ACC0=BBE_MULNX16J(B0,V0);
            ACC0=BBE_MOVNX40_FROMC40(BBE_RADDNX40C(ACC0));
            ACC0=BBE_ADDNX40(ACC0,Z1);
            Z1=ACC0;

            ACC0=BBE_RNDADJNX40(ACC0, rnd14); 
            Z0=BBE_PACKVNX40(ACC0,rnd14);            

            BBE_SPNX16_I(Z0,pZst,0);

            step_m=0;
            XT_MOVEQZ(step_m, 1, count);//step_m= (count==0)? 1:0;
            pZst =(     xb_vecNx16*)XT_ADDX4( step_m, (uintptr_t)pZst);
                       
            step_m=step_m*40;
            shift=BBE_MOVVSA32(step_m);
            Z1=BBE_SRLNX40(Z1, shift);  

            XT_MOVEQZ(step_m, (SB-2*M), count);//step_m= (count==0)? 1:0;
            pBld =(const xb_vecNx16*)XT_ADDX2( step_m, (uintptr_t)pBld);
            pBst =(      xb_vecNx16*)XT_ADDX2( step_m, (uintptr_t)pBst );  
            
            XT_MOVEQZ(step_m, (m*BBE_SIMD_WIDTH), count);//step_m= (count==0)? 1:0;
            nbytes=XT_ADDX2(-BBE_SIMD_WIDTH,nbytes);
            nbytes=XT_ADDX2(step_m,nbytes);         
        }
            
        nbytes=M*4;

        pBld=(const xb_vecNx16 * )(B);
        pBst=(      xb_vecNx16 * )(B);

        pZld = (const xb_vecNx16*)Z;
        pV   = (const xb_vecNx16*)V;

        al_V=BBE_LA_PP(pV);
        count=0;
        for (l=0; l<L*m; l++)
        { 
            count=BBE_ADDMOD16U(count,modulo); //count=(count+1)&m;
            
            al_Bld=BBE_LA_PP(pBld);
            BBE_LAVNX16_XP(B0, al_Bld, pBld, nbytes); 
            BBE_LAVNX16_XP(V0, al_V, pV, nbytes);
            BBE_LPNX16_XP(Z0, pZld, 0);
            Z0=BBE_REPNX16C(Z0, 0);

            ACC0=BBE_MULRNX16(B0,_1Q14,rnd14);
            BBE_MULSNX16C(ACC0,Z0,V0);
            B0=BBE_PACKVNX40(ACC0,rnd14);

            al_Bst=BBE_ZALIGN();
            BBE_SAVNX16_XP(B0, al_Bst, pBst, nbytes);            
            BBE_SAPOS_FP( al_Bst, pBst); 
                       
            step_m=0;
            XT_MOVEQZ(step_m, 1, count);//step_m= (count==0)? 1:0;
            pZld =(const xb_vecNx16*)XT_ADDX4( step_m, (uintptr_t)pZld);
            
            XT_MOVEQZ(step_m, (SB-2*M), count);//step_m= (count==0)? 1:0;
            pBld =(const xb_vecNx16*)XT_ADDX2( step_m, (uintptr_t)pBld);
            pBst =(      xb_vecNx16*)XT_ADDX2( step_m, (uintptr_t)pBst );  
            
            XT_MOVEQZ(step_m, (m*BBE_SIMD_WIDTH), count);//step_m= (count==0)? 1:0;
            nbytes=XT_ADDX2(-BBE_SIMD_WIDTH,nbytes);
            nbytes=XT_ADDX2(step_m,nbytes);            
        }
        return;
    }

    /* Generic case: for P>1 */
    /* apply Housholder vector to the matrix B */
    if (P>1)
    {
        int nbytes=P*4;
        int next_row=-1, test_bytes=-1;
        int step_m=0;
        int count_0,modulo_0;

        p=((P+7)/8);
        count_0=0;
        modulo_0=(p*0x010000)+1;

        pBld=(const xb_vecNx16 * )(B);
        pZst=(      xb_vecNx16 * )(Z);
        pV  =(const xb_vecNx16 * )(V);
        al_V=BBE_LA_PP(pV);

        for (l=0; l<L*p; l++)
        { 
            count_0=BBE_ADDMOD16U(count_0,modulo_0); //count=(count+1)&m; 
             
            next_row=nbytes;
            test_bytes=XT_ADDX2( -BBE_SIMD_WIDTH, nbytes);
            XT_MOVGEZ(next_row, 32, test_bytes);//count_bytes= (nbytes>=0)? 32:nbytes; 
            next_row=4*P-next_row;

            BBE_LPNX16_IP(V0, pV, 4);
            V0= BBE_REPNX16C(V0, 0);
            al_Bld=BBE_LA_PP(pBld);
            BBE_LAVNX16_XP(B0, al_Bld, pBld, nbytes); 
            ACC0=BBE_MULRNX16J(B0,V0,rnd14); 
            pBld=(const xb_vecNx16*)XT_ADD( next_row, (uintptr_t)pBld);
            for(m=0; m<M-1; m++)
            {            
                BBE_LPNX16_IP(V0, pV, 4);
                V0= BBE_REPNX16C(V0, 0);
                al_Bld=BBE_LA_PP(pBld);
                BBE_LAVNX16_XP(B0, al_Bld, pBld, nbytes); 
                BBE_MULANX16J(ACC0,B0,V0);
                pBld=(const xb_vecNx16*)XT_ADD( next_row, (uintptr_t)pBld);
            }                
            Z0=BBE_PACKVNX40(ACC0,rnd14);
            BBE_SVNX16_XP(Z0,pZst,2*BBE_SIMD_WIDTH);

            step_m=0;
            XT_MOVEQZ(step_m, (p*4), count_0);//step_m= (count_0==0)? 1:0;

            nbytes=XT_ADDX2( -BBE_SIMD_WIDTH, nbytes);
            nbytes=XT_ADDX8( step_m, nbytes);

            XT_MOVEQZ(step_m, (SB-2*P), count_0);//step_m= (count_0==0)? 1:0;
            pBld=(const xb_vecNx16*)XT_ADD( (-next_row-4*P*(M-1)), (uintptr_t)pBld);//next column            
            pBld=(const xb_vecNx16*)XT_ADDX2( step_m, (uintptr_t)pBld );//next matrix

            step_m=-M;
            XT_MOVEQZ(step_m, 0, count_0);//step_m= (count_0==0)? 0:-M;           
            pV=(const xb_vecNx16*)XT_ADDX4( step_m, (uintptr_t)pV);             
        }
    
        nbytes=P*4;
        count_0=0;

        pBld=(const xb_vecNx16 *)(B);
        pBst=(      xb_vecNx16 *)(B);
        pZld=(const xb_vecNx16 *)(Z);
        pV  =(const xb_vecNx16 *)(V);

        for (l=0; l<L*p; l++)
        { 
            count_0=BBE_ADDMOD16U(count_0,modulo_0); //count=(count+1)&m; 
             
            next_row=nbytes;
            test_bytes=XT_ADDX2( -BBE_SIMD_WIDTH, nbytes);
            XT_MOVGEZ(next_row, 32, test_bytes);//count_bytes= (nbytes>=0)? 32:nbytes; 
            next_row=4*P-next_row;
                
            BBE_LVNX16_XP(Z0, pZld, 2*BBE_SIMD_WIDTH);
            for(m=0; m<M; m++)
            {            
                BBE_LPNX16_XP(V0, pV, 4);
                V0= BBE_REPNX16C(V0, 0);                
                al_Bld=BBE_LA_PP(pBld);
                BBE_LAVNX16_XP(B0, al_Bld, pBld, nbytes);   

                ACC0=BBE_MULRNX16(B0,_1Q14,rnd14);
                BBE_MULSNX16C(ACC0,Z0,V0);
                B0=BBE_PACKVNX40(ACC0,rnd14);
                            
                al_Bst=BBE_ZALIGN();
                BBE_SAVNX16_XP(B0, al_Bst, pBst, nbytes);            
                BBE_SAPOS_FP( al_Bst, pBst); 

                pBld=(const xb_vecNx16*)XT_ADD( next_row, (uintptr_t)pBld);
                pBst=(      xb_vecNx16*)XT_ADD( next_row, (uintptr_t)pBst);             
            }  
            step_m=0;
            XT_MOVEQZ(step_m, (p*4), count_0);//step_m= (count_0==0)? 1:0;

            nbytes=XT_ADDX2( -BBE_SIMD_WIDTH, nbytes);           
            nbytes=XT_ADDX8( step_m, nbytes); 
    
            XT_MOVEQZ(step_m, (SB-2*P), count_0);//step_m= (count_0==0)? 1:0;
            pBld=(const xb_vecNx16*)XT_ADD( (-next_row-4*P*(M-1)), (uintptr_t)pBld);//next column            
            pBld=(const xb_vecNx16*)XT_ADDX2( step_m, (uintptr_t)pBld );//next matrix
                        
            pBst=(      xb_vecNx16*)XT_ADD( (-next_row-4*P*(M-1)), (uintptr_t)pBst);//next column            
            pBst=(      xb_vecNx16*)XT_ADDX2( step_m, (uintptr_t)pBst );//next matrix

            step_m=-M;
            XT_MOVEQZ(step_m, 0, count_0);//step_m= (count_0==0)? 0:-M;           
            pV=(const xb_vecNx16*)XT_ADDX4( step_m, (uintptr_t)pV);
        }
    }
}

/*-------------------------------------------------------------------------
Update right side of equations for QR process for block ordered matrices.
Matrix sizes SB,SV are selected as usual for complex block ordered matrix 
sequencies, i.e. total size is rounded up to the closest bigger multiple of 
BBE_SIMD_WIDTH/2==8 elements or, if it is lesser, to the closest bigger 
multiple of degree of 2. 
SB=size(M*P)
SV=size(((2*M-N+1)*N/2+N)*L)
Scratch size in bytes is defined by cqr_calc_qbmxnn_getScratchSize(M,N,P,L)
functions

Input:
 M, N, P      dimensional parameters
 L            Number of matrices
Input/output:
 B[L][SB]     On input it is the sequence of L complex matrices B. 
              At the end of the process, matrices Z replace input
              matrices A. In a case of non-square matrices (N!=M), 
              only N*P elements of each output matrix will be valid.
Input:
 V[SV]        Sequence of L Housholder rotation vectors 

Restrictions:
1. B, V, pScr must not overlap
2. All matrices and the scratch must not overlap and must be aligned 
   on 32-byte boundary 
3. Number of matrices L must be positive
4. M and N must be positive multiples of 4
5. N <= M
---------------------------------------------------------------------------*/
void cqr_calc_qbmxnxpn (void *pScr,
                          complex_fract16* _B,
                    const complex_fract16* _V,
                    int M, int N, int P,
                    int L)
{
        int16_t* B=(      int16_t*)_B;
  const int16_t* V=(const int16_t*)_V;

/*
Reference code:
% compute Q'B matrix
% input:
% V  - sequence of Housholder vectors  [(2*M-N+1)*N/2,1]
% Fi - common rotation diagonal matrix [Nx1]
% R  - upper triangle decomposition
function [B] = cqr_calcQB(B,V,Fi)
[M, P] = size(B); 
[N, t] = size(Fi);
Z=zeros(M,P);
im=1;
for m=1:N
    v=V(im:im+M-m);
    im=im+(M-m+1);
    Bm=B(m:end,:);
    Bm=(Bm-2*v*v'*Bm);
    B(m:end,:)=Bm;
end
B=diag([Fi;ones(M-N,1)])'*B;
*/
    int16_t* Z=(int16_t*)pScr;
    int m;
    int SB=2*getSpace(M*P);
    const int16_t* pV;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B   ,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V   ,2*BBE_SIMD_WIDTH);
    NASSERT(L>0);
    NASSERT(N%4==0 && N>0);
    NASSERT(M%4==0 && M>0);
    if (P<1) return;

    {
        int sz=(int)cqr_calc_qbmxnxpn_getScratchSize(M, N, P, L);
        xb_vecNx16 * pWR = (xb_vecNx16 *)pScr;
        xb_vecNx16 veczero = BBE_ZERONX16();
        for (m = 0; m < sz; m += 2 * BBE_SIMD_WIDTH)
        {
            BBE_SVNX16_IP(veczero, pWR, 2 * BBE_SIMD_WIDTH);
        }
    }

    /* scale down B by 1 bit right */
    cqrnScaleB(B,SB*L);

    /* apply Housholder vectors */
    for (pV=V,m=0; m<N; m++)
    {
        cqrnUpdateB(Z,B+2*m*P,pV,(M-m),N,P,SB,L); 
        pV+=2*(M-m)*L;
    }

    /* Rotate matrix B by diagonal matrix Fi' */
    if (P>1)    cqrnRotateB(B,V+(2*M-N+1)*N*L,N,P,SB,L);
    else        cqrnRotateB1(pScr,B,V+(2*M-N+1)*N*L,N,SB,L);
} /* cqr_calc_qbmxnxpn() */

size_t cqr_calc_qbmxnxpn_getScratchSize (int M, int N,int P,int L)
{
    size_t Zsize,Bsize;
    NASSERT(L>0);
    NASSERT(N%4==0 && N>0);
    NASSERT(M%4==0 && M>0);
    (void)M;
    if (P==1) Zsize = L*2*sizeof(int16_t);
    else      Zsize = L*2*sizeof(int16_t)*((P+7)/8)*BBE_SIMD_WIDTH/2;
    Bsize = N*2*sizeof(int16_t)*BBE_SIMD_WIDTH/2;
    return XT_MAX(Zsize,Bsize);
} /* cqr_calc_qbmxnxpn_getScratchSize() */
#else
DISCARD_FUN(void,cqr_calc_qbmxnxpn,(void *pScr,
                          complex_fract16* B,
                    const complex_fract16* V,
                    int M, int N,int P,
                    int L))
size_t cqr_calc_qbmxnxpn_getScratchSize(int M, int N,int P,int L)
{
    (void)M;(void)N;(void)P;(void)L;
    return 0;
}
#endif
