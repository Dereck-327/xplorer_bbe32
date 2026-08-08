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
    BBE32 code for Cholesky backward recursion, block format
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


// get allocated space per one matrix
static int getSpace(int S)
{
    int m;
    m=30-XT_NSA(S);
    m=XT_MIN(m,(LOG2_BBE_SIMD_WIDTH-1));
    // round up to the  next multiple of 8 or lesser degree of 2
    S=(((S-1)>>m)+1)<<m;
    return S;
}

/*----------------------------------------------------------------------------------
   reversing R matrices for easier readings by rows (diagonal elements are omitted):
   original R    transformed R
   0 1 3 6 a     d 8 c 4 7 b 1 3 6 a
     2 4 7 b
       5 8 c
         9 d
           e

   Input:
   R[L][SR]        L input matrices
   Rt[L*N*(N-1)]   stream of L trasposed matrices
----------------------------------------------------------------------------------*/
void cholnTransformR(int16_t* Rt,const int16_t* R,int N,int L)
#if 0
{
    int SR=2*getSpace((N*(N+1))>>1);
    int n,m,l;
    const int16_t* pR;
    for(l=0; l<L; l++)
    {
        for (n=0; n<N; n++)
        {
            pR=R+l*SR+(N-n)*(N-n+3)-2;
            for (m=0; m<n; m++)
            {
                Rt[0]=pR[0];
                Rt[1]=pR[1];
                pR+=2*(N-n+m+1);
                Rt+=2;
            }
        }
    }
}
#else
{
    const int SR= 2*getSpace((N*(N+1))>>1);
    int n,m,l;
    const int16_t* restrict pR;
    const int16_t* restrict pR2= R;
    int delta, delta2, dpR;
    const int delta3= 2*(N+1)*2;
    xb_vecNx16 vTmp;
    xb_vecNx16 * pRt;
    valign v_Rt;
    __Pragma("loop_count min=1")
    for(l=0; l<L; l++)
    {
        delta2= delta3;
        dpR= N;
        __Pragma("loop_count min=4")
        for (n=0; n<N; n++)
        {
            pR=(const int16_t*)XT_ADDX2((dpR)*(dpR+3)-2,(uintptr_t)pR2);
            delta= delta2;
            for (m=0; m<n; m++)
            {
                BBE_LPNX16_XP(vTmp,pR,delta);
                BBE_SPNX16_IP(vTmp,Rt,2*2);
                delta= XT_ADDI_N(delta, 2*2);
            }
            delta2= XT_ADDI(delta2, -2*2);
            dpR= XT_ADDI(dpR,-1);
        }
        pR2=(const int16_t*)XT_ADDX2(SR,(uintptr_t)pR2);
    }
    vTmp = BBE_ZERONX16();
    pRt = (xb_vecNx16 *)Rt;
    v_Rt = BBE_ZALIGN();
    BBE_SANX16_IP(vTmp, v_Rt, pRt);
    BBE_SAPOS_FP(v_Rt, pRt);
}
#endif

#endif
