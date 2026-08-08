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
#include "matinvgjnxnnf_common.h"
#include "common.h"

#if HAVE_VFPU

/*----------------------------------------------
compute reciprocal of main diagonal and replace 
main diagonal with ones

Input/output:
Z[L][SA]  L matrices (only each SA element is 
          updated)
Output:
recip[L]  reciprocals of main diagonal elements

----------------------------------------------*/
void matinvgjnxnnf_recipPivot(float32_t* recip, float32_t *z, int L, int SA)
#if 0
{
    int l;
        for (l=0; l<L; l++)
        {
            /* process pivot row */
            recip[l]=1.0f/z[SA*l]; z[SA*l]=1.f;
        }
}
#else
{
          xtfloat* restrict pRecip;
    const xtfloat* restrict pZkk0;
    const xtfloat* restrict pZkk1;
          xtfloat* restrict pZwr;
    pRecip=(      xtfloat*)recip;
    pZkk0 =(const xtfloat*)z;
    pZkk1 =(const xtfloat*)z;
    pZwr  =(      xtfloat*)z;
    int l;
    for (l=0; l<L; l++)
    {
        xb_vecN_2xf32 zkk0,zkk1,t;
        BBE_LSN_2XF32_XP(zkk0,pZkk0,SA*sizeof(float32_t));
        BBE_LSN_2XF32_XP(zkk1,pZkk1,SA*sizeof(float32_t));
        {
            xb_vecN_2xf32 vt1, vt2, vt3;
            vt1=BBE_RECIP0N_2XF32 (zkk0);
            vt2=BBE_CONSTN_2XF32 (1);
            BBE_MULSN_2XF32 (vt2, zkk0, vt1);
            BBE_MULANN_2XF32(vt1, vt1, vt2);
            vt3=BBE_CONSTN_2XF32 (1);
            BBE_MULSN_2XF32 (vt3, zkk1, vt1);
            BBE_MULANN_2XF32 (vt1, vt1, vt3);
            t=vt1;
        }
        BBE_SSN_2XF32_IP(t,pRecip,sizeof(float32_t));
    }
    __Pragma("no_reorder")
    for (l=0; l<L; l++)
    {
        BBE_SSN_2XF32_XP(BBE_CONSTN_2XF32(1),pZwr,SA*sizeof(float32_t));
    }
}
#endif

#endif
