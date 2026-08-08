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
    Complex Matrix Gauss-Jordan inversion, floating point complex data, block 
    format
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#include "NatureDSP_types.h"
#include "cmatinvgjnxnnf_common.h"
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
void cmatinvgjnxnnf_recipPivot(complex_float* recip, complex_float *z, int L, int SA)
#if 0
{
    int l;
    for (l=0; l<L; l++)
    {
        /* process pivot row */
        recip[l]=recipc(z[SA*l]); z[SA*l]=makecomplexf(1.f,0.f);
    }
}
#else
{
          xtcomplexfloat* restrict pRecip;
    const xtcomplexfloat* restrict pZkk0;
    const xtcomplexfloat* restrict pZkk1;
          xtcomplexfloat* restrict pZwr;
    pRecip=(      xtcomplexfloat*)recip;
    pZkk0 =(const xtcomplexfloat*)z;
    pZkk1 =(const xtcomplexfloat*)z;
    pZwr  =(      xtcomplexfloat*)z;
    int l;
    for (l=0; l<L; l++)
    {
        xb_vecN_4xcf32 zkk0,zkk1,t;
        BBE_LSN_4XCF32_XP(zkk0,pZkk0,SA*sizeof(complex_float));
        BBE_LSN_4XCF32_XP(zkk1,pZkk1,SA*sizeof(complex_float));
        { 
            xb_vecN_2xf32 a, c;
            xb_vecN_2xf32 t1, t2, t5, t6, t7, t8;
            c=BBE_MOVN_2XF32_FROMN_4XCF32(zkk0);
            t6=BBE_MULN_2XF32 (c, c);
            t7=BBE_MULMN_2XF32 (c, c, 0, 3);
            t8=BBE_ADDN_2XF32 (t6, t7);
            t1=BBE_RECIP0N_2XF32 (t8);
            t2=BBE_CONSTN_2XF32 (1);
            BBE_MULSN_2XF32  (t2, t8, t1);
            BBE_MULANN_2XF32 (t1, t1, t2);
            t5=BBE_CONSTN_2XF32 (1);
            BBE_MULSN_2XF32  (t5, t8, t1);
            BBE_MULANN_2XF32 (t1, t1, t5);
            c=BBE_MOVN_2XF32_FROMN_4XCF32(zkk1);
            a=BBE_MULMN_2XF32 (t1, c, 2, 12);
            t=BBE_MOVN_4XCF32_FROMN_2XF32(a);
        }
        //t=BBE_RECIPN_4XCF32(zkk0);
        BBE_SSN_4XCF32_IP(t,pRecip,sizeof(complex_float));
    }
    __Pragma("no_reorder")
    for (l=0; l<L; l++)
    {
        BBE_SSN_4XCF32_XP(1.f,pZwr,SA*sizeof(complex_float));
    }
}
#endif

#endif
