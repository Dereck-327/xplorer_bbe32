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
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
  NatureDSP_Baseband library. Matrix Operations
    Complex Matrix-Matrix/Matrix-Vector Multiply
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"

#if !(HAVE_MULPC && HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, cmatmul3x3n,(complex_fract16 * restrict z, 
            const complex_fract16 * restrict x, 
            const complex_fract16 * restrict y, 
            int L, int Q))
#else
/*
value of 1 enables use another algorithm to calculate z[8].
It works faster, but without rounding
*/
/*-------------------------------------------------------------------------
Complex Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of complex matrices or vectors. Both the block order and 
streaming order are allowed for input/output matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand complex matrices
y[L*Sy]     Sequence of right-hand complex matrices
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices 
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of complex result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, 3x3*3x3->3x3, Sx=16, Sy=16, Sz=16
   Restrictions:
     None
*/
void cmatmul3x3n(complex_fract16 * restrict z,
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y,
             int L, int Q )
{
    static
    const uint16_t ALIGN(32) aSel[BBE_SIMD_WIDTH]= {	0x0404,	0x0505,	0x0a04,	0x0b05,	0x1004,	0x1105,	0xff0a,	0xff0b,	0xff0a,	0xff0b,	0xff0a,	0xff0b,	0xff10,	0xff11,	0xff10,	0xff11};

    const xb_vecNx16 * restrict pXrd = (const xb_vecNx16 *)x;
    const xb_vecNx16 * restrict pYrd = (const xb_vecNx16 *)y;
    xb_vecNx16 * restrict pZwr = (xb_vecNx16 *)z;
    
    int l;
    xb_vecNx16 vTmp;
    xb_vecNx16 vX0, vX1, vY0, vY1, vZ;
    xb_vecNx16 vXsel, vYsel;
    vselN vsX2;
    vselN vsY_last_sel;
    xb_c40 wvR;
    xb_vecNx40 wvZ;
    const vsaN  q= BBE_MOVVSA32(Q);

    /* check restrictions */
    NASSERT_ALIGN(x,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y,(2*BBE_SIMD_WIDTH));
    NASSERT_ALIGN(z,(2*BBE_SIMD_WIDTH));
    NASSERT(Q>=0 && Q<=16);
    if (L <= 0) return;

    vTmp= BBE_LVNX16_I((const xb_vecNx16*)aSel,0);
    vsX2= BBE_MOVVSELNX16(vTmp,0);
    vsY_last_sel= BBE_MOVVSELNX16(vTmp,8);

    __Pragma( "loop_count min=1" )
    for ( l=0; l<L; l++ )
    {
        /* load entire X matrix */
        BBE_LVNX16_IP(vX0,pXrd,2*BBE_SIMD_WIDTH);
        BBE_LPNX16_IP(vX1,pXrd,2*BBE_SIMD_WIDTH);
        /* load entire Y matrix */
        BBE_LVNX16_IP(vY0,pYrd,2*BBE_SIMD_WIDTH);
        BBE_LPNX16_IP(vY1,pYrd,2*BBE_SIMD_WIDTH);

        vXsel= BBE_SHFLNX16I(vX0,BBE_SHFLI_MMC3X3X3X3_OFFSET_M1_0_STEP_1);
        vYsel= BBE_SHFLNX16I(vY0,BBE_SHFLI_MMC3X3X3X3_OFFSET_M2_0_STEP_1);

        wvZ= BBE_MULRNX16C(vXsel, vYsel, q);

        vXsel= BBE_SHFLNX16I(vX0,BBE_SHFLI_MMC3X3X3X3_OFFSET_M1_0_STEP_2);
        vYsel= BBE_SHFLNX16I(vY0,BBE_SHFLI_MMC3X3X3X3_OFFSET_M2_0_STEP_2);
        BBE_MULANX16C(wvZ, vXsel, vYsel);

        vXsel= BBE_SELNX16(vX1,vX0,vsX2);
        vYsel= BBE_SELNX16I(vY1,vY0,BBE_SELI_MMC3X3X3X3_OFFSET_M2_0_STEP_3);
        BBE_MULANX16C(wvZ, vXsel, vYsel);

        /* pack and save res 8 elements */
        vZ = BBE_PACKVNX40(wvZ, q);
        BBE_SVNX16_IP(vZ, pZwr, 2*BBE_SIMD_WIDTH);
        
        /* calculate 8-th element */
        vXsel= BBE_SELNX16I(vX1, vX0,BBE_SELI_ROTATE_LEFT_4);
        vYsel= BBE_SELNX16(vY1, vY0, vsY_last_sel);
        wvZ= BBE_MULNX16C(vXsel, vYsel);

        wvR= BBE_RADDNX40C(wvZ);
        wvZ= BBE_MOVNX40_FROMC40(wvR);

        /* pack and save res 8 elements */
        vZ = BBE_PACKVNX40(wvZ, q);
        BBE_SVNX16_IP(vZ, pZwr, 2*BBE_SIMD_WIDTH);
    }
} /* cmatmul3x3n() */
#endif
