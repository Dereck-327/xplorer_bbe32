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
  NatureDSP_Baseband library. Vector Mathematics
    Dual Peak Search
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_vector.h"

/*-------------------------------------------------------------------------
Dual Peak Search 

Description: These functions retrieve the maximum (minimum) and next-to-
maximum (next-to-minimum) values of vector elements. They output both the 
peak values and their indices. 

Representation:
vmax,vmin    16-bit signed fixed-point format
             Special values for vmax and vmin (see note 1) are -32768 and
             32767, respectively
vmaxf,vminf  IEEE-754 Std. single precision floating-point format
             Special values for vmaxf and vminf are -HUGE_VALF and +HUGE_VALF,
             respectively

Notes:
1. Each kind of dual-peak search function reserves a special value (see
   above) to maintain internal invariants during the search process. However,
   it is still legal for the input vector to contain special values among
   other data. If this is the case, then the following peculiarities should
   be considered:
     A) Input vector elements equal to the special value are ignored.
     B) If the total number K of input vector elements distinct from the special
        value is less than 2 (i.e. K==0 or K==1), then the last 2-K entries of output
        vector idx[2] are assigned zero, and the last 2-K elements of output vector
        m[2] are assigned the special value.
2. If the peak value is encountered more than once in the input data vector, then
   it will be reported twice in m[2], and idx[2] will contain indices of the first
   two occurencies of the peak value.
3. For floating-point functions NAN values in vectors are ignored. If, however,
   the input vector contains NANs only, then the entries of output vector idx[2] are
   assigned zero, and the elements of output vector m[2] are assigned the special value.
Parameters:
Input:
x[N]     Input data vector
N        Length of input data vector
Output:
m[2]     2 peak values is descending (vmax) or ascending (vmin) order
idx[2]   Indices of 2 peak elements; optional

Restrictions:
x        Aligned on 32-byte boundary
x,m,idx  Must not overlap
N        Multiple of 16 (vmax, vmin) or 8 (maxf, xminf)
-------------------------------------------------------------------------*/

void vmin(int16_t * restrict m,
    int16_t * restrict idx,
    const int16_t * restrict x, int N)
{
#ifdef BBE_SETDUALMAX

    static const uint16_t ALIGN(32) FlagsV[BBE_SIMD_WIDTH] = { 0x0003, 0x000c, 0x0030, 0x00c0, 0x0300, 0x0c00, 0x3000, 0xc000, 0x0000, 0x000c, 0x0030, 0x00c0, 0x0300, 0x0c00, 0x3000, 0xc000 };
    const xb_vecNx16 * restrict pFlagsV = (const xb_vecNx16 *)FlagsV;
    xb_vecNx16 flagsv;

    int n;
    int32_t ALIGN(32) PeakVal[BBE_SIMD_WIDTH];
    int32_t ALIGN(32) PeakInd[BBE_SIMD_WIDTH];
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    const xb_vecNx16 * pPeakVal = (const xb_vecNx16 *)(PeakVal);
    const xb_vecNx16 * pPeakInd = (const xb_vecNx16 *)(PeakInd);
    xb_vecNx16 x0, x1;
    xb_vecNx40 xw, peakv, idxv, zm, id, z1;
    vboolN bLane01, bLane23, bLane45, bLane67, _b1, i0;
    vboolN bLane89, bLaneAB, bLaneCD, bLaneEF;
    int Peak0, Peak1, Ind0, Ind1;

    if (N <= 0)
    {
        if (m)
        {
            m[0] = (int16_t)0x7fff;
            m[1] = (int16_t)0x7fff;
        }
        if (idx)
        {
            idx[0] = 0;
            idx[1] = 0;
        }
        return;
    }
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(N > 0 && N % BBE_SIMD_WIDTH == 0);

    // init registers BBE_MAX and BBE_MAX2 by INTMIN value
    // init registers BBE_MAXIDX and BBE_MAXIDX2 by zero
    BBE_SETDUALMAX(0xffff8001);

    for (n = 0; n < (N >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        BBE_LVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH);

        xw = BBE_UNPKSNX16(x0); // to avoid overflow at inversing values 32767 and -32768 use wide vector
        xw = BBE_NEGNX40(xw);   // to find minimum value we will find maximum value in array of -x[i]

        BBE_DUALMAXWNX32(xw, 0); // update registers: BBE_MAX, BBE_MAX2, BBE_MAXIDX and BBE_MAXIDX2
    }

    // init bool flags to get register BBE_MAXIDX
    flagsv = BBE_LVNX16_I(pFlagsV, 0);
    bLane01 = BBE_EXTRBN(flagsv, 0);
    bLane23 = BBE_EXTRBN(flagsv, 1);
    bLane45 = BBE_EXTRBN(flagsv, 2);
    bLane67 = BBE_EXTRBN(flagsv, 3);
    bLane89 = BBE_EXTRBN(flagsv, 4);
    bLaneAB = BBE_EXTRBN(flagsv, 5);
    bLaneCD = BBE_EXTRBN(flagsv, 6);
    bLaneEF = BBE_EXTRBN(flagsv, 7);

    _b1 = BBE_EXTRBN(flagsv, 8);
    _b1 = BBE_NOTB(_b1);              // set "true" to all lanes

    //zm = BBE_MOVWINT40(N - 1);            //  replicate (N-1)
    //z1 = BBE_MOVWINT40(0xffff8001);       // replicate "MININT"

    // get registers BBE_MAX and BBE_MAX_IDX
    PeakVal[0] = (int32_t)RUR_BBE_MAX_0();
    PeakVal[1] = (int32_t)RUR_BBE_MAX_1();
    PeakVal[2] = (int32_t)RUR_BBE_MAX_2();
    PeakVal[3] = (int32_t)RUR_BBE_MAX_3();
    PeakVal[4] = (int32_t)RUR_BBE_MAX_4();
    PeakVal[5] = (int32_t)RUR_BBE_MAX_5();
    PeakVal[6] = (int32_t)RUR_BBE_MAX_6();
    PeakVal[7] = (int32_t)RUR_BBE_MAX_7();
    PeakInd[0] = (int32_t)BBE_SELMAXIDX(bLane01, 0);
    PeakInd[1] = (int32_t)BBE_SELMAXIDX(bLane23, 0);
    PeakInd[2] = (int32_t)BBE_SELMAXIDX(bLane45, 0);
    PeakInd[3] = (int32_t)BBE_SELMAXIDX(bLane67, 0);
    PeakInd[4] = (int32_t)BBE_SELMAXIDX(bLane89, 0);
    PeakInd[5] = (int32_t)BBE_SELMAXIDX(bLaneAB, 0);
    PeakInd[6] = (int32_t)BBE_SELMAXIDX(bLaneCD, 0);
    PeakInd[7] = (int32_t)BBE_SELMAXIDX(bLaneEF, 0);

    // move state register BBE_MAX2 to BBE_MAX
    // move state register BBE_MAXIDX2 to BBE_MAXIDX
    BBE_MOVDUALMAXT(_b1);

    // get registers BBE_MAX2 and BBE_MAX_IDX2
    PeakVal[8] = (int32_t)RUR_BBE_MAX_0();
    PeakVal[9] = (int32_t)RUR_BBE_MAX_1();
    PeakVal[10] = (int32_t)RUR_BBE_MAX_2();
    PeakVal[11] = (int32_t)RUR_BBE_MAX_3();
    PeakVal[12] = (int32_t)RUR_BBE_MAX_4();
    PeakVal[13] = (int32_t)RUR_BBE_MAX_5();
    PeakVal[14] = (int32_t)RUR_BBE_MAX_6();
    PeakVal[15] = (int32_t)RUR_BBE_MAX_7();
    PeakInd[8] = (int32_t)BBE_SELMAXIDX(bLane01, 0);
    PeakInd[9] = (int32_t)BBE_SELMAXIDX(bLane23, 0);
    PeakInd[10] = (int32_t)BBE_SELMAXIDX(bLane45, 0);
    PeakInd[11] = (int32_t)BBE_SELMAXIDX(bLane67, 0);
    PeakInd[12] = (int32_t)BBE_SELMAXIDX(bLane89, 0);
    PeakInd[13] = (int32_t)BBE_SELMAXIDX(bLaneAB, 0);
    PeakInd[14] = (int32_t)BBE_SELMAXIDX(bLaneCD, 0);
    PeakInd[15] = (int32_t)BBE_SELMAXIDX(bLaneEF, 0);


    BBE_LVNX16_XP(x0, pPeakVal, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(x1, pPeakVal, -2 * BBE_SIMD_WIDTH);
    peakv = BBE_MOVSWV(x1, x0);
    BBE_LVNX16_XP(x0, pPeakInd, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(x1, pPeakInd, -2 * BBE_SIMD_WIDTH);
    idxv = BBE_MOVSWV(x1, x0);

    xw = BBE_RMAXNX40_STEP0(peakv); // first minimum (negative maximum) value
    xw = BBE_RMAXNX40_STEP1(xw);
    Peak0 = BBE_MOVAW32(xw);
    xw = BBE_REPNX40(xw, 0);
    i0 = BBE_EQNX40(xw, peakv);
    if (m)
    {
        m[0] = (-Peak0); // inverse "-x[i]" to get correct result
    }

    zm = N - 1;                               //  replicate (N-1)
    id = BBE_MOVNX40T(idxv, zm, i0);
    id = BBE_RMINNX40_STEP0(id); //first index
    id = BBE_RMINNX40_STEP1(id);
    Ind0 = BBE_MOVAW32(id);
    id = BBE_REPNX40(id, 0);
    i0 = BBE_EQNX16(id, idxv);
    if (idx)
    {
        idx[0] = (Ind0);
    }

    z1 = BBE_MOVWA40(0x000000ff, 0xffff8001); // replicate "MININT"
    peakv = BBE_MOVNX40T(z1, peakv, i0);
    xw = BBE_RMAXNX40_STEP0(peakv); // second minimum (negative maximum) value
    xw = BBE_RMAXNX40_STEP1(xw);
    Peak1 = BBE_MOVAW32(xw);
    xw = BBE_REPNX40(xw, 0);
    i0 = BBE_EQNX40(xw, peakv);
    if (m)
    {
        m[1] = (-Peak1);
    }

    BBE_LVNX16_XP(x0, pPeakInd, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(x1, pPeakInd, -2 * BBE_SIMD_WIDTH);
    idxv = BBE_MOVSWV(x1, x0);
    zm = N - 1;                               //  replicate (N-1)
    id = BBE_MOVNX40T(idxv, zm, i0);
    id = BBE_RMINNX40_STEP0(id); // second index
    id = BBE_RMINNX40_STEP1(id);
    Ind1 = BBE_MOVAW32(id);
    if (idx)
    {
        idx[1] = (Ind1);
    }

#else // OLD (bitexact code but it don't use instructions "dual max")

    int n;
    unsigned idn0, idn1;
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;

    xb_vecNx16  x0, m0, x1, m1, mx;
    xb_vecNx16  z0, z1, z2, zm;
    xb_vecNx16  id0, id1, tmp, id, id01;
    vboolN      i0, i1, i2, i_equ, i_ind;
    xb_int16    a, b, q0, q1;

    if (N <= 0)
    {
        if (m)
        {
            m[0] = (int16_t)0x7fff;
            m[1] = (int16_t)0x7fff;
        }
        if (idx)
        {
            idx[0] = 0;
            idx[1] = 0;
        }
        return;
    }
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(N>0 && N % BBE_SIMD_WIDTH == 0);

    z0 = 0;
    z1 = BBE_MOVVINX16(BBE_MOVVI_INT16_MAXINT);
    zm = (N - 1);              //  replicate (N-1)
    id = BBE_SEQNX16();
    z2 = BBE_MOVVINT16(BBE_SIMD_WIDTH);
    m0 = m1 = z1;
    id0 = id1 = 0;
    idn0 = idn1 = 0;
    for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        BBE_LVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_BMINNX16(i0, mx, x0, m0);
        x1 = BBE_MOVNX16T(m0, x0, i0);

        BBE_MOVIDXNX16F(tmp, id0, idn0, i0);
        BBE_MOVIDXNX16T(id0, id0, idn1, i0);

        BBE_BMINNX16(i1, m1, x1, m1);
        id1 = BBE_MOVNX16T(tmp, id1, i1);

        // n*BBE_SIMD_WIDTH+0,....n*BBE_SIMD_WIDTH+BBE_SIMD_WIDTH-1
        m0 = mx;
    }
    a = BBE_RMINNX16(m0); //first maximum value
    x0 = BBE_MOVNX16_FROM16(a);
    x0 = BBE_REPNX16(x0, 0);
    i0 = BBE_EQNX16(x0, m0);
    id = BBE_MOVNX16T(id0, zm, i0);
    q0 = BBE_RMINNX16(id); //first index

    x0 = BBE_MOVNX16_FROM16(q0);
    x0 = BBE_REPNX16(x0, 0);
    i0 = BBE_EQNX16(x0, id0);

    tmp = BBE_MOVNX16T(z1, m0, i0);
    i_equ = BBE_EQNX16(tmp, m1);         // flags: "min values is equivalent"
    BBE_BMINNX16(i_ind, id01, id0, id1); // flags: "min index of min value"
    i_equ = BBE_ANDB(i_equ, i_ind);      // flags: "choise min value with min index from id0 (else from id1)"

    BBE_BMINNX16(i2, m1, tmp, m1);
    i2 = BBE_ORB(i2, i_equ);             // use all flags
    b = BBE_RMINNX16(m1); //second maximum value 
    x1 = BBE_MOVNX16_FROM16(b);
    x1 = BBE_REPNX16(x1, 0);
    i1 = BBE_EQNX16(x1, m1);
    id1 = BBE_MOVNX16T(id0, id1, i2);
    id1 = BBE_MOVNX16T(id1, zm, i1);
    q1 = BBE_RMINNX16(id1); //second index

    if (m)
    {
        m[0] = (a);
        m[1] = (b);
    }
    if (idx)
    {
        idx[0] = (q0);
        idx[1] = (q1);
    }

#endif
} /* vmin() */
