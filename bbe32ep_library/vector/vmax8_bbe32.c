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
    Search for 8 Top Value Elements
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
Search for 8 Top Value Elements

Description: Functions search for 8 top values in the input vector. They
output both the peak values and respective indices.

Representation:
vmax8   16-bit signed fixed-point format
vmax8f  IEEE-754 Std. single precision floating-point format

Notes:
1. The function requires read/write access to the input vector x[N], because 
   it performs some data manipulations on it. Before returning, it always 
   restores original data.
2. Both top-8 search functions reserve a special value to maintain internal
   invariants during the search process. For vmax8 this special value is -32768,
   and for vmax8f it is -HUGE_VALF. It is still legal for the input vector to
   contain special values among other data. If this is the case, then the
   following peculiarities should be considered:
     A) Input vector elements equal to the special value are ignored.
     B) If the total number K of input vector elements distinct from the special
        value is less than 8 (i.e. K<8), then the last 8-K entries of output vector
        idx[8] are assigned zero, and the last 8-K elements of output vector m[8]
        are assigned the special value.
3. All elements of the input vector x[N] holding the same top value will be
   reported in search results, unless the total limit of 8 peak elements is
   exhausted.
4. For floating-point functions NAN values in vectors are ignored. If, however,
   the input vector contains NANs only, then the entries of output vector idx[8] are
   assigned zero, and the elements of output vector m[8] are assigned the special value.

Example (vmax8):
     Input:   N:    16
              x[N]: 1,-32768,1,2,3,-32768,-32768,-32768,-32768,-32768,
                    -32768,-32768,-32768,-32768,-32768,-32768
     Note that the number of input value distinct from -32768 is less than 8!
     Also note the duplicated value at positions 0 and 2.
     Results: m[8]: 3,2,1,1,-32768,-32768,-32768,-32768
            idx[8]: 4,3,0,2,0,0,0,0

Parameters:
Input:
x[N]     Input data vector
N        Length of input data vector
Output:
m[8]     8 maximum values in descending order
idx[8]   Indices of 8 maximum values in the input vector

Restrictions:
x        Aligned on 32-byte boundary
x,m,idx  Must not overlap
N        Multiple of 16 (vmax8) or 8 (vmax8f)
-------------------------------------------------------------------------*/

#ifndef BBE_SETDUALMAX

DISCARD_FUN(void, vmax8, (int16_t * restrict m,
    int16_t * restrict idx,
    int16_t * restrict x, int N))
#else

void vmax8 ( int16_t * restrict m, 
             int16_t * restrict idx, 
             int16_t * restrict x, int N )
{
    static const int16_t ALIGN(16) FlagsV[BBE_SIMD_WIDTH] = { 0x0003, 0x000c, 0x0030, 0x00c0, 0x0300, 0x0c00, 0x3000, (int16_t)0xc000, 0x0000, 0x000c, 0x0030, 0x00c0, 0x0300, 0x0c00, 0x3000, (int16_t)0xc000 };
    const xb_vecNx16 * restrict pFlagsV = (const xb_vecNx16 *)FlagsV;
    xb_vecNx16 flagsv, tmp;

    int n, k;
    int16_t ALIGN(16) PeakVal[BBE_SIMD_WIDTH];
    int16_t ALIGN(16) PeakInd[BBE_SIMD_WIDTH];
    const xb_vecNx16 * restrict pX;
    xb_vecNx16 x0, peakv, idxv, zm, id, z1;
    vboolN bLane01, bLane23, bLane45, bLane67, _b1, i0;
    vboolN bLane89, bLaneAB, bLaneCD, bLaneEF;
    xb_int16    a, b, q0, q1;

    if (N <= 0)
    {
        for (n = 0; n<8; n++)
        {
            if (m) m[n] = MIN_INT16;
            if (idx) idx[n] = 0;
        }
        return;
    }
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(N>0 && N % BBE_SIMD_WIDTH == 0);

    // init bool flags to get register BBE_MAXIDX
    flagsv = BBE_LVNX16_I(pFlagsV, 0);
    bLane01 = BBE_EXTRBN(flagsv, 0);
    bLane23 = BBE_EXTRBN(flagsv, 1);
    bLane45 = BBE_EXTRBN(flagsv, 2);
    bLane67 = BBE_EXTRBN(flagsv, 3);
    tmp = BBE_SELNX16I(flagsv, flagsv, BBE_SELI_ROTATE_RIGHT_4);
    bLane89 = BBE_EXTRBN(tmp, 0);
    bLaneAB = BBE_EXTRBN(tmp, 1);
    bLaneCD = BBE_EXTRBN(tmp, 2);
    bLaneEF = BBE_EXTRBN(tmp, 3);

    _b1 = BBE_EXTRBN(flagsv, 8);
    _b1 = BBE_NOTB(_b1);              // set "true" to all lanes

    //zm = BBE_MOVVINT16(N - 1);            //  replicate (N-1)
    zm = N - 1;                           //  replicate (N-1)
    z1 = BBE_MOVVINX16(BBE_MOVVI_Q15_M1); // replicate MININT

    for (k = 0; k < 8; k += 2)
    {
        // init registers BBE_MAX and BBE_MAX2 by INTMIN value
        // init registers BBE_MAXIDX and BBE_MAXIDX2 by zero
        BBE_SETDUALMAX(0x80000000);

        pX = (const xb_vecNx16 *)x;
        for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
        {
            BBE_LVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
            BBE_DUALMAXNX16(x0); // update registers: BBE_MAX, BBE_MAX2, BBE_MAXIDX and BBE_MAXIDX2
        }

        // get registers BBE_MAX and BBE_MAX_IDX
        PeakVal[0] = (int16_t)RUR_BBE_MAX_0();
        PeakVal[1] = (int16_t)RUR_BBE_MAX_1();
        PeakVal[2] = (int16_t)RUR_BBE_MAX_2();
        PeakVal[3] = (int16_t)RUR_BBE_MAX_3();
        PeakVal[4] = (int16_t)RUR_BBE_MAX_4();
        PeakVal[5] = (int16_t)RUR_BBE_MAX_5();
        PeakVal[6] = (int16_t)RUR_BBE_MAX_6();
        PeakVal[7] = (int16_t)RUR_BBE_MAX_7();
        PeakInd[0] = (int16_t)BBE_SELMAXIDX(bLane01, 0);
        PeakInd[1] = (int16_t)BBE_SELMAXIDX(bLane23, 0);
        PeakInd[2] = (int16_t)BBE_SELMAXIDX(bLane45, 0);
        PeakInd[3] = (int16_t)BBE_SELMAXIDX(bLane67, 0);
        PeakInd[4] = (int16_t)BBE_SELMAXIDX(bLane89, 0);
        PeakInd[5] = (int16_t)BBE_SELMAXIDX(bLaneAB, 0);
        PeakInd[6] = (int16_t)BBE_SELMAXIDX(bLaneCD, 0);
        PeakInd[7] = (int16_t)BBE_SELMAXIDX(bLaneEF, 0);

        // move state register BBE_MAX2 to BBE_MAX
        // move state register BBE_MAXIDX2 to BBE_MAXIDX
        BBE_MOVDUALMAXT(_b1);

        // get registers BBE_MAX2 and BBE_MAX_IDX2
        PeakVal[ 8] = (int16_t)RUR_BBE_MAX_0();
        PeakVal[ 9] = (int16_t)RUR_BBE_MAX_1();
        PeakVal[10] = (int16_t)RUR_BBE_MAX_2();
        PeakVal[11] = (int16_t)RUR_BBE_MAX_3();
        PeakVal[12] = (int16_t)RUR_BBE_MAX_4();
        PeakVal[13] = (int16_t)RUR_BBE_MAX_5();
        PeakVal[14] = (int16_t)RUR_BBE_MAX_6();
        PeakVal[15] = (int16_t)RUR_BBE_MAX_7();
        PeakInd[ 8] = (int16_t)BBE_SELMAXIDX(bLane01, 0);
        PeakInd[ 9] = (int16_t)BBE_SELMAXIDX(bLane23, 0);
        PeakInd[10] = (int16_t)BBE_SELMAXIDX(bLane45, 0);
        PeakInd[11] = (int16_t)BBE_SELMAXIDX(bLane67, 0);
        PeakInd[12] = (int16_t)BBE_SELMAXIDX(bLane89, 0);
        PeakInd[13] = (int16_t)BBE_SELMAXIDX(bLaneAB, 0);
        PeakInd[14] = (int16_t)BBE_SELMAXIDX(bLaneCD, 0);
        PeakInd[15] = (int16_t)BBE_SELMAXIDX(bLaneEF, 0);

        peakv = BBE_LVNX16_I((xb_vecNx16 *)PeakVal, 0);
        idxv = BBE_LVNX16_I((xb_vecNx16 *)PeakInd, 0);

        a = BBE_RMAXNX16(peakv); //first maximum value
        x0 = BBE_MOVNX16_FROM16(a);
        x0 = BBE_REPNX16(x0, 0);
        i0 = BBE_EQNX16(x0, peakv);
        id = BBE_MOVNX16T(idxv, zm, i0);
        q0 = BBE_RMINNX16(id); //first index

        x0 = BBE_MOVNX16_FROM16(q0);
        x0 = BBE_REPNX16(x0, 0);
        i0 = BBE_EQNX16(x0, idxv);

        peakv = BBE_MOVNX16T(z1, peakv, i0);
        b = BBE_RMAXNX16(peakv); //second maximum value
        x0 = BBE_MOVNX16_FROM16(b);
        x0 = BBE_REPNX16(x0, 0);
        i0 = BBE_EQNX16(x0, peakv);
        id = BBE_MOVNX16T(idxv, zm, i0);
        q1 = BBE_RMINNX16(id); //second index

        if (m)
        {
            m[k] = (a);
            m[k + 1] = (b);
        }
        if (idx)
        {
            idx[k] = (q0);
            idx[k + 1] = (q1);
        }

        NASSERT((int)q0 >= 0 && (int)q0 < N);
        NASSERT((int)q1 >= 0 && (int)q1 < N);

        x[(int)q0] = MIN_INT16;
        x[(int)q1] = MIN_INT16;
    }

    // restore source array
    for (k = 7; k >= 0; k--)
    {
        x[idx[k]] = m[k];
    }
} /* vmax8() */

#endif
