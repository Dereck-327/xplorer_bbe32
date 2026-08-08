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

void vmax ( int16_t * restrict m, 
            int16_t * restrict idx, 
      const int16_t * restrict x, int N )
{
#ifdef BBE_SETDUALMAX

    static const int16_t ALIGN(16) FlagsV[BBE_SIMD_WIDTH] = { 0x0003, 0x000c, 0x0030, 0x00c0, 0x0300, 0x0c00, 0x3000, (int16_t)0xc000, 0x0000, 0x000c, 0x0030, 0x00c0, 0x0300, 0x0c00, 0x3000, (int16_t)0xc000 };
    const xb_vecNx16 * restrict pFlagsV = (const xb_vecNx16 *)FlagsV;
    xb_vecNx16 flagsv, tmp;

    int n;
    int16_t ALIGN(16) PeakVal[BBE_SIMD_WIDTH];
    int16_t ALIGN(16) PeakInd[BBE_SIMD_WIDTH];
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    xb_vecNx16 x0, peakv, idxv, zm, id, z1;
    vboolN bLane01, bLane23, bLane45, bLane67, _b1, i0;
    vboolN bLane89, bLaneAB, bLaneCD, bLaneEF;
    xb_int16    a, b, q0, q1;

    if (N <= 0)
    {
        if (m)
        {
            m[0] = (int16_t)0x8000;
            m[1] = (int16_t)0x8000;
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

    // init registers BBE_MAX and BBE_MAX2 by INTMIN value
    // init registers BBE_MAXIDX and BBE_MAXIDX2 by zero
    BBE_SETDUALMAX(0x80000000);

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
        m[0] = (a);
        m[1] = (b);
    }
    if (idx)
    {
        idx[0] = (q0);
        idx[1] = (q1);
    }

#elif 0 // OLD: This code is short (as Xtensa Xplorer sample) but is not bitexact because "dual max" instructions get MIN LANE, not MIN INDEX.
    //      For example, in test 22 (2 equivalent maximum values) instruction BBE_SELMAXIDX() returns incorrect index of maximum value. 

    int n;
    unsigned int max0, max1, idx0, idx1;
    const xb_vecNx16 * restrict pX = (const xb_vecNx16 *)x;
    xb_vecNx16 x0;
    vboolN bMax0, bMax1, bMax12;

    //TMP: TO DEBUG
    /*
    unsigned int MAX_0, MAX_1, MAX_2, MAX_3;
    unsigned int IDX_0, IDX_1;
    unsigned int MAX2_0, MAX2_1, MAX2_2, MAX2_3;
    unsigned int IDX2_0, IDX2_1;
    unsigned int BbeIdx;
    */

    if (N <= 0)
    {
        if (m)
        {
            m[0] = (int16_t)0x8000;
            m[1] = (int16_t)0x8000;
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

    // init registers BBE_MAX and BBE_MAX2 by INTMIN value
    // init registers BBE_MAXIDX and BBE_MAXIDX2 by zero
    BBE_SETDUALMAX(0x80000000);

    //TMP: TO DEBUG
    /*
    MAX_0 = RUR_BBE_MAX_0(); MAX_1 = RUR_BBE_MAX_1(); MAX_2 = RUR_BBE_MAX_2(); MAX_3 = RUR_BBE_MAX_3();
    IDX_0 = RUR_BBE_MAXIDX_0(); IDX_1 = RUR_BBE_MAXIDX_1();
    MAX2_0 = RUR_BBE_MAX2_0(); MAX2_1 = RUR_BBE_MAX2_1(); MAX2_2 = RUR_BBE_MAX2_2(); MAX2_3 = RUR_BBE_MAX2_3();
    IDX2_0 = RUR_BBE_MAXIDX2_0(); IDX2_1 = RUR_BBE_MAXIDX2_1();
    BbeIdx = RUR_BBE_IDX();
    */

    for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        BBE_LVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH);

        BBE_DUALMAXNX16(x0); // update registers: BBE_MAX, BBE_MAX2, BBE_MAXIDX and BBE_MAXIDX2
    }

    //TMP: TO DEBUG
    /*
    MAX_0 = RUR_BBE_MAX_0(); MAX_1 = RUR_BBE_MAX_1(); MAX_2 = RUR_BBE_MAX_2(); MAX_3 = RUR_BBE_MAX_3();
    IDX_0 = RUR_BBE_MAXIDX_0(); IDX_1 = RUR_BBE_MAXIDX_1();
    MAX2_0 = RUR_BBE_MAX2_0(); MAX2_1 = RUR_BBE_MAX2_1(); MAX2_2 = RUR_BBE_MAX2_2(); MAX2_3 = RUR_BBE_MAX2_3();
    IDX2_0 = RUR_BBE_MAXIDX2_0(); IDX2_1 = RUR_BBE_MAXIDX2_1();
    BbeIdx = RUR_BBE_IDX();
    */

    // get first max value from register BBE_MAX and flags in which lane is placed max value
    BBE_RBDUALMAXR(max0, bMax0);

    // get index of first max value from register  BBE_MAXIDX
    idx0 = BBE_SELMAXIDX(bMax0, 0 /* was 1 */);

    // get flags:  BBE_MAX2i > BBE_MAXi
    bMax12 = BBE_GTMAXNX16();

    bMax12 = BBE_ORB(bMax12, bMax0); // flags: "BBE_MAX2i > BBE_MAXi" exclude lane of max0 value (in this lane flag is 1)

    // merge state registers BBE_MAX and BBE_MAX2 into state register BBE_MAX
    // merge state registers BBE_MAXIDX and BBE_MAXIDX2 into state register BBE_MAXIDX
    BBE_MOVDUALMAXT(bMax12);

    //TMP: TO DEBUG
    /*
    MAX_0 = RUR_BBE_MAX_0(); MAX_1 = RUR_BBE_MAX_1(); MAX_2 = RUR_BBE_MAX_2(); MAX_3 = RUR_BBE_MAX_3();
    IDX_0 = RUR_BBE_MAXIDX_0(); IDX_1 = RUR_BBE_MAXIDX_1();
    MAX2_0 = RUR_BBE_MAX2_0(); MAX2_1 = RUR_BBE_MAX2_1(); MAX2_2 = RUR_BBE_MAX2_2(); MAX2_3 = RUR_BBE_MAX2_3();
    IDX2_0 = RUR_BBE_MAXIDX2_0(); IDX2_1 = RUR_BBE_MAXIDX2_1();
    BbeIdx = RUR_BBE_IDX();
    */

    // get second max value from register BBE_MAX and flags in which lane is placed max value
    BBE_RBDUALMAXR(max1, bMax1);

    // get index of second max value from register  BBE_MAXIDX
    idx1 = BBE_SELMAXIDX(bMax1, 0 /* was 1 */);

    if (m)
    {
        m[0] = (max0);
        m[1] = (max1);
    }
    if (idx)
    {
        idx[0] = ((int16_t)max0 == (int16_t)0x8000) ? (0) : (idx0);
        idx[1] = ((int16_t)max1 == (int16_t)0x8000) ? (0) : (idx1);
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
            m[0] = (int16_t)0x8000;
            m[1] = (int16_t)0x8000;
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
    z1 = BBE_MOVVINX16(BBE_MOVVI_Q15_M1); // replicate MININT
    zm = (N - 1);              //  replicate (N-1)
    id = BBE_SEQNX16();                   // sequence: 0,1,2, ... 7
    z2 = BBE_MOVVINT16(BBE_SIMD_WIDTH);   // replicate count of lanes
    m0 = m1 = z1;
    id0 = id1 = 0;
    idn0 = idn1 = 0;
    for (n = 0; n<(N >> LOG2_BBE_SIMD_WIDTH); n++)
    {
        BBE_LVNX16_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
        BBE_BMAXNX16(i0, mx, x0, m0);
        x1 = BBE_MOVNX16T(m0, x0, i0);

        BBE_MOVIDXNX16F(tmp, id0, idn1, i0);
        BBE_MOVIDXNX16T(id0, id0, idn0, i0);

        BBE_BMAXNX16(i1, m1, x1, m1);
        id1 = BBE_MOVNX16T(tmp, id1, i1);

        // n*BBE_SIMD_WIDTH+0,....n*BBE_SIMD_WIDTH+BBE_SIMD_WIDTH-1
        m0 = mx;

    }
    a = BBE_RMAXNX16(m0); //first maximum value
    x0 = BBE_MOVNX16_FROM16(a);
    x0 = BBE_REPNX16(x0, 0);
    i0 = BBE_EQNX16(x0, m0);
    id = BBE_MOVNX16T(id0, zm, i0);
    q0 = BBE_RMINNX16(id); //first index

    x0 = BBE_MOVNX16_FROM16(q0);
    x0 = BBE_REPNX16(x0, 0);
    i0 = BBE_EQNX16(x0, id0);

    tmp = BBE_MOVNX16T(z1, m0, i0);
    i_equ = BBE_EQNX16(tmp, m1);         // flags: "max values is equivalent"
    BBE_BMINNX16(i_ind, id01, id0, id1); // flags: "min index of max value"
    i_equ = BBE_ANDB(i_equ, i_ind);      // flags: "choise max value with min index from id0 (else from id1)"
    BBE_BMAXNX16(i2, m1, tmp, m1);
    i2 = BBE_ORB(i2, i_equ);             // use all flags
    b = BBE_RMAXNX16(m1); //second maximum value
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
} /* vmax() */
