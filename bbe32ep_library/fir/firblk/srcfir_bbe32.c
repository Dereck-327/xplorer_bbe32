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
    NatureDSP_Baseband library. FIR filters and Related Functions
    Block Complex FIR filter w/ Real Symmetric Impulse Response
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fir.h"
#include "NatureDSP_Baseband_id.h"
#include "srcfir_common.h"

/*-------------------------------------------------------------------------
Block Complex FIR filter w/ Real Symmetric Impulse Response

Passes complex input data through a direct-form FIR filter with real 
coefficients and symmetric impulse response. The filter calculates N output 
samples using (M+1)/2 coefficients and requires last M+N-1 samples in the delay 
line.

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in complex samples) depends on FIR order M and is defined by srcfir_algDelay(M)

Precision: 16-bit data, 16-bit coefficients, 16-bit outputs, everything in Q15

Parameters:
Input:
objmem      Allocated memory block
h[(M+1)/2]  Filter coefficients (half the number of filter taps); h[0] is to
            be multiplied by the newest sample,  Q15
N           Length of sample block
M           Length of filter
x[N]        Input samples, Q15
Output:
y[N]        Output samples, Q15

Restrictions:
x,y         Must not overlap
x,y         Aligned on 32-byte boundary
M           17,33 or a multiple of 16
N           Multiple of 8
M>0

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=16,17,32 and 33.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, srcfir_init returns NULL handle.
---------------------------------------------------------------------------*/

/* Instance pointer validation number. */
#define SRCFIR_MAGIC 0x2e361b4f

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

/* Object allocation */
size_t srcfir_alloc ( int M )
{
    int delLen;

    NASSERT(M>0);
    NASSERT(M == 17 || M == 33 || (M >= 16 && !(M & 15)));

    delLen = 2 * 2 * (M + BBE_SIMD_WIDTH / 2);

    if (M>32 && M % 16 == 0)
    {
        return (ALIGNED_SIZE(sizeof(srcfir_t), 4)
            + // Delay line
            ALIGNED_SIZE(2 * delLen * sizeof(int16_t), (2 * BBE_SIMD_WIDTH))
            + // Coefficients
            ALIGNED_SIZE((M + 1) / 2 * sizeof(int16_t), (2 * BBE_SIMD_WIDTH)));
    }
    // for M==16,17,32,33
    return ALIGNED_SIZE(sizeof(srcfir_t), 4) +
        ALIGNED_SIZE(      2 * M * sizeof(int16_t), (2 * BBE_SIMD_WIDTH)) +
        ALIGNED_SIZE((M + 1) / 2 * sizeof(int16_t), (2 * BBE_SIMD_WIDTH));
} /* srcfir_alloc() */

/* Object initialization */
srcfir_handle_t srcfir_init ( void * mem, int M, const int16_t * h )
{
    srcfir_ptr_t srcfir;
    void *       ptr;
    int16_t *    delLine;
    int          delLen;
    int16_t *    coef;
    int m;

    NASSERT(mem && h);
    NASSERT(M == 17 || M == 33 || (M >= 16 && !(M & 15)));

    if (M>32 && !(M & 15))
    {
        //
        // Introduce a delay for some combinations of M/N.
        //

        delLen = 2 * 2 * (M + BBE_SIMD_WIDTH / 2);

        //
        // Partition the memory block.
        //

        ptr = mem;
        srcfir = (srcfir_ptr_t)ALIGNED_ADDR(ptr, 4);
        ptr = srcfir + 1;
        delLine = (int16_t*)ALIGNED_ADDR(ptr, 32);
        ptr = delLine + 2 * delLen;
        coef = (int16_t*)ALIGNED_ADDR(ptr, 32);
        ptr = coef + (M + 1) / 2;

        NASSERT((int8_t*)ptr - (int8_t*)mem <= (int)srcfir_alloc(M));

        //
        // Copy the filter coefficients, zero the delay line.
        //
        for (m = 0; m<(M >> 1); m++)
        {
            coef[m] = h[m];
        }
        if (M & 1)
        {
            coef[m + 0] = h[m]; coef[m + 1] = 0;
        }

        for (m = 0; m < 2 * delLen; m++)
        {
            delLine[m] = 0;
        }     

        //
        // Initialize the filter instance.
        //
        srcfir->magic = SRCFIR_MAGIC;
        srcfir->M = M;
        srcfir->coef = coef;
        srcfir->delayLine = delLine;
        srcfir->p0 = delLine;
        srcfir->p1 = delLine + delLen;

        srcfir->process = srcfir_process_16m;
    }
    else
    {
        //----------------------------------------------------
        // for M=16,17,32,33
        //----------------------------------------------------
        // Partition the memory block.

        ptr = mem;
        srcfir = (srcfir_ptr_t)ALIGNED_ADDR(ptr, 4);
        ptr = srcfir + 1;
        coef = (int16_t*)ALIGNED_ADDR(ptr, 32);
        ptr = coef + (M + 1) / 2;
        delLine = (int16_t*)ALIGNED_ADDR(ptr, 32);
        ptr = delLine + 2 * M;

        NASSERT((int8_t*)ptr - (int8_t*)mem <= (int)srcfir_alloc(M));

        srcfir->magic = SRCFIR_MAGIC;
        srcfir->M = M;
        srcfir->coef = coef;
        srcfir->delayLine = delLine;

        // Copy the filter coefficients, zero the delay line.
        {
            xb_vecNx16 t;
            const xb_vecNx16 *ph = (const xb_vecNx16 *)h;
            xb_vecNx16 *d = (xb_vecNx16 *)srcfir->delayLine;
            valign vh;
            int i;

            vh = BBE_LA_PP(ph);
            BBE_LAVNX16_XP(t, vh, ph, ((M + 1) >> 1) * sizeof(int16_t));
            BBE_SVNX16_I(t, (xb_vecNx16 *)srcfir->coef, 0);
            BBE_LAVNX16_XP(t, vh, ph, ((M + 1) >> 1) * sizeof(int16_t));
            BBE_SVNX16_I(t, (xb_vecNx16 *)srcfir->coef, 2 * BBE_SIMD_WIDTH);
            if (M == 17)
            {
                ((int16_t *)srcfir->coef)[8] = h[8];
            }
            else if (M == 33)
            {
                ((int16_t *)srcfir->coef)[16] = h[16];
            }
            t = 0;
            for (i = 0; i<2 * (M - 1 + BBE_SIMD_WIDTH / 2 - 1) / BBE_SIMD_WIDTH; i++)
                BBE_SVNX16_IP(t, d, 2 * BBE_SIMD_WIDTH);
        }

        // set process handler
        switch (M)
        {
        case 16:    srcfir->process = srcfir_process_16; break;
        case 17:    srcfir->process = srcfir_process_17; break;
        case 32:    srcfir->process = srcfir_process_32; break;
        case 33:    srcfir->process = srcfir_process_33; break;
        }
    }

    return NatureDSP_Baseband_isPresent(srcfir->process) ? (srcfir) : NULL;
} /* srcfir_init() */

/* Update the delay line and compute filter output */
void srcfir_process ( srcfir_handle_t handle, 
                      complex_fract16 * restrict y, const complex_fract16 * restrict x , int N)
{
    srcfir_t *srcfir = (srcfir_t *)handle;
    if (N <= 0) return;
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(srcfir && srcfir->magic == SRCFIR_MAGIC && y && x);
    srcfir->process(srcfir, srcfir->delayLine, y, x, srcfir->coef, srcfir->M, N);
} /* srcfir_process() */

/* Return the algorithmic delay, in samples. */
int srcfir_algDelay( int M )
{
  return (0);

} /* srcfir_algDelay() */
