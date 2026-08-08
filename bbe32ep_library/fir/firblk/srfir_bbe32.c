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
    Block Real FIR filter w/ Symmetric Impulse Response
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
#include "srfir_common.h"

/*-------------------------------------------------------------------------
Block Real FIR filter w/ Symmetric Impulse Response

Passes real input data through a direct-form FIR filter with real coefficients
and symmetric impulse response. The filter calculates N output samples using
(M+1)/2 coefficients and requires last M+N-1 samples in the delay line.

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in samples) depends on FIR order M and is defined by srfir_algDelay(M)

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
M           2,4,8 or a multiple of 16
N           Multiple of 16  
M>0

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=2,4,8,16 and 32.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, srfir_init returns NULL handle.
---------------------------------------------------------------------------*/

/* Instance pointer validation number. */
#define SRFIR_MAGIC 0xbacf0d84

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

/* Object allocation */
size_t srfir_alloc(int M)
{
    int delLen;

    NASSERT(M > 0 && (!(M & 15) || (M == 2) || (M == 4) || (M == 8) || (M == 16)));

    delLen = 2 * 2 * (M + BBE_SIMD_WIDTH);
    return (ALIGNED_SIZE(sizeof(srfir_t), 4)
        + // Delay line
        ALIGNED_SIZE(2 * delLen * sizeof(int16_t), 32)
        + // Coefficients
        ALIGNED_SIZE((M + 1) / 2 * sizeof(int16_t), 32));
} /* srfir_alloc() */

/* Object initialization */
srfir_handle_t srfir_init ( void * mem, int M, const int16_t * h )
{
    srfir_ptr_t srfir;
    void *      ptr;
    int16_t *   delLine;
    int         delLen;
    int16_t *   coef;
    int m;

    NASSERT(mem && h);
    NASSERT(M > 0 && (!(M & 15) || (M == 2) || (M == 4) || (M == 8) || (M == 16)));

    //
    // Introduce a delay for some combinations of M
    //

    delLen = 2 * 2 * (M + BBE_SIMD_WIDTH);

    //
    // Partition the memory block
    //

    ptr = mem;
    srfir = (srfir_ptr_t)ALIGNED_ADDR(ptr, 4);
    ptr = srfir + 1;
    delLine = (int16_t*)ALIGNED_ADDR(ptr, 32);
    ptr = delLine + delLen;
    coef = (int16_t*)ALIGNED_ADDR(ptr, 32);
    ptr = coef + (M + 1) / 2;

    NASSERT((int8_t*)ptr - (int8_t*)mem <= (int)srfir_alloc(M));

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

    for (m = 0; m < delLen; m++)
    {
        delLine[m] = 0;
    }

    //
    // Initialize the filter instance.
    //

    srfir->magic = SRFIR_MAGIC;
    srfir->M = M;
    srfir->coef = coef;
    srfir->delayLine = delLine;
    srfir->p0 = delLine;
    srfir->p1 = delLine + 2 * (M + BBE_SIMD_WIDTH);

    srfir->procFxn = ((M == 2)  ? srfir_process_2  :
                      (M == 4)  ? srfir_process_4  :
                      (M == 8)  ? srfir_process_8  :
                      (M == 16) ? srfir_process_16 : 
                      (M == 32) ? srfir_process_32 : srfir_process_16m);

    return NatureDSP_Baseband_isPresent(srfir->procFxn) ? (srfir) : NULL;
} /* srfir_init() */

/* Update the delay line and compute filter output */
void srfir_process ( srfir_handle_t  handle, 
                    int16_t * restrict y, const int16_t * restrict x , int N)
{
    srfir_ptr_t srfir = (srfir_ptr_t)handle;
    if (N <= 0) return;
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT(srfir != NULL && srfir->magic == SRFIR_MAGIC);
    srfir->procFxn(srfir, srfir->delayLine, y, x, srfir->coef, srfir->M, N);
} /* srfir_process() */

/* Return the algorithmic delay, in samples. */
int srfir_algDelay ( int M )
{
  return (0);

} /* srfir_algDelay() */
