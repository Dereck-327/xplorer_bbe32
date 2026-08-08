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
    Block Complex FIR Filter with Real Coefficients
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
#include "rcfirf_common.h"

/*-------------------------------------------------------------------------
Block Complex FIR Filter with Real Coefficients

Computes a complex FIR filter (direct-form) using real IR stored in vector h.
The complex data input is stored in vector x. The filter output result is
stored in vector y. The filter calculates N output samples using M coefficients
and requires last M+N-1 samples in the delay line.

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in complex samples) depends on FIR order M and defined by rcfir[f]_algDelay(M)

Representation:
rcfir   16-bit signed fixed-point format
        Filter coefficients are Q15
        Number of fractional bits for input/output samples is user-difined
rcfirf  IEEE-754 Std. single precision floating-point format for filter 
        coefficients and input/output samples

Parameters:
Input:
objmem  Allocated memory block
h[M]    Filter coefficients; h[0] is to be multiplied by the newest sample
N       Length of sample block
M       Length of filter
x[N]    Complex input samples
Output:
y[N]    Complex output samples

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 8 (rcfir) or 4 (rcfirf)
M       2,4,8 or a positive multiple of 16

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=2,4,8 and 16.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, rcfir[f]_init returnd NULL handle.
-------------------------------------------------------------------------*/

/* Instance pointer validation number. */
#define RCFIRF_MAGIC 0xdb5afd8d

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

/* Filter instance structure. */
typedef struct
{
    uint32_t          magic;   /* Instance pointer validation number */
    int               M;       /* Filter length                      */
    const float32_t * h;       /* Filter coefficients                */
          float32_t * d;       /* Delay line for samples             */
    proc_fxn_t        procFxn; /* Filter processing function         */
}
rcfirf_t;

/* Object allocation */
size_t rcfirf_alloc(int M)
{
    int M1 = (M + 3)&~3;

    NASSERT(M == 2 || M == 4 || M == 8 || M == 16 || (M>16 && !(M & 3)));

    return (ALIGNED_SIZE(sizeof(rcfirf_t), 4)
            + // Delay line
            ALIGNED_SIZE(M1 * sizeof(complex_float) + 4 * BBE_SIMD_WIDTH, 32)
            + // Filter coefficients
            ALIGNED_SIZE(M1 * sizeof(float32_t), 32));
} // rcfirf_alloc()

/* Object initialization */
rcfirf_handle_t rcfirf_init(void * objmem, int M, const float32_t * restrict h)
{
    rcfirf_t* rcfir;
    void * ptr;
    int M1 = (M + 3)&~3;
    float32_t * restrict pd;
    float32_t * restrict ph;
    int m;

    NASSERT(objmem && h);
    NASSERT(M == 2 || M == 4 || M == 8 || M == 16 || (M>16 && !(M & 3)));

    /* Partition the memory block */
    ptr = objmem;
    rcfir = (rcfirf_t*)ALIGNED_ADDR(ptr, 4);
    ptr = rcfir + 1;
    ph = (float32_t*)ALIGNED_ADDR(ptr, 32);
    ptr = ph + M;
    pd = (float32_t*)ALIGNED_ADDR(ptr, 32);
    ptr = pd + 2 * M + 4 * BBE_SIMD_WIDTH / sizeof(float32_t);
    ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)rcfirf_alloc(M));

    /* Copy filter coefficients in reverted order and zero the delay line. */
    for (m = 0; m < M; m++) ph[m] = h[M - 1 - m];
    for (; m < M1; m++) ph[m] = 0.f;  /* pad remaining coefficients with zeroes */
    for (m = 0; m < 2 * M1 + BBE_SIMD_WIDTH; m++) pd[m] = 0.f;

    rcfir->magic = RCFIRF_MAGIC;
    rcfir->M = M;
    rcfir->h = ph;
    rcfir->d = pd;
    rcfir->procFxn = ((M == 2) ? rcfirf_process_2 :
                      (M == 4) ? rcfirf_process_4 :
                      (M == 8) ? rcfirf_process_8 :
                      (M == 16)? rcfirf_process_16:
                                 rcfirf_process_8m);

    return NatureDSP_Baseband_isPresent(rcfir->procFxn) ? (rcfir) : NULL;
} // rcfirf_init()

/* Update the delay line and compute filter output */
void rcfirf_process(rcfirf_handle_t _rcfir, complex_float * restrict y, const complex_float * restrict x, int N)
{
    int M;
    rcfirf_t * rcfir = (rcfirf_t *)_rcfir;
    NASSERT(rcfir && rcfir->magic == RCFIRF_MAGIC && y && x);

    if (N <= 0) return;
    M = rcfir->M;

    NASSERT(!(N & 3));
    NASSERT(M == 2 || M == 4 || M == 8 || M == 16 || (M>16 && !(M & 3)));
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT(rcfir->procFxn);

    rcfir->procFxn(y, x, rcfir->h, rcfir->d, M, N);
} //rcfirf_process()

/* Return the algorithmic delay, in samples. */
int rcfirf_algDelay ( int M )
{
    return (0);
} // rcfirf_algDelay()
