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
    Block Real FIR Filter
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
#include "bkfirf_common.h"

/*-------------------------------------------------------------------------
Block Real FIR Filter

Computes a real FIR filter (direct-form) using IR stored in vector h. The
real data input is stored in vector x. The filter output result is stored
in vector y. The filter calculates N output samples using M coefficients
and requires last M+N-1 samples in the delay line.

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in samples) depends on FIR order M and is defined by bkfir[f]_algDelay(M).

Representation:
bkfir   16-bit signed fixed-point format
        Filter coefficients are Q15
        Number of fractional bits for input/output samples is user-difined
bkfirf  IEEE-754 Std. single precision floating-point format for filter 
        coefficients and input/output samples

Parameters:
Input:
objmem  Allocated memory block
h[M]    Filter coefficients; h[0] is to be multiplied by the newest sample
N       Length of sample block
M       Length of filter
x[N]    Input samples
Output:
y[N]    Output samples

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 16 (bkfir) or 8 (bkfirf)
M       2,4,8 or a positive multiple of 16

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=2,4,8 and 16.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, bkfir[f]_init returns NULL handle.
-------------------------------------------------------------------------*/

/* Instance pointer validation number. */
#define BKFIRF_MAGIC 0x6543face

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
bkfirf_t;

/* Calculate the memory block size for an FIR filter with given attributes. */
size_t bkfirf_alloc(int M)
{
    int M1 = (M + 7)&~7;

    NASSERT(((M >= 16) && !(M & 7)) || M == 2 || M == 4 || M == 8);

    return (ALIGNED_SIZE(sizeof(bkfirf_t), 4)
            + // Delay line
            ALIGNED_SIZE(M1 * sizeof(float32_t) + 4 * BBE_SIMD_WIDTH, 32)
            + // Filter coefficients
            ALIGNED_SIZE(M1 * sizeof(float32_t), 32));
} // bkfirf_alloc()

/* Initialize the filter structure. The delay line is zeroed. */
bkfirf_handle_t bkfirf_init(void * objmem, int M, const float32_t * h)
{
    bkfirf_t* bkfir;
    void * ptr;
    int M1 = (M + 7)&~7;
    float32_t * restrict pd;
    float32_t * restrict ph;
    int m;

    NASSERT(objmem && h);
    NASSERT(((M >= 16) && !(M & 7)) || M == 2 || M == 4 || M == 8);

    /* Partition the memory block */
    ptr = objmem;
    bkfir = (bkfirf_t*)ALIGNED_ADDR(ptr, 4);
    ptr = bkfir + 1;
    ph = (float32_t*)ALIGNED_ADDR(ptr, 32);
    ptr = ph + M;
    pd = (float32_t*)ALIGNED_ADDR(ptr, 32);
    ptr = pd + M + 4 * BBE_SIMD_WIDTH / sizeof(float32_t);
    ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)bkfirf_alloc(M));

    /* Copy filter coefficients in reverted order and zero the delay line. */
    for (m = 0; m < M; m++) ph[m] = h[M-1-m];
    for (; m < M1; m++) ph[m] = 0.f;  /* pad remaining coefficients with zeroes */
    for (m = 0; m < M1 + BBE_SIMD_WIDTH; m++) pd[m] = 0.f;

    bkfir->magic = BKFIRF_MAGIC;
    bkfir->M = M;
    bkfir->h = ph;
    bkfir->d = pd;
    bkfir->procFxn = ((M == 2) ? bkfirf_process_2 :
                      (M == 4) ? bkfirf_process_4 :
                      (M == 8) ? bkfirf_process_8 :
                      (M == 16)? bkfirf_process_16:
                                 bkfirf_process_8m);
    
    return NatureDSP_Baseband_isPresent(bkfir->procFxn) ? (bkfir) : NULL;
} // bkfirf_init()

/* process block of samples */
void bkfirf_process(bkfirf_handle_t _bkfir, float32_t * restrict  y, const float32_t * restrict  x, int N)
{
    int M;
    bkfirf_t * bkfir = (bkfirf_t *)_bkfir;
    NASSERT(bkfir && bkfir->magic == BKFIRF_MAGIC && y && x);

    if (N <= 0) return;
    M = bkfir->M;

    NASSERT(!(N & 7));
    NASSERT(((M >= 16) && !(M & 7)) || M == 2 || M == 4 || M == 8);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT(bkfir->procFxn);

    bkfir->procFxn(bkfir->d, y, x, bkfir->h, M, N);
} // bkfirf_process()

/* Return the algorithmic delay, in samples. */
int bkfirf_algDelay(int M) 
{ 
    return 0;
} // bkfirf_algDelay()
