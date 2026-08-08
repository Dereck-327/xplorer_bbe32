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
    Decimating Block Complex FIR Filter
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
#include "firdecf_common.h"

/*-------------------------------------------------------------------------
Decimating Block Complex FIR Filter

Computes a complex FIR filter (direct-form) with decimation using real IR 
stored in vector h. The complex data input is stored in vector x. The filter
output result is stored in vector y. The filter calculates N output samples
using M coefficients and requires last D*N+M-1 samples in the delay line.

NOTE:
To avoid aliasing, the IR should be synthesized in such a way that filter pass
band is limited by input sample frequency divided by 2*D.

Representation:
firdec   16-bit signed fixed-point format
         Filter coefficients are Q15
         Number of fractional bits for input/output samples is user-difined
firdecf  IEEE-754 Std. single precision floating-point format for filter 
         coefficients and input/output samples

Parameters:
Input:
D        Decimation factor
N        Length of output sample block
M        Length of filter
h[M]     Filter coefficients; h[0] is to be multiplied by the newest 
         sample
x[N*D]   Input complex samples
Output:
y[N]     Output complex samples

Restrictions:
x,y      Must not overlap
x,y      Aligned on 32-byte boundary
N        Multiple of 8 (firdec) or 4 (firdecf)
M        2,4,8 or a positive multiple of 16 for D=2,3,4; or 
         a positive multiple of 16 for D>4
D>1

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for
filter lengths M=2,4,8,16 and 32 and decimation factors D=2,3 and 4, in
any combination.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, firdec[f]_init returns NULL handle.
-------------------------------------------------------------------------*/

/* Instance pointer validation number. */
#define FIRDECF_MAGIC 0xfcb66012

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

/* Decimator instance structure. */
typedef struct
{
    uint32_t          magic;   /* Instance pointer validation number */
    int               M;       /* Filter length                      */
    int               D;       /* Decimation factor                  */
    const float32_t * h;       /* Filter coefficients                */
          float32_t * d;       /* Delay line for samples             */
    proc_fxdxn_t      procFxn; /* Filter processing function         */
}
firdecf_t;

static const tFirFxdxns * getLayout(tFilterLayout * pFltr, int M, int D)
{
    const tFirFxdxns* pFxdxn;
    // Select the filter processing function, delay line length and coefficients number.
    pFxdxn = ((D == 2) ? (
        (M == 2)  ? &firdecf_2d_2m :
        (M == 4)  ? &firdecf_2d_4m :
        (M == 8)  ? &firdecf_2d_8m :
        (M == 16) ? &firdecf_2d_16m:
                    &firdecf_2d_xm):
              (D == 3) ? (
        (M == 2)  ? &firdecf_3d_2m :
        (M == 4)  ? &firdecf_3d_4m :
        (M == 8)  ? &firdecf_3d_8m :
        (M == 16) ? &firdecf_3d_16m:
                    &firdecf_3d_xm):
              (D == 4) ? (
        (M == 2)  ? &firdecf_4d_2m :
        (M == 4)  ? &firdecf_4d_4m :
        (M == 8)  ? &firdecf_4d_8m :
        (M == 16) ? &firdecf_4d_16m:
                    &firdecf_4d_xm):
        &firdecf_xd_xm
        );
    pFxdxn->initAllocFxdxn->fnalloc(pFltr, M, D);
    return pFxdxn;
};

/* Calculate the memory block size for a decimator with given attributes. */
size_t firdecf_alloc(int D, int M)
{
    tFilterLayout fltr;
    NASSERT(((M == 2 || M == 4 || M == 8 || M == 16) && (D == 2 || D == 3 || D == 4)) || (D >= 2 && M % 16 == 0));

    getLayout(&fltr, M, D);
    return (ALIGNED_SIZE(sizeof(firdecf_t), 4) +
        ALIGNED_SIZE(fltr.delLength * 2 * sizeof(float32_t), 2 * BBE_SIMD_WIDTH) +
        ALIGNED_SIZE(fltr.coefNum       * sizeof(float32_t), 2 * BBE_SIMD_WIDTH));
} // firdecf_alloc()

/* Initialize the decimator structure. The delay line is zeroed. */
firdecf_handle_t firdecf_init(void * objmem, int D, int M, const float32_t * restrict h)
{
    const tFirFxdxns * pFxdxn;
    tFilterLayout fltr;
    firdecf_t * firdec;
    void * ptr;
    float32_t * restrict pd;
    float32_t * restrict ph;
    int m;

    NASSERT(objmem && h);
    NASSERT(((M == 2 || M == 4 || M == 8 || M == 16) && (D == 2 || D == 3 || D == 4)) || (D >= 2 && M % 16 == 0));

    pFxdxn = getLayout(&fltr, M, D);

    /* Partition the memory block */
    ptr = objmem;
    firdec = (firdecf_t*)ALIGNED_ADDR(ptr, 4);
    ptr = firdec + 1;
    ph = (float32_t*)ALIGNED_ADDR(ptr, 32);
    ptr = ph + fltr.coefNum;
    pd = (float32_t*)ALIGNED_ADDR(ptr, 32);
    ptr = pd + 2 * fltr.delLength;
    ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)firdecf_alloc(D, M));

    /* Copy filter coefficients in reverted order and zero the delay line. */
    pFxdxn->initAllocFxdxn->fninit(ph, h, M, D);
    for (m = 0; m < 2 * fltr.delLength; m++) pd[m] = 0.f;

    firdec->magic = FIRDECF_MAGIC;
    firdec->M = fltr.coefNum;
    firdec->D = D;
    firdec->h = ph;
    firdec->d = pd;
    firdec->procFxn = pFxdxn->procFxdxn;

    return NatureDSP_Baseband_isPresent(firdec->procFxn) ? (firdec) : NULL;
} // firdecf_init()

/* process block of samples */
void firdecf_process(firdecf_handle_t handle, complex_float * restrict y, const complex_float * restrict x, int N)
{
    int M, D;
    firdecf_t * firdec = (firdecf_t*)handle;
    NASSERT(firdec && firdec->magic == FIRDECF_MAGIC && y && x);
    M = firdec->M;
    D = firdec->D;

    if (N <= 0) return;

    NASSERT(!(N & 3));
    //NASSERT(((M == 2 || M == 4 || M == 8 || M == 16) && (D == 2 || D == 3 || D == 4)) || (D >= 2 && M % 16 == 0));
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(firdec->d);
    NASSERT_ALIGN32(firdec->h);
    NASSERT(firdec->procFxn);

    firdec->procFxn(y, x, firdec->h, firdec->d, M, N, D);
} // firdecf_process()
