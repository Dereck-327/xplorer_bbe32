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
    Interpolating Block Complex FIR Filter
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
#include "firinterpf_common.h"

/*-------------------------------------------------------------------------
Interpolating Block Complex FIR Filter

Computes a complex FIR filter (direct-form) with interpolation using real
IR stored in vector h. The complex data input is stored in vector x. The
filter output result is stored in vector y. The filter calculates N*D complex
output samples using M*D coefficients and requires last N+M-1 samples in the
delay line.

Representation:
firinterp   16-bit signed fixed-point format
            Filter coefficients are Q15
            Number of fractional bits for input/output samples is user-difined
firinterpf  IEEE-754 Std. single precision floating-point format for filter 
            coefficients and input/output samples

Parameters:
Input:
D           Interpolation ratio 
N           Length of input sample block
M           Length of subfilter. Total length of filter is M*D
h[M*D]      Filter coefficients; h[0] is to be multiplied by the newest 
            sample,Q15
x[N]        Input complex samples
Output:
y[N*D]      Output complex samples

Restrictions:
x,y         Must not overlap
x,y         Aligned on 32-byte boundary
N           Multiple of 8 (firinterp) or 4 (firinterpf)
M           2,4,8 or a positive multiple of 16 for D=2,3,4,6,12; or 
            a positive multiple of 8 for other D
D>1

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for
subfilter lengths M=2,4,8,16 and 32 and interpolation factors D=2,3 and 4,
in any combination.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, firinterp[f]_init returns NULL handle.
-------------------------------------------------------------------------*/

/* Instance pointer validation number. */
#define FIRINTERPF_MAGIC 0x514c3ea8

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
firinterpf_t;

static const tFirFxdxns * getLayout(tFilterLayout * pFltr, int M, int D)
{
    const tFirFxdxns* pFxdxn;
    // Select the filter processing function, delay line length and coefficients number.
    pFxdxn = ((D == 2) ? (
        (M == 2)  ? &interpf_2d_2m :
        (M == 4)  ? &interpf_2d_4m :
        (M == 8)  ? &interpf_2d_8m :
        (M == 16) ? &interpf_2d_16m:
                    &interpf_2d_xm):
              (D == 3) ? (
        (M == 2)  ? &interpf_3d_2m :
        (M == 4)  ? &interpf_3d_4m :
        (M == 8)  ? &interpf_3d_8m :
        (M == 16) ? &interpf_3d_16m:
                    &interpf_3d_xm):
              (D == 4) ? (
        (M == 2)  ? &interpf_4d_2m :
        (M == 4)  ? &interpf_4d_4m :
        (M == 8)  ? &interpf_4d_8m :
        (M == 16) ? &interpf_4d_16m:
                    &interpf_4d_xm):
              (D == 6) ? (
        (M == 2)  ? &interpf_6d_2m :
        (M == 4)  ? &interpf_6d_4m :
        (M == 8)  ? &interpf_6d_8m :
        (M == 16) ? &interpf_6d_16m:
                    &interpf_6d_xm):
              (D == 12) ? (
        (M == 2)  ? &interpf_12d_2m :
        (M == 4)  ? &interpf_12d_4m :
        (M == 8)  ? &interpf_12d_8m :
        (M == 16) ? &interpf_12d_16m:
                    &interpf_12d_xm):
        &interpf_xd_xm
        );
    pFxdxn->initAllocFxdxn->fnalloc(pFltr, M, D);
    return pFxdxn;
};

/* Calculate the memory block size for an interpolator with given attributes. */
size_t firinterpf_alloc(int D, int M)
{
	size_t sz;
    tFilterLayout fltr;
    NASSERT(((M == 2 || M == 4 || M == 8 || M == 16) && (D == 2 || D == 3 || D == 4 || D == 6 || D == 12)) || (D >= 2 && M % 8 == 0));

    getLayout(&fltr, M, D);
    sz=(ALIGNED_SIZE(sizeof(firinterpf_t), 4) +
        ALIGNED_SIZE(fltr.delLength * 2 * sizeof(float32_t), 2 * BBE_SIMD_WIDTH) +
        ALIGNED_SIZE(fltr.coefNum       * sizeof(float32_t), 2 * BBE_SIMD_WIDTH));
	return sz;
} // firinterpf_alloc()

/* Initialize the interpolator structure. The delay line is zeroed. */
firinterpf_handle_t firinterpf_init(void * objmem, int D, int M, const float32_t * restrict h)
{
    const tFirFxdxns * pFxdxn;
    tFilterLayout fltr;
    firinterpf_t * firinterp;
    void * ptr;
    float32_t * restrict pd;
    float32_t * restrict ph;
    int m;

    NASSERT(objmem && h);
    NASSERT(((M == 2 || M == 4 || M == 8 || M == 16) && (D == 2 || D == 3 || D == 4 || D == 6 || D == 12)) || (D >= 2 && M % 8 == 0));

    pFxdxn = getLayout(&fltr, M, D);

    /* Partition the memory block */
    ptr = objmem;
    firinterp = (firinterpf_t*)ALIGNED_ADDR(ptr, 4);
    ptr = firinterp + 1;
    ph = (float32_t*)ALIGNED_ADDR(ptr, 32);
    ptr = ph + fltr.coefNum;
    pd = (float32_t*)ALIGNED_ADDR(ptr, 32);
    ptr = pd + 2 * fltr.delLength;
    ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)firinterpf_alloc(D, M));

    /* Copy filter coefficients in reverted order and zero the delay line. */
    pFxdxn->initAllocFxdxn->fninit(ph, h, M, D);
    for (m = 0; m < 2 * fltr.delLength; m++) pd[m] = 0.f;

    firinterp->magic = FIRINTERPF_MAGIC;
    firinterp->M = M;
    firinterp->D = D;
    firinterp->h = ph;
    firinterp->d = pd;
    firinterp->procFxn = pFxdxn->procFxdxn;

    return NatureDSP_Baseband_isPresent(firinterp->procFxn) ? (firinterp) : NULL;
} // firinterpf_init()

/* process block of samples */
void firinterpf_process(firinterpf_handle_t handle, complex_float * restrict y, const complex_float * restrict x, int N)
{
    int M, D;
    firinterpf_t * firinterp = (firinterpf_t*)handle;
    NASSERT(firinterp && firinterp->magic == FIRINTERPF_MAGIC && y && x);
    M = firinterp->M;
    D = firinterp->D;

    if (N <= 0) return;

    NASSERT(!(N & 3));
    NASSERT(((M == 2 || M == 4 || M == 8 || M == 16) && (D == 2 || D == 3 || D == 4 || D == 6 || D == 12)) || (D >= 2 && M % 8 == 0));
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    NASSERT_ALIGN32(firinterp->d);
    NASSERT_ALIGN32(firinterp->h);
    NASSERT(firinterp->procFxn);

    firinterp->procFxn(y, x, firinterp->h, firinterp->d, M, N, D);
} // firinterpf_process()
