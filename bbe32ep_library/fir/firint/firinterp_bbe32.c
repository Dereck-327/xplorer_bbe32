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
#include "firinterp_common.h"

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
#define FIRINTERP_MAGIC 0x514c3ea8

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

static const tFirFxdxns* getLayout(tFilterLayout* pFltr, int M, int D)
{
    const tFirFxdxns* pFxdxn;
    //
    // Select the filter processing function, delay line length and 
    // coefficients number.
    pFxdxn = ((D == 2) ? (
        (M == 2) ? &interp_2d_2_8n :
        (M == 4) ? &interp_2d_4_8n :
        (M == 8) ? &interp_2d_8_8n :
        (M == 16) ? &interp_2d_16_8n :
        (M == 32) ? &interp_2d_32_8n :
        &interp_2d_mx_8n
        ) :

        (D == 3) ? (
        (M == 2) ? &interp_3d_2_8n :
        (M == 4) ? &interp_3d_4_8n :
        (M == 8) ? &interp_3d_8_8n :
        (M == 16) ? &interp_3d_16_8n :
        (M == 32) ? &interp_3d_32_8n :
        &interp_3d_mx_8n
        ) :
        (D == 4) ? (
        (M == 2) ? &interp_4d_2_8n :
        (M == 4) ? &interp_4d_4_8n :
        (M == 8) ? &interp_4d_8_8n :
        (M == 16) ? &interp_4d_16_8n :
        (M == 32) ? &interp_4d_32_8n :
        &interp_4d_mx_8n
        ) :
        (D == 6) ? (
        (M == 2) ? &interp_6d_2_8n :
        (M == 4) ? &interp_6d_4_8n :
        (M == 8) ? &interp_6d_8_8n :
        (M == 16) ? &interp_6d_16_8n :
        (M == 32) ? &interp_6d_32_8n :
        &interp_6d_mx_8n
        ) :
        (D == 12) ? (
        (M == 2) ? &interp_12d_2_8n :
        (M == 4) ? &interp_12d_4_8n :
        (M == 8) ? &interp_12d_8_8n :
        (M == 16) ? &interp_12d_16_8n :
        (M == 32) ? &interp_12d_32_8n :
        &interp_12d_mx_8n
        ) :
        (
        &interp_dx_mx_8n
        )
        );
    pFxdxn->initAllocFxdxn->fnalloc(pFltr, M, D);
    return pFxdxn;
};

/* Object allocation */
size_t firinterp_alloc (int D, int M)
{
    NASSERT(M>1);
    NASSERT(D>1);
    tFilterLayout fltr;
    getLayout(&fltr, M, D);
    return (ALIGNED_SIZE(sizeof(firinterp_), 4)
        + // Delay line
        ALIGNED_SIZE(2 * fltr.delLength, 2 * BBE_SIMD_WIDTH)
        + // Filter coefficients
        ALIGNED_SIZE(2 * fltr.coefNum, 2 * BBE_SIMD_WIDTH));
} /* firinterp_alloc() */

/* Object initialization */
firinterp_handle_t firinterp_init ( void * objmem, int D, int M, const int16_t * h )
{
    const tFirFxdxns  *   pFxdxn;
    firinterp_        *   firinterp;
    void              *   ptr;
    int16_t           *   delLine;
    int16_t           *   coef;
    int m;
    tFilterLayout fltr;
    NASSERT(objmem != NULL);

    //
    // Select the filter processing function, delay line length and 
    // coefficients number.
    //
    pFxdxn = getLayout(&fltr, M, D);

    //
    // Partition the memory block
    //

    ptr = objmem;
    firinterp = (firinterp_ *)ALIGNED_ADDR(ptr, 4);
    ptr = firinterp + 1;
    delLine = (int16_t*)ALIGNED_ADDR(ptr, 2 * BBE_SIMD_WIDTH);
    ptr = delLine + fltr.delLength;
    coef = (int16_t*)ALIGNED_ADDR(ptr, 2 * BBE_SIMD_WIDTH);
    ptr = coef + fltr.coefNum;

    ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)firinterp_alloc(D, M));

    pFxdxn->initAllocFxdxn->fninit(coef, h, M, D);

    //
    // Zero the delay line.
    //

    for (m = 0; m<fltr.delLength; m++)
    {
        delLine[m] = 0;
    }
    //
    // Initialize the filter instance.
    //

    firinterp->magic = FIRINTERP_MAGIC;
    firinterp->M = M;
    firinterp->D = D;
    firinterp->coef = coef;
    firinterp->procFxdxn = pFxdxn->procFxdxn;
    firinterp->delayLine = delLine;

    return NatureDSP_Baseband_isPresent(firinterp->procFxdxn) ? ((firinterp_handle_t)firinterp) : NULL;
} /* firinterp_init() */

/* Update the delay line and compute filter output */
void firinterp_process ( firinterp_handle_t handle, 
                         complex_fract16 * restrict y, const complex_fract16 * restrict x , int N)
{
    firinterp_ *firinterp = (firinterp_*)handle;
    if (N <= 0) return;
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);
    firinterp->procFxdxn(handle, (int16_t *)y, (int16_t *)x, firinterp->coef, firinterp->delayLine, firinterp->M, N, firinterp->D);
} /* firinterp_process() */
