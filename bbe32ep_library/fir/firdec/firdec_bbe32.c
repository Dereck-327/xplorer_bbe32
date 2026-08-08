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
#include "firdec_common.h"

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
#define FIRDEC_MAGIC 0xfcb66012

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

#define sz_i16   sizeof(int16_t)

typedef struct
{
    uint32_t   magic;
    int        M, D;
    int16_t *  coef;
    int16_t *  delayLine;
    proc_fxdxn_t proc_fxn;
}
firdec_;

static const tFirFxdxns* getLayout(tFilterLayout* pFltr, int M, int D)
{
    const tFirFxdxns* pFxdxn;
    //
    // Select the filter processing function, delay line length and 
    // coefficients number.
    pFxdxn = ((D == 2) ? (
        (M == 2) ? &firdec_2d_2_8n :
        (M == 4) ? &firdec_2d_4_8n :
        (M == 8) ? &firdec_2d_8_8n :
        (M == 16) ? &firdec_2d_16_8n :
        (M == 32) ? &firdec_2d_32_8n :
        &firdec_2d_x_8n
        ) :
        (D == 3) ? (
        (M == 2) ? &firdec_3d_2_8n :
        (M == 4) ? &firdec_3d_4_8n :
        (M == 8) ? &firdec_3d_8_8n :
        (M == 16) ? &firdec_3d_16_8n :
        (M == 32) ? &firdec_3d_32_8n :
        &firdec_3d_x_8n
        ) :
        (D == 4) ? (
        (M == 2) ? &firdec_4d_2_8n :
        (M == 4) ? &firdec_4d_4_8n :
        (M == 8) ? &firdec_4d_8_8n :
        (M == 16) ? &firdec_4d_16_8n :
        (M == 32) ? &firdec_4d_32_8n :
        (M % (BBE_SIMD_WIDTH * 2)) ? &firdec_4d_xodd_8n : &firdec_4d_xeven_8n
        ) :
        (&firdec_dx_8n)
        );
    pFxdxn->initAllocFxdxn->fnalloc(pFltr, M, D);
    return pFxdxn;
};

/* Object allocation */
size_t firdec_alloc ( int D, int M )
{
    tFilterLayout fltr;
    NASSERT(((M == 2 || M == 4 || M == 8 || M == 16 || M == 32) && (D == 2 || D == 3 || D == 4)) ||
        (D >= 2 && M % 16 == 0));

    getLayout(&fltr, M, D);
    return (ALIGNED_SIZE(sizeof(firdec_), 4) +
        ALIGNED_SIZE(fltr.delLength * 2 * sz_i16, 2 * BBE_SIMD_WIDTH) +
        ALIGNED_SIZE(fltr.coefNum       * sz_i16, 2 * BBE_SIMD_WIDTH));
} /* firdec_alloc() */

/* Object initialization */
firdec_handle_t firdec_init ( void * objmem,  int D, int M, const int16_t * h )
{
    const tFirFxdxns  *   pFxdxn;
    firdec_ *  fir;
    int16_t *  delayLine;
    int16_t *  coef;
    void    *   ptr;
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
    fir = (firdec_*)ALIGNED_ADDR(ptr, 4);
    ptr = fir + 1;
    delayLine = (int16_t*)ALIGNED_ADDR(ptr, 2 * BBE_SIMD_WIDTH);
    ptr = delayLine + 2 * fltr.delLength;
    coef = (int16_t*)ALIGNED_ADDR(ptr, 2 * BBE_SIMD_WIDTH);
    ptr = coef + fltr.coefNum;

    NASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)firdec_alloc(D, M));
    
    pFxdxn->initAllocFxdxn->fninit(coef, h, M, D);

    //
    // Zero the delay line.
    //
    for (m = 0; m<2 * fltr.delLength; m++)
    {
        delayLine[m] = 0;
    }

    //
    // Initialize the filter instance.
    //
    fir->magic = FIRDEC_MAGIC;
    fir->M = fltr.coefNum;
    fir->D = D;
    fir->coef = coef;
    fir->delayLine = delayLine;
    fir->proc_fxn = pFxdxn->procFxdxn;

    return NatureDSP_Baseband_isPresent(fir->proc_fxn) ? (fir) : NULL;
} /* firdec_init() */

/* Update the delay line and compute filter output */
void firdec_process ( firdec_handle_t handle, 
                      complex_fract16 * restrict y, const complex_fract16 * restrict x, int N)
{
    firdec_ * fir = (firdec_*)handle;

    if (N <= 0) return;

    NASSERT(fir && fir->magic == FIRDEC_MAGIC);

    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);

    NASSERT_ALIGN32(fir->delayLine);
    NASSERT_ALIGN32(fir->coef);

    (*fir->proc_fxn)((int16_t *)y, (int16_t *)x, fir->coef, fir->delayLine, fir->M, N, fir->D);
} /* firdec_process() */
