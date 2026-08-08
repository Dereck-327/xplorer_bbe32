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
#include "rcfir_common.h"

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
#define RCFIR_MAGIC 0xdb5afd8d

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

#define sz_i16   sizeof(int16_t)

typedef struct
{
    uint32_t        magic;     // Instance pointer validation number
    int             M;         // Number of filter coefficients
    const int16_t * coef;      // Filter coefficients
    proc_fxn_t      procFxn;   // Filter processing function
    int16_t *       delayLine; // Delay line for samples
} rcfir_;

/* Object allocation */
size_t rcfir_alloc ( int M )
{
  NASSERT( M==2 || M==4 || M==8 || M==16 ||( M>16 && !(M&15) ) );
  
  return ( ALIGNED_SIZE( sizeof( rcfir_ ), 4 )
           + // Delay line
           ALIGNED_SIZE( ( 2*M + 2*BBE_SIMD_WIDTH )*sz_i16, 32 )
           + // Filter coefficients
           ALIGNED_SIZE( M*sz_i16, 32 ) );
} /* rcfir_alloc() */

/* Object initialization */
rcfir_handle_t rcfir_init ( void * objmem, int M, const int16_t * restrict h )
{
    rcfir_  *   rcfir;
    void    *   ptr;
    int16_t *   delLine;
    int16_t *   coef;

    int m;

    NASSERT(objmem && h);
    NASSERT(M == 2 || M == 4 || M == 8 || M == 16 || (M>16 && !(M & 15)));

    rcfir = NULL;

    if (M>16)
    {
        //
        // Partition the memory block
        //
        ptr     = objmem;
        rcfir   = (rcfir_ *)ALIGNED_ADDR(ptr, 4);
        ptr     = rcfir + 1;
        delLine = (int16_t*)ALIGNED_ADDR(ptr, 32);
        ptr     = delLine + 2 * M + 2 * BBE_SIMD_WIDTH;
        coef    = (int16_t *)ALIGNED_ADDR(ptr, 32);
        ptr     = coef + M;

        ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)rcfir_alloc(M));

        //
        // Copy the filter coefficients in reverted order and zero the delay line.
        //

        for (m = 0; m<M / 2; m++)
        {
            coef[2 * m + 0] = h[M - 1 - 2 * m - 1];
            coef[2 * m + 1] = h[M - 1 - 2 * m - 0];
        }

        for (m = 0; m<2 * M + 2 * BBE_SIMD_WIDTH; m++)
        {
            delLine[m] = 0;
        }

        //
        // Initialize the filter instance.
        //

        rcfir->magic = RCFIR_MAGIC;
        rcfir->M = M;
        rcfir->coef = coef;
        rcfir->delayLine = delLine;
        rcfir->procFxn = rcfir_process_16m;
    }
    else
    {
        ptr = objmem;
        rcfir = (rcfir_ *)ALIGNED_ADDR(ptr, 4);
        ptr = rcfir + 1;
        coef = (int16_t *)ALIGNED_ADDR(ptr, 32);
        ptr = coef + M;
        delLine = (int16_t*)ALIGNED_ADDR(ptr, 32);
        ptr = delLine + 2 * M + 2 * BBE_SIMD_WIDTH;
        
        ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)rcfir_alloc(M));

        rcfir->magic = RCFIR_MAGIC;
        rcfir->M = M;
        rcfir->coef = coef;
        rcfir->delayLine = delLine;

        //
        // Copy the filter coefficients in reverted order and zero the delay line.
        //
        {
            valign vh;
            const xb_vecNx16 *ph;
            xb_vecNx16 t;

            // clean delay line
            t = BBE_ZERONX16();
            BBE_SVNX16_I(t, (xb_vecNx16*)delLine, 0);
            if (M == 16)
            {
                BBE_SVNX16_I(t, (xb_vecNx16*)delLine, 2 * BBE_SIMD_WIDTH);
            }

            // copy IR
            ph = (const xb_vecNx16 *)h;
            vh = BBE_LA_PP(ph);
            BBE_LAVNX16_XP(t, vh, ph, sizeof(int16_t)*M);
            BBE_SVNX16_I(t, (xb_vecNx16*)coef, 0);
        }

        switch (M)
        {
            case 2:  rcfir->procFxn = rcfir_process_2;   break;
            case 4:  rcfir->procFxn = rcfir_process_4;   break;
            case 8:  rcfir->procFxn = rcfir_process_8;   break;
            case 16: rcfir->procFxn = rcfir_process_16;  break;
            default: rcfir->procFxn = NULL;
        }
    }
    return NatureDSP_Baseband_isPresent(rcfir->procFxn) ? (rcfir) : NULL;
} /* rcfir_init() */

/* Update the delay line and compute filter output */
void rcfir_process ( rcfir_handle_t  handle, 
                     complex_fract16 * restrict y, const complex_fract16 * restrict x, int N )
{
    int M;
    rcfir_ *rcfir = (rcfir_*)handle;
    NASSERT(rcfir && rcfir->magic == RCFIR_MAGIC && y && x);

    if (N <= 0) return;
    M = rcfir->M;

    NASSERT(!(N & 7));
    NASSERT(M == 2 || M == 4 || M == 8 || M == 16 || (M>16 && !(M & 15)));
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);

    NASSERT(rcfir->procFxn);

    rcfir->procFxn(rcfir->delayLine, y, x, rcfir->coef, rcfir->M, N);
} /* rcfir_process() */

/* Return the algorithmic delay, in samples. */
int rcfir_algDelay ( int M )
{
  return (0);

} /* rcfir_algDelay() */
