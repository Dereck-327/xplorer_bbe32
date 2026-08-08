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
    Block Complex FIR Filter
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
#include "cxfir_common.h"

/*-------------------------------------------------------------------------
Block Complex FIR Filter

Computes a complex FIR filter (direct-form) using complex IR stored in 
vector h. The complex data input is stored in vector x. The filter output
result is stored in vector y. The filter calculates N output samples using
M coefficients and requires last M+N-1 samples in the delay line. 

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in samples) depends on FIR order M and is defined by cxfir[f]_algDelay(M).

Representation:
cxfir   16-bit signed fixed-point format
        Filter coefficients are Q15
        Number of fractional bits for input/output samples is user-difined
cxfirf  IEEE-754 Std. single precision floating-point format for filter
        coefficients and input/output samples

Parameters:
Input:
objmem  Allocated memory block
h[2*M]  Filter coefficients; h[0]+j*h[1] is to be multiplied by the newest
        sample,
N       Length of sample block
M       Length of filter
x[N]    Complex input samples
Output:
y[N]    Complex output samples

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 8 (cxfir) or 4 (cxfirf)
M       2,4,8 or a positive multiple of 16

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=2,4,8 and 16.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, cxfir[f]_init returns NULL handle.
-------------------------------------------------------------------------*/

/* Instance pointer validation number. */
#define CXFIR_MAGIC 0xac289df1

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
} cxfir_;

/* Object allocation */
size_t cxfir_alloc ( int M )
{
    NASSERT(M == 2 || M == 4 || M == 8 || M == 16 || (M>16 && !(M & 15)));

    if (!(M & 7))
    {
        return (ALIGNED_SIZE(sizeof(cxfir_), 4)
            + // Delay line
            ALIGNED_SIZE((2 * M + 2 * BBE_SIMD_WIDTH)*sz_i16, 32)
            + // Filter coefficients
            ALIGNED_SIZE(2 * M*sz_i16, 32));
    }
    else
    {
        return (ALIGNED_SIZE(sizeof(cxfir_), 4)
            + // Delay line
            ALIGNED_SIZE(2 * (M + 8 - 1)*sz_i16, 32)
            + // Filter coefficients
            ALIGNED_SIZE(2 * M*sz_i16, 32));
    }
} /* cxfir_alloc() */

/* Object initialization */
cxfir_handle_t cxfir_init(void * objmem, int M, const complex_fract16 * restrict h)
{
  cxfir_  *   cxfir;
  void    *   ptr;
  int16_t *   delLine;
  int16_t *   coef;

  int m;

  NASSERT( objmem && h );
  NASSERT( M==2 || M==4 || M==8 || M==16 || ( M>16 && !(M&15) ) );

  cxfir = NULL;

  if (M>16 && !(M&15))
  {
    //
    // Partition the memory block
    //
    ptr     = objmem;
    cxfir   = (cxfir_*)ALIGNED_ADDR( ptr, 4 );
    ptr     = cxfir + 1;
    delLine = (int16_t*)ALIGNED_ADDR( ptr, 32 );
    ptr     = delLine + 2 * M + 2*BBE_SIMD_WIDTH;
    coef    = (int16_t*)ALIGNED_ADDR( ptr, 32 );
    ptr     = coef + 2 * M;

    ASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)cxfir_alloc( M ) );

    //
    // Copy the filter coefficients in reverted order and zero the delay line.
    //

    for ( m=0; m<M; m++ )
    {
        coef[2 * m + 0] = h[M - 1 - m].s.re;
        coef[2 * m + 1] = h[M - 1 - m].s.im;
    }

    for ( m=0; m<2*M + 2*BBE_SIMD_WIDTH; m++ )
    {
      delLine[m] = 0;
    }

    //
    // Initialize the filter instance.
    //

    cxfir->magic     = CXFIR_MAGIC;
    cxfir->M         = M;
    cxfir->coef      = coef;
    cxfir->delayLine = delLine;
    cxfir->procFxn   = cxfir_process_8m;
  }
  else
  {

    ptr     = objmem;
    cxfir   = (cxfir_*)ALIGNED_ADDR( ptr, 4 );
    ptr     = cxfir + 1;
    delLine = (int16_t*)ALIGNED_ADDR( ptr, 32 );
    ptr     = delLine + 2*(M + 8 - 1);
    coef    = (int16_t*)ALIGNED_ADDR( ptr, 32 );
    ptr     = coef + 2*M;
    ASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)cxfir_alloc( M ) );
    
    // Copy the filter coefficients in reverted order and zero the delay line.
    for ( m=0; m<M; m++ )
    {
        coef[2 * m + 0] = h[M - 1 - m].s.re;
        coef[2 * m + 1] = h[M - 1 - m].s.im;
    }

    for ( m=0; m<M+8-1; m++ )
    {
      delLine[2*m+0] = 0;
      delLine[2*m+1] = 0;
    }

    switch (M)
    {
        case 2:  cxfir->procFxn = cxfir_process_2;  break;
        case 4:  cxfir->procFxn = cxfir_process_4;  break;
        case 8:  cxfir->procFxn = cxfir_process_8;  break;
        case 16: cxfir->procFxn = cxfir_process_16; break;
        default: cxfir->procFxn = NULL;
    }
    //
    // Initialize the filter instance.
    //
    cxfir->magic     = CXFIR_MAGIC;
    cxfir->M         = M;
    cxfir->coef      = coef;
    cxfir->delayLine = delLine;
  }
  return NatureDSP_Baseband_isPresent(cxfir->procFxn) ? (cxfir):NULL;
} /* cxfir_init() */

/* Update the delay line and compute filter output */
void cxfir_process ( cxfir_handle_t  handle, 
                     complex_fract16 * restrict y, const complex_fract16 * restrict x, int N)
{
    int M;
    cxfir_ *cxfir = (cxfir_*)handle;
    NASSERT(cxfir && cxfir->magic == CXFIR_MAGIC && y && x);

    if (N <= 0) return;
    M = cxfir->M;

    NASSERT(!(N & 7));
    NASSERT(M == 2 || M == 4 || M == 8 || M == 16 || (M>16 && !(M & 15)));
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);

    NASSERT(cxfir->procFxn);

    cxfir->procFxn(y, x, cxfir->coef, cxfir->delayLine, M, N);
} /* cxfir_process() */

/* Return the algorithmic delay, in samples. */
int cxfir_algDelay ( int M )
{
  return (0);

} /* cxfir_algDelay() */
