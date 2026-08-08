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
#include "bkfir_common.h"

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
#define BKFIR_MAGIC 0x6543face

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

#define sz_i16   sizeof(int16_t)

/* Filter instance structure. */
typedef struct
{
    uint32_t        magic;     // Instance pointer validation number
    int             M;         // Number of filter coefficients
    const int16_t * coef;      // Filter coefficients
    proc_fxn_t      procFxn;   // Filter processing function
    int16_t *       delayLine; // Delay line for samples
} bkfir_t, *bkfir_ptr_t;

/* Object allocation */
size_t bkfir_alloc ( int M )
{
  NASSERT(((M >= 16) && !(M & 15)) || M == 2 || M == 4 || M == 8);

  return ( ALIGNED_SIZE( sizeof( bkfir_t ), 4 )
           + // Delay line
           ALIGNED_SIZE( ( M + 2*BBE_SIMD_WIDTH )*sz_i16, 32 )
           + // Filter coefficients
           ALIGNED_SIZE( M*sz_i16, 32 ) );
} /* bkfir_alloc() */

/* Object initialization */
bkfir_handle_t bkfir_init ( void * objmem, int M, const int16_t * restrict h )
{
  bkfir_ptr_t bkfir;
  void *      ptr;
  int16_t *   delLine;
  int16_t *   coef;
  int m;

  NASSERT( objmem && h );
  NASSERT(((M >= 16) && !(M & 15)) || M == 2 || M == 4 || M == 8);

  //
  // Partition the memory block
  //
  ptr     = objmem;
  bkfir   = (bkfir_ptr_t)ALIGNED_ADDR( ptr, 4 );
  ptr     = bkfir + 1;
  delLine = (int16_t*)ALIGNED_ADDR( ptr, 32 );
  ptr     = delLine + M + 2*BBE_SIMD_WIDTH;
  coef    = (int16_t *)ALIGNED_ADDR( ptr, 32 );
  ptr     = coef + M;
  ASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)bkfir_alloc( M ) );

  //
  // Copy the filter coefficients in reverted order and zero the delay line.
  //
  for ( m=0; m<M/2; m++ )
  {
    coef[2*m+0] = h[M-1-2*m-1];
    coef[2*m+1] = h[M-1-2*m-0];
  }
  for ( m=0; m<(M+2*BBE_SIMD_WIDTH); m++ )
  {
    delLine[m] = 0;
  }

  //
  // Initialize the filter instance.
  //
  bkfir->magic     = BKFIR_MAGIC;
  bkfir->M         = M;
  bkfir->coef      = coef;
  bkfir->delayLine = delLine;
  bkfir->procFxn =  ((M ==  2 )? bkfir_process_2  :
                     (M ==  4 )? bkfir_process_4  :
                     (M ==  8 )? bkfir_process_8  :
                    ((M == 16 )? bkfir_process_16 : bkfir_process_16m));
 
  return NatureDSP_Baseband_isPresent(bkfir->procFxn) ? (bkfir):NULL;
} /* bkfir_init() */

/* Update the delay line and compute filter output */
void bkfir_process ( bkfir_handle_t  handle, 
                     int16_t * restrict y, const int16_t * restrict x , int N)
{
    int M;
    bkfir_ptr_t bkfir = (bkfir_ptr_t)handle;
    NASSERT(bkfir && bkfir->magic == BKFIR_MAGIC && y && x);

    if (N <= 0) return;
    M = bkfir->M;

    NASSERT(!(N & 15));
    NASSERT(((M >= 16) && !(M & 15)) || M == 2 || M == 4 || M == 8);
    NASSERT_ALIGN32(y);
    NASSERT_ALIGN32(x);

    NASSERT(bkfir->procFxn);

    bkfir->procFxn(bkfir->delayLine, y, x, bkfir->coef, M, N);
} /* bkfir_process() */

/* Return the algorithmic delay, in samples. */
int bkfir_algDelay ( int M )
{
  return (0);

} /* bkfir_algDelay() */
