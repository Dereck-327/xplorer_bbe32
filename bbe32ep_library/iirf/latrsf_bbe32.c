/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
  NatureDSP_Baseband library. IIR part
    Lattice real block IIR, floating point, streaming version
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/*-------------------------------------------------------------------------
  Streaming Lattice Block Real IIR

  Passes a few streams of real data through an autoregressive lattice IIR
  filter with real reflection coefficients. Please refer to the NatureDSP
  Baseband Library Reference for details on the filter structure.

  Representation:
  latrs      16-bit signed fixed-point format for input/output samples 
             and reflection coefficients
  latrsf     IEEE-754 Std. single precision floating-point format for
             input/output data and reflection coefficients

  Implementation of lattice filter for fixed-point data format has two flavors:
  Low_Noise  Performs 32x16-bit multiplications in the feedback branch to attain
             a precise output sample; feedforward products that update the delay
             line are still produced by 16x16-bit multiplications, but the full
             32-bit result is preserved for each delay element to be used on
             future iterations
  Fast       Makes use of native 16x16-bit multipliers and allocates 16 bits for
             a delay element

  Note:
    Streaming lattice block real IIR filter may by applied to complex
    data stored in streaming format, which is equivalent to passing a
    complex signal through a complex IIR with real coefficients. The number
    of real data streams L as it is seen by the filter implementation must
    be two times the number of actual complex streams.

  Methods:
    latrs[f]_alloc()   - calculate memory block size for a filter instance
    latrs[f]_init()    - initialize the filter instance
    latrs[f]_process() - pass a block of input samples to the filter and return
                         filter output samples
  Restrictions:
    - input/output data arrays must be aligned on 32-byte boundary
    - filter order cannot exceed 8, M <= 8
    - number of streams must be a multiple of 16 (latrs) or 8 (latrsf)
---------------------------------------------------------------------------*/  

#include <string.h>
/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Baseband_iir.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Filter instance and processing functions. */
#include "latrsf_common.h"

#if !(HAVE_VFPU)
DISCARD_FUN(void, latrsf_process,( latrsf_handle_t      _latrs,
                                   float32_t * restrict r,
                             const float32_t *          x,
                                   int                  N ))
/*---------------------------------------------------------------------------
   Calculate the memory block size for a filter with given attributes.
   Input:
     L     Number of streams
     M     Filter order
     opt   Fixed-point implementation variant selector, LATRS_OPT_LOW_NOISE
           or LATRS_OPT_FAST
   Output:
     Returns required memory block size, in bytes, or zero if failed
   Restrictions:
     M<=8  Filter order cannot exceed 8 
     L     Must be a multiple of 16 (latrs) or 8 (latrsf)
---------------------------------------------------------------------------*/
size_t latrsf_alloc( int L, int M )
{
    return (0);
}
/*---------------------------------------------------------------------------
   Initialize the filter instance. Delay elements are zeroed.
   Input:
     objmem   Memory block of latrs_alloc(L,M,opt) or latrsf_alloc(L,M) bytes
     L        Number of streams
     M        Filter order
     opt      Fixed-point implementation variant selector, LATRS_OPT_LOW_NOISE
              or LATRS_OPT_FAST
     coef[M]  Filter coefficients; coef[0] matches the newest delayed sample.
              Fixed-point data format for latrs_init() is Q15
     gain     Total gain coefficient to be applied to the input. Fixed-point
              data format for latrs_init() is Q15
   Output:
     Returns the filter handle, or zero if failed.
   Restrictions:
     M<=8     Filter order must not exceed 8
     L        Must be a multiple of 16 (latrs) or 8 (latrsf)
---------------------------------------------------------------------------*/
latrsf_handle_t latrsf_init( void      * objmem,
                             int         L,
                             int         M,
                       const float32_t * coef,
                             float32_t   gain )
{
    return NULL;
}
#else

/* Instance pointer validation number. */
#define MAGIC     0x5BB52112

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
        ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
        ( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

/*---------------------------------------------------------------------------
   Calculate the memory block size for a filter with given attributes.
   Input:
     L     Number of streams
     M     Filter order
     opt   Fixed-point implementation variant selector, LATRS_OPT_LOW_NOISE
           or LATRS_OPT_FAST
   Output:
     Returns required memory block size, in bytes, or zero if failed
   Restrictions:
     M<=8  Filter order cannot exceed 8 
     L     Must be a multiple of 16 (latrs) or 8 (latrsf)
---------------------------------------------------------------------------*/
size_t latrsf_alloc( int L, int M )
{
  int M_;

  ASSERT( (L%8)==0 );
  if ( M<=0 || L<=0 )
  {
    return 0;
  }

  M_ = (M+BBE_SIMD_WIDTH/2-1) & ~(BBE_SIMD_WIDTH/2-1);

  return ( ALIGNED_SIZE( sizeof( latrsf_t ), sizeof(int) )
           + // Lattice delay line
           ALIGNED_SIZE( L*M*sizeof(float32_t), BBE_SIMD_WIDTH/2*sizeof(float32_t) )
           + // Filter coefficients
           ALIGNED_SIZE( M_*sizeof(float32_t), BBE_SIMD_WIDTH/2*sizeof(float32_t) ) );
} /* latrsf_alloc() */

/*---------------------------------------------------------------------------
   Initialize the filter instance. Delay elements are zeroed.
   Input:
     objmem   Memory block of latrs_alloc(L,M,opt) or latrsf_alloc(L,M) bytes
     L        Number of streams
     M        Filter order
     opt      Fixed-point implementation variant selector, LATRS_OPT_LOW_NOISE
              or LATRS_OPT_FAST
     coef[M]  Filter coefficients; coef[0] matches the newest delayed sample.
              Fixed-point data format for latrs_init() is Q15
     gain     Total gain coefficient to be applied to the input. Fixed-point
              data format for latrs_init() is Q15
   Output:
     Returns the filter handle, or zero if failed.
   Restrictions:
     M<=8     Filter order must not exceed 8
     L        Must be a multiple of 16 (latrs) or 8 (latrsf)
---------------------------------------------------------------------------*/
latrsf_handle_t latrsf_init( void      * objmem,
                             int         L,
                             int         M,
                       const float32_t * coef,
                             float32_t   gain )
{
  latrsf_ptr_t  latrs;
  void        * ptr;
  float32_t   * delLine;
  float32_t   * cf;
  int m, M_;

  //
  // Validate the arguments.
  //
  
  ASSERT((L%8)==0 && objmem && coef );
  if ( M <= 0 || L <= 0 )
  {
    return NULL;
  }
  
  M_ = (M+BBE_SIMD_WIDTH/2-1) & ~(BBE_SIMD_WIDTH/2-1);

  //
  // Partition the memory block.
  //

  ptr     = objmem;
  latrs   = (latrsf_ptr_t)ALIGNED_ADDR( ptr, sizeof(int) );
  ptr     = latrs + 1;
  delLine = (float32_t*)ALIGNED_ADDR( ptr, BBE_SIMD_WIDTH/2*sizeof(float32_t) );
  ptr     = delLine + L*M;
  cf      = (float32_t*)ALIGNED_ADDR( ptr, BBE_SIMD_WIDTH/2*sizeof(float32_t) );
  ptr     = cf + M_;

  NASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)latrsf_alloc( L, M ) );

  //
  // Copy and pad the filter coefficients, zero the delay line.
  //

  for ( m=0; m<M; m++ )
  {
    cf[m] = coef[m];
  }
  for ( m=M; m<M_; m++ )
  {
    cf[m] = 0.0f;
  }

  memset( delLine, 0, L*M*sizeof(float32_t) );

  //
  // Initialize the filter instance.
  //

  memset( latrs, 0, sizeof(*latrs) );

  latrs->magic   = MAGIC;
  latrs->L       = L;
  latrs->M       = M;
  latrs->delLine = delLine;
  latrs->coef    = cf;
  latrs->gain    = gain;
  {
        typedef latrsf_processFxns* fxns_ptr;
        static const fxns_ptr fxns[]=
        {
            latrsf_process1,
            latrsf_process2,
            latrsf_process3,
            latrsf_process4,
            latrsf_process5,
            latrsf_process6,
            latrsf_process7,
            latrsf_process8,
            latrsf_processX
        };
        M = (M>9) ? 9 : M;
        latrs->fxns= fxns[M-1];
  }

  return (latrs);
} /* latrsf_init() */

/*---------------------------------------------------------------------------
   Pass a block of input samples to the filter and return filter output
   samples
   Input:
     _latrs,_latrsf  Filter handle
     N               Input/output signal chunk size, in samples (per stream)
     x[N*L]          Input samples
   Output:
     r[N*L]          Output samples
   Restrictions:
     x,r             Must be aligned on 32-byte boundary
---------------------------------------------------------------------------*/
void latrsf_process( latrsf_handle_t      _latrs,
                     float32_t * restrict r,
               const float32_t *          x,
                     int                  N )
{
  latrsf_ptr_t latrs = (latrsf_ptr_t)_latrs;

  NASSERT( latrs && latrs->magic == MAGIC );
  NASSERT( r && x );
  NASSERT_ALIGN32( r );
  NASSERT_ALIGN32( x );
  NASSERT( latrs->fxns );
  NASSERT( latrs->M > 0 );
  
  if(N>0) latrs->fxns(latrs,r,x,N);
} /* latrsf_process() */

#endif /* if !(HAVE_VFPU) */
