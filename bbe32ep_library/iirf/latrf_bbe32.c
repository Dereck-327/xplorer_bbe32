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
    Lattice real block IIR, floating point
    C code optimized for BBE32
    Integrit, 2006-2016
*/

/*-------------------------------------------------------------------------
  Lattice Real Block IIR

  Passes real input signal through an autoregressive lattice IIR filter with
  real reflection coefficients. Please refer to the NatureDSP Baseband Library
  Reference for details on the filter structure.

  Representation:
  latr       16-bit signed fixed-point format for input/output samples 
             and reflection coefficients
  latrf      IEEE-754 Std. single precision floating-point format for
             input/output data and reflection coefficients

  Implementation of lattice filter for fixed-point data format has two flavors:
  Low_Noise  Performs 32x16-bit multiplications in the feedback branch to attain
             a precise output sample; feedforward products that update the delay
             line are still produced by 16x16-bit multiplications, but the full
             32-bit result is preserved for each delay element to be used on
             future iterations
  Fast       Makes use of native 16x16-bit multipliers and allocates 16 bits for
             a delay element

  Methods:
    latr[f]_alloc()   - calculate memory block size for a filter instance
    latr[f]_init()    - initialize the filter instance
    latr[f]_process() - pass a block of input samples to the filter and return
                        filter output samples
  Restrictions:
    - filter order must not exceed 8: M <= 8
---------------------------------------------------------------------------*/

#include <string.h>
/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Baseband_iir.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Filter instance and processing functions. */
#include "latrf_common.h"

#if !(HAVE_VFPU)
DISCARD_FUN(void, latrf_process,( latrf_handle_t       _latr,
                                  float32_t * restrict r,
                            const float32_t *          x,
                                  int                  N ))
/*-------------------------------------------------------------------------
   Calculate the memory block size for a filter with given attributes.
   Input:
     M     Filter order
     opt   Fixed-point implementation variant selector, LATR_OPT_LOW_NOISE
           or LATR_OPT_FAST
   Output:
     Returns required memory block size, in bytes, or zero if failed
   Restrictions:
     M<=8  Filter order must not exceed 8
-------------------------------------------------------------------------*/
size_t latrf_alloc( int M )
{
    return (0);
}
/*-------------------------------------------------------------------------
   Initialize the filter instance. Delay elements are zeroed.
   Input:
     objmem   Memory block of latr_alloc(M,opt) or latrf_alloc(M) bytes
     M        filter order
     opt      Fixed-point implementation variant selector, LATR_OPT_LOW_NOISE
              or LATR_OPT_FAST
     coef[M]  Filter coefficients; coef[0] matches the newest delayed sample.
              Fixed-point data format for latr_init() is Q15
     gain     Total gain coefficient to be applied to the input. Fixed-point
              data format for latr_init() is Q15
   Output:
     Returns the filter handle, or zero if failed.
   Restrictions:
     M<=8     Filter order must not exceed 8
-------------------------------------------------------------------------*/
latrf_handle_t latrf_init( void      * objmem,
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
        (void*)( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

/*-------------------------------------------------------------------------
   Calculate the memory block size for a filter with given attributes.
   Input:
     M     Filter order
     opt   Fixed-point implementation variant selector, LATR_OPT_LOW_NOISE
           or LATR_OPT_FAST
   Output:
     Returns required memory block size, in bytes, or zero if failed
   Restrictions:
     M<=8  Filter order must not exceed 8
-------------------------------------------------------------------------*/
size_t latrf_alloc( int M )
{
  if ( M <= 0 )
  {
    return 0;
  }

  M = (M+BBE_SIMD_WIDTH/2-1) & ~(BBE_SIMD_WIDTH/2-1);

  return ( ALIGNED_SIZE( sizeof( latrf_t ), sizeof(int) )
           + // Lattice delay line
           ALIGNED_SIZE( M*sizeof(float32_t), BBE_SIMD_WIDTH/2*sizeof(float32_t) )
           + // Filter coefficients
           ALIGNED_SIZE( M*sizeof(float32_t), BBE_SIMD_WIDTH/2*sizeof(float32_t) ) );
} /* latrf_alloc() */

/*-------------------------------------------------------------------------
   Initialize the filter instance. Delay elements are zeroed.
   Input:
     objmem   Memory block of latr_alloc(M,opt) or latrf_alloc(M) bytes
     M        filter order
     opt      Fixed-point implementation variant selector, LATR_OPT_LOW_NOISE
              or LATR_OPT_FAST
     coef[M]  Filter coefficients; coef[0] matches the newest delayed sample.
              Fixed-point data format for latr_init() is Q15
     gain     Total gain coefficient to be applied to the input. Fixed-point
              data format for latr_init() is Q15
   Output:
     Returns the filter handle, or zero if failed.
   Restrictions:
     M<=8     Filter order must not exceed 8
-------------------------------------------------------------------------*/
latrf_handle_t latrf_init( void      * objmem,
                           int         M,
                     const float32_t * coef,
                           float32_t   gain )
{
  latrf_ptr_t  latr;
  void       * ptr;
  float32_t  * delLine;
  float32_t  * cf;
  int m, M_;

  //
  // Validate the arguments.
  //
  
  ASSERT( objmem && coef );
  if ( M <= 0 )
  {
    return NULL;
  }

  //
  // Partition the memory block.
  //
  
  M_ = (M+BBE_SIMD_WIDTH/2-1) & ~(BBE_SIMD_WIDTH/2-1);

  ptr     = objmem;
  latr    = (latrf_ptr_t)ALIGNED_ADDR( ptr, sizeof(int) );
  ptr     = latr + 1;
  delLine = (float32_t*)ALIGNED_ADDR( ptr, BBE_SIMD_WIDTH/2*sizeof(float32_t) );
  ptr     = delLine + M_;
  cf      = (float32_t*)ALIGNED_ADDR( ptr, BBE_SIMD_WIDTH/2*sizeof(float32_t) );
  ptr     = cf + M_;

  NASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)latrf_alloc( M ) );

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

  memset( delLine, 0, M_*sizeof(float32_t) );

  //
  // Initialize the filter instance.
  //

  memset( latr, 0, sizeof(*latr) );

  latr->magic   = MAGIC;
  latr->M       = M;
  latr->delLine = delLine;
  latr->coef    = cf;
  latr->gain    = gain;
  {
        typedef latrf_processFxns* fxns_ptr;
        static const fxns_ptr fxns[]=
        {
            latrf_process2,
            latrf_process4,
            latrf_process6,
            latrf_process8,
            latrf_processX
        };
        M = (M>9) ? 9 : M;
        latr->fxns= fxns[(M-1)>>1];
  }

  return (latr);
} /* latrf_init() */

/*-------------------------------------------------------------------------
   Pass a block of input samples to the filter and return filter output
   samples
   Input:
     _latr,_latrf  Filter handle
     N             Input/output signal chunk size, in samples
     x[N]          Input samples
   Output:
     r[N]          Output samples
   Restrictions:
     None 
-------------------------------------------------------------------------*/
void latrf_process( latrf_handle_t       _latr,
                    float32_t * restrict r,
              const float32_t *          x,
                    int                  N )
{
  latrf_ptr_t latr = (latrf_ptr_t)_latr;

  NASSERT( latr && latr->magic == MAGIC );
  NASSERT( r && x );
  NASSERT( latr->fxns );
  NASSERT( latr->M>0 );
  
  if(N>0) latr->fxns(latr,r,x,N);
} /* latrf_process() */

#endif /* if !(HAVE_VFPU) */
