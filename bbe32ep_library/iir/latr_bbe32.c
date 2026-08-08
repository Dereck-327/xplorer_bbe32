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
    Lattice real block IIR
    C code optimized for BBE32
    IntegrIT, 2006-2017
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
/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_iir.h"

/* Instance pointer validation number. */
#define MAGIC     0x5BB5210F

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
      ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      ( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

/* Lattice filter processing function. */
typedef void (proc_fxn_t)( int16_t * restrict r,
                           int16_t * restrict d,
                     const int16_t *          x,
                     const int16_t *          coef,
                           int16_t            gain,
                           int N );
typedef proc_fxn_t *  proc_fxn_ptr_t;

/* Filter instance structure. */
typedef struct tag_latr_t
{
  int              magic;   // Instance pointer validation number
  int              opt;     // LATR_OPT_LOW_NOISE | LATR_OPT_FAST
  proc_fxn_ptr_t   procFxn; // Filter processing function
  int16_t        * delLine; // Lattice delay line
  int16_t        * coef;    // Filter coefficients
  int16_t          gain;    // Total gain
} latr_t, *latr_ptr_t;

/* Lattice real block IIR processing functions, Fast fixed-point implementation. */
proc_fxn_t latr_sp_proc2;
proc_fxn_t latr_sp_proc4;
proc_fxn_t latr_sp_proc6;
proc_fxn_t latr_sp_proc8;
static const proc_fxn_ptr_t latr_sp_proc_tbl[4] = 
{
    latr_sp_proc2,
    latr_sp_proc4,
    latr_sp_proc6,
    latr_sp_proc8
};

/* Lattice real block IIR processing functions, Low Noise fixed-point implementation. */
proc_fxn_t latr_dp_proc2;
proc_fxn_t latr_dp_proc4;
proc_fxn_t latr_dp_proc6;
proc_fxn_t latr_dp_proc8;
static const proc_fxn_ptr_t latr_dp_proc_tbl[4] = 
{
    latr_dp_proc2,
    latr_dp_proc4,
    latr_dp_proc6,
    latr_dp_proc8
};

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

size_t latr_alloc ( int M, int opt )
{
  ASSERT( M<=8 );
  if ( M <= 0 )
  {
    return 0;
  }

  switch ( opt )
  {
    case LATR_OPT_FAST:

        return ( ALIGNED_SIZE( sizeof( latr_t ), sizeof(int) )
                 + // Lattice delay line
                 ALIGNED_SIZE( BBE_SIMD_WIDTH*sizeof(int16_t), sizeof(int16_t)*BBE_SIMD_WIDTH )
                 + // Filter coefficients
                 ALIGNED_SIZE( BBE_SIMD_WIDTH*sizeof(int16_t), sizeof(int16_t)*BBE_SIMD_WIDTH ) );

    case LATR_OPT_LOW_NOISE:

        return ( ALIGNED_SIZE( sizeof( latr_t ), sizeof(int) )
                 + // Lattice delay line
                 ALIGNED_SIZE( BBE_SIMD_WIDTH*sizeof(int32_t), sizeof(int16_t)*BBE_SIMD_WIDTH )
                 + // Filter coefficients
                 ALIGNED_SIZE( BBE_SIMD_WIDTH*sizeof(int16_t), sizeof(int16_t)*BBE_SIMD_WIDTH ) );
    default:

        return 0;
  }
} /* latr_alloc() */

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

latr_handle_t latr_init ( void * objmem, int M, int opt,
                         const int16_t * coef, int16_t gain )
{
    
  latr_ptr_t       latr;
  void           * ptr;
  int16_t        * delLine;
  int              delSize;
  int16_t        * cf;
  proc_fxn_ptr_t   procFxn;
  int m, M_;

  //
  // Validate the arguments.
  //
  
  ASSERT( M<=8 && objmem && coef );
  if ( M <= 0 )
  {
    return NULL;
  }

  M_ = M;
  M = (M+1) & ~1;

  //
  // Select the taps/coefs number and processing function.
  //

  switch ( opt )
  {
  case LATR_OPT_FAST:

    procFxn = latr_sp_proc_tbl[(M>>1)-1];
    delSize = BBE_SIMD_WIDTH;
    break;

  case LATR_OPT_LOW_NOISE:
      
    procFxn = latr_dp_proc_tbl[(M>>1)-1];
    delSize = 2*BBE_SIMD_WIDTH;
    break;

  default:
    return NULL;
  }

  //
  // Partition the memory block.
  //

  ptr     = objmem;
  latr    = (latr_ptr_t)ALIGNED_ADDR( ptr, sizeof(int) );
  ptr     = latr + 1;
  delLine = (int16_t*)ALIGNED_ADDR( ptr, sizeof(int16_t)*BBE_SIMD_WIDTH );
  ptr     = delLine + delSize;
  cf      = (int16_t*)ALIGNED_ADDR( ptr, sizeof(int16_t)*BBE_SIMD_WIDTH );
  ptr     = cf + BBE_SIMD_WIDTH;

  NASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)latr_alloc( M, opt ) );

  //
  // Copy and pad the filter coefficients, zero the delay line.
  //

  for ( m=0; m<M_; m++ )
  {
    cf[m] = coef[m];
  }
  for ( m=M_; m<BBE_SIMD_WIDTH; m++ )
  {
    cf[m] = 0;
  }

  memset( delLine, 0, delSize*sizeof(int16_t) );

  //
  // Initialize the filter instance.
  //

  memset( latr, 0, sizeof(*latr) );

  latr->magic   = MAGIC;
  latr->opt     = opt;
  latr->procFxn = procFxn;
  latr->delLine = delLine;
  latr->coef    = cf;
  latr->gain    = gain;

  return (latr);

} /* latr_init() */

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

void latr_process ( latr_handle_t      _latr,
                   int16_t * restrict r,
             const int16_t *          x,
                   int                N )
{
  latr_ptr_t latr = (latr_ptr_t)_latr;

  NASSERT( latr && latr->magic == MAGIC );
  NASSERT( r && x );

  //
  // Run the processing function.
  //

  latr->procFxn( r,
                 latr->delLine,
                 x,
                 latr->coef,
                 latr->gain,
                 N );
} /* latr_process() */
