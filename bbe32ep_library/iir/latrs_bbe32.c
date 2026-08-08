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
    Lattice real block IIR, streaming version
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
/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_iir.h"
/* Processing functions declarations. */
#include "latrs_common.h"

/* Instance pointer validation number. */
#define MAGIC     0x4587BFD1

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
      ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      ( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

/* Filter instance structure. */
typedef struct tag_latr_t
{
  int              magic;   // Instance pointer validation number
  int              opt;     // LATRS_OPT_LOW_NOISE | LATRS_OPT_FAST
  int              L;       // Number of streams
  proc_fxn_ptr_t   procFxn; // Filter processing function
  int16_t        * delLine; // Lattice delay line
  int16_t        * coef;    // Filter coefficients
  int16_t          gain;    // Total gain

} latrs_t, *latrs_ptr_t;

static const proc_fxn_ptr_t sp_fxns[]=
{ NULL,
  latrs_sp_proc1,
  latrs_sp_proc2,
  latrs_sp_proc3,
  latrs_sp_proc4,
  latrs_sp_proc5,
  latrs_sp_proc6,
  latrs_sp_proc7,
  latrs_sp_proc8
};
static const proc_fxn_ptr_t dp_fxns[]=
{ NULL,
  latrs_dp_proc1,
  latrs_dp_proc2,
  latrs_dp_proc3,
  latrs_dp_proc4,
  latrs_dp_proc5,
  latrs_dp_proc6,
  latrs_dp_proc7,
  latrs_dp_proc8
};

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

size_t latrs_alloc ( int L, int M, int opt )
{

  ASSERT( M<=8 && (L%8)==0 );
  if ( M<=0 || L<=0 )
  {
    return 0;
  }

  switch ( opt )
  {
    case LATRS_OPT_FAST:

        return ( ALIGNED_SIZE( sizeof( latrs_t ), sizeof(int) )
                 + // Lattice delay line
                 ALIGNED_SIZE( M*L*sizeof(int16_t), sizeof(int16_t)*BBE_SIMD_WIDTH )
                 + // Filter coefficients
                 ALIGNED_SIZE( M*sizeof(int16_t), sizeof(int16_t)*2 ) );

    case LATRS_OPT_LOW_NOISE:

        return ( ALIGNED_SIZE( sizeof( latrs_t ), sizeof(int) )
                 + // Lattice delay line
                 ALIGNED_SIZE( M*L*sizeof(int32_t), sizeof(int16_t)*BBE_SIMD_WIDTH )
                 + // Filter coefficients
                 ALIGNED_SIZE( M*sizeof(int16_t), sizeof(int16_t)*2 ) );
    default:

        return 0;
  }


} /* latrs_alloc() */

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

latrs_handle_t latrs_init ( void * objmem, int L, int M, int opt,
                           const int16_t * coef, int16_t gain )
{
    
  latrs_ptr_t      latrs;
  void           * ptr;
  proc_fxn_ptr_t   procFxn;
  int16_t        * delLine;
  int              delSize;
  int16_t        * cf;

  int m;

  //
  // Validate the arguments.
  //
  
  ASSERT( M<=8 && (L%16)==0 && objmem && coef );
  if ( M <= 0 || L <= 0 )
  {
    return NULL;
  }

  //
  // Select the taps/coefs number and processing function.
  //

  switch ( opt )
  {
  case LATRS_OPT_FAST:

    delSize = M*L;
    procFxn = sp_fxns[M];
    break;

  case LATRS_OPT_LOW_NOISE:

    delSize = M*L*2;
    procFxn = dp_fxns[M];
    break;

  default:
    return NULL;
  }

  //
  // Partition the memory block.
  //

  ptr     = objmem;
  latrs   = (latrs_ptr_t)ALIGNED_ADDR( ptr, sizeof(int) );
  ptr     = latrs + 1;
  delLine = (int16_t*)ALIGNED_ADDR( ptr, sizeof(int16_t)*BBE_SIMD_WIDTH );
  ptr     = delLine + delSize;
  cf      = (int16_t*)ALIGNED_ADDR( ptr, sizeof(int16_t)*2 );
  ptr     = cf + M;

  NASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)latrs_alloc( L, M, opt ) );

  //
  // Copy and pad the filter coefficients, zero the delay line.
  //

  for ( m=0; m<M; m++ )
  {
    cf[m] = coef[m];
  }

  memset( delLine, 0, delSize*sizeof(int16_t) );

  //
  // Initialize the filter instance.
  //

  memset( latrs, 0, sizeof(*latrs) );

  latrs->magic   = MAGIC;
  latrs->opt     = opt;
  latrs->L       = L;
  latrs->delLine = delLine;
  latrs->coef    = cf;
  latrs->procFxn = procFxn;
  latrs->gain    = gain;

  return (latrs);

} /* latrs_init() */

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

void latrs_process ( latrs_handle_t     _latrs,
                     int16_t * restrict r,
               const int16_t *          x,
                     int                N )
{
  latrs_ptr_t latrs = (latrs_ptr_t)_latrs;

  ASSERT( latrs && latrs->magic == MAGIC && latrs->procFxn );
  ASSERT( r && x );
  NASSERT_ALIGN32( r );
  NASSERT_ALIGN32( x );

  ( *latrs->procFxn )( r,
                       latrs->delLine,
                       x,
                       latrs->coef,
                       latrs->gain,
                       N,
                       latrs->L );

} /* latrs_process() */
