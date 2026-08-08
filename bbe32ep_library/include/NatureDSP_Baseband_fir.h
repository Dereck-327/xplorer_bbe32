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
 * NatureDSP_Baseband Library API
 * FIR Filters and Related Functions
 */

#ifndef __NATUREDSP_BASEBAND_FIR_H
#define __NATUREDSP_BASEBAND_FIR_H

#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
FIR filters and Related Functions
bkfir        Block FIR filter, real data, real coefficients
cxfir        Block FIR filter, complex data, complex coefficients
rcfir        Block FIR filter, complex data, real coefficients
srfir        Block FIR filter, real data, real symmetric impulse response
srcfir       Block FIR filter, complex data, real symmetric impulse response
firdec       Decimating block FIR filter, complex data, real coefficients
firinterp    Interpolating block complex FIR filter
fir_convol   Circular/linear convolution for complex data
fir_rconvol  Circular/linear convolution for real data
fir_xcorr    Circular/linear correlation for complex data
fir_rxcorr   Circular/linear correlation for real data
fir_acorr    Autocorrelation for complex data
fir_racorr   Autocorrelation for real data
despread     Despreading
fir_blms     Blockwise Adaptive LMS Algorithm for Complex Data
===========================================================================*/

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
/* Filter instance structure is defined internally. */
typedef void * bkfir_handle_t;
typedef void * bkfirf_handle_t;
/* Object allocation */
size_t bkfir_alloc ( int M );
size_t bkfirf_alloc( int M );
/* Object initialization */
bkfir_handle_t  bkfir_init ( void * objmem, int M, const int16_t   * restrict h );
bkfirf_handle_t bkfirf_init( void * objmem, int M, const float32_t * restrict h );
/* Update the delay line and compute filter output */
void bkfir_process ( bkfir_handle_t  handle, int16_t   * restrict y, const int16_t   * restrict x , int N );
void bkfirf_process( bkfirf_handle_t handle, float32_t * restrict y, const float32_t * restrict x , int N );
/* Return the algorithmic delay, in samples. */
int bkfir_algDelay ( int M );
int bkfirf_algDelay( int M );

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
/* Filter instance structure is defined internally. */
typedef void * cxfir_handle_t;
typedef void * cxfirf_handle_t;
/* Object allocation */
size_t cxfir_alloc ( int M );
size_t cxfirf_alloc( int M );
/* Object initialization */
cxfir_handle_t  cxfir_init ( void * objmem, int M, const complex_fract16 * restrict h );
cxfirf_handle_t cxfirf_init( void * objmem, int M, const complex_float   * restrict h );
/* Update the delay line and compute filter output */
void cxfir_process ( cxfir_handle_t  handle, complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );
void cxfirf_process( cxfirf_handle_t handle, complex_float   * restrict y, const complex_float   * restrict x, int N );
/* Return the algorithmic delay, in samples. */
int cxfir_algDelay ( int M );
int cxfirf_algDelay( int M );

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
/* Filter instance structure is defined internally. */
typedef void * rcfir_handle_t;
typedef void * rcfirf_handle_t;
/* Object allocation */
size_t rcfir_alloc ( int M );
size_t rcfirf_alloc( int M );
/* Object initialization */
rcfir_handle_t rcfir_init ( void * objmem, int M, const int16_t   * restrict h );
rcfir_handle_t rcfirf_init( void * objmem, int M, const float32_t * restrict h );
/* Update the delay line and compute filter output */
void rcfir_process ( rcfir_handle_t  handle, complex_fract16* restrict y, const complex_fract16 * restrict x, int N );
void rcfirf_process( rcfirf_handle_t handle, complex_float  * restrict y, const complex_float   * restrict x, int N );
/* Return the algorithmic delay, in samples. */
int rcfir_algDelay ( int M );
int rcfirf_algDelay( int M );

/*-------------------------------------------------------------------------
Block Real FIR filter w/ Symmetric Impulse Response

Passes real input data through a direct-form FIR filter with real coefficients
and symmetric impulse response. The filter calculates N output samples using
(M+1)/2 coefficients and requires last M+N-1 samples in the delay line.

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in samples) depends on FIR order M and is defined by srfir_algDelay(M)

Precision: 16-bit data, 16-bit coefficients, 16-bit outputs, everything in Q15

Parameters:
Input:
objmem      Allocated memory block
h[(M+1)/2]  Filter coefficients (half the number of filter taps); h[0] is to
            be multiplied by the newest sample,  Q15
N           Length of sample block
M           Length of filter
x[N]        Input samples, Q15
Output:
y[N]        Output samples, Q15

Restrictions:
x,y         Must not overlap
x,y         Aligned on 32-byte boundary
M           2,4,8 or a multiple of 16
N           Multiple of 16  
M>0

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=2,4,8,16 and 32.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, srfir_init returns NULL handle.
---------------------------------------------------------------------------*/
/* Filter instance structure is defined internally. */
typedef void * srfir_handle_t;
/* Object allocation */
size_t srfir_alloc(  int M );
/* Object initialization */
srfir_handle_t srfir_init( void * mem, int M, const int16_t * h );
/* Update the delay line and compute filter output */
void srfir_process( srfir_handle_t  handle, 
                    int16_t * restrict y, const int16_t * restrict x , int N);
/* Return the algorithmic delay, in samples. */
int srfir_algDelay( int M );

/*-------------------------------------------------------------------------
Block Complex FIR filter w/ Real Symmetric Impulse Response

Passes complex input data through a direct-form FIR filter with real 
coefficients and symmetric impulse response. The filter calculates N output 
samples using (M+1)/2 coefficients and requires last M+N-1 samples in the delay 
line.

IMPORTANT NOTE:
Due to the performance reasons, implementation may introduce additional
algorithmic delay (group delay) to the output data. The amount of this delay 
(in complex samples) depends on FIR order M and is defined by srcfir_algDelay(M)

Precision: 16-bit data, 16-bit coefficients, 16-bit outputs, everything in Q15

Parameters:
Input:
objmem      Allocated memory block
h[(M+1)/2]  Filter coefficients (half the number of filter taps); h[0] is to
            be multiplied by the newest sample,  Q15
N           Length of sample block
M           Length of filter
x[N]        Input samples, Q15
Output:
y[N]        Output samples, Q15

Restrictions:
x,y         Must not overlap
x,y         Aligned on 32-byte boundary
M           17,33 or a multiple of 16
N           Multiple of 8
M>0

Note on performance:
Most efficient operation (maximal MACs per cycle count) is achieved for filter
lengths M=16,17,32 and 33.

Note on availability:
Depending on available ISA options, some combinations of filter parameters
may not be supported. In that case, srcfir_init returns NULL handle.
---------------------------------------------------------------------------*/
/* Filter instance structure is defined internally. */
typedef void * srcfir_handle_t;
/* Object allocation */
size_t srcfir_alloc( int M );
/* Object initialization */
srcfir_handle_t srcfir_init( void * mem, int M, const int16_t * h );
/* Update the delay line and compute filter output */
void srcfir_process( srcfir_handle_t handle, 
                     complex_fract16 * restrict y, const complex_fract16 * restrict x , int N);
/* Return the algorithmic delay, in samples. */
int srcfir_algDelay( int M );

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
/* Filter instance structure is defined internally. */
typedef void * firdec_handle_t ;
typedef void * firdecf_handle_t ;
/* Object allocation */
size_t firdec_alloc ( int D, int M );
size_t firdecf_alloc( int D, int M );
/* Object initialization */
firdec_handle_t  firdec_init ( void * objmem,  int D, int M, const int16_t   * h ); 
firdecf_handle_t firdecf_init( void * objmem,  int D, int M, const float32_t * h ); 
/* Update the delay line and compute filter output */
void firdec_process ( firdec_handle_t handle,  complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );
void firdecf_process( firdec_handle_t handle,  complex_float   * restrict y, const complex_float   * restrict x, int N );

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
/* Filter instance structure is defined internally. */
typedef void * firinterp_handle_t;
typedef void * firinterpf_handle_t;
/* Object allocation */
size_t firinterp_alloc ( int D, int M );
size_t firinterpf_alloc( int D, int M );
/* Object initialization */
firinterp_handle_t firinterp_init  ( void * objmem, int D, int M, const int16_t   * h );
firinterpf_handle_t firinterpf_init( void * objmem, int D, int M, const float32_t * h );
/* Update the delay line and compute filter output */
void firinterp_process ( firinterp_handle_t  handle, complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );
void firinterpf_process( firinterpf_handle_t handle, complex_float   * restrict y, const complex_float   * restrict x, int N );

/*-------------------------------------------------------------------------
Circular/Linear Convolution for Complex Data

Compute circular or linear convolution between complex-valued vectors x (of
length N) and y (of length M) resulting in vector r of N (circular) or
N+M-1 (linear) complex values.

MATLAB code for circular convolution:
  for n=1:N
    r(n) = 2*sum(x(1+mod(n-(1:M),N)).*y(1:M));
  end

MATLAB code for linear convolution:
  for n=1:N+M-1
    ix = max(n-N+1,1):min(n,M);
    r(n) = 2*sum(x(n+1-ix).*y(ix));
  end

Representation:
fir_convol   Signed fixed-point format
             Input data are 16-bit Q15, output data are 32-bit Q31
fir_convolf  IEEE-754 Std. single precision floating-point format for
             input/output data

Parameters:
Input:
x[N]         Left-hand complex data sequence
y[M]         Right-hand complex data sequence
Output:
r[]          Complex output data. Size is N for circular convolution, 
             or (N+M-1) for linear convolution
Restrictions:
x,y,r        Must not overlap
x,y,r        Must be aligned on 32-byte boundary
N,M          Multiples of 8 (fir_convol_circ, fir_convol_lin) or 4 
             (fir_convolf_circ, fir_convolf_lin)
N>=M         N must be greater than or equal to M
-------------------------------------------------------------------------*/
void fir_convol_circ ( complex_fract32   * restrict r,
                 const complex_fract16   * restrict x,
                 const complex_fract16   * restrict y,
                 int N, int M );
void fir_convolf_circ( complex_float * restrict r,
                 const complex_float * restrict x,
                 const complex_float * restrict y,
                 int N, int M );
void fir_convol_lin  ( complex_fract32   * restrict r,
                 const complex_fract16   * restrict x,
                 const complex_fract16   * restrict y,
                 int N, int M );
void fir_convolf_lin ( complex_float * restrict r,
                 const complex_float * restrict x,
                 const complex_float * restrict y,
                 int N, int M );

/*-------------------------------------------------------------------------
Circular/Linear Convolution for Real Data

Compute circular or linear convolution between real-valued vectors x (of
length N) and y (of length M) resulting in vector r of N (circular) or
N+M-1 (linear) real values.

MATLAB code for circular convolution:
  for n=1:N
    r(n) = 2*sum(x(1+mod(n-(1:M),N)).*y(1:M));
  end

MATLAB code for linear convolution:
  for n=1:N+M-1
    ix = max(n-N+1,1):min(n,M);
    r(n) = 2*sum(x(n+1-ix).*y(ix));
  end

Representation: IEEE-754 Std. single precision floating-point format for
             input/output data

Parameters:
Input:
x[N]    Left-hand data sequence
y[M]    Right-hand data sequence
Output:
r[]     Output data. Size is N for circular convolution, or N+M-1 
        for linear convolution
Restrictions:
x,y,r   Must not overlap
x,y,r   Must be aligned on 32-byte boundary
N,M     Multiples of 8
N>=M    N must be greater than or equal to M
-------------------------------------------------------------------------*/
void fir_rconvolf_circ( float32_t * restrict r,
                  const float32_t * restrict x,
                  const float32_t * restrict y,
                  int N, int M );
void fir_rconvolf_lin ( float32_t * restrict r,
                  const float32_t * restrict x,
                  const float32_t * restrict y,
                  int N, int M );

/*-------------------------------------------------------------------------
Circular/Linear Cross-Correlation for Complex Data

Estimates the circular/linear cross-correlation between complex-valued
vectors x (of length N) and y (of length M) resulting in vector r. It
is similar to convolution, but y is read in opposite direction.

MATLAB code for circular cross-correlation:
  for n=1:N
    r(n) = 2*sum(x(1+mod(n-1+(0:M-1),N)).*conj(y(1:M)));
  end

MATLAB code for linear cross-correlation:
  for n=1:M+N-1
    ix = max(n-M+1,1):min(n,N);
    r(n) = 2*sum(x(ix).*conj(y(M-(n-ix))));
  end;

Representation:
fir_xcorr   Signed fixed-point format
            Input data are 16-bit Q15, output data are 32-bit Q31
fir_xcorrf  IEEE-754 Std. single precision floating-point format for
            input/output data

Parameters:
Input:
x[N]        Left-hand complex data sequence
y[M]        Right-hand complex data sequence
N           Length of x
M           Length of y
Output:
r[]         Complex output data. Size is N for circular cross-correlation,
            or (N+M-1) for linear cross-correlation
Restrictions:
x,y,r       Must not overlap
x,y,r       Must be aligned on 32-byte boundary
N,M         Multiples of 8 (fir_xcorr_circ, fir_xcorr_lin) or 4 
            (fir_xcorrf_circ, fir_xcorrf_lin)
N>=M        N must be greater than or equal to M

NOTES:
1.fir_xcorr[f]_lin returns the same output as MATLAB xcorr but omits first 
N-M zeros
-------------------------------------------------------------------------*/
void fir_xcorr_circ ( complex_fract32   * restrict r,
                const complex_fract16   * restrict x,
                const complex_fract16   * restrict y,
                int N, int M );
void fir_xcorrf_circ( complex_float * restrict r,
                const complex_float * restrict x,
                const complex_float * restrict y,
                int N, int M );
void fir_xcorr_lin  ( complex_fract32   * restrict r,
                const complex_fract16   * restrict x,
                const complex_fract16   * restrict y,
                int N, int M );
void fir_xcorrf_lin ( complex_float * restrict r,
                const complex_float * restrict x,
                const complex_float * restrict y,
                int N, int M );

/*-------------------------------------------------------------------------
Circular/Linear Cross-Correlation for Real Data

Estimates the circular/linear cross-correlation between real-valued vectors
x (of length N) and y (of length M) resulting in vector r. It is similar to
convolution, but y is read in opposite direction.

MATLAB code for circular cross-correlation:
  for n=1:N
    r(n) = 2*sum(x(1+mod(n-1+(0:M-1),N)).*y(1:M));
  end

MATLAB code for linear cross-correlation:
  for n=1:M+N-1
    ix = max(n-M+1,1):min(n,N);
    r(n) = 2*sum(x(ix).*y(M-(n-ix)));
  end;

Representation: IEEE-754 Std. single precision floating-point format for
             input/output data

Parameters:
Input:
x[N]         Left-hand data sequence
y[M]         Right-hand data sequence
N            Length of x
M            Length of y
Output:
r[]          Output data. Size is N for circular cross-correlation, or N+M-1 for
             linear cross-correlation
Restrictions:
x,y,r        Must not overlap
x,y,r        Must be aligned on 32-byte boundary
N,M          Multiples of 8
N>=M         N must be greater than or equal to M

NOTES:
1. fir_rxcorrf_lin returns the same output as MATLAB xcorr but omits first 
   N-M zeros
-------------------------------------------------------------------------*/
void fir_rxcorrf_circ( float32_t * restrict r,
                 const float32_t * restrict x,
                 const float32_t * restrict y,
                 int N, int M );
void fir_rxcorrf_lin ( float32_t * restrict r,
                 const float32_t * restrict x,
                 const float32_t * restrict y,
                 int N, int M );

/*-------------------------------------------------------------------------
Autocorrelation for a Ñomplex Data Vector

Estimates the auto-correlation of complex-valued vector x, positive side
only. Returns autocorrelation of length N. For an input vector of N complex
samples x[0..N-1] the computation follows the MATLAB code given below:

  for n = 1:N
    r(n) = sum(x(n:N).*conj(x(1:(N+1-n))));
  end

Representation:
fir_acorr   Signed fixed-point format
            Input data are 16-bit Q15, output data are 32-bit Q31
fir_acorrf  IEEE-754 Std. single precision floating-point format for
            input/output data

Parameters:
Input:
x[N]        Complex input data
N           Length of x
Output:
r[N]        Complex output data

Restrictions:
x,r         Must not overlap
x,r         Must be aligned on 32-byte boundary
N           Multiple of 8 (fir_acorr) or 4 (fir_acorrf)
-------------------------------------------------------------------------*/
void fir_acorr ( complex_fract32   * restrict r,
           const complex_fract16   * restrict x,
           int N );
void fir_acorrf( complex_float * restrict r,
           const complex_float * restrict x,
           int N );

/*-------------------------------------------------------------------------
Autocorrelation for a Real Data Vector

Estimates the auto-correlation of real-valued vector x, positive side 
only. Returns autocorrelation of length N. For an input vector of N 
samples x[0..N-1] the computation follows the MATLAB code given below:

  for n = 1:N
    r(n) = sum(x(n:N).*x(1:(N+1-n)));
  end

Representation: IEEE-754 Std. single precision floating-point format for
             input/output data

Parameters:
Input:
x[N]         Input data
N            Length of x
Output:
r[N]         Output data

Restrictions:
x,r          Must not overlap
x,r          Must be aligned on 32-byte boundary
N            Multiple of 8 (fir_racorrf)
-------------------------------------------------------------------------*/
void fir_racorrf( float32_t * restrict r,
            const float32_t * restrict x,
            int N );

/*-------------------------------------------------------------------------
Despreading

Perform a special kind of correlation utilized in Direct-Sequence Spread Spectrum
(DSSS) modulation technique, with Spreading Factor (SF) 4 or 8. Despreading
functions multiply the input sequence of 16-bit complex samples x[] by a coded
pseudo-noise sequence of 2-bit QPSK symbols pn[], then reduce the product sequence
by adding together 4 or 8 contiguous products. Resulting sums are shifted to the
right by rsh bit positions, saturated to -32768 or 32767 and stored to the output
array y[].

For the coded sequence multiplicand pn[], there are two variants of 2-bit QPSK
codeset identified with the qpsk_type argument:

         | QPKS symbol for | QPKS symbol for 
   Dibit | qpsk_type == 0  | qpsk_type == 1
  -------+-----------------+-----------------
     00  |      1 + j      |       1
     01  |     -1 + j      |      -j
     10  |      1 - j      |      -1
     11  |     -1 - j      |       j

First 8 consecutive dibits of the coded sequence are concatenated in a 16-bit
word pn[0], with the very first dibit going into 2 LSBs of pn[0] and the last
dibit - into 2 MSBs. The next 8 dibits are stored in pn[1], and so on.

Parameters:
  Input:
  x[N]          Input signal, 16-bit complex samples
  pn[N/8]       Pseudo-noise coded sequence of 2-bit QPSK symbols
  qpsk_type     QPSK codeset selection, 0 or 1
  rsh           Right shift amount for reduced product sums, 0..31
  Output:
  y[N/SF]       Output signal, 16-bit complex samples. Fixed point position
                for the output signal is Qx-rsh, where Qx is the fixed point
                position for the input signal.
Restrictions:
  x[],y[],pn[]  Must not overlap
  x[],y[],pn[]  Must be aligned on 32-byte boundary
  N             Multiple of 8
 -------------------------------------------------------------------------*/
void despread4( complex_fract16 * restrict y,
          const complex_fract16 * restrict x,
          const int16_t * restrict pn,
          int N, int qpsk_type, int rsh );

void despread8( complex_fract16 * restrict y,
          const complex_fract16 * restrict x,
          const int16_t * restrict pn,
          int N, int qpsk_type, int rsh );

/*-------------------------------------------------------------------------
Blockwise Adaptive LMS Algorithm for Complex Data

Blockwise LMS algorithm performs filtering of complex input samples, 
computation of error over a block of reference samples and makes blockwise
update of IR to minimize the error output.
Algorithm includes FIR filtering, calculation of correlation between the 
error output and reference signal and IR taps update based on that 
correlation.
NOTES: 
1.  For N=1 this algorithm is equivalent to standard LMS, however a bigger 
    block size reduces the computational overhead keeping adaptation 
    properties.
2.  Right selection of N depends on the change rate of impulse response. 
    However, on static or slow varying channels convergence rate depends 
    on selected mu and M, but not on N.
3.  Computation of filter output is done using only higher 16-bit words of 
    IR coefficients (in fact converting them to Q15). However, when it 
    performs coefficient update, it uses full Q30 accuracy. 
4.  Functions use BBE32EP 40-bit accumulators for all computational steps. 
    Thus, the user needs to avoid very long blocks to prevent overflows 
    during IR update step. Typically, N should not exceed 128 to guarantee 
    the 40-bit range.

Parameters:
Temporary:
pScr          Scratch memory area of FIR_BLMS_SCRATCH(M,N) bytes
Input:
h[M]          Complex impulse response, Q30
r[N]          Reference (near end) complex data vector, Q15. First in time 
              value is in r[0].
x[(N+M-1)]    Input (far end) complex data vector, Q15. First in time value 
              is in x[0].
norm          Normalization factor: power of signal multiplied by N, Q31
mu            Adaptation coefficient (LMS step), Q15
N             Length of data block
M             Length of h
Output
e[N]          Estimated error, Q15
h[M]          Updated impulse response, Q30

Restrictions:
pScr,e,h,r,x  Must not overlap
pScr,e,h,r,x  Aligned on 32-byte boundary
N             Multiple of 8 
M             Multiple of 8 (applies to fir_blms() function only)
-------------------------------------------------------------------------*/
void fir_blms  ( void * pScr, complex_fract16 * restrict e, complex_fract32 * restrict h,
                 const complex_fract16 * restrict r,
                 const complex_fract16 * restrict x,
                 int32_t norm, int16_t mu,
                 int N, int M
                 );
void fir_blms4 ( void * pScr, complex_fract16 * restrict e, complex_fract32 * restrict h,
                 const complex_fract16 * restrict r,
                 const complex_fract16 * restrict x,
                 int32_t norm, int16_t mu,
                 int N
                 );
void fir_blms8 ( void * pScr, complex_fract16 * restrict e, complex_fract32 * restrict h,
                 const complex_fract16 * restrict r,
                 const complex_fract16 * restrict x,
                 int32_t norm, int16_t mu,
                 int N
                 );
void fir_blms16( void * pScr, complex_fract16 * restrict e, complex_fract32 * restrict h,
                 const complex_fract16 * restrict r,
                 const complex_fract16 * restrict x,
                 int32_t norm, int16_t mu,
                 int N
                 );
void fir_blms32( void * pScr, complex_fract16 * restrict e, complex_fract32 * restrict h,
                 const complex_fract16 * restrict r,
                 const complex_fract16 * restrict x,
                 int32_t norm, int16_t mu,
                 int N
                 );

#define FIR_BLMS_SCRATCH(M,N)  (4*4*(M)+4*4*(N)+8*(M))

#ifdef __cplusplus
};
#endif

#endif /* __NATUREDSP_BASEBAND_FIR_H */
