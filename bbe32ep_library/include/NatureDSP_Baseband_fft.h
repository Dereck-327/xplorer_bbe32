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
 * FFT Routines
 */

#ifndef __NATUREDSP_BASEBAND_FFT_H
#define __NATUREDSP_BASEBAND_FFT_H

#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
Fixed point FFT routines:
cfftas      Radix-2 forward FFT on complex data, auto scaling
cifftas     Radix-2 inverse FFT on complex data, auto scaling
tfft        Radix-2 forward FFT on complex data with reduced twiddle table
tifft       Radix-2 inverse FFT on complex data with reduced twiddle table
bcfft       Blockwise radix-2 forward FFT on complex data, no data scaling
bcifft      Blockwise radix-2 inverse FFT on complex data, no data scaling
cnfft       Mixed radix forward FFT on complex data, auto scaling
cinfft      Mixed radix inverse FFT on complex data, auto scaling
bcnfft      Blockwise mixed radix forward FFT on complex data, no data scaling
bcinfft     Blockwise mixed radix inverse FFT on complex data, no data scaling
rfft        Radix-2 fixed point forward FFT on real data, auto scaling
rifft       Radix-2 fixed point inverse FFT forming real data, auto scaling

Floating point FFT routines:
cfftf      FFT on complex data
cifftf     inverse FFT on complex data
tfftf      FFT on complex data with reduced twiddle table
tifft      inverse FFT on complex data with reduced twiddle table
bcfft      Blockwise forward FFT on complex data
bcifftf    Blockwise inverse FFT on complex data
brfft      Blockwise forward FFT on real data
brifftf    Blockwise inverse FFT forming real data
rfftf      forward FFT on real data
rifftf     inverse FFT forming real data

===========================================================================*/

/*-------------------------------------------------------------------------
Radix-2 forward FFT on complex data, auto scaling

Description: These functions make forward FFT on complex data of power of 2
sizes: N=2^n, n=4..15. Functions with _norm suffix expect input data to be
normalized, i.e. the minimum number of redundant sign bits over x[]
(a.k.a the common block exponent) should be zero. Neglecting to normalize
data leads to significant loss in transform quality. On the contrary, regular
variants with no _norm suffix allow for non-zero common block exponent, but
they appear slightly slower due to internal data normalization.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[N]        Complex input signal
    bexp        Common block exponent, that is the minimum number of redundant
                sign bits over input data x[]
  Output:            
    y[N]        Output spectrum samples
  Returned value:
                Total shift amount applied throughout the transform to scale
                the data, with positive numbers corresponding to the right
                shift. _norm-suffixed functions return strictly positive
                values, while for regular variants the total shift amount is
                bi-directional.
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int cfftas16         ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cfftas16_norm    ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cfftas32         ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cfftas32_norm    ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cfftas64         ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cfftas64_norm    ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cfftas128        ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cfftas128_norm   ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cfftas256        ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cfftas256_norm   ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cfftas512        ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cfftas512_norm   ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cfftas1024       ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cfftas1024_norm  ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cfftas2048       ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cfftas2048_norm  ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cfftas4096       ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cfftas4096_norm  ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cfftas8192       ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cfftas8192_norm  ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cfftas16384      ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cfftas16384_norm ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cfftas32768      ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cfftas32768_norm ( complex_fract16 * restrict y, complex_fract16 * restrict x           );

/*-------------------------------------------------------------------------
Radix-2 inverse FFT on complex data, auto scaling

Description: These functions make inverse FFT on complex data of power of 2
sizes: N=2^n, n=4..15. Functions with _norm suffix expect input data to be
normalized, i.e. the minimum number of redundant sign bits over x[]
(a.k.a the common block exponent) should be zero. Neglecting to normalize
data leads to significant loss in transform quality. On the contrary, regular
variants with no _norm suffix allow for non-zero common block exponent, but
they appear slightly slower due to internal data normalization.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[N]        Input spectrum samples
    bexp        Common block exponent, that is the minimum number of redundant
                sign bits over input data x[]
  Output:            
    y[N]        Complex output signal
  Returned value:
                Total shift amount applied throughout the transform to scale
                the data, with positive numbers corresponding to the right
                shift. _norm-suffixed functions return strictly positive
                values, while for regular variants the total shift amount is
                bi-directional.
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int cifftas16         ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cifftas16_norm    ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cifftas32         ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cifftas32_norm    ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cifftas64         ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cifftas64_norm    ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cifftas128        ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cifftas128_norm   ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cifftas256        ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cifftas256_norm   ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cifftas512        ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cifftas512_norm   ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cifftas1024       ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cifftas1024_norm  ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cifftas2048       ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cifftas2048_norm  ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cifftas4096       ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cifftas4096_norm  ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cifftas8192       ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cifftas8192_norm  ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cifftas16384      ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cifftas16384_norm ( complex_fract16 * restrict y, complex_fract16 * restrict x           );
int cifftas32768      ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cifftas32768_norm ( complex_fract16 * restrict y, complex_fract16 * restrict x           );

/*-------------------------------------------------------------------------
Radix-2 forward FFT on complex data with reduced twiddle table, auto scaling

Description: These functions make forward FFT on complex data of power of 2
sizes: N=2^n, n=13..15. As opposed to regular FFT routines, these use
smaller twiddle factor tables but work a bit slower.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[N  ]      Complex input signal
    bexp        Common block exponent, that is the minimum number of redundant
                sign bits over input data x[]
  Output:            
    y[N]      Output spectrum samples
  Returned value:
                Total shift amount applied throughout the transform to scale
                the data. Total shift is bi-directional, with positive numbers
                corresponding to the right shift.
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int tfft8192  ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int tfft16384 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int tfft32768 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );

/*-------------------------------------------------------------------------
Radix-2 inverse FFT on complex data with reduced twiddle table, auto scaling

Description: These functions make inverse FFT on complex data of power of 2
sizes: N=2^n, n=13..15. As opposed to regular FFT routines, these use
smaller twiddle factor tables but work a bit slower.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[N  ]      Input spectrum samples
    bexp        Common block exponent, that is the minimum number of redundant
                sign bits over input data x[]
  Output:            
    y[N]        Complex output signal
  Returned value:
                Total shift amount applied throughout the transform to scale
                the data. Total shift is bi-directional, with positive numbers
                corresponding to the right shift.
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int tifft8192  ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int tifft16384 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int tifft32768 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );

/*-------------------------------------------------------------------------
Blockwise radix-2 forward FFT on complex data, no data scaling

Description: These functions make forward FFT on L blocks, each of N=2^n
complex samples, where n=4..7. It is user's responsibility to pre-scale input
data in such a way that FFT calculation overflows are avoided.

Precision: 16-bit input, 16-bit output
Scaling  : none

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[L][N]   Complex input signal
  Output:          
    y[L][N]   Output spectrum samples
  Returned value:
                None
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
void bcfft16 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcfft32 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcfft64 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcfft128( complex_fract16 * restrict y, complex_fract16 * restrict x, int L );

/*-------------------------------------------------------------------------
Blockwise radix-2 inverse FFT on complex data, no data scaling

Description: These functions make inverse FFT on L blocks, each of N=2^n
complex samples, where n=4..7. It is user's responsibility to pre-scale input
data in such a way that FFT calculation overflows are avoided.

Precision: 16-bit input, 16-bit output
Scaling  : none

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[L][N]   Complex input signal
  Output:          
    y[L][N]   Output spectrum samples
  Returned value:
                None
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
void bcifft16 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcifft32 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcifft64 ( complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcifft128( complex_fract16 * restrict y, complex_fract16 * restrict x, int L );

/*-------------------------------------------------------------------------
Mixed radix forward FFT on complex data, auto scaling
  
Description: These functions make forward FFT on complex data of the following
sizes: N = 12,24,36,48,60,72,96,108,120,144,180,192,216,240,288,300,324,360,
384,432,480,540,576,600,648,720,768,864,900,960,972,1080,1152,1200,1536.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

  Parameters:
  Temprorary:
    pScr[]      Scratch memory area of CNFFT_SCRATCH_SIZE(N) bytes
  Input:
    S           Required input/output buffer size may exceed actual data size:
                S >= N. Use CNFFT_BUF_SIZE(N) macro to determine the minimum
                buffer size expressed in complex 16-bit elements
    x[S]        N complex samples of input signal 
    bexp        Common block exponent, that is the minimum number of redundant
                sign bits over input data x[]
  Output:            
    y[S]        N complex samples of output spectrum
  Returned value:
                Total shift amount applied throughout the transform to scale
                the data. Total shift is bi-directional, with positive numbers
                corresponding to the right shift.
Restrictions:
  x,y,pScr      Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int cnfft12   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft24   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft36   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft48   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft60   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft72   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft96   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft108  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft120  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft144  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft180  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft192  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft216  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft240  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft288  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft300  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft324  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft360  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft384  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft432  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft480  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft540  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft576  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft600  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft648  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft720  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft768  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft864  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft900  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft960  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft972  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft1080 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft1152 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft1200 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cnfft1536 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );

// Size of input/output buffers for cnfft functions, in complex_fract16 elements
#define CNFFT_BUF_SIZE(N)  ((      \
         ((N)==12  ?   32  :0) +\
         ((N)==96  ?   192  :0) +\
         ((N)==192  ?   384  :0) +\
         ((N)==288  ?   576  :0) +\
         ((N)==480  ?   960  :0) +\
         ((N)==384  ?   768  :0) +\
         ((N)==576  ?   1152  :0) +\
         ((N)==768  ?   1536  :0) +\
         ((N)==864  ?   1728  :0) +\
         ((N)==960  ?   1920  :0) +\
         ((N)==1152  ?   2304  :0) +\
         ((N)==1200  ?   2432  :0) +\
         ((N)==1536  ?   3072  :0) +\
         ((N)==24  ?   64  :0) +\
         ((N)==36  ?   128  :0) +\
         ((N)==60  ?   192  :0) +\
         ((N)==72  ?   192  :0) +\
         ((N)==108  ?   320  :0) +\
         ((N)==324  ?   896  :0) +\
         ((N)==972  ?   2624  :0) +\
         ((N)==120  ?   288  :0) +\
         ((N)==180  ?   480  :0) +\
         ((N)==216  ?   480  :0) +\
         ((N)==300  ?   768  :0) +\
         ((N)==360  ?   768  :0) +\
         ((N)==540  ?   1344  :0) +\
         ((N)==600  ?   1248  :0) +\
         ((N)==648  ?   1344  :0) +\
         ((N)==900  ?   2208  :0) +\
         ((N)==1080  ?   2208  :0) +\
         ((N)==48  ?   128  :0) +\
         ((N)==144  ?   320  :0) +\
         ((N)==240  ?   512  :0) +\
         ((N)==432  ?   896  :0) +\
         ((N)==720  ?   1472  :0) \
         )/2)

// Scratch area size, in bytes.
#define CNFFT_SCRATCH_SIZE(N) ((\
         ((N)==108  ?   320  :0) + \
         ((N)==972  ?   2624  :0) +\
         ((N)==180  ?   480  :0) +\
         ((N)==768  ?   1536  :0) +\
         ((N)==300  ?   768  :0) +\
         ((N)==648  ?   1344  :0))*2)

/*-------------------------------------------------------------------------
Mixed radix inverse FFT on complex data, auto scaling
  
Description: These functions make inverse FFT on complex data of the following
sizes: N = 12,24,36,48,60,72,96,108,120,144,180,192,216,240,288,300,324,360,
384,432,480,540,576,600,648,720,768,864,900,960,972,1080,1152,1200,1536.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

  Parameters:
  Temprorary:
    pScr[]      Scratch memory area of CINFFT_SCRATCH_SIZE(N) bytes
  Input:            
    S           Required input/output buffer size may exceed actual data size:
                S >= N. Use CINFFT_BUF_SIZE(N) macro to determine the minimum
                buffer size expressed in complex 16-bit elements
    x[S]        N complex samples of input spectrum
    bexp        Common block exponent, that is the minimum number of redundant
                sign bits over input data x[]
  Output:            
    y[S]        N complex samples of output signal
  Returned value:
                Total shift amount applied throughout the transform to scale
                the data. Total shift is bi-directional, with positive numbers
                corresponding to the right shift.
Restrictions:
  x,y,pScr      Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int cinfft12   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft24   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft36   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft48   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft60   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft72   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft96   ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft108  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft120  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft144  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft180  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft192  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft216  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft240  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft288  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft300  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft324  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft360  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft384  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft432  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft480  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft540  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft576  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft600  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft648  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft720  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft768  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft864  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft900  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft960  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft972  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft1080 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft1152 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft1200 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );
int cinfft1536 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int bexp );

// Size of input/output buffers for cnfft functions, in complex_fract16 elements
#define CINFFT_BUF_SIZE(_N_)  CNFFT_BUF_SIZE(_N_)

// Scratch area size, in bytes
#define CINFFT_SCRATCH_SIZE(_N_) CNFFT_SCRATCH_SIZE(_N_)

/*-------------------------------------------------------------------------
Blockwise mixed radix forward FFT on complex data, no data scaling
  
Description: These functions make forward FFT on complex data of the following
sizes: N = 12,24,36,48,60,72,96,108,120. It is user's responsibility to pre-scale input
data in such a way that FFT calculation overflows are avoided.

Precision: 16-bit input, 16-bit output
Scaling  : none

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

  Parameters:
  Temprorary:
    pScr[]      Scratch memory area of BCNFFT_SCRATCH_SIZE(N) bytes
  Input:
    S           Required input/output buffer size may exceed actual data size:
                S >= N. Use BCNFFT_BUF_SIZE(N) macro to determine the minimum
                buffer size expressed in complex 16-bit elements
    x[L][S]     Complex input signal
  Output:            
    y[L][S]     Output spectrum samples
  Returned value:
                None
Restrictions:
  x,y,pScr      Must not overlap and must be aligned on 32-byte boundary
  L>0           The number of blocks must be positive
-------------------------------------------------------------------------*/
void bcnfft12  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcnfft24  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcnfft36  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcnfft48  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcnfft60  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcnfft72  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcnfft96  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcnfft108 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcnfft120 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );

// Size of input/output buffers for cnfft functions, in complex_fract16 elements
#define BCNFFT_BUF_SIZE(N) ((        \
         ((N)==12  ?   32  :0) +\
         ((N)==96  ?   192   :0) +\
         ((N)==24  ?   48  :0) +\
         ((N)==36  ?   80  :0) +\
         ((N)==48  ?   96  :0) +\
         ((N)==60  ?   160  :0) +\
         ((N)==72  ?   192  :0) +\
         ((N)==108  ?   288  :0) +\
         ((N)==120  ?   240  :0) )/2)
         
// Scratch area size, in bytes
#define BCNFFT_SCRATCH_SIZE(N,L) (  \
         ( (N)==108 ? (L)*288*2 : 0 ) )

/*-------------------------------------------------------------------------
Blockwise mixed radix inverse FFT on complex data, no data scaling
  
Description: These functions make inverse FFT on complex data of the following
sizes: N = 12,24,36,48,60,72,96,108,120. It is user's responsibility to pre-scale input
data in such a way that FFT calculation overflows are avoided.

Precision: 16-bit input, 16-bit output
Scaling  : none

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

  Parameters:
  Temprorary:
    pScr[]      Scratch memory area of BCINFFT_SCRATCH_SIZE(N) bytes
  Input:
    S           Required input/output buffer size may exceed actual data size:
                S >= 2*N. Use BCINFFT_BUF_SIZE(N) macro to determine the minimum
                buffer size expressed in complex 16-bit elements
    x[L][S]     Complex input signal
  Output:            
    y[L][S]     Output spectrum samples
  Returned value:
                None
Restrictions:
  x,y,pScr      Must not overlap and must be aligned on 32-byte boundary
  L>0           The number of blocks must be positive
-------------------------------------------------------------------------*/
void bcinfft12  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcinfft24  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcinfft36  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcinfft48  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcinfft60  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcinfft72  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcinfft96  ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcinfft108 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );
void bcinfft120 ( void * restrict pScr, complex_fract16 * restrict y, complex_fract16 * restrict x, int L );

// Size of input/output buffers for cnfft functions, in complex_fract16 elements
#define BCINFFT_BUF_SIZE(_N_)  BCNFFT_BUF_SIZE(_N_)

// Scratch area size, in bytes
#define BCINFFT_SCRATCH_SIZE(_N_,L) BCNFFT_SCRATCH_SIZE(_N_,L)

/*-------------------------------------------------------------------------
Radix-2 forward FFT on real data, auto scaling

Description: These functions make FFT on real data of length N=2^n, n=4..15.
The algorithm exploits the symmetry properties of the FFT: first, a complex
FFT of half the original size is applied to input data, then the resulting
spectrum undergoes a postprocessing procedure which results in complex spectrum
of real input data.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:      
    x[N]          Real input signal
    bexp          Common block exponent, that is the minimum number of redundant
                  sign bits over input data x[N]
  Output:      
    y[(N/2+1)]    Output spectrum samples. 
  Returned value:
                  Total shift amount applied throughout the transform to scale
                  the data. Total shift is bi-directional, with positive numbers
                  corresponding to the right shift.
Restrictions:
  x,y             Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int rfft16    ( complex_fract16 * restrict y, int16_t * restrict x, int bexp );
int rfft32    ( complex_fract16 * restrict y, int16_t * restrict x, int bexp );
int rfft64    ( complex_fract16 * restrict y, int16_t * restrict x, int bexp );
int rfft128   ( complex_fract16 * restrict y, int16_t * restrict x, int bexp );
int rfft256   ( complex_fract16 * restrict y, int16_t * restrict x, int bexp );
int rfft512   ( complex_fract16 * restrict y, int16_t * restrict x, int bexp );
int rfft1024  ( complex_fract16 * restrict y, int16_t * restrict x, int bexp );
int rfft2048  ( complex_fract16 * restrict y, int16_t * restrict x, int bexp );
int rfft4096  ( complex_fract16 * restrict y, int16_t * restrict x, int bexp );
int rfft8192  ( complex_fract16 * restrict y, int16_t * restrict x, int bexp );
int rfft16384 ( complex_fract16 * restrict y, int16_t * restrict x, int bexp );
int rfft32768 ( complex_fract16 * restrict y, int16_t * restrict x, int bexp );

/*-------------------------------------------------------------------------
Radix-2 inverse FFT forming real data, auto scaling

Description: These functions make inverse FFT forming real data of length
N=2^n, n=4..15. Algorithm exploits the symmetry properties of the FFT:
the input spectrum is modified in such a way that a complex-valued inverse FFT
of half the original size that is applied to the transformed spectrum actually
results in real data.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:
    x[(N/2+1)]    Input spectrum samples
    bexp          Common block exponent, that is the minimum number of redundant
                  sign bits over input data x[]
  Output:
    y[N]          Real output signal
  Returned value:
                  Total shift amount applied throughout the transform to scale
                  the data. Total shift is bi-directional, with positive numbers
                  corresponding to the right shift.
Restrictions:
  x,y             Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int rifft16    ( int16_t * restrict y, complex_fract16 * restrict x, int bexp );
int rifft32    ( int16_t * restrict y, complex_fract16 * restrict x, int bexp );
int rifft64    ( int16_t * restrict y, complex_fract16 * restrict x, int bexp );
int rifft128   ( int16_t * restrict y, complex_fract16 * restrict x, int bexp );
int rifft256   ( int16_t * restrict y, complex_fract16 * restrict x, int bexp );
int rifft512   ( int16_t * restrict y, complex_fract16 * restrict x, int bexp );
int rifft1024  ( int16_t * restrict y, complex_fract16 * restrict x, int bexp );
int rifft2048  ( int16_t * restrict y, complex_fract16 * restrict x, int bexp );
int rifft4096  ( int16_t * restrict y, complex_fract16 * restrict x, int bexp );
int rifft8192  ( int16_t * restrict y, complex_fract16 * restrict x, int bexp );
int rifft16384 ( int16_t * restrict y, complex_fract16 * restrict x, int bexp );
int rifft32768 ( int16_t * restrict y, complex_fract16 * restrict x, int bexp );

/*-------------------------------------------------------------------------
Radix-2 floating point forward FFT on complex data

Description: These functions make forward FFT on complex data of power of 2
sizes: N=2^n, n=4..15. 

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[N]        Complex input signal
  Output:            
    y[N]        Output spectrum samples
  Returned value:
                zero
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int cfft16f         ( complex_float * restrict y, complex_float * restrict x );
int cfft32f         ( complex_float * restrict y, complex_float * restrict x );
int cfft64f         ( complex_float * restrict y, complex_float * restrict x );
int cfft128f        ( complex_float * restrict y, complex_float * restrict x );
int cfft256f        ( complex_float * restrict y, complex_float * restrict x );
int cfft512f        ( complex_float * restrict y, complex_float * restrict x );
int cfft1024f       ( complex_float * restrict y, complex_float * restrict x );
int cfft2048f       ( complex_float * restrict y, complex_float * restrict x );
int cfft4096f       ( complex_float * restrict y, complex_float * restrict x );
int cfft8192f       ( complex_float * restrict y, complex_float * restrict x );
int cfft16384f      ( complex_float * restrict y, complex_float * restrict x );
int cfft32768f      ( complex_float * restrict y, complex_float * restrict x );

/*-------------------------------------------------------------------------
Radix-2 floating point inverse FFT on complex data

Description: These functions make inverse FFT on complex data of power of 2
sizes: N=2^n, n=4..15. 

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[N]        Input spectrum samples
  Output:            
    y[N]        Complex output signal
  Returned value:
                zero
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int cifft16f         ( complex_float * restrict y, complex_float * restrict x);
int cifft32f         ( complex_float * restrict y, complex_float * restrict x);
int cifft64f         ( complex_float * restrict y, complex_float * restrict x);
int cifft128f        ( complex_float * restrict y, complex_float * restrict x);
int cifft256f        ( complex_float * restrict y, complex_float * restrict x);
int cifft512f        ( complex_float * restrict y, complex_float * restrict x);
int cifft1024f       ( complex_float * restrict y, complex_float * restrict x);
int cifft2048f       ( complex_float * restrict y, complex_float * restrict x);
int cifft4096f       ( complex_float * restrict y, complex_float * restrict x);
int cifft8192f       ( complex_float * restrict y, complex_float * restrict x);
int cifft16384f      ( complex_float * restrict y, complex_float * restrict x);
int cifft32768f      ( complex_float * restrict y, complex_float * restrict x);

/*-------------------------------------------------------------------------
Radix-2 floating point forward FFT on complex data with reduced twiddle table

Description: These functions make forward FFT on complex data of power of 2
sizes: N=2^n, n=13..15. As opposed to regular FFT routines, these use
smaller twiddle factor tables but work a bit slower.

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[N  ]      Complex input signal
  Output:            
    y[N]      Output spectrum samples
  Returned value:
                zero
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int tfft8192f  ( complex_float * restrict y, complex_float * restrict x );
int tfft16384f ( complex_float * restrict y, complex_float * restrict x );
int tfft32768f ( complex_float * restrict y, complex_float * restrict x );

/*-------------------------------------------------------------------------
Radix-2 floating point inverse FFT on complex data with reduced twiddle table

Description: These functions make inverse FFT on complex data of power of 2
sizes: N=2^n, n=13..15. As opposed to regular FFT routines, these use
smaller twiddle factor tables but work a bit slower.

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[N  ]      Input spectrum samples
  Output:            
    y[N]        Complex output signal
  Returned value:
                zero
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int tifft8192f  ( complex_float * restrict y, complex_float * restrict x );
int tifft16384f ( complex_float * restrict y, complex_float * restrict x );
int tifft32768f ( complex_float * restrict y, complex_float * restrict x );

/*-------------------------------------------------------------------------
Blockwise radix-2 floating point forward FFT on complex data

Description: These functions make forward FFT on L blocks, each of N=2^n
complex samples, where n=4..7. 

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[L][N]   Complex input signal
  Output:          
    y[L][N]   Output spectrum samples
  Returned value:
                None
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
void bcfft16f ( complex_float * restrict y, complex_float * restrict x, int L );
void bcfft32f ( complex_float * restrict y, complex_float * restrict x, int L );
void bcfft64f ( complex_float * restrict y, complex_float * restrict x, int L );
void bcfft128f( complex_float * restrict y, complex_float * restrict x, int L );

/*-------------------------------------------------------------------------
Blockwise radix-2 floating point inverse FFT on complex data

Description: These functions make inverse FFT on L blocks, each of N=2^n
complex samples, where n=4..7

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[L][N]   Complex input spectrum
  Output:          
    y[L][N]   Output complex data
  Returned value:
                None
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
void bcifft16f ( complex_float * restrict y, complex_float * restrict x, int L );
void bcifft32f ( complex_float * restrict y, complex_float * restrict x, int L );
void bcifft64f ( complex_float * restrict y, complex_float * restrict x, int L );
void bcifft128f( complex_float * restrict y, complex_float * restrict x, int L );

/*-------------------------------------------------------------------------
Blockwise radix-2 floating point forward FFT on real data

Description: These functions make forward real FFT on L blocks, each of N=2^n
complex samples, where n=4..7. 
FFT implementation for real signal exploits the symmetry properties of the 
Fourier Transform: first, a blockwise complex FFT of half the original size 
is applied to input data, and then the resulting spectrum undergoes a 
conversion procedure which results in complex spectrum of real input data.
NOTES:
1. Bit-reversing permutation is done here. 
2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
   the call. 
3. As the complex spectrum of a real FFT is conjugate symmetric about the 
   midpoint, the brfftf functions generate only the first (N/2)+1 points 
   of the FFT, with the first point being the DC component, and the last 
   point – the Nyquist frequency component. For real signal these two 
   spectral components have zero imaginary part, thus they can be packed 
   in a single 'complex' value. Finally, the the m-th row of the output 
   matrix will contain the following values:
    - DC component of the signal in y[m][0].re
    - Nyquist frequency component in y[m][0].im
    - first half of the complex spectrum in y[m][1]...y[m][(N/2)-1]
   These data are used to reconstruct N values that are stored into the 
   m-th row of real output matrix y[m][N]. The reconstruction is 
   accomplished through a special conversion of input spectrum, such 
   that subsequent invocation of blockwise complex inverse FFT of size 
   N/2 actually produces the desired real signal.

Representation: floating point


Parameters:
  Input:            
    x[L][N]   Real input signal
  Output:          
    y[L][N/2] Half of output spectrum, see note above
  Returned value:
                None
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
void brfft16f ( complex_float * restrict y, float32_t * restrict x, int L );
void brfft32f ( complex_float * restrict y, float32_t * restrict x, int L );
void brfft64f ( complex_float * restrict y, float32_t * restrict x, int L );
void brfft128f( complex_float * restrict y, float32_t * restrict x, int L );

/*-------------------------------------------------------------------------
Blockwise inverse floating point FFT forming real data

These functions make inverse FFT forming real data on L blocks.
FFT implementation for real signal exploits the symmetry properties of the 
Fourier Transform: first, a blockwise complex FFT of half the original size
is applied to input data, and then the resulting spectrum undergoes a 
conversion procedure which results in complex spectrum of real input data.
NOTES:
1.  Bit-reversing permutation is done here.
2.  IFFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after
    the call.
3.  Spectrum of real signal is conjugate-symmetric about the midpoint, thus 
    real inverse FFT functions require only the first half of spectrum and 
    Nyquist frequency bin (altogether N/2+1 complex points) to reconstruct 
    the time domain signal. Moreover, the DC component and Nyquist 
    frequency component have zero imaginary part, thus they can be packed 
    in a single 'complex' value in input data. So, the m-th row of the input 
    matrix shall contain the following values:
    - DC component of the signal in x[m][0].re
    - Nyquist frequency component in x[m][0].im
    - first half of the complex spectrum in x[m][1]...x[m][(N/2)-1]
These data are used to reconstruct N values that are stored into the m-th 
row of real output matrix y[m][N]. The reconstruction is accomplished 
through a special conversion of input spectrum, such that subsequent 
invocation of blockwise complex inverse FFT of size N/2 actually produces 
the desired real signal.

Representation: floating point

Parameters:
  Input:            
    x[L][N/2] half of input spectrum, see note above
  Output:          
    y[L][N]   Output real samples
  Returned value:
                None
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
void brifft16f ( float32_t * restrict y, complex_float * restrict x, int L );
void brifft32f ( float32_t * restrict y, complex_float * restrict x, int L );
void brifft64f ( float32_t * restrict y, complex_float * restrict x, int L );
void brifft128f( float32_t * restrict y, complex_float * restrict x, int L );

/*-------------------------------------------------------------------------
Radix-2 forward floating point FFT on real data

Description: These functions make FFT on real data of length N=2^n, n=4..15.
The algorithm exploits the symmetry properties of the FFT: first, a complex
FFT of half the original size is applied to input data, then the resulting
spectrum undergoes a postprocessing procedure which results in complex spectrum
of real input data.

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:      
    x[N]          Real input signal
  Output:      
    y[(N/2+1)]    Output spectrum samples. 
  Returned value:
                  zero
Restrictions:
  x,y             Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int rfft16f    ( complex_float * restrict y, float32_t * restrict x);
int rfft32f    ( complex_float * restrict y, float32_t * restrict x);
int rfft64f    ( complex_float * restrict y, float32_t * restrict x);
int rfft128f   ( complex_float * restrict y, float32_t * restrict x);
int rfft256f   ( complex_float * restrict y, float32_t * restrict x);
int rfft512f   ( complex_float * restrict y, float32_t * restrict x);
int rfft1024f  ( complex_float * restrict y, float32_t * restrict x);
int rfft2048f  ( complex_float * restrict y, float32_t * restrict x);
int rfft4096f  ( complex_float * restrict y, float32_t * restrict x);
int rfft8192f  ( complex_float * restrict y, float32_t * restrict x);
int rfft16384f ( complex_float * restrict y, float32_t * restrict x);
int rfft32768f ( complex_float * restrict y, float32_t * restrict x);

/*-------------------------------------------------------------------------
Radix-2 inverse floating point FFT forming real data

Description: These functions make inverse FFT forming real data of length
N=2^n, n=4..15. Algorithm exploits the symmetry properties of the FFT:
the input spectrum is modified in such a way that a complex-valued inverse FFT
of half the original size that is applied to the transformed spectrum actually
results in real data.

Representation: floating point

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:
    x[(N/2+1)]    Input spectrum samples
  Output:
    y[N]          Real output signal
  Returned value:
                  zero
Restrictions:
  x,y             Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
int rifft16f    ( float32_t * restrict y, complex_float * restrict x );
int rifft32f    ( float32_t * restrict y, complex_float * restrict x );
int rifft64f    ( float32_t * restrict y, complex_float * restrict x );
int rifft128f   ( float32_t * restrict y, complex_float * restrict x );
int rifft256f   ( float32_t * restrict y, complex_float * restrict x );
int rifft512f   ( float32_t * restrict y, complex_float * restrict x );
int rifft1024f  ( float32_t * restrict y, complex_float * restrict x );
int rifft2048f  ( float32_t * restrict y, complex_float * restrict x );
int rifft4096f  ( float32_t * restrict y, complex_float * restrict x );
int rifft8192f  ( float32_t * restrict y, complex_float * restrict x );
int rifft16384f ( float32_t * restrict y, complex_float * restrict x );
int rifft32768f ( float32_t * restrict y, complex_float * restrict x );
#ifdef __cplusplus
};
#endif

#endif /* __NATUREDSP_BASEBAND_FFT_H */
