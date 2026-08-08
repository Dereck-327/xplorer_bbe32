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
 * Complex Math Functions
 */

#ifndef __NATUREDSP_BASEBAND_COMPLEX_H
#define __NATUREDSP_BASEBAND_COMPLEX_H

#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
Complex Math functions
vdividecx,cdiv      Complex division (vector, scalar)
vpolar,spolar       Polar to Cartesian Conversion (vector, scalar)
vcabs,scabs         Magnitude of complex number (vector, scalar)
varg,sarg           Phase angle of complex number (vector, scalar)
vconj,sconj         Conjugate of complex number (vector, scalar)
cbexp,cexp          Complex exponential (vector, scalar)
cbexplin,cexplin    Complex exponential with linearly evolving phase (vector, scalar)
vslope              Correction of vector slope

For applications where computing speed is critical, the library defines a set 
of 'fast' vectorized floating point functions which use simpler 
approximations, but still keep reasonable level of accuracy. These 'fast' 
versions differ from normal variants as following:
  - accuracy is worse by 1 ULP comparing with non-fast versions
  - domain is narrower (for example, for trigonometric functions it is 
    +/-256pi, not +/-32768pi)
  - functions do not handle inputs which are outside of domain, do not define 
    specific return values for such cases (i.e. for negative inputs for 
    logarithm) or just may return results with worse accuracy
  - usually, these functions do not handle infinities
  - functions do not make special handling of signaling NaNs; all NaNs are 
    treated as inputs outside the domain and return value is not defined
  - functions do not touch errno and may generate some exceptions if input is 
    outside domain
===========================================================================*/

/*-------------------------------------------------------------------------
Division of Complex Numbers

Representation:
vdividecx         16-bit signed fixed-point format (Q15)
vdividecxf,cdivf  IEEE-754 Std. single precision floating-point format

Fixed-point function returns the fractional and exponential portion of real
and imaginary components of the division result. Fixed-point format for the
fractional part is 16-bit Q(15-exp), where exp denotes the respective 
exponential value. Full division result can be restored in 48-bit Q31 format
by sign extending the fractional part to 64 bits and shifting it to the left
by 16+exp bit positions

Special cases:
      x   |    y    |  Result |  Extra Conditions    
  --------|---------|---------|---------------------
    +/-0  |  +/-0   |   NaN   |
     x    | +/-inf  |    0    | x is a finite number          (floating-point functions)
   +/-inf |   y     | +/-inf  | y is a finite number (see *)
   +/-inf | +/-inf  |   NaN   |
  --------|---------|---------|---------------------
     0    |   0     |   not   |                               (fixed-point functions)
          |         | defined |

* - quotent of infinite number to finite number may have NaN in its real or 
    imaginary part depending on signs of input arguments

Input domain for 'fast' version vfastdividecxf is limited by usage of direct 
formula without normalization of inputs, so input domain is:
|real(x)*real(y)|<Inf
|imag(x)*imag(y)|<Inf
|real(y)*real(y)|<Inf
|imag(y)*imag(y)|<Inf
1.1755e-038 < |real(x)*real(y) + imag(x)*imag(y)| < Inf
1.1755e-038 < |real(y)*real(y) + imag(y)*imag(y)| < Inf
The output value is not defined outside of this range or accuracy is degraded.

Accuracy:
vdividecx         1.5 LSB of the fractional part
vdividecxf,cdivf  2 ULP
vfastdividecxf    2 ULP

Input:
x[N]      Input vector of dividends
y[N]      Input vector of divisors
N         Length of vectors
Output:
vdividecx:
fract[N]  Fractional part of quotients, Q(15-exp); if non-zero, then
            8192<=|fract|<32768
exp[2*N]  Exponential part of quotients, -14..16. Exponential values for
          real and imaginary components are interleaved, with exponential
          parts of real components going first (at even indices)
vdividecxf,cdivf
z[N]      Quotients
Restrictions:
x,y,z,fract,exp   Aligned on 32-byte boundary
x,y,z,fract,exp   Must not overlap
N                 Multiple of 8 (vdividecx) or 4 (vdividecxf, vfastdividecxf)
---------------------------------------------------------------------------*/

void vdividecx ( complex_fract16 * restrict fract, 
                 int16_t         * restrict exp, 
           const complex_fract16 * restrict x, 
           const complex_fract16 * restrict y,
           int N );
void vdividecxf( complex_float   * restrict z,
           const complex_float   * restrict x,
           const complex_float   * restrict y,
           int N );
void vfastdividecxf( complex_float   * restrict z,
           const complex_float   * restrict x,
           const complex_float   * restrict y,
           int N );

complex_float cdivf( complex_float x, complex_float y );


/*-------------------------------------------------------------------------
Polar to Cartesian Conversion 

Description: These functions convert input pairs of magnitude and phase
values into Cartesian coordiantes of respective points on the complex plane. 

Representation:
vpolar           16-bit signed fixed-point format
                 Input phase data phase[] are Q15 angular values
                 normalized by pi.
                 Input magnitude data x[] are of abn arbitrary fixed-
                 point format Qx.
                 Fixed-point format for resulting coordinates on the 
                 complex plane z[] is Qx-sh.
vpolarf,spolarf  IEEE-754 Std. single precision floating-point format
                 Input phase data phase[] are in radians.

Accuracy:
9.2e-5 (3 in Q15) - for vpolar
2 ULP - vpolarf, spolarf
3 ULP - vfastpolarf

Note for vpolarf,spolarf:
For the floating-point functions, input phase data should belong to the
range [-102940.0, 102940.0], otherwise the conversion result is (0,0).

Input domain for 'fast' version vfastpolarf()
|a|<804.2477
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
r[N]       Magnitude data
a[N]       Phase data
sh         Right bit shift amount, 0..15 (vpolar)
N          Length of vectors
Output:
z[N]       Points on the complex plane

Restrictions:
z,x,phase  Aligned on 32-byte boundary
z,x,phase  Must not overlap
N          Multiple of 16 (vpolar) or 8 (vpolarf,vfastpolarf)
-------------------------------------------------------------------------*/
void vpolar  ( complex_fract16 * restrict z, 
               const int16_t   * restrict r, 
               const int16_t   * restrict a, 
               int sh, int N );
void vpolarf ( complex_float   * restrict z, 
               const float32_t * restrict r, 
               const float32_t * restrict a, 
               int N );
void vfastpolarf ( complex_float   * restrict z, 
               const float32_t * restrict r, 
               const float32_t * restrict a, 
               int N );
complex_float spolarf( float32_t r, float32_t a );

/*-------------------------------------------------------------------------
Cartesian To Polar Conversion

Description: These functions convert Cartesian coordinates of points on
the complex plane to the polar system (magnitude and phase).

Representation: IEEE-754 Std. single precision floating-point format

Input domain for 'fast' version vfastcartesianf():
|real(x)|>1.1755e-038
|imag(x)|>1.1755e-038
1.1755e-038 < |real(x)*real(x)+ imag(x)*imag(x)| < Inf 
The output value is not defined outside of this range or accuracy is degraded.

Accuracy: 
2 ULP for vcartesianf(),scartesianf()
3 ULP for vfastcartesianf()

Parameters:
Input:
x[N]       Cartesian coordinates
N          Length of vectors
Output:
z[N]       Magnitude data
phase[N]   Phase data

Restrictions:
z,phase,x  Aligned on 32-byte boundary
z,phase,x  Must not overlap
N          Multiple of 8
-------------------------------------------------------------------------*/
void vcartesianf( float32_t     * restrict z,
                  float32_t     * restrict phase,
            const complex_float * restrict x,
            int N );
void vfastcartesianf(float32_t     * restrict z,
                  float32_t     * restrict phase,
            const complex_float * restrict x,
            int N );
float32_t scartesianf( complex_float x, float32_t * phase );

/*-------------------------------------------------------------------------
Magnitude Of Complex Number

Description: These functions multiply complex input number by its conjugate
and take the square root of the result.

Representation:
vcabs          16-bit signed fixed-point format
               Number of fractional bits for input data Qx is user-defined.
               Fixed-point format for output data is Qz = Qx+sh/2, where sh
               is the shift control argument. If a resulting value is too
               large to be represented in Qz format, then it is saturated
               by 32767.
vcabsf,scabsf  IEEE-754 Std. single precision floating-point format


Input domain for 'fast' version vfastcabsf():
1.1755e-038 < |real(x)*real(x)+ imag(x)*imag(x)| < Inf 
The output value is not defined outside of this range or accuracy is degraded.

Accuracy:
1 LSB - vcabs
2 ULP - vcabsf,scabsf
3 ULP - vfastcabsf

Parameters:
Input:
x[N]  Complex numbers
sh    (vcabs only) bi-directional shift control, an even integer from
      the range [-2,30]. Positive value corresponds to a left shift
N     Length of vectors
Output:
z[N]  Magnitudes

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 16 (vcabs) or 8 (vcabsf,vfastcabsf)
-------------------------------------------------------------------------*/
void vcabs  ( int16_t   * restrict z, const complex_fract16 * restrict x, int sh, int N );
void vcabsf ( float32_t * restrict z, const complex_float   * restrict x, int N );
void vfastcabsf ( float32_t * restrict z, const complex_float   * restrict x, int N );

float32_t scabsf( complex_float x );

/*-------------------------------------------------------------------------
Conjugate Of Complex Number

Description: These functions negate the imaginary part of complex input 
numbers.

Representation: IEEE-754 Std. single precision floating-point format

Accuracy: Exact

Parameters:
Input:
x[N]  Complex numbers
N     Length of vectors
Output:
z[N]  Conjugated input numbers

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 4
-------------------------------------------------------------------------*/
void vconjf( complex_float * restrict z, const complex_float * restrict x, int N );

complex_float sconjf( complex_float x );

/*-------------------------------------------------------------------------
Phase Angle Of Complex Number

Description: These functions compute phase angles of input complex numbers
and return the results in radians: [-pi,pi].

Representation: IEEE-754 Std. single precision floating-point format

Accuracy: 
2 ULP - vargf, sargf
3 ULP - vfastargf

Input domain for 'fast' version vfastargf:
1.1755e-038<|real(x)|<Inf
1.1755e-038<|imag(x)|<Inf
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
x[N]  Complex numbers
N     Length of vectors
Output:
z[N]  Phase angles

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
void vargf    ( float32_t * restrict z, const complex_float * restrict x, int N );
void vfastargf( float32_t * restrict z, const complex_float * restrict x, int N );

float32_t sargf( complex_float x );

/*-------------------------------------------------------------------------
Complex Exponential

Description: These functions compute complex number lying on the unit 
circle from input phase data.

Representation:
cbexp         16-bit signed fixed-point format
              Input data phase[] are Q15 angular values normalized by pi.
              Output data z[] format is Q(15-sh), where sh is the shift
              control argument
cbexpf,sexpf  IEEE-754 Std. single precision floating-point format
              Input phase data phase[] are in radians. See the Note.
              
Accuracy:
3.7e-4 (12 in Q15) for cbexp
2 ULP for cbexpf, sexpf
3 ULP for cbfastexpf

Note for (non-fast versions):
For the floating-point functions, input phase data should belong to the
range [-102940.0, 102940.0], otherwise the respective result is 0+0j.

Input domain for 'fast' version cfastbexpf()
|phase|<804.2477
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
phase[N]  Phase values
sh        (cbexp only) left bit shift amount, [0..15]
N         Length of vectors
Output:
z[N]      Results

Restrictions:
z,phase   Aligned on 32-byte boundary
z,phase   Must not overlap
N         Multiple of 16 (cbexp) or 8 (cbexpf,cbfastexpf)
-------------------------------------------------------------------------*/
void cbexp      ( complex_fract16 * restrict z, const int16_t   * restrict phase, int sh, int N );
void cbexpf     ( complex_float   * restrict z, const float32_t * restrict phase, int N );
void cbfastexpf ( complex_float   * restrict z, const float32_t * restrict phase, int N );

complex_float sexpf ( float32_t phase );

/*-------------------------------------------------------------------------
Complex Exponential with Linearly Evolving Phase

Description: These functions compute a series of complex numbers lying on
the unit circle, with the phase angle being incremented for each successive
complex number.

Representation:
cbexplin,           16-bit signed fixed-point format
cbexplin_fast       Initial phase and phase increment (phase0 and delta) are 
                    Q15 angular values normalized by pi.
                    Output data z[] format is Q15.
cbexplinf,          IEEE-754 Std. single precision floating-point format
cexplinf            Initial phase and phase increment (phase0 and delta) are
                    in radians. See Note 3.

Accuracy:
cbexplin            Absolute error for re/im components does not exceed
                    3.7e-4 (12 in Q15)
cbexplin_fast       Absolute error for re/im components for n-th result is
                    ~(floor(n/8)+1)*3.7e-4, n=0..N-1
cbexplinf,cexplinf  2 ULP 
cbfastexplinf       3 ULP

Notes:
1. cbexplin() function computes each complex value from a full-precision
   linearly evolving angle, thus the function is able to generate long data
   sequences with no risk of phase errors accumulation.
2. cbexplin_fast() function is approximately twice faster but less accurate 
   because it performs successive rotations by small angle on each iteration.
   Avoid using this function for long sequences or tiny phase increments.
3. For the non-fast floating-point functions (cbexplinf, cexplinf) , the 
   phase (initial and updated) should belong to the range [-102940.0, 
   102940.0], otherwise the respective resultis 0+0j.
4. Scalar function cexplinf uses the initial phase angle (phase0) to compute 
   a single complex number, and returns the once incremented phase value. 

Input domain for 'fast' version (cbfastexplinf)
|phase|<804.2477
for all initial and updated phases 
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
phase0  Initial phase angle (phase of the first complex number in a series)
delta   Phase increment step
N       Length of output vector
Output:
z[N]    Complex numbers
Returned Value:
        Updated phase, phase0+N*delta

Restrictions:
z       Aligned on 32-byte boundary
z       Must not overlap
N       Multiple of 16 (cbexplln, cbexplln_fast) or 8 (cbexplinf,cbfastexplinf)
-------------------------------------------------------------------------*/
int16_t   cbexplin      ( complex_fract16 * restrict z,  int16_t   phase0, int16_t   delta, int N );
int16_t   cbexplin_fast ( complex_fract16 * restrict z,  int16_t   phase0, int16_t   delta, int N );
float32_t cbexplinf     ( complex_float   * restrict z,  float32_t phase0, float32_t delta, int N );
float32_t cbfastexplinf ( complex_float   * restrict z,  float32_t phase0, float32_t delta, int N );

float32_t cexplinf ( complex_float * restrict z, float32_t phase0, float32_t delta );

/*-------------------------------------------------------------------------
Correction of Vector Slope

Description: Function adjusts the phase of complex values from input vector
by a sequence of linearly evolving angle values.

Representation: 16-bit signed fixed-point

Accuracy: absolute error for re/im component does not exceed 4.0e-4 (13 in Q15).

Note:
Phase adjustment of k-th complex input value is performed through the 
multiplication by a conjugated rotating factor lying on the unit circle. This
factor is calculated from the full-precision phase angle: phase0+k*delta, thus
the function is able to process long data sequences with no risk of phase 
errors accumulation.

Parameters:
Input:
phase0   Initial phase divided by pi, Q15
delta    Phase increment divided by pi, Q15
N        Size of input/output array
Input/Output:
x[N]     Input/output complex values, Q15; real and imaginary parts are 
         Interleaved with the real part stored at even indices
Returned Value:
         Updated phase, phase0+N*delta

Restrictions:
x        Aligned on 32-byte boundary
N        Multiple of 16
-------------------------------------------------------------------------*/
int16_t vslope ( complex_fract16 * restrict x, int16_t phase0, int16_t delta, int N );

#ifdef __cplusplus
};
#endif

#endif /* __NATUREDSP_BASEBAND_COMPLEX_H */
