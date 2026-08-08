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
 * Math Functions
 */

#ifndef __NATUREDSP_BASEBAND_MATH_H
#define __NATUREDSP_BASEBAND_MATH_H

#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
Math functions
vrecip,srecip    Reciprocal (vector, scalar)
vrsqrt,srsqrt    Reciprocal square root (vector, scalar)
vsqrt,ssqrt      Square root (vector, scalar)
vfastrecip16     Fast Vector reciprocal
vfastrsqrt       Fast vector rsqrt
vdivide,sdivide  Division (vector, scalar)
vfmodf,sfmodf    Modulus (vector, scalar)
vlogn,slogn      Natural logarithm (vector, scalar)
vlog2,slog2      Base-2 logarithm (vector, scalar)
vlog10,slog10    Base 10 logarithm (vector, scalar)
valogn,salogn    Natural antilogarithm
valog10,salog10  Base-10 antilogarithm
vldexpf,sldexpf  Modify the Exponent of a Floating-Point Number
vpowf,spowf      Raise To a Power
vsine,ssine      Sine (vector, scalar)
vcos,scos        Cosine (vector, scalar)
vtan,stan        Tangent (vector, scalar)
vcot,scot        Cotangent (vector, scalar)
vasin,sasin      Arcsine (vector, scalar)
vatan2,satan2    Full arctangent (vector, scalar)
vatan16,satan16  Arctangent (vector, scalar)
vtanh,stanh      Hyperbolic tangent (vector, scalar)
vsinh,ssinh      Hyperbolic sine (vector, scalar)
vcosh,scosh      Hyperbolic cosine (vector, scalar)
vint2float,
sint2float       Integer to floating value conversion (vector, scalar)
vfloat2int,
sfloat2int       Floating To Integer Value Conversion (vector, scalar)
fix2hf           Fixed to Floating Point Conversion (vector, scalar)
hf2fix           Floating to Fixed Point Conversion (vector, scalar)
vfloorf,sfloorf  Floating-point floor (vector, scalar)
vceilf,sceilf    Floating-point ceil (vector, scalar)
vabsf,sabsf      Absolute value (vector, scalar)
vcopysign,
scopysign        Copy sign (vector, scalar)
vcountones16
vcountones32 
countones16 
countones32      Count one bits in a word (vector, scalar)
vclip,sclip      Clipping (vector, scalar)
vavg,savg        Average of two arguments (vector, scalar)


Vector mathematic functions include vectorized mathematics, trigonometric and 
transcendental functions. Typically, they require non-overlapping arguments. 
Most functions require aligned inputs/outputs and number of data which is a 
multiple of 16 real 16-bit data or 8 for real floating point data.
A number of library functions supersede standard floating-point mathematical 
functions defined in <math.h>, as listed below:
acosf(), asinf(), atanf(), atan2f(), cosf(), sinf(), tanf(), cotf(), coshf(), 
sinhf(), tanhf(), floorf(), ceilf(), fmodf(), ldexpf(), logf(), log10f(), 
expf(), alog10f(), powf(),copysignf(), sqrtf()
All these functions conform to ISO/IEC 9899 standard (commonly referred to 
as C99) in respect to function semantics, parameters and return value 
specification, however have slightly different name, so might be used 
together with functions from standard C-library supported by the core. These 
floating-point mathematical functions and their vectorized counterparts 
handle error conditions in special Aforementioned functions follow the next 
basic rules:
  - Each function executes as if it were a single operation, and may 
    generate any of 'invalid', 'overflow' or 'divide-by-zero' floating-
    point exceptions only to reflect the result of that operation.
  - A domain error occurs if input argument(s) fall out of the function 
    domain as defined in function specification. In such a case, the 
    function assigns EDOM to the integer expression errno , raises the 
    'invalid' floating-point exception, and returns a quiet NaN.
  - NaN as an input argument is a special kind of domain error. Namely, the 
    integer expression errno acquires EDOM and returned value is a quiet NaN, 
    but the function raises the 'invalid' floating-point exception only if 
    the input argument is a signaling NaN.
  - A pole error occurs if the mathematical function has an exact infinite 
    result as the finite input argument(s) are approached in the limit (for 
    example, cotf(0.0f)). The specification of each function lists all 
    potential pole errors. On a pole error, the function assigns ERANGE to 
    the integer expression errno, raises the 'divide-by-zero' floating-point 
    exception and returns the properly signed infinity value.
  - A floating-point result overflows if the magnitude of the mathematical 
    result is finite but so large that the target floating-point type cannot 
    represent the mathematical result without extraordinary round-off error 
    (for example, expf(100.0f)). If a function detects a floating-point 
    result overflow, it assigns ERANGE to the integer expression errno, 
    raises the 'overflow' floating-point exception and returns the properly 
    signed infinity value.
Vectorized functions obey the same rules of error condition treatment, with 
the following clarifications:
  - When a vectorized floating-point function is applied to data vector(s), 
    it raises the same set of floating-point exceptions as if the 
    corresponding scalar function was sequentially invoked for each slice of 
    input data that comprise the data vector(s).
  - ERANGE and EDOM values are mutually exclusive in respect to evaluation 
    result of the integer expression errno. By design, the vectorized 
    functions give the EDOM value precedence over the ERANGE value. That is, 
    if a vectorized function detects both a domain error and a pole error 
    (or a result overflow) for distinct slices of input data vector(s), then 
    the expression errno evaluates to EDOM after the function returns. This 
    behavior does not depend on the order the errors actually occurred.
Accuracy for each floating point function is defined specifically for each 
one but, as general rule, all of them provide accuracy not worse than 2 ULP in 
specified domain.
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
Reciprocal

Description: Evaluate the reciprocal of input value x and store result
to y: y = 1/x.

Representation:
vrecip16,srecip16  16-bit signed fixed-point format
vrecipf            IEEE-754 Std. single precision floating-point format

Fixed-point routines compute reciprocals for Q15 input data, and return the
fractional and exponential parts of the result. Since the reciprocal of 
a 16-bit Q15 is at least 1.0 in magnitude, functions return fractional part
frac in Q(15-exp) format, where exp is the exponential part of the respective
result. Full result can be restored in 32-bit Q15 format by sign extending
the fractional part to 32 bits and shifting it to the left by exp bit positions.
Scalar fixed-point function returns packed 32-bit result, where exponential part
resides in 16 MSBs and fractional part is located in 16 LSBs.

Special cases:
      x    |  Result |  Extra Conditions    
  ---------|---------|---------------------
   +/-Inf  | +/-0    | vrecipf,srecipf
   +/-0    | +/-Inf  |  
  ---------|---------|---------------------
     0     |   not   | vrecip16,srecip16
           | defined |

Input domain for vfastrecipf:
|x|>2.94e-39, |x|<Inf
The output value is not defined outside of this range.

Accuracy:
vrecip16,srecip16  1 LSB of the fractional part
vrecipf,srecipf    1 ULP 
vfastrecipf        2 ULP 

Parameters:
Input:
x[N]           Input data vector
N              Length of vectors
Output:
vrecip16,srecip16
frac[N]        Fractional part of reciprocals, Q(15-exp)
exp[N]         Exponent of reciprocals (1...16)
vrecip:
y[N]           Reciprocals

Restrictions:
y,x,fract,exp  Aligned on 32-byte boundary
y,x,fract,exp  Must not overlap
N              Multiple of 16 (vrecip16), 8 (vrecipf,vfastrecipf)
-------------------------------------------------------------------------*/
void vrecip16 ( int16_t   * restrict fract, 
                int16_t   * restrict exp, 
          const int16_t   * restrict x, 
          int N);
void vrecipf  ( float32_t * restrict y,
          const float32_t * restrict x,
          int N );
void vfastrecipf  
              ( float32_t * restrict y, 
          const float32_t * restrict x, 
          int N );

uint32_t  srecip16 ( int16_t   x );
float32_t srecipf  ( float32_t x );

/*-------------------------------------------------------------------------
Reciprocal Square Root

Description: Evaluate the reciprocal square root of input value x and store
result to y: y = 1/x^0.5.

Representation:
vrsqrt,srsqrt  32-bit signed fixed-point format
               Number of fractional bits for input data Qx is user-
               defined, provided that it is even.
               Fixed-point format Qy for output data is Qy = 31-Qx/2.
               For example, if Qx == 30 then Qy = 31-30/2 = 16
vrsrtf         IEEE-754 Std. single precision floating-point format

Accuracy:
vrsqrt,srsqrt  1.1e-4 (worst case relative error)
vrsrtf,srsqrtf 2 ULP
vfastrsrtf     3 ULP

Notes (for non-fast versions)::
1. Fixed-point functions return MIN_INT32 (0x80000000) for a negative or
   zero input value.
2. Floating-point reciprocal square root conforms to IEEE-754 Std rSqrt 
   operation in respect to signaling error conditions by means of floating-
   point exceptions.
3. Floating-point reciprocal square root limits the range of allowable
   input values, as follows:
   A) If x<0, functions raise the "invalid" floating-point exception,
      assign EDOM to errno and set output value y to NaN.
   B) If x==+/-0, functions set output value y to +/-HUGE_VALF, raise the
      "divide by zero" floating-point exception, and assign ERANGE to errno.

Input domain for vfastrsqrtf():
x>=+0 && x<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vrsqrt) or 8 (vrsqrtf,vfastrsqrtf)
-------------------------------------------------------------------------*/
void vrsqrt  ( int32_t   * restrict y,
         const int32_t   * restrict x,
         int N );
void vrsqrtf ( float32_t * restrict y,
         const float32_t * restrict x,
         int N );
void vfastrsqrtf ( float32_t * restrict y,
         const float32_t * restrict x,
         int N );


int32_t   srsqrt  ( int32_t   x );
float32_t srsqrtf ( float32_t x );

/*-------------------------------------------------------------------------
Square Root

Description: These functions compute the positive square root of input
vector elements.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
vsqrtf, ssqrtf  2 ULP 
vfastsqrtf      3 ULP

Notes:
1. Square root functions conform to ANSI C requirements on standard math library
   functions in respect to treatment of errno and floating-point exceptions.
2. For a negative input value, functions raise the "invalid" floating-point
   exception, assign EDOM to errno and set the respective output value z to NaN.
3. Negative zero (-0) is considered as a valid input, the result is also -0.

Input domain for vfastrsqrtf():
x>=+0 && x<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
void vsqrtf ( float32_t * restrict y,
        const float32_t * restrict x,
        int N );
void vfastsqrtf ( float32_t * restrict y,
        const float32_t * restrict x,
        int N );
float32_t ssqrtf ( float32_t x );

/*-------------------------------------------------------------------------
Reciprocal for Pseudo-Floating Point Format

Description: calculates y=1/x; input and output values are represented
by 7-bit signed exponent and 16-bit signed mantissa, Q(15+exp).
If the magnitude of result is so large that the result cannot be
represented in the pseudo-floating point format, then the function
returns 32767/-32768,-64 (+/-1.8e19) depending on the sign bit of the
input mantissa.

Relative accuracy: 7.7e-3 (~7 bits of precision) in worst case

Parameters:
Input:
xmant[N]                Mantissa of input values, -32768..32767
xexp[N]                 Exponent of input values, -64..63
Output:
ymant[N]                Mantissa of output values, -32768..32767
yexp[N]                 Exponent of output values, -64..63

Restrictions:
xmant,ymant,xexp,yexp   Aligned on 32-byte boundary
xmant,ymant,xexp,yexp   Must not overlap
N                       Multiple of 16
-------------------------------------------------------------------------*/
void vfastrecip16 ( int16_t * restrict ymant,
                    int16_t * restrict yexp,
              const int16_t * restrict xmant,
              const int16_t * restrict xexp,
              int N );

/*-------------------------------------------------------------------------
Inverse Square Root for Pseudo-Floating Point Format

Description: calculates y=1/sqrt(x); input and output values are represented
by 16-bit signed mantissa and 7-bit signed exponent, Q(15+exp). Returns 
-32768,0 (-1.0) for input values with negative mantissa. For zero mantissa
the result is 32767,-64 (~1.8e19).

Relative accuracy: 2.0e-3 (~9 bits of precision) in worst case

Parameters:
Input:
xmant[N]                Mantissa of input values, -32768..32767
xexp[N]                 Exponent of input values, -64..62
Output:
ymant[N]                Mantissa of output values, -32768,16384..32767
yexp[N]                 Exponent of output values, -64,-40..31

Restrictions:
xmant,ymant,xexp,yexp   Aligned on 32-byte boundary
xmant,ymant,xexp,yexp   Must not overlap
N                       Multiple of 16
-------------------------------------------------------------------------*/
void vfastrsqrt ( int16_t * restrict ymant,
                  int16_t * restrict yexp,
            const int16_t * restrict xmant,
            const int16_t * restrict xexp,
            int N );

/*-------------------------------------------------------------------------
Division

Representation:
vdivide,sdivide    16-bit signed fixed-point format (Q15)
vdividef,sdividef  IEEE-754 Std. single precision floating-point format

Fixed-point function returns the fractional and exponential portion of the 
division result. Fixed-point format for the fractional part is 16-bit 
Q(15-exp), where exp denotes the respective exponential value. Full division 
result can be restored in 48-bit Q31 format by sign extending the fractional 
part to 64 bits and shifting it to the left by 16+exp bit positions
Scalar fixed-point function returns packed 32-bit result, where exponential part
resides in 16 MSBs and fractional part is located in 16 LSBs.

Special cases:
      x   |    y    |  Result |  Extra Conditions    
  --------|---------|---------|---------------------
    +/-0  |  +/-0   |   NaN   |
     x    | +/-inf  |    0    | x is a finite number  (floating-point functions)
   +/-inf |   y     | +/-inf  | y is a finite number 
   +/-inf | +/-inf  |   NaN   |
  --------|---------|---------|---------------------
     0    |   0     |   not   |                       (fixed-point functions)
          |         | defined |

Accuracy:
vdivide,sdivide   1 LSB of the fractional part
vdividef,dividef  1 ULP

Parameters:
Input:
x[N]      Input vector of dividends
y[N]      Input vector of divisors
N         Length of vectors
Output:
vdivide,sdivide
fract[N]  Fractional part of quotients, Q(15-exp); if non-zero, then
            8192<=|fract|<32768
exp[N]    Exponential part of quotients, -14..16
vdividef,sdividef
z[N]      Quotients

Restrictions:
z,x,y,fract,exp   Aligned on 32-byte boundary
z,x,y,fract,exp   Must not overlap
N                 Multiple of 16 (vdivide) or 8 (vdividef)
-------------------------------------------------------------------------*/
void vdivide  ( int16_t * restrict fract, 
                int16_t * restrict exp, 
          const int16_t * restrict x, 
          const int16_t * restrict y,
          int N );
void vdividef ( float32_t * restrict z, 
          const float32_t * restrict x, 
          const float32_t * restrict y,
          int N );

uint32_t  sdivide ( int16_t   x, int16_t   y );
float32_t sdividef( float32_t x, float32_t y );

/*-------------------------------------------------------------------------
Modulus

Description: these functions compute the floating-point remainder that
results from dividing the first argument by the second argument. The result
is less than the second argument and has the same sign as the first argument.

Data format: IEEE-754 Std. single precision floating-point

Special cases:
    x    |   y    | Result | Extra Conditions
  -------+--------+--------+-------------------
   +/-0  |   y    |  +/0   | y!=0
    x    | +/-Inf |   x    | x != +/-Inf
  +/-Inf |   y    |  NaN   | for any y
    x    | +/-0   |  NaN   | for any x
    x    |   y    |  NaN   | |x/y|>=2^24


Accuracy: 2 ULP

Notes:
1. Modulus functions conform to ANSI C requirements on standard math library
   functions in respect to treatment of errno and floating-point exceptions.
2. Modulus functions limit the range of allowable input values, as follows:
   A) If |x/y|>=2^24, then the respective result z is set to NaN
   B) If x==+/-Inf and/or y==+/-0, functions set output value z to NaN, raise
      the "invalid" floating-point exception, and assign EDOM to errno.

Parameters:
Input:
x[N]    Input values
y[N]    Modulus values
N       Length of input/output vectors
Output:
z[N]    Results

Restrictions:
z,x,y   Aligned on 32-byte boundary
z,x,y   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
void vfmodf ( float32_t * restrict z, 
        const float32_t * restrict x,
        const float32_t * restrict y,
        int N );

float32_t sfmodf ( float32_t x, float32_t y );

/*-------------------------------------------------------------------------
Logarithms

Description: These function compute base-2, base-10 or natural logarithm of
input data.

Representation:
vlog2,vlogn,vlog10,     Signed fixed-point format
slog2,slogn,slog10      Input data are 32-bit Q16.15, results are 16-bit Q4.11.
                        Here are a few examples for the base-2 logarithm:

                          Function | Input Data Q16.15 <real> | Result Q4.11 <real>
                        -----------+--------------------------+--------------------
                           slog2   | 65536 <2.0>              | 2048 <1.0>
                           slog2   | 2147483647 <65535.99997> | 32767 <15.9995>
                           slog2   | 1 <3.052e-5>             | -30720 <-15.0>
                        -----------+--------------------------+--------------------
vlog2f,vlognf,vlog10f,  IEEE-754 Std. single precision floating-point format
slog2f,slognf,slog10f

Accuracy:
1 LSB for the fixed-point functions
2 ULP for the floating-point functions

Notes:
1. Fixed-point Functions return -32768 for a negative or zero input.
2. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
3. Floating point functions limit the range of allowable input values:
   A) If x<0, the result is set to NaN, errno is assigned the value EDOM, and
      "invalid" floating-point exception is raised
   B) If x==0, the result is set to minus infinity, errno is assigned the value 
      ERANGE, and "divide-by-zero" floating-point exception is raised

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vlog2,vlogn,vlog10) or 8 (vlog2f,vlognf,vlog10f)
-------------------------------------------------------------------------*/
void vlog2   ( int16_t   * restrict y, const int32_t   * restrict x, int N );
void vlogn   ( int16_t   * restrict y, const int32_t   * restrict x, int N );
void vlog10  ( int16_t   * restrict y, const int32_t   * restrict x, int N );
void vlog2f  ( float32_t * restrict y, const float32_t * restrict x, int N );
void vlognf  ( float32_t * restrict y, const float32_t * restrict x, int N );
void vlog10f ( float32_t * restrict y, const float32_t * restrict x, int N );

int16_t   slog2   ( int32_t   x );
int16_t   slogn   ( int32_t   x );
int16_t   slog10  ( int32_t   x );
float32_t slog2f  ( float32_t x );
float32_t slognf  ( float32_t x );
float32_t slog10f ( float32_t x );

/*-------------------------------------------------------------------------
Antilogarithm and Exponential

Description: These functions compute base-10 or natural antilogarithm of
input data.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 2 ULP.

Note:
Antilogarithm functions conform to ANSI C requirements on standard math 
library functions in respect to treatment of errno and floating-point
exceptions.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
void valog10f ( float32_t * restrict y, const float32_t * restrict x, int N );
void valognf  ( float32_t * restrict y, const float32_t * restrict x, int N );

float32_t salog10f ( float32_t x );
float32_t salognf  ( float32_t x );

/*-------------------------------------------------------------------------
Modify the Exponent of a Floating-Point Number

Description: These functions multiply input values by 2^n, where n is an
exponent adjustment term. If the result overflows, functions return 
HUGE_VALF with the proper sign. If, on the contrary, the result underflows,
functions return zero with proper sign.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: exact

Note:
Exponent modification functions conform to ANSI C requirements on standard
math library functions in respect to treatment of errno and floating-point
exceptions.

Parameters:
Input:
x[N]    Input data
n[N]    Exponent adjustment terms
N       Length of input/output data vectors
Output:
y[N]    Results

Restrictions:
y,x,n   Aligned on 32-byte boundary
y,x,n   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
void vldexpf ( float32_t * restrict y, 
         const float32_t * restrict x, 
         const int32_t   * restrict n, 
         int N );

float32_t sldexpf ( float32_t x, int32_t n );

/*-------------------------------------------------------------------------
Raise To a Power

Description: These functions compute the value of the first argument x 
raised to the power of the second argument y.

Data format: IEEE-754 Std. single precision floating-point.

Special cases:
      x   |   y    | Result |  Extra Conditions    
  --------+--------+--------+---------------------
    +/-0  | y      | +/-inf | odd y<0
    +/-0  | y      | +inf   | even y<0
    +/-0  | y      | +/-0   | odd y>0
    +/-0  | y      | 0      | even y>0
    +/-1  | +/-inf | 1      | 
    1     | y      | 1      | any y including NaN
    x     | +/-0   | 1      | any x including NaN
    x     | y      | NaN    | finite x<0 and finite non-integer y (see note 2)
    x     | -inf   | +inf   | |x|<1
    x     | -inf   | 0      | |x|>1
    x     | +inf   | 0      | |x|<1
    x     | +inf   | +inf   | |x|>1
    -inf  | y      | -0     | y an odd integer <0
    -inf  | y      | 0      | y<0 and not an odd integer
    -inf  | y      | -inf   | y an odd integer >0
    -inf  | y      | +inf   | y>0 and not an odd integer
    +inf  | y      | 0      | y<0
    +inf  | y      | +inf   | y>0

Accuracy: 2 ULP under condition that |y|<=100

Notes:
1. Raise to a power functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-point
   exceptions.
2. If x<0 is finite, y is finite and not an integer value, then the respective
   result z is set to NaN, errno is assigned the value EDOM, and the "invalid"
   floating-point exception is raised.

Parameters:
Input:
x[N]    Input data
y[N]    Power values
N       Length of input/output data vectors
Output:
z[N]    Results

Restrictions:
z,x,y   Aligned on 32-byte boundary
z,x,y   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
void vpowf ( float32_t * restrict z, 
       const float32_t * restrict x, 
       const float32_t * restrict y, 
       int N );

float32_t spowf ( float32_t x, float32_t y );

/*-------------------------------------------------------------------------
Sine/Cosine 

Description: These functions compute sine or cosine of input data

Representation:
vsine,vcos,    16-bit signed fixed-point format Q15 for input/output data
ssine,scos     It is assumed that input angular values are normalized by pi.
               That is, Fixed-point functions actually compute sin(pi*x) or
               cos(pi*x)
vsinef,vcosf,  IEEE-754 Std. single precision floating-point format for
ssinef,scosf   input/output data. Input data are treated as angular values
               specified in radians. Floating-point functions limit the
               rangw of allowable input values, see note 2.

Accuracy:
2 LSB - vsine(),vcos(), ssine(), scos(),
2 ULP - vsinef(), vcosf(), ssinef(), scosf(),
3 ULP - vfastsinef(),vfastcosf()

Notes for non-fast versions:
1. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
2. Floating-point functions require that input value belongs to the 
   closed range [-102940.0,102940.0], otherwise the respective result
   is NaN.

Input domain for 'fast' versions vfastsinef(),vfastcosf()
|x|<804.2477
The output value is not defined outside of this range or accuracy is degraded

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vsine,vcos) or 8 (vsinef,vcosf,vfastsinef,vfastcosf)
-------------------------------------------------------------------------*/
void vsine      ( int16_t   * restrict y, const int16_t   * restrict x, int N );
void vcos       ( int16_t   * restrict y, const int16_t   * restrict x, int N );
void vsinef     ( float32_t * restrict y, const float32_t * restrict x, int N );
void vcosf      ( float32_t * restrict y, const float32_t * restrict x, int N );
void vfastsinef ( float32_t * restrict y, const float32_t * restrict x, int N );
void vfastcosf  ( float32_t * restrict y, const float32_t * restrict x, int N );

int16_t   ssine  ( int16_t   x );
int16_t   scos   ( int16_t   x );
float32_t ssinef ( float32_t x );
float32_t scosf  ( float32_t x );

/*-------------------------------------------------------------------------
Tangent 

Description: These functions compute tangent of input data

Representation:
vtan,stan    Signed fixed-point format
             Input data are 16-bit Q15 angular values normalized by pi,
             i.e. fixed-point functions actually compute tan(pi*x).
             Output data are 32-bit Q16.15 values.
vtanf,stanf  IEEE-754 Std. single precision floating-point format for
             input/output data. Input data are treated as angular values
             specified in radians. Floating-point functions limit the
             rangw of allowable input values, see note 3.

Accuracy:
For the fixed-point functions, accuracy depends on the input value x,
as shown in the table below:
   Range of |x|    | Absolute error  | Relative error
-------------------+-----------------+----------------
 [-pi/4; pi/4]     |     1 (Q15)     |
 [pi/4; 7pi/16]    |    15 (Q15)     |    4.6e-4
 [7pi/16; 31pi/64] |   242 (Q15)     |    1.5e-3
-------------------+-----------------+----------------
2 ULP for vtanf(),stanf()
3 ULP for vfasttanf()

Notes for non-fast versions:
1. Fixed-point function result is not defined if input value x is
   +/-pi/2 (+/-8192 in Q15 normalized by pi).
2. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
3. Floating-point functions require that input value belongs to the 
   closed range [-9099.0,9099.0], otherwise the respective result is NaN.

Input domain for 'fast' version vfasttanf():
|x|<804.2477
The output value is not defined outside of this range or accuracy is 
degraded

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 16 (vtan) or 8 (vtanf,vfasttanf)
-------------------------------------------------------------------------*/
void vtan      ( int32_t   * restrict y, const int16_t   * restrict x, int N );
void vtanf     ( float32_t * restrict y, const float32_t * restrict x, int N );
void vfasttanf ( float32_t * restrict y, const float32_t * restrict x, int N );

int32_t   stan  ( int16_t   x );
float32_t stanf ( float32_t x );

/*-------------------------------------------------------------------------
Cotangent 

Description: These functions compute cotangent of input data.

Data format: IEEE-754 Std. single precision floating-point.
Input data are treated as angular values specified in radians.

Accuracy: 
2 ULP - vcotf(), scotf()
3 ULP - vfastcotf()

Notes for non-fast versions:
1. Cotangent functions conform to ANSI C requirements on standard math
   library functions in respect to treatment of errno and floating-point
   exceptions.
2. Cotangent functions functions require that input value belongs to the 
   closed range [-9099.0,9099.0], otherwise the respective result is NaN.

Input domain for 'fast' version vfastcotf():
|x|<804.2477, x!=0
The output value is not defined outside of this range or accuracy is 
degraded

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
void vcotf     ( float32_t * restrict y, const float32_t * restrict x, int N );
void vfastcotf ( float32_t * restrict y, const float32_t * restrict x, int N );

float32_t scotf ( float32_t x );

/*-------------------------------------------------------------------------
Arcsine

Description: These functions compute arcsine of input data. Functions output
is in radians.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
2 ULP - vasinf(), sasinf()
3 ULP - vfastasinf()

Notes for non-fast versions:
1. Arcsine functions conform to ANSI C requirements on standard math
   library functions in respect to treatment of errno and floating-point
   exceptions.
2. Input values should belong to [-1,1], otherwise functions raise the
   "invalid" floating-point exception, assign EDOM to errno and set the
   respective output value to NaN.

Input domain for 'fast' version vfastasinf():
|x|<=1
The output value is not defined outside of this range

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
void vasinf     ( float32_t * restrict y, const float32_t * restrict x, int N );
void vfastasinf ( float32_t * restrict y, const float32_t * restrict x, int N );

float32_t sasinf ( float32_t x );

/*-------------------------------------------------------------------------
Arccosine

Description: These functions compute arccosine of input data. Functions output
is in radians.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy:
2 ULP - vacosf(), sacosf()
3 ULP - vfastacosf()

Notes:
1. Arccosine functions conform to ANSI C requirements on standard math
   library functions in respect to treatment of errno and floating-point
   exceptions.
3. Input values should belong to [-1,1], otherwise functions raise the
   "invalid" floating-point exception, assign EDOM to errno and set the
   respective output value to NaN.

Input domain for 'fast' version vfastacosf():
|x|<=1
The output value is not defined outside of this range

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
void vacosf     ( float32_t * restrict y, const float32_t * restrict x, int N );
void vfastacosf ( float32_t * restrict y, const float32_t * restrict x, int N );

float32_t sacosf ( float32_t x );

/*-------------------------------------------------------------------------
Full Arctangent 

Description: These functionx calculate four quadrant arctangent of complex 
argument x and return arg(x)/pi in Q15 format. Q15 is also assumed for
argument. If real and imaginary components are zero, the result is also zero.

Absolute phase accuracy: 1 LSB 

Parameters:
Input:
x[2*N]  Complex input data. Real and imaginary parts are interleaved with
        real parts stored at even indices.
N       Length of vectors
Output:
z[N]    Result

Restrictions:
z,x     Aligned on 32-byte boundary
z,x     Must not overlap
N       Multiple of 16
-------------------------------------------------------------------------*/
void vatan2_16 ( int16_t * restrict z, const int16_t * restrict x,  int N );

int16_t satan2_16 ( int16_t re, int16_t im );

/*-------------------------------------------------------------------------
Full Arctangent (Floating-Point)

Description: These functions calculate full-quadrant arctangent of 
ratio y/x and output results in radians.

Data format: IEEE-754 Std. single precision floating-point

Special cases:
     y    |   x   |  Result   |  Extra Conditions    
  --------+-------+-----------+---------------------
   +/-0   | -0    | +/-pi     |
   +/-0   | +0    | +/-0      |
   +/-0   |  x    | +/-pi     | x<0
   +/-0   |  x    | +/-0      | x>0
   y      | +/-0  | -pi/2     | y<0
   y      | +/-0  |  pi/2     | y>0
   +/-y   | -inf  | +/-pi     | finite y>0
   +/-y   | +inf  | +/-0      | finite y>0
   +/-inf | x     | +/-pi/2   | finite x
   +/-inf | -inf  | +/-3*pi/4 | 
   +/-inf | +inf  | +/-pi/4   |

Notes for non-fast versions:
1. Full arctangent functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
2. Full arctangent functions assign EDOM to errno whenever x==0 and y==0.
   The resulting value depends on signs of x and y, see the Special Cases
   above.

Accuracy: 
2 ULP for vatan2f(),satan2f()
3 ULP for vfastatan2f()

Input domain for 'fast' version:
1.1755e-038 < |x| < Inf
1.1755e-038 < |y| < Inf
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
y[N]    Numerator values
x[N]    Denominator values
N       Length of input/output vectors
Output:
z[N]    Results

Restrictions:
z,y,x   Aligned on 32-byte boundary
z,y,x   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
void vatan2f ( float32_t * restrict z,
         const float32_t * restrict y,
         const float32_t * restrict x,
         int N );
void vfastatan2f ( float32_t * restrict z,
         const float32_t * restrict y,
         const float32_t * restrict x,
         int N );

float32_t satan2f ( float32_t y, float32_t x );

/*-------------------------------------------------------------------------
Arctangent 

Description: These functions compute the principal value of arctangent.

Representation:
vatan16,satan16  16-bit signed fixed-point format
                 Input data are Q15. Functions compute atan(x)/(pi/4)
                 and output results in Q15 format.
vatanf,satanf    IEEE-754 Std. single precision floating-point format
vfastatanf       Functions compute atan(x) and output results in radians

Special cases:
    Input | Result 
   -------+--------
    +inf  |  pi/2  (floating-point functions)
    -inf  | -pi/2  

Accuracy:
1 LSB for fixed point fixed point functions
1 ULP for vatanf(), satanf()
2 ULP for vfastatanf()

Notes:
1. These functions are much faster than full-quadrant arctangent atan2,
   so they are preferable when the full phase is not required.
2. Floating-point functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.

Input domain for 'fast' version vfastatanf():
|x|<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
z[N]  Results

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 16 (vatan16) or 8 (vatanf,vfastatanf)
-------------------------------------------------------------------------*/
void vatan16 ( int16_t   * restrict z, 
         const int16_t   * restrict x, 
         int N );
void vatanf  ( float32_t * restrict z, 
         const float32_t * restrict x, 
         int N );
void vfastatanf  ( float32_t * restrict z, 
         const float32_t * restrict x, 
         int N );
int16_t   satan16 ( int16_t   x );
float32_t satanf  ( float32_t x );

/*-------------------------------------------------------------------------
Hyperbolic Tangent

Description: These functions compute hyperbolic tangent of input data

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
2 ULP for vtanhf(), stanhf()
3 ULP for vfasttanhf()

Note for non-fast version
Hyperbolic tangent functions conform to ANSI C requirements on standard
math library functions in respect to treatment of errno and floating-
point exceptions.

Input domain for 'fast' version vfasttanhf():
|x|<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
void vtanhf     ( float32_t * restrict y, const float32_t * restrict x, int N );
void vfasttanhf ( float32_t * restrict y, const float32_t * restrict x, int N );

float32_t stanhf ( float32_t x );

/*-------------------------------------------------------------------------
Hyperbolic Sine

Description: These functions compute hyperbolic sine of input data

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
2 ULP for vsinhf(), ssinhf()
3 ULP for vfastsinhf()

Notes for non-fast versions:
1. Hyperbolic sine functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
2. Due to limited dynamic range of single precision floating-point format,
   hyperbolic sine result for an input value x such that |x|>89.41599 is
   sign(x)*HUGE_VALF.

Input domain for 'fast' version vfastsinhf():
|x|<89.41599
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
void vsinhf     ( float32_t * restrict y, const float32_t * restrict x, int N );
void vfastsinhf ( float32_t * restrict y, const float32_t * restrict x, int N );

float32_t ssinhf ( float32_t x );

/*-------------------------------------------------------------------------
Hyperbolic Cosine

Description: These functions compute hyperbolic cosine of input data

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
2 ULP for vcoshf(), scoshf()
3 ULP for vfastcoshf()

Notes for non-fast versions:
1. Hyperbolic cosine functions conform to ANSI C requirements on standard
   math library functions in respect to treatment of errno and floating-
   point exceptions.
2. Due to limited dynamic range of single precision floating-point format,
   hyperbolic cosine result for an input value x such that |x|>89.41599 is
   HUGE_VALF.

Input domain for 'fast' version vfastcoshf():
|x|<89.41599
The output value is not defined outside of this range.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
void vcoshf     ( float32_t * restrict y, const float32_t * restrict x, int N );
void vfastcoshf ( float32_t * restrict y, const float32_t * restrict x, int N );

float32_t scoshf ( float32_t x );

/*-------------------------------------------------------------------------
Integer To Floating Value Conversion

Description: These functions convert integer input values to floating 
values and scale them by 2^t.

Data format: Signed 32-bit integer on input, 
             IEEE-754 Std. single precision floating-point on output.

Parameters:
Input:
x[N]  Input integers
N     Length of input/output data vectors
t     Scale factor
Output:
y[N]  Conversion results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
t     Must belong to [-126,126]
-------------------------------------------------------------------------*/
void vint2float ( float32_t * restrict y, 
            const int32_t   * restrict x,
            int t, int N );

float32_t sint2float ( int32_t x, int t );

/*-------------------------------------------------------------------------
Floating To Integer Value Conversion

Description: These functions scale input floating input values by 2^-t and
convert them to integers with saturation.

Data format: IEEE-754 Std. single precision floating-point on input,
             signed 32-bit integer on output
Parameters:
Input:
x[N]  Input floating values
N     Length of input/output data vectors
t     Scale factor
Output:
y[N]  Conversion results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
t     Must belong to [-126,126]
-------------------------------------------------------------------------*/
void vfloat2int ( int32_t   * restrict y, 
            const float32_t * restrict x, 
            int t, int N );

int32_t sfloat2int ( float32_t x, int t );

/*-------------------------------------------------------------------------
Fixed to Floating-Point Conversion

Description: Function converts 16-bit signed fixed-point data with 
user-defined point position to 16-bit half-precision floating-point
format (IEEE 754-2008 binary16). 

Notes:
1. IEEE 754-2008 binary16 encoding format provides 1 sign bit, 10 trailing
   significand bits and 5 exponent bits.
2. If absolute input value |x|*2^-q exceeds 65504, the result is saturated to
   +/-Inf (0x7C00 or 0xFC00). In particular, if point position is less than
   or equal to -16 (q<=-16), any nonzero fixed-point value is converted to
   +/-Inf.
3. If absolute input value belongs to the closed range of [2^-24,2^-15], then
   the conversion results in a subnormal half-precision floating-point number
   as defined by IEEE 754-2008.
4. If input value x*2^-q belongs to the open range of (-2^-24,2^-24), then it
   is converted to zero of original sign (0x0000 or 0x8000). Specifically, any
   fixed-point value is converted to +/-0 when the point position is greater
   than or equal to 40 (q>=40).
5. If input value cannot be exactly represented in half-precision floating-
   point format, then it is rounded to the next representable value toward zero.

Parameters:
Input:
x[N]    Input fixed-point data, Q(q)
q       Point position for input data
N       Size of input/output arrays
Output:
y[N]    Output half-precision floating-point values

Restrictions:
y,x     Aligned on 32-byte boundary
y,x     Must not overlap
N       Multiple of 16
-------------------------------------------------------------------------*/
void fix2hf ( int16_t * restrict y, 
        const int16_t * restrict x, 
        int q, int N );

/*-------------------------------------------------------------------------
Floating to Fixed-Point Conversion

Description: Function converts 16-bit half-precision floating-point data
(IEEE 754-2008 binary16 format) to 16-bit signed fixed-point format with
user-defined point position.

Notes:
1. +/-Inf is converted to 32767 or -32768, respectively.
2. NaN values (zero exponent field, non-zero significand field) are treated as
   +/-Inf depending on the sign bit.
3. If a finite input value is too large in its magnitude to be represented
   in the fixed-point format with user-specified point position, then
   it is saturated to -32768 or 32767 depending on the input sign.
4. If the input value cannot be exactly represented in fixed-point format
   with user-specified point position, then it is rounded to the next
   representable value toward -Inf.
5. The sensible range for the point position q is [-15,39]. That is, for q<=-16
   the result of conversion is always zero or -1 (depending on the input sign),
   unless the input value is +/-Inf or NaN. For q>=39, any non-zero floating-
   point value is converted to either -32768 or 32767, depending on the input sign.
6. Given an arbitrary floating-point value, the optimum point position for the
   floating-to-fixed conversion is q = 14-e, where e is the unbiased exponent
   of the floating-point value. For a normal positive floating-point value, the
   optimum point position provides a normalized fixed-point result. For a
   normal negative value on input, the number of redundant sign bits for the
   conversion result is at most 1.

Parameters:
Input:
x[N]    Input half-precision floating-point values 
q       Point position for output data
N       Size of input/output arrays
Output:
y[N]    Output fixed-point data, Q(q)

Restrictions:
y,x     Aligned on 32-byte boundary
y,x     Must not overlap
N       Multiple of 16
-------------------------------------------------------------------------*/
void hf2fix ( int16_t * restrict y, 
        const int16_t * restrict x, 
        int q, int N );

/*-------------------------------------------------------------------------
Floating-Point Floor

Description: These functions return the largest integral value that is not
greater than an input value.

Data format: IEEE-754 Std. single precision floating-point.

Note:
Floor functions conform to ANSI C requirements on standard math library
functions in respect to treatment of errno and floating-point exceptions.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
void vfloorf ( float32_t * restrict y, const float32_t * restrict x, int N );

float32_t sfloorf ( float32_t x );

/*-------------------------------------------------------------------------
Floating-Point Ceil

Description: These functions return the smallest integral value that is not
less than an input value.

Data format: IEEE-754 Std. single precision floating-point.

Note:
Ceil functions conform to ANSI C requirements on standard math library
functions in respect to treatment of errno and floating-point exceptions.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
void vceilf ( float32_t * restrict y, const float32_t * restrict x, int N );

float32_t sceilf ( float32_t x );

/*-------------------------------------------------------------------------
Absolute Value

Description: These functions return absolute value for each input value.

Data format: IEEE-754 Std. single precision floating-point.

Notes:
In terms of IEEE-754 Std, this is a quiet-computational operation which
treats floating-point numbers and NaNs alike, and does not raise any 
floating-point exceptions.

Parameters:
Input:
x[N]  Input data
N     Length of input/output data vectors
Output:
y[N]  Results

Restrictions:
y,x   Aligned on 32-byte boundary
y,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
void vabsf ( float32_t * restrict y, const float32_t * restrict x, int N );

float32_t sabsf ( float32_t x );

/*-------------------------------------------------------------------------
Copy Sign

Description: These functions copy the sign bit of the second argument to the
first argument and return the result.

Data format: IEEE-754 Std. single precision floating-point.

Note:
In terms of IEEE-754 Std, this is a quiet-computational operation which
treats floating-point numbers and NaNs alike, and does not raise any 
floating-point exceptions.

Parameters:
Input:
x[N]    Input data
y[N]    Sign data
N       Length of input/output data vectors
Output:
z[N]    Results

Restrictions:
z,x,y   Aligned on 32-byte boundary
z,x,y   Must not overlap
N       Multiple of 8
-------------------------------------------------------------------------*/
void vcopysignf ( float32_t * restrict z, 
            const float32_t * restrict x, 
            const float32_t * restrict y, 
            int N );

float32_t scopysignf ( float32_t x , float32_t y);

/*-------------------------------------------------------------------------
Clipping

Description: These functions limit the absolute value of inputs. If 
magnitude of input value is less than the second argument then it is left
unchanged. Otherwise it is replaced by absolute value of the second argument
or its negation, depending on the sign of the input value.

Data format: IEEE-754 Std. single precision floating-point.

Parameters:
Input:
x[N]  Input data
y     Limiting value
N     Length of input/output data vectors
Output:
z[N]  Results

Restrictions:
z,x   Aligned on 32-byte boundary
z,x   Must not overlap
N     Multiple of 8
-------------------------------------------------------------------------*/
void vclipf ( float32_t * restrict z, 
        const float32_t * restrict x, 
              float32_t            y,
        int N );

float32_t sclipf ( float32_t x, float32_t y );

/*-------------------------------------------------------------------------
Average of Two Values

Description: These functions compute the average of two arguments.

Data format: IEEE-754 Std. single precision floating-point.

Accuracy: 
2 ULP

Input domain for 'fast' version vfastavgf():
|x+y|<Inf
The output value is not defined outside of this range.

Parameters:
Input:
x[N]   Input data
y[N]   Input data
N      Length of input/output data vectors
Output:
z[N]   Results

Restrictions:
z,x,y  Aligned on 32-byte boundary
z,x,y  Must not overlap
N      Multiple of 8
-------------------------------------------------------------------------*/
void vavgf ( float32_t * restrict z, 
       const float32_t * restrict x, 
       const float32_t * restrict y,
       int N );
void vfastavgf ( float32_t * restrict z, 
       const float32_t * restrict x, 
       const float32_t * restrict y,
       int N );
float32_t savgf ( float32_t x, float32_t y );

/*-------------------------------------------------------------------------
Count One Bits in a Word

Description: Functions count the number of one bits in the number or each number of a vector.

Data format: 16-bit/32-bit fixed-point format

Accuracy: exact

Parameters:
Input:
x[N]   Input data, 16-bit/32-bit
N      Length of input/output data vectors
Output:
z[N]   Results, 16-bit/32-bit

Restrictions:
z,x,y  Aligned on 32-byte boundary
z,x,y  Must not overlap
N      Multiple of 16 (vcountones16) or 8 (vcountones32)
-------------------------------------------------------------------------*/
void vcountones16(int16_t *z, const int16_t   *x, int N);
void vcountones32(int32_t *z, const int32_t   *x, int N);

int countones16(int16_t x);
int countones32(int32_t x);

#ifdef __cplusplus
};
#endif

#endif /* __NATUREDSP_BASEBAND_MATH_H */
