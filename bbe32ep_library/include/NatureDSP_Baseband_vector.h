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
 * Vector Operations
 */

#ifndef __NATUREDSP_BASEBAND_VECTOR_H
#define __NATUREDSP_BASEBAND_VECTOR_H

#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
Vector Operations
rvdot           Dot product of real vectors
cvdot           Dot product of complex vectors
rvadd           Real vectors sum
cvadd           Complex vectors sum
vpower          Sum of squares of a vector
vmag            Square root of sum of squares 
vnorm           Vector normalization
vbexp,sbexp     Common block exponent (vector, scalar)
vbexp_fast      Common exponent (fast vector variant)
rvmean,
cvmean          Mean of Vector Elements
vmax,vmin       Dual Peak Search
vmax8           Search for 8 Top Values
vthreshold      Find Values above Threshold
rmelements,
cmelements      Move the Elements at Given Indices
rrelements, 
crelements      Remove the Elements at Given Indices (remove holes)
aunwrap         Angle Unwrapping
===========================================================================*/

/*-------------------------------------------------------------------------
Dot Product of Real Vectors

Description: These routines take two real vectors and calculate their dot 
product.

Representation:
rvdot   Signed fixed-point format
        Input vectors x and y are 16-bit signed data of arbitrary formats
        Qx and Qy. Dot product is computed in a 40-bit accumulator, which
        is then rounded, shifted to the right by rsh bit positions and
        saturated to form a 32-bit result with Qx+Qy-rsh fractional bits.
rvdotf  IEEE-754 Std. single precision floating-point format for
        input vectors and dot product result

Parameters:
Input:
x[N]    Input vector
y[N]    Input vector
rsh     Right shift amount (rvdot)
N       Length of vectors
Returned Value:
Dot product result

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 16 (rvdot) or 8 (rvdotf)
rsh>=0  Right shift amount must be non-negative
-------------------------------------------------------------------------*/
int32_t   rvdot  ( const int16_t   * restrict x,
                   const int16_t   * restrict y,
                   int rsh, 
                   int N );
float32_t rvdotf ( const float32_t * restrict x,
                   const float32_t * restrict y,
                   int N );

/*-------------------------------------------------------------------------
Dot Product of Complex Vectors

Description: These routines take two complex vectors and calculate their
dot product.

Representation:
cvdot   Signed fixed-point format
        Input vectors x and y are comprised of 32-bit complex data with
        16-bit signed real and imaginary components of arbitrary formats
        Qx and Qy. Complex dot product is computed in a pair of 40-bit 
        accumulators, which are then rounded, shifted to the right by rsh
        bit positions and saturated to form a 64-bit complex result with
        Qx+Qy-rsh fractional bits in 32-bit real and imaginary components.
cvdotf  IEEE-754 Std. single precision floating-point format for real/imaginary
        components of 64-bit input data and dot product result

Parameters:
Input:
x[N]    Input vector
y[N]    Input vector
rsh     Right shift amount (cvdot)
N       Length of vectors
Returned Value:
Dot product result

Restrictions:
x,y     Must not overlap
x,y     Aligned on 32-byte boundary
N       Multiple of 8 (cvdot) or 4 (cvdotf)
rsh>=0  Right shift amount must be non-negative
-------------------------------------------------------------------------*/
complex_fract32 cvdot  ( const complex_fract16 * restrict x,
                         const complex_fract16 * restrict y,
                         int rsh,
                         int N );
complex_float   cvdotf ( const complex_float   * restrict x,
                         const complex_float   * restrict y,
                         int N );

/*-------------------------------------------------------------------------
Real Vectors Sum

Description: These routines perform pairwise summation of real vectors.

Representation:
rvadd   Signed fixed-point format. 16-bit inputs, 16-bit saturated results
rvaddf  IEEE-754 Std. single precision floating-point format for input/output
        data

Parameters:
Input:
x[N]   Input vector
y[N]   Input vector
N      Length of vectors 
Output:
z[N]   Sum of input vectirs

Restrictions:
z,x,y  Must not overlap
z,x,y  Aligned on 32-byte boundary
N      Multiple of 16 (rvadd) or 8 (rvaddf)
-------------------------------------------------------------------------*/
void rvadd  ( int16_t   * restrict z,
        const int16_t   * restrict x,
        const int16_t   * restrict y,
        int N );
void rvaddf ( float32_t * restrict z,
        const float32_t * restrict x,
        const float32_t * restrict y,
        int N );

/*-------------------------------------------------------------------------
Complex Vectors Sum

Description: These routines perform pairwise summation of complex vectors.

Representation:
cvadd   Signed fixed-point format. Input vectors are comprised of 32-bit
        complex elements with 16-bit real and imaginary components. Output
        vector elements are 32-bit complex values with 16-bit saturated
        real and imaginary parts.
cvaddf  IEEE-754 Std. single precision floating-point format for real/imaginary
        components of 64-bit input/output data

Parameters:
Input:
x[2*N]  Input complex data,Q15
y[2*N]  Input complex data,Q15
N       Length of vectors
Output:
z[2*N]  Output complex data,Q15

Restrictions:
z,x,y   Must not overlap
x,y,z   Aligned on 32-byte boundary
N       Multiple of 8 (cvadd) or 4 (cvaddf) 
-------------------------------------------------------------------------*/
void cvadd  ( complex_fract16 * restrict z,
        const complex_fract16 * restrict x,
        const complex_fract16 * restrict y,
        int N );
void cvaddf ( complex_float   * restrict z,
        const complex_float   * restrict x,
        const complex_float   * restrict y,
        int N );

/*-------------------------------------------------------------------------
Sum of Squares of a Vector

Description: These routines compute the power of a vector.

Representation:
vpower   Signed fixed-point format
         Input vector elements are 16-bit signed data of arbitrary format
         Qx. Sum of squared values is computed in a 40-bit accumulator,
         which is then rounded, shifted to the right by rsh bit positions and
         saturated to form a 32-bit result with 2*Qx-rsh fractional bits.
vpowerf  IEEE-754 Std. single precision floating-point format for the
         input vector and the result

Parameters:
Input:
x[N]     Input vector
rsh      Right shift amount (vpower)
N        Length of input vector 
Returned Value:
Sum of squares over the input vector

Restrictions:
x        Aligned on 32-byte boundary
N        Multiple of 16 (vpower) or 8 (vpowerf)
rsh>=0   Right shift amount must be non-negative
-------------------------------------------------------------------------*/
int32_t   vpower  ( const int16_t   * restrict x,
                    int rsh,
                    int N );
float32_t vpowerf ( const float32_t * restrict x,
                    int N );

/*-------------------------------------------------------------------------
Square Root of Sum of Squares 

Description: These routines compute the magnitude of a vector

Representation:
vmag    Signed fixed-point format
        Input vector elements are 16-bit signed data of arbitrary format
        Qx. Sum of squared values is computed in a 40-bit accumulator,
        which is then shifted to the right by rsh bit positions. Square
        root approximation of the shifted sum is formatted as a 16-bit
        non-negative value with Qx-rsh/2 fractional bits.
vmagf   IEEE-754 Std. single precision floating-point format for
        input vector and magnitude estimation result

Parameters:
Input:
x[N]    Input data
rsh     Right shift amount (vmag)
N       Length of vector 
Returned Value:
Square root of the sum of squares over the input vector

Restrictions:
x       Aligned on 32-byte boundary
N       Multiple of 16 (vmag) or 8 (vmagf)
rsh>=0  Right shift amount must be non-negative
-------------------------------------------------------------------------*/
int16_t   vmag  ( const int16_t * restrict x,
                  int rsh,
                  int N );
float32_t vmagf ( const float32_t * restrict x,
                  int N );

/*-------------------------------------------------------------------------
Vector Normalization

Description: This function scale all data in the input vector x by a power
of 2 factor and returns results in the output vector y. Most often it is used
in conjunction with common block exponent function vbexp() to normalise vector
data, i.e. to shift all elements of a vector to the left such that the minimum
number of redundant sign bits over the output vector is zero. It can, however,
be utilised for shifting data by arbitrary number of bit positions: for positive
shift amount it is a saturating left shift, for negative shift amount it is a
sign-extending right shift.

Representation:
vnorm  16-bit signed fixed-point data

Parameters:
x[N]   Input data
t      Shift amount, [-16..16]
N      Length of input/output vectors
Output:  
y[N]   Output data

Restrictions:
x,y    Aligned on 32-byte boundary
x,y    Must not overlap
N      Multiple of 16
-------------------------------------------------------------------------*/
void vnorm  ( int16_t   * restrict y,
        const int16_t   * restrict x,
        int t,
        int N );

/*-------------------------------------------------------------------------
Common Block Exponent

Description: These functions compute base-2 exponent adjustment term needed
to normalize data in the input vector. Exact meaning of normalization depends
on the data format (see below).

Representation: 
vbexp,vbexp_fast        16-bit signed fixed-point format
                        For each input value functions count the number of
                        redundant sign bits (as if the value was loaded in
                        a 32-bit register) and return the minimum result over
                        the input data vector.
sbexp                   32-bit signed fixed-point format
                        Count the number of redundant sign bits and return 
                        the result.
vbexpf,sbexpf           IEEE-754 Std. single precision floating-point format
                        For each finite input value x, functions estimate the
                        integer E(x), such that 0.5 <= |x|*2^E(x) < 1 and
                        -128 <= E(x) <= 148. The minimum value of E(x) over
                        input data vector is the result.

Special cases:
   x    |  Result |    Extra Conditions    
--------+---------|---------------------------
0       |    0    |
+/-Inf  | -129    | floating-point functions
NaN     |    0    |
--------|---------|---------------------------
0       |   31    |
-32768  |   16    | fixed-point functions
32767   |   16    |

Parameters:
Input:
x[N]    Input data
N       Length of data vector
Returned Value: exponent adjustment term, or zero if N<=0
Restrictions:
vbexp(), vbexpf():
  No restrictions
vbexp_fast():
  x     Aligned on 32-byte boundary
  N     Multiple of 16
-------------------------------------------------------------------------*/
int vbexp      ( const int16_t   * restrict x, int N );
int vbexp_fast ( const int16_t   * restrict x, int N );
int vbexpf     ( const float32_t * restrict x, int N );

int sbexp  ( int32_t   x );
int sbexpf ( float32_t x );

/*-------------------------------------------------------------------------
Mean of Vector Elements

Description: Compute the mean value over all elements of given vector (real 
or complex).

Representation:
rvmean,cvmean    16-bit signed fixed-point format
rvmeanf,cvmeanf  IEEE-754 Std. single precision floating-point format

Parameters:
Input:
x[N]    Input vector
N       Length of input vector, in real or complex samples
Output:
m[1]    Mean value, real or complex

Restrictions:
x       Aligned on 32-byte boundary
x,m     Must not overlap
N       Must be a multiple of either:
          4 (cvmeanf), or
          8 (cvmean, rvmeanf), or
          16 (rvmean)
-------------------------------------------------------------------------*/
// Process real data:
void rvmean ( int16_t   * m, const int16_t   * restrict x, int N );
void rvmeanf( float32_t * m, const float32_t * restrict x, int N );
// Process complex data:
void cvmean ( complex_fract16 * m, const complex_fract16 * restrict x, int N );
void cvmeanf( complex_float   * m, const complex_float   * restrict x, int N );

/*-------------------------------------------------------------------------
Dual Peak Search 

Description: These functions retrieve the maximum (minimum) and next-to-
maximum (next-to-minimum) values of vector elements. They output both the 
peak values and their indices. 

Representation:
vmax,vmin    16-bit signed fixed-point format
             Special values for vmax and vmin (see note 1) are -32768 and
             32767, respectively
vmaxf,vminf  IEEE-754 Std. single precision floating-point format
             Special values for vmaxf and vminf are -HUGE_VALF and +HUGE_VALF,
             respectively

Notes:
1. Each kind of dual-peak search function reserves a special value (see
   above) to maintain internal invariants during the search process. However,
   it is still legal for the input vector to contain special values among
   other data. If this is the case, then the following peculiarities should
   be considered:
     A) Input vector elements equal to the special value are ignored.
     B) If the total number K of input vector elements distinct from the special
        value is less than 2 (i.e. K==0 or K==1), then the last 2-K entries of output
        vector idx[2] are assigned zero, and the last 2-K elements of output vector
        m[2] are assigned the special value.
2. If the peak value is encountered more than once in the input data vector, then
   it will be reported twice in m[2], and idx[2] will contain indices of the first
   two occurencies of the peak value.
3. For floating-point functions NAN values in vectors are ignored. If, however,
   the input vector contains NANs only, then the entries of output vector idx[2] are
   assigned zero, and the elements of output vector m[2] are assigned the special value.
Parameters:
Input:
x[N]     Input data vector
N        Length of input data vector
Output:
m[2]     2 peak values is descending (vmax) or ascending (vmin) order
idx[2]   Indices of 2 peak elements; optional

Restrictions:
x        Aligned on 32-byte boundary
x,m,idx  Must not overlap
N        Multiple of 16 (vmax, vmin) or 8 (maxf, xminf)
-------------------------------------------------------------------------*/
void vmax  ( int16_t   * restrict m, 
             int16_t   * restrict idx, 
       const int16_t   * restrict x, int N );
void vmaxf ( float32_t * restrict m, 
             int16_t   * restrict idx, 
      const  float32_t * restrict x, int N );
void vmin  ( int16_t   * restrict m, 
             int16_t   * restrict idx, 
       const int16_t   * restrict x, int N );
void vminf ( float32_t * restrict m, 
             int16_t   * restrict idx, 
       const float32_t * restrict x, int N );

/*-------------------------------------------------------------------------
Search for 8 Top Value Elements

Description: Functions search for 8 top values in the input vector. They
output both the peak values and respective indices.

Representation:
vmax8   16-bit signed fixed-point format
vmax8f  IEEE-754 Std. single precision floating-point format

Notes:
1. The function requires read/write access to the input vector x[N], because 
   it performs some data manipulations on it. Before returning, it always 
   restores original data.
2. Both top-8 search functions reserve a special value to maintain internal
   invariants during the search process. For vmax8 this special value is -32768,
   and for vmax8f it is -HUGE_VALF. It is still legal for the input vector to
   contain special values among other data. If this is the case, then the
   following peculiarities should be considered:
     A) Input vector elements equal to the special value are ignored.
     B) If the total number K of input vector elements distinct from the special
        value is less than 8 (i.e. K<8), then the last 8-K entries of output vector
        idx[8] are assigned zero, and the last 8-K elements of output vector m[8]
        are assigned the special value.
3. All elements of the input vector x[N] holding the same top value will be
   reported in search results, unless the total limit of 8 peak elements is
   exhausted.
4. For floating-point functions NAN values in vectors are ignored. If, however,
   the input vector contains NANs only, then the entries of output vector idx[8] are
   assigned zero, and the elements of output vector m[8] are assigned the special value.

Example (vmax8):
     Input:   N:    16
              x[N]: 1,-32768,1,2,3,-32768,-32768,-32768,-32768,-32768,
                    -32768,-32768,-32768,-32768,-32768,-32768
     Note that the number of input value distinct from -32768 is less than 8!
     Also note the duplicated value at positions 0 and 2.
     Results: m[8]: 3,2,1,1,-32768,-32768,-32768,-32768
            idx[8]: 4,3,0,2,0,0,0,0

Parameters:
Input:
x[N]     Input data vector
N        Length of input data vector
Output:
m[8]     8 maximum values in descending order
idx[8]   Indices of 8 maximum values in the input vector

Restrictions:
x        Aligned on 32-byte boundary
x,m,idx  Must not overlap
N        Multiple of 16 (vmax8) or 8 (vmax8f)
-------------------------------------------------------------------------*/
void vmax8  ( int16_t   * restrict m, 
              int16_t   * restrict idx, 
              int16_t   * restrict x, int N );
void vmax8f ( float32_t * restrict m, 
              int16_t   * restrict idx, 
              float32_t * restrict x, int N );

/*-------------------------------------------------------------------------
Find Values Above Threshold 

Description: function collects indices of all input vector elements that hold
a value greater than the designated threshold.

Representation:
vthreshold   16-bit signed fixed-point format
vthresholdf  IEEE-754 Std. single precision floating-point format

Parameters:
Input:
thr       Threshold value
x[N]      Input data
N         Size of input array
Output:
idx[N+1]  Indices of elements with a value greater the threshold
NOTE: 
extra cell should be reserved in the end of array
Returns:  Number of found elements

Restrictions:
x,idx     Aligned on 32-byte boundary
x,idx     Must not overlap
N         Multiple of 16 (vthreshold) or 8 (vthresholdf)
-------------------------------------------------------------------------*/
int vthreshold  ( int16_t   * restrict idx, 
                  int16_t              thr,
            const int16_t   * restrict x,
            int N );
int vthresholdf ( int16_t   * restrict idx, 
                  float32_t            thr,
            const float32_t * restrict x,
            int N );

/*-------------------------------------------------------------------------
Move the Elements at Given Indices 

Description: <r|c>melements<M> Functions retrieve each M-th real or complex
element of the input vector x (x[0],x[M],x[2*M],...) and stores it to the
output vector y.

Representation: 16-bit fixed-point data

Parameters:
Input:
x[N*M]  Input data vector
N       Size of output data vector
Output:
y[N]    Output data vector

Restrictions:
x,y     Aligned on 32-byte boundary
x,y     Must not overlap
N       Multiple of 16 for real data, or a multiple of 8 for complex data
-------------------------------------------------------------------------*/
void rmelements2  ( int16_t * restrict y, const int16_t * restrict x, int N );
void rmelements3  ( int16_t * restrict y, const int16_t * restrict x, int N );
void rmelements4  ( int16_t * restrict y, const int16_t * restrict x, int N );
void rmelements6  ( int16_t * restrict y, const int16_t * restrict x, int N );
void rmelements12 ( int16_t * restrict y, const int16_t * restrict x, int N );
void rmelements24 ( int16_t * restrict y, const int16_t * restrict x, int N );

void cmelements2  ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );
void cmelements3  ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );
void cmelements4  ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );
void cmelements6  ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );
void cmelements12 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );
void cmelements24 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );

/*-------------------------------------------------------------------------
Remove the Elements at Given Indices 

Description: <r|c>relements<M> copies real or complex elements from the input
vector x to the output vector y, with each M-th element (x[0],x[M],x[2M],...)
being discarded.

Representation: 16-bit fixed-point data

Parameters:
Input:
x[N*M]      Input data vector
N           Size of input data vector divided by M
Output:
y[N*(M-1)]  Output data vector

Restrictions:
x,y         Aligned on 32-byte boundary
x,y         Must not overlap
N           Multiple of 16 for real data, or a multiple of 8 for complex data
-------------------------------------------------------------------------*/
void rrelements2  ( int16_t * restrict y, const int16_t * restrict x, int N );
void rrelements3  ( int16_t * restrict y, const int16_t * restrict x, int N );
void rrelements4  ( int16_t * restrict y, const int16_t * restrict x, int N );
void rrelements6  ( int16_t * restrict y, const int16_t * restrict x, int N );
void rrelements12 ( int16_t * restrict y, const int16_t * restrict x, int N );
void rrelements24 ( int16_t * restrict y, const int16_t * restrict x, int N );

void crelements2  ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );
void crelements3  ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );
void crelements4  ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );
void crelements6  ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );
void crelements12 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );
void crelements24 ( complex_fract16 * restrict y, const complex_fract16 * restrict x, int N );

/*-------------------------------------------------------------------------
Angle Unwrapping

Description: These Functions correct phase angles to produce smooth phase 
plots. They adjust phase angles from the input vector by adding a multiple
of 2pi when the absolute difference between two consecutive elements of the
input vector is greater than or equal to the jump tolerance. Tolerance values
smaller than pi are equivalent to the tolerance value of pi.

Representation:
aunwrap   16-bit signed fixed-point format
          It is assumed that all real-world angles are normalized by pi when 
          converted to a fixed-point format. For example, -pi represented
          in Q14 format would be -1*2^14 == -16384.
          Fixed-point format of input phase angles x[N] is Q15. 
          To accommodate the output data format to possible 2pi cycles, the
          number of fractional bits for unwrapped angles y[N] is specified
          through the q argument.
          Fixed-point format for the jump tolerance value tol is Q14.
aunwrapf  IEEE-754 Std. single precision floating-point format for input/output
          data

Parameters:
Input:
x[N]      Original phase angles
tol       Jump tolerance. Typical value is pi (16384 in Q14).
q         Fixed-point position for output data, 0..15 (aunwrap)
N         Size of input/output arrays
Output:
y[N]      Unwrapped phase angles

Restrictions:
x,y       Aligned on 32-byte boundary
x,y       Must not overlap
N         Multiple of 16 (aunwrap) or 8 (aunwrapf)
-------------------------------------------------------------------------*/
void aunwrap  ( int16_t * restrict y, const int16_t * restrict x, int16_t tol, int q, int N );
void aunwrapf ( float32_t * restrict y, const float32_t * restrict x, float32_t tol, int N );

#ifdef __cplusplus
};
#endif

#endif /* __NATUREDSP_BASEBAND_VECTOR_H */
