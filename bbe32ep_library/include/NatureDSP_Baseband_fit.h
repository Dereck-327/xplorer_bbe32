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
 * Fitting and Interpolation Routines
 */

#ifndef __NATUREDSP_BASEBAND_FIT_H
#define __NATUREDSP_BASEBAND_FIT_H

#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
Fitting and Interpolation Routines
pfit         Polynomial fitting and interpolation for fixed-point data
[r]proot     Find roots of a polynomial with complex/real coefficients
vinterp      Vector interpolation
interp_tap   2/4 tap interpolators
===========================================================================*/

/*-------------------------------------------------------------------------
Polynomial Fitting and Interpolation

Description: the pfit functions fit (in least squares sense) a degree N 
polynomial to input data sampled at a points grid of length M, and use
that polynomial to interpolate data at arbitrary query points. Namely,
the pfitN_grid() functions compute the Vandermonde matrix for the sample
points grid and perform the Cholesky decomposition of that matrix, the
pfitN_process() functions calculate the least squares solution for the
polynomial coefficients. Finally, the pfitN_eval() functions evaluate
the polynomial at query points.

Please refer to the NatureDSP Baseband Library Reference for full details
on these functions.

Representation:
pfit_gridN,      16-bit fixed-point data. Parameter specifications denote
pfit_processN,   fixed-point format for various data items
pfit_evalN       
pfitf_gridN,     IEEE-754 Std single precision floating-point data
pfitf_processN,  
pfitf_evalN     

Note:
Number of fractional bits specidied for various input/output arguments below apply
for the fixed-point variant

Parameters:
Input:
N                     Degree of polynomial, 1..6
M                     Number of sample points
P                     Number of query points
maxIter               Number of least squares solution enhancement iterations. Right 
                      choice depends on required accuracy, the ad-hoc value is (N+1)/2
x[M]                  Sample points grid, Q15 or floating point
y[M]                  Sampled data values, Q15 or floating point
xi[P]                 Query points, Q15 or floating point
M'=(M+7)&(~7), N'=8   for floating point API
M'=(M+15)&(~15),N'=16 for fixed-point API

Intermediate:
V[M'*8]               Vandermonde matrix, Q15 or floating point
R[N'*8]               Upper triangular Cholesky factor of matrix V, Q11 or floating point
Output:
yi[P]                 Data values interpolated at query points, Q15 or floating point
p[N+1]                Polynomial coefficients, Q8.23 or floating point
Temporary:
pScr                  Scratch memory area. To determine the scratch area size required by
                      a function pfitN_<fun>, use the respective helper function 
                      pfit_<fun>_getScratchSize(M,N)

Restrictions:
M>N                   The number of sample points must exceed the degree of polynomial
x,y,xi,yi,V,R,p,pScr  Must not overlap
V,R,pScr              Aligned on 32-byte boundary
---------------------------------------------------------------------------*/
void pfit_grid1
(
  void * restrict          pScr, /* (temp) Scratch memory area                          */
  int16_t * restrict       V,    /* (out ) Vandermonde matrix                           */
  int16_t * restrict       R,    /* (out ) Upper triangular Cholesky factor of matrix V */
  const int16_t * restrict x,    /* (in  ) Sample points grid                           */
  int M                          /* (in  ) Number of sample points                      */
);
void pfit_grid2 ( void * restrict pScr, int16_t * restrict  V, int16_t * restrict R, const int16_t * restrict x, int M );
void pfit_grid3 ( void * restrict pScr, int16_t * restrict  V, int16_t * restrict R, const int16_t * restrict x, int M );
void pfit_grid4 ( void * restrict pScr, int16_t * restrict  V, int16_t * restrict R, const int16_t * restrict x, int M );
void pfit_grid5 ( void * restrict pScr, int16_t * restrict  V, int16_t * restrict R, const int16_t * restrict x, int M );
void pfit_grid6 ( void * restrict pScr, int16_t * restrict  V, int16_t * restrict R, const int16_t * restrict x, int M );

void pfitf_grid1( void * restrict pScr, float32_t * restrict  V, float32_t * restrict R, const float32_t * restrict x, int M );
void pfitf_grid2( void * restrict pScr, float32_t * restrict  V, float32_t * restrict R, const float32_t * restrict x, int M );
void pfitf_grid3( void * restrict pScr, float32_t * restrict  V, float32_t * restrict R, const float32_t * restrict x, int M );
void pfitf_grid4( void * restrict pScr, float32_t * restrict  V, float32_t * restrict R, const float32_t * restrict x, int M );
void pfitf_grid5( void * restrict pScr, float32_t * restrict  V, float32_t * restrict R, const float32_t * restrict x, int M );
void pfitf_grid6( void * restrict pScr, float32_t * restrict  V, float32_t * restrict R, const float32_t * restrict x, int M );


/* Return the scratch area size, in bytes. */
size_t pfit_grid_getScratchSize( int M, int N );
size_t pfitf_grid_getScratchSize( int M, int N );

void pfit_process1
(
  void * restrict          pScr, /* (temp) Scratch memory area         */
  int32_t * restrict       p,    /* (out ) Polynomial coefficients     */
  const int16_t * restrict V,    /* (in  ) Vandermonde matrix          */
  const int16_t * restrict R,    /* (in  ) Cholesky factor of matrix V */
  const int16_t * restrict y,    /* (in  ) Sampled data values         */
  int M,                         /* (in  ) Number of sample points     */
  int maxIter                    /* (in  ) Number of iterations        */
);
void pfit_process2 ( void * restrict pScr, int32_t * restrict p, const int16_t * restrict V, const int16_t * restrict R, const int16_t * restrict y, int M, int maxIter );
void pfit_process3 ( void * restrict pScr, int32_t * restrict p, const int16_t * restrict V, const int16_t * restrict R, const int16_t * restrict y, int M, int maxIter );
void pfit_process4 ( void * restrict pScr, int32_t * restrict p, const int16_t * restrict V, const int16_t * restrict R, const int16_t * restrict y, int M, int maxIter );
void pfit_process5 ( void * restrict pScr, int32_t * restrict p, const int16_t * restrict V, const int16_t * restrict R, const int16_t * restrict y, int M, int maxIter );
void pfit_process6 ( void * restrict pScr, int32_t * restrict p, const int16_t * restrict V, const int16_t * restrict R, const int16_t * restrict y, int M, int maxIter );

void pfitf_process1( void * restrict pScr, float32_t * restrict p, const float32_t * restrict V, const float32_t * restrict R, const float32_t * restrict y, int M, int maxIter );
void pfitf_process2( void * restrict pScr, float32_t * restrict p, const float32_t * restrict V, const float32_t * restrict R, const float32_t * restrict y, int M, int maxIter );
void pfitf_process3( void * restrict pScr, float32_t * restrict p, const float32_t * restrict V, const float32_t * restrict R, const float32_t * restrict y, int M, int maxIter );
void pfitf_process4( void * restrict pScr, float32_t * restrict p, const float32_t * restrict V, const float32_t * restrict R, const float32_t * restrict y, int M, int maxIter );
void pfitf_process5( void * restrict pScr, float32_t * restrict p, const float32_t * restrict V, const float32_t * restrict R, const float32_t * restrict y, int M, int maxIter );
void pfitf_process6( void * restrict pScr, float32_t * restrict p, const float32_t * restrict V, const float32_t * restrict R, const float32_t * restrict y, int M, int maxIter );

/* Return the scratch area size, in bytes. */
size_t pfit_process_getScratchSize(int M, int N);
size_t pfitf_process_getScratchSize(int M, int N);

void pfit_eval1
(
  int16_t * restrict       yi, /* (out) Data values interpolated at query points */
  const int16_t * restrict xi, /* (in ) Query points                             */
  const int32_t * restrict p,  /* (in ) Polynomial coefficients                  */
  int P                        /* (in ) Number of query points                   */
);
void pfit_eval2 ( int16_t   * restrict yi, const int16_t   * restrict xi, const int32_t   * restrict p, int P );
void pfit_eval3 ( int16_t   * restrict yi, const int16_t   * restrict xi, const int32_t   * restrict p, int P );
void pfit_eval4 ( int16_t   * restrict yi, const int16_t   * restrict xi, const int32_t   * restrict p, int P );
void pfit_eval5 ( int16_t   * restrict yi, const int16_t   * restrict xi, const int32_t   * restrict p, int P );
void pfit_eval6 ( int16_t   * restrict yi, const int16_t   * restrict xi, const int32_t   * restrict p, int P );

void pfitf_eval1( float32_t * restrict yi, const float32_t * restrict xi, const float32_t * restrict p, int P );
void pfitf_eval2( float32_t * restrict yi, const float32_t * restrict xi, const float32_t * restrict p, int P );
void pfitf_eval3( float32_t * restrict yi, const float32_t * restrict xi, const float32_t * restrict p, int P );
void pfitf_eval4( float32_t * restrict yi, const float32_t * restrict xi, const float32_t * restrict p, int P );
void pfitf_eval5( float32_t * restrict yi, const float32_t * restrict xi, const float32_t * restrict p, int P );
void pfitf_eval6( float32_t * restrict yi, const float32_t * restrict xi, const float32_t * restrict p, int P );

/*---------------------------------------------------------------------------
Find Roots of a Polynomial with Complex/Real Coefficients

Given input array c[N+1] containing N+1 complex/real coefficients of a 
degree N polynomial, find its N roots (possibly with repetitions) and
store results to the output array r[N].

Representation: IEEE-754 Std single precision floating-point data

Parameters:
Temporary:
  pScr      Scratch area. Required size (in bytes) is defined by functions 
            [r]prootf_getScratchSize(N)
Input:
  N         Polynomial degree
  c[N+1]    Polynomial coefficients in descending powers order, i.e. 
            p(x) = c[0]*x^N + ... + c[N-1]*x + c[N]
Output:
  r[N]      Roots of the polynomial

Restrictions:
  pScr,r,c  Must not overlap and must be aligned on 32-byte boundary
  N>0       Polynomial degree must be positive
  c[0]!=0   The leading coefficient must be non-zero
---------------------------------------------------------------------------*/
/* Complex polynomial coefficients */
void prootf (
          void * pScr, 
          complex_float * restrict r,
    const complex_float * restrict c,
    int N );
/* Real polynomial coefficients */
void rprootf ( 
          void * pScr, 
          complex_float * restrict r,
    const float32_t     * restrict c,
    int N );

/* Return the scratch area size, in bytes. */
size_t prootf_getScratchSize  ( int N );
size_t rprootf_getScratchSize ( int N );

/*-------------------------------------------------------------------------
Vector Interpolation

For real-valued target function: given array y[M] of the target function
sampled at M control points x[M], find values yi[N] of the target function
at points xi[N].

For complex-valued target function: given array y[2*M] of the target function
sampled at M real control points x[M], find values yi[2*N] of the target
function at points xi[N].

There are two options for the interpolation method:
  - for each of xi[N], locate the nearest control point (the neighbour) and
    return the corresponding value of the target function
  - perform the linear interpolation between the two closest control points

Representation: 16-bit signed fixed-point

Methods:
  vinterp_neib()       - interpolate a real-valued function by locating
                         the nearest control point
  vinterp_cneib()      - interpolate a complex-valued function by locating
                         the nearest control point
  vinterp_lin()        - perform the linear interpolation for a 
                         real-valued target function
  vinterp_clin()       - perform the linear interpolation for a 
                         complex-valued target function
---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------
Find values of the target function at points xi[N] by locating the
nearest control point.

Representation: 16-bit signed fixed-point

Parameters:
Input:
  C          1 for real data (vinterp_neib), 2 for complex (vinterp_cneib)
  xi[N]      Desired points
  y[C*M]     Values of the target function sampled at control points
  x[M]       Control points
Output:
  yi[C*N]    Values of the target function at points xi[N]
Temporary:
  scr[]      Scratch memory area. To determine the scratch area size required by
             a function vinterp_<fun>, use the respective helper function 
             vinterp_<fun>_getScratchSize(M,N)

Restrictions:
  scr,x[M],y[C*M],xi[N],yi[C*N]  Must not overlap
  x[M],xi[N]                     Must be sorted in the ascending order
  scr,x[M],xi[N],yi[C*N]         Must be aligned on 32-byte boundary 
  y[C*M]                         Must be aligned on (C*2)-byte boundary
  M>=1                           The minimum number of control points is 1
  N                              Must be a multiple of 16
---------------------------------------------------------------------------*/
void vinterp_neib  ( void *    restrict scr,
                     int16_t * restrict yi,
               const int16_t *          xi,
               const int16_t *          y,
               const int16_t *          x,
               int M, int N );
void vinterp_cneib ( void *    restrict scr,
                     int16_t * restrict yi,
               const int16_t *          xi,
               const int16_t *          y,
               const int16_t *          x,
               int M, int N );

/* Return the scratch area size, in bytes. */
size_t vinterp_neib_getScratchSize  ( int M, int N );
size_t vinterp_cneib_getScratchSize ( int M, int N );

/*---------------------------------------------------------------------------
Find values of the target function at points xi[N] by linear
interpolation between the two closest control points.

Representation: 16-bit signed fixed-point

Parameters:
Input:
  C          1 for real data (vinterp_lin), 2 for complex (vinterp_clin)
  scr[]      Scratch area, VINTERP_[C]LIN_SCRATCH(M,N) bytes
  xi[N]      Desired points
  y[C*M]     Values of the target function sampled at control points
  x[M]       Control points
Output:
  yi[C*N]    Values of the target function at points xi[N]
Temporary:
  scr[]      Scratch memory area. To determine the scratch area size required by
             a function vinterp_<fun>, use the respective helper function 
             vinterp_<fun>_getScratchSize(M,N)

Restrictions:
  scr,x[M],y[C*M],xi[N],yi[C*N]  Must not overlap
  x[M],xi[N]                     Must be sorted in the ascending order
  x[M]                           Must be distinct
  scr,x[M],y[C*M],xi[N],yi[C*N]  Must be aligned on 32-byte boundary 
  M>=2                           The minimum number of control points is 2
  N                              Must be a multiple of 16

Note:
  Interpolation result is not defined for those desired points that do
  not belong to the closed interval between the left-most and right-most
  control points: xi<x[0] || xi>x[M-1].
---------------------------------------------------------------------------*/
void vinterp_lin  ( void *    restrict scr,
                    int16_t * restrict yi,
              const int16_t *          xi,
              const int16_t *          y,
              const int16_t *          x,
              int M, int N );
void vinterp_clin ( void *    restrict scr,
                    int16_t * restrict yi,
              const int16_t *          xi,
              const int16_t *          y,
              const int16_t *          x,
              int M, int N );

/* Return the scratch area size, in bytes. */
size_t vinterp_lin_getScratchSize  ( int M, int N );
size_t vinterp_clin_getScratchSize ( int M, int N );

/*-------------------------------------------------------------------------
2/4 Tap Interpolators

Representation: 16-bit signed fixed-point

Functions calculate 
y[k]=taps[0]*x1[k]+taps[1]*x2[k];
or 
y[k]=taps[0]*x1[k]+taps[1]*x2[k]+taps[2]*x3[k]+taps[3]*x4[k];
for all elements of input vectors
where taps[] is real vector and x<1..4>[] are complex vectors

Input:
N             size of complex vectors
x1[2*N]       Input complex vectors of length N, Q15
x2[2*N] 
x3[2*N] 
x4[2*N] 
taps[2 or 4]  Taps, Q15
Output:
y[2*N]        Output complex vector of length N, Q15

Return value:
none

Restrictions:
x1,x2,x3,x4,y       Must be aligned on 32-byte boundary
x1,x2,x3,x4,y,taps  Must not overlap
N                   Must be a multiple of 8
-------------------------------------------------------------------------*/
void interp_2tap ( int16_t * restrict y,
             const int16_t * restrict x1,
             const int16_t * restrict x2,
             const int16_t * restrict taps,
             int N );
void interp_4tap ( int16_t * restrict y,
             const int16_t * restrict x1,
             const int16_t * restrict x2,
             const int16_t * restrict x3,
             const int16_t * restrict x4,
             const int16_t * restrict taps,
             int N );

#ifdef __cplusplus
};
#endif

#endif /* __NATUREDSP_BASEBAND_FIT_H */
