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
  NatureDSP_Baseband library. Vector Operations
    Angle Unwrapping
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Vector Operations. */
#include "NatureDSP_Baseband_vector.h"

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
#if !HAVE_VFPU
DISCARD_FUN(void,aunwrapf,( float32_t * restrict y, const float32_t * restrict x, float32_t tol, int N ))
#else
void aunwrapf ( float32_t * restrict y, const float32_t * restrict x, float32_t tol, int N )
{ 
	// MATLAB code:
	//  function y=aunwrap(x,thr)
	//  sz=size(x);
	//  N = length(x);
	//  x=reshape(x,N,1);
	//  % Unwrap phase angles. Algorithm minimizes the incremental phase variation
	//  % by constraining it to the range [-pi,pi]
	//  dx = diff(x,1,1);           % Incremental phase variations
	//  dxs = mod(dx+pi,2*pi) - pi;      % Equivalent phase variations in [-pi,pi)
	//  dx_corr = dxs - dx;         % Incremental phase corrections
	//  dx_corr(abs(dx)<thr) = 0;   % Ignore correction when incr. variation
	//                              % is < CUTOFF
	//  
	//  % Integrate corrections and add to x to produce smoothed phase values
	//  x(2:N) = x(2:N) + cumsum(dx_corr);
	//  y=reshape(x,sz);
	int n;
	vboolN_2 bfirst;
	xb_vecN_2xf32 y0, x0, x1, dx,dxs,dx_corr;
	const xb_vecN_2xf32* restrict pX;
	const xb_vecN_2xf32* restrict pYrd;
		  xb_vecN_2xf32* restrict pY;
	static const union { uint32_t u; float32_t f; } pi = { 0x40490fdb };
	static const union { uint32_t u; float32_t f; } inv2pi = { 0x3e22f983 };

	NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
	NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
	NASSERT(N % (BBE_SIMD_WIDTH / 2) == 0);

	if (N <= 0) return;
	tol=XT_ABS_S(tol);
	pX = (const xb_vecN_2xf32*)x;
	pY = (      xb_vecN_2xf32*)y;
	/* if tolerance is NaN just copy input to the output */
	if( !(tol==tol))
	{
		for (n = 0; n < (N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
		{
			BBE_LVN_2XF32_IP(x0, pX, 2 * BBE_SIMD_WIDTH);
			BBE_SVN_2XF32_IP(x0, pY, 2 * BBE_SIMD_WIDTH);
		}
		return;
	}
	x1 = BBE_LVN_2XF32_I(pX,0);
	x1 = BBE_REPN_2XF32(x1, 0);
	for (n = 0; n < (N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
	{
		vboolN_2 small;
		x0 = x1;
		BBE_LVN_2XF32_IP(x1, pX, 2 * BBE_SIMD_WIDTH);
		x0 = BBE_SELN_2XF32I(x1,x0,  BBE_SELI_ROTATE_LEFT_2);
		dx = BBE_SUBN_2XF32(x1, x0);
		dxs = BBE_CONSTN_2XF32(3); 
		BBE_MULAN_2XF32(dxs, dx, inv2pi.f);
    #if 1
    dx_corr = dx;
    BBE_MULSN_2XF32(dx_corr,BBE_FIFLOORN_2XF32(dxs), 2.f*pi.f);
    dx_corr = BBE_SUBN_2XF32(dx_corr, dx);
    #else
    dx_corr =  BBE_MULN_2XF32(BBE_FIFLOORN_2XF32(dxs), -2.f*pi.f);
    #endif
		small = BBE_OLTN_2XF32(BBE_ABSN_2XF32(dx), tol);
		dx_corr = BBE_MOVN_2XF32T(BBE_CONSTN_2XF32(0), dx_corr, small);
		BBE_SVN_2XF32_IP(dx_corr, pY, 2 * BBE_SIMD_WIDTH);
	}

#if 1
	__Pragma("no_reorder")
	pX   = (const xb_vecN_2xf32*)x;
	pYrd = (const xb_vecN_2xf32*)y;
	pY   = (      xb_vecN_2xf32*)y;
	y0 = BBE_CONSTN_2XF32(0);
	bfirst=BBE_LTRN_2I(1);
	for (n = 0; n < (N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
	{
		xb_vecN_2xf32 sum;
		BBE_LVN_2XF32F_IP(sum, pYrd, 2 * BBE_SIMD_WIDTH,bfirst); bfirst=BBE_LTRN_2I(0);
		BBE_LVN_2XF32_IP(x0, pX  , 2 * BBE_SIMD_WIDTH);
		sum = BBE_ADDN_2XF32(sum, BBE_SELN_2XF32I(sum, BBE_CONSTN_2XF32(0), BBE_SELI_ROTATE_LEFT_2));
		sum = BBE_ADDN_2XF32(sum, BBE_SELN_2XF32I(sum, BBE_CONSTN_2XF32(0), BBE_SELI_ROTATE_LEFT_4));
		sum = BBE_ADDN_2XF32(sum, BBE_SELN_2XF32I(sum, BBE_CONSTN_2XF32(0), BBE_SELI_ROTATE_LEFT_8));
		y0 = BBE_ADDN_2XF32(sum, y0);
		x0 = BBE_ADDN_2XF32(x0, y0);
		BBE_SVN_2XF32_IP(x0, pY, 2 * BBE_SIMD_WIDTH);
		y0 = BBE_REPN_2XF32(y0, (BBE_SIMD_WIDTH/2 - 1));
	}
#else
	__Pragma("no_reorder")
	pYrd = (const xb_vecN_2xf32*)y;
	pY   = (      xb_vecN_2xf32*)y;
	for (n = 0; n < (N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
	{
		xb_vecN_2xf32 sum;
		BBE_LVN_2XF32_IP(sum, pYrd, 2 * BBE_SIMD_WIDTH);
		sum = BBE_ADDN_2XF32(sum, BBE_SELN_2XF32I(sum, BBE_CONSTN_2XF32(0), BBE_SELI_ROTATE_LEFT_2));
		sum = BBE_ADDN_2XF32(sum, BBE_SELN_2XF32I(sum, BBE_CONSTN_2XF32(0), BBE_SELI_ROTATE_LEFT_4));
		sum = BBE_ADDN_2XF32(sum, BBE_SELN_2XF32I(sum, BBE_CONSTN_2XF32(0), BBE_SELI_ROTATE_LEFT_8));
		BBE_SVN_2XF32_IP(sum, pY, 2 * BBE_SIMD_WIDTH);
	}
	pX   = (const xb_vecN_2xf32*)x;
	y0 = BBE_CONSTN_2XF32(0);
	pYrd = (const xb_vecN_2xf32*)y;
	pY   = (      xb_vecN_2xf32*)y;
	__Pragma("no_reorder")
	for (n = 0; n < (N >> (LOG2_BBE_SIMD_WIDTH - 1)); n++)
	{
		xb_vecN_2xf32 sum;
		BBE_LVN_2XF32_IP(sum, pYrd, 2 * BBE_SIMD_WIDTH);
		BBE_LVN_2XF32_IP(x0, pX  , 2 * BBE_SIMD_WIDTH);
		y0 = BBE_ADDN_2XF32(sum, y0);
		x0 = BBE_ADDN_2XF32(x0, y0);
		BBE_SVN_2XF32_IP(x0, pY, 2 * BBE_SIMD_WIDTH);
		y0 = BBE_REPN_2XF32(y0, (BBE_SIMD_WIDTH/2 - 1));
	}
#endif
}
#endif
