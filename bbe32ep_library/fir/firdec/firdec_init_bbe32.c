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
    NatureDSP_Baseband library. FIR filters and Related Functions
    Decimating Block Complex FIR Filter
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
#include "firdec_common.h"

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


// generic allocation
static void _alloc_gen(tFilterLayout* pFltr,int M, int D)
{
    pFltr->delLength = (M + (BBE_SIMD_WIDTH / 2 - 1))&~(BBE_SIMD_WIDTH / 2 - 1);
    pFltr->coefNum   = M;
}

// simple initialization
static void _init_gen(int16_t* restrict coef, const int16_t* restrict h, int M, int D)
{
    int m,d,coefNum;
    coefNum=((M+D-1)/D)*D;
    for ( d=0; d<D; d++ )
    {
      for ( m=0; m<coefNum/D; m++ )
      {
        int ix = coefNum-1 - m*D - d;

        coef[d*coefNum/D+m] = ( ix < M ? h[ix] : 0 );
      }
    }
}

//for M=2,D=2
static void _alloc_d2_m2(tFilterLayout* pFltr,int M,int D)
{
    pFltr->delLength = M+(-M&((BBE_SIMD_WIDTH/2)-1));
    pFltr->coefNum   = BBE_SIMD_WIDTH;
}

static void _init_d2_m2 (int16_t* restrict coef, const int16_t* restrict h, int M, int D)
{
  int m;
  for ( m=0; m<BBE_SIMD_WIDTH/4; m++ )
  {
    coef[m*4+0] = h[1];
    coef[m*4+1] = 0;
    coef[m*4+2] = h[0];
    coef[m*4+3] = 0;
  }
}

// general D
static void _alloc_d2_mx(tFilterLayout* pFltr,int M, int D)
{
    pFltr->delLength = M + 16*D + (-M&((16*2)-1));
    pFltr->coefNum   = M;
}

static void _init_d2_mx (int16_t* restrict coef, const int16_t* restrict h, int M, int D)
{
    int s,d,m;
    for (s = 0; s<(M / BBE_SIMD_WIDTH); s++)
    {
      for ( d=0; d<D; d++ )
      {
          for (m = 0; m<BBE_SIMD_WIDTH / D; m++)
        {
              int ix = M - 1 - m*D - d - BBE_SIMD_WIDTH*s;

              coef[d*BBE_SIMD_WIDTH / D + m + s*BBE_SIMD_WIDTH] = (ix < M ? h[ix] : 0);
        }
      }
    }
}

/*
    For D=3
*/
static void _init_d3_m8 (int16_t* restrict coef, const int16_t* restrict h, int M, int D)
{
    coef[0] = h[7]; // Bank 0, ofs 0
    coef[1] = h[4]; // Bank 0, ofs 1
    coef[2] = h[6]; // Bank 1, ofs 0
    coef[3] = h[3]; // Bank 1, ofs 1
    coef[4] = h[5]; // Bank 2, ofs 0
    coef[5] = h[2]; // Bank 2, ofs 1
    coef[6] = h[1]; // Bank 0, ofs 2
    coef[7] = h[0]; // Bank 1, ofs 2
}

#if 1
// for M=16,D=3
static void _alloc_d3_m16(tFilterLayout* pFltr,int M, int D)
{
    pFltr->delLength = 3*BBE_SIMD_WIDTH/2;
    pFltr->coefNum   = 18;
}

static void _alloc_d3_m32(tFilterLayout* pFltr,int M, int D)
{
    pFltr->delLength = 6*BBE_SIMD_WIDTH/2;
    pFltr->coefNum   = 36;
}

static void _init_d3_m32(int16_t* restrict coef, const int16_t* restrict h, int M, int D)
{
    int m,d,coefNum;
    coefNum=36;
    for ( d=0; d<D; d++ )
    {
      for ( m=0; m<coefNum/D; m++ )
      {
        int ix = coefNum-1 - m*D - d;

        coef[d*coefNum/D+m] = ( ix < M ? h[ix] : 0 );
      }
    }
}
#endif

static void _alloc_d3_mx(tFilterLayout* pFltr,int M,int D)
{
    pFltr->delLength = M + 16 * D + (-M&((BBE_SIMD_WIDTH * 2) - 1));
    pFltr->coefNum = 3 * ((M / 3 + ((BBE_SIMD_WIDTH / 2) - 1))&(~((BBE_SIMD_WIDTH / 2) - 1)));
}

static void _init_d3_mx(int16_t* restrict coef, const int16_t* restrict h, int M, int D)
{
  int s,d,m;
  int coefNum = 3 * ((M / 3 + ((BBE_SIMD_WIDTH / 2) - 1))&(~((BBE_SIMD_WIDTH / 2) - 1)));
  int n = 24;

  NASSERT( D==3 );

  /*
  -------+-----+-----+--//--+-----+-----+-----+-----+-----+-----+-----+-----+
  bank 2 | M-2 |     |      |  22 |  19 |  16 |  13 |  10 |  7  |  4  |   1 |
  -------+-----+-----+--//--+-----+-----+-----+-----+-----+-----+-----+-----+
  bank 1 | M-1 |     |      |  23 |  20 |  17 |  14 |  11 |  8  |  5  |   2 |
  -------+-----+-----+--//--+-----+-----+-----+-----+-----+-----+-----+-----+
  */

  for ( s=0; s<(coefNum/n); s++)
  {
    for ( d=0; d<2; d++ )
    {
      for ( m=0; m<n/3; m++ )
      {
        int ix = coefNum -1 - m*D - d - n*s;
        coef[d*n / 3 + m + s*BBE_SIMD_WIDTH] = (ix < M ? h[ix] : 0);
      }
    }
  }

  /*
  -------+-----+-----+--//--+-----+-----+-----+-----+-----+-----+-----+-----+
  bank 0 | M-3 |     |      |  21 |  18 |  15 |  12 |  9  |  6  |  3  |  0  |
  -------+-----+-----+--//--+-----+-----+-----+-----+-----+-----+-----+-----+
  */
  for ( s=0; s<(coefNum/n); s++)
  {    
      for ( m=0; m<n/3; m++ )
      {
        int ix = coefNum -1 - m*D - d - n*s;
        coef[2 * coefNum / 3 + s*(BBE_SIMD_WIDTH/2)+m] = (ix < M ? h[ix] : 0);
      }
  }
}

static void _alloc_dx(tFilterLayout* pFltr,int M,int D)
{
    //pFltr->delLength = ( D*N + M ); // N<-present
    pFltr->delLength = ( D*128 + M );  // NSamples == 64
    pFltr->coefNum   = M;
}
static void _init_dx (int16_t* restrict coef, const int16_t* restrict h, int M, int D)
{
    int m;
    for ( m=0; m<M; m++ ) coef[m] =h[M-1 - m];
}

const tFirDecAlloc firdec_alloc_gen={_alloc_gen,_init_gen};
const tFirDecAlloc firdec_alloc_d2_m2={_alloc_d2_m2,_init_d2_m2};
const tFirDecAlloc firdec_alloc_d2_mx={_alloc_d2_mx,_init_d2_mx};
const tFirDecAlloc firdec_alloc_d3_m16={_alloc_d3_m16,_init_gen};     // not used
const tFirDecAlloc firdec_alloc_d3_m32={_alloc_d3_m32,_init_d3_m32};  // not used
const tFirDecAlloc firdec_alloc_d3_mx={_alloc_d3_mx,_init_d3_mx};
const tFirDecAlloc firdec_alloc_d3_m8 ={_alloc_gen,_init_d3_m8};
const tFirDecAlloc firdec_alloc_d3_m4 ={_alloc_gen,_init_dx};
const tFirDecAlloc firdec_alloc_dx    ={_alloc_dx,_init_dx};
