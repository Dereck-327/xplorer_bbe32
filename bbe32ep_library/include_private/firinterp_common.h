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

#ifndef FIRINTERP_COMMON_H__
#define FIRINTERP_COMMON_H__

#include "NatureDSP_types.h"
#include "common.h"

typedef void (*proc_fxdxn_t)(  void* handle,  
                                 int16_t *  restrict  y,
                           const int16_t *  restrict  x,
                           const int16_t *  restrict  coef,
                                 int16_t *  restrict  delayLine,
                                      int   M,
                                      int   N,
                                      int   D
                          );

typedef struct
{
  int coefNum;
  int delLength;
}
tFilterLayout;

typedef struct
{
    void (*fnalloc)(tFilterLayout* pFltr,int M,int D);
    void (*fninit )(int16_t* restrict coef, const int16_t* restrict h, int M, int D);
}
tFirInitAlloc_int;


typedef struct
{
  const tFirInitAlloc_int *   initAllocFxdxn;
  proc_fxdxn_t                   procFxdxn;
}
tFirFxdxns;

typedef struct
{
  int32_t         magic;      // Instance pointer validation number
  int             M;          // Length of subfilter
  int             D;          // Interpolation ratio
  proc_fxdxn_t    procFxdxn;  // Filter processing function
  int16_t *       coef;       // Filter coefficients, aligned
  int16_t *       delayLine;  // Delay line for samples, aligned  
}
firinterp_;

extern const tFirFxdxns interp_2d_2_8n  ;
extern const tFirFxdxns interp_2d_4_8n  ;
extern const tFirFxdxns interp_2d_8_8n  ;
extern const tFirFxdxns interp_2d_16_8n ;
extern const tFirFxdxns interp_2d_32_8n ;
extern const tFirFxdxns interp_2d_mx_8n ;
extern const tFirFxdxns interp_3d_2_8n  ;
extern const tFirFxdxns interp_3d_4_8n  ;
extern const tFirFxdxns interp_3d_8_8n  ;
extern const tFirFxdxns interp_3d_16_8n ;
extern const tFirFxdxns interp_3d_32_8n ;
extern const tFirFxdxns interp_3d_mx_8n ;
extern const tFirFxdxns interp_4d_2_8n  ;
extern const tFirFxdxns interp_4d_4_8n  ;
extern const tFirFxdxns interp_4d_8_8n  ;
extern const tFirFxdxns interp_4d_16_8n ;
extern const tFirFxdxns interp_4d_32_8n ;
extern const tFirFxdxns interp_4d_mx_8n ;
extern const tFirFxdxns interp_6d_2_8n  ;
extern const tFirFxdxns interp_6d_4_8n  ;
extern const tFirFxdxns interp_6d_8_8n  ;
extern const tFirFxdxns interp_6d_16_8n ;
extern const tFirFxdxns interp_6d_32_8n ;
extern const tFirFxdxns interp_6d_mx_8n ;
extern const tFirFxdxns interp_12d_2_8n ;
extern const tFirFxdxns interp_12d_4_8n ;
extern const tFirFxdxns interp_12d_8_8n ;
extern const tFirFxdxns interp_12d_16_8n;
extern const tFirFxdxns interp_12d_32_8n;
extern const tFirFxdxns interp_12d_mx_8n;
extern const tFirFxdxns interp_dx_mx_8n ;

extern const tFirInitAlloc_int firinterp_gen;
extern const tFirInitAlloc_int firinterp_dx;
extern const tFirInitAlloc_int firinterp_dx_m_32;
extern const tFirInitAlloc_int firinterp_d_6_12_mx;

/* interleaving functions */
void firinterp_12d_intlv(int16_t *  y,int   N);
void firinterp_6d_intlv (int16_t *  y,int   N);

#endif


