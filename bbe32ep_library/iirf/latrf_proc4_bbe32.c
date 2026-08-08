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
    Lattice real block IIR, floating point
    C code optimized for BBE32
    Integrit, 2006-2016
*/

/* M=3...4 */

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Filter instance and processing functions. */
#include "latrf_common.h"

#if HAVE_VFPU

/*-------------------------------------------------------------------------
   Pass a block of input samples to the filter 
   Input:
     latr   Filter handle
     N      Input/output signal chunk size, in samples
     x[N]   Input samples
   Output:
     r[N]   Output samples
   Restrictions:
     None 
-------------------------------------------------------------------------*/
void latrf_process4( latrf_t   *          latr,
                     float32_t * restrict r,
               const float32_t *          x,
                     int                  N )
{
          float32_t * restrict delLine;
    const float32_t * restrict coef;

    float32_t d0,d1,d2,d3,c0,c1,c2,c3;
    float32_t t0;
    float32_t scale;

    int n;

    NASSERT((latr->M==3)||(latr->M==4));
    delLine = latr->delLine;
    coef    = latr->coef;
    scale   = latr->gain;
    d0=delLine[0];
    d1=delLine[1];
    d2=delLine[2];
    d3=delLine[3];
    c0=coef[0];
    c1=coef[1];
    c2=coef[2];
    c3=coef[3];

    for ( n=0; n<N; n++ )
    {
        t0  = x[n]*scale;
        t0 -= d3 * c3;
        t0 -= d2 * c2;
        d3  = d2 + t0 * c2 ;
        t0 -= d1 * c1;
        d2  = d1 + t0 * c1 ;
        t0 -= d0 * c0 ;
        d1  = d0 + t0 * c0 ;
        d0  = t0;
        r[n] = t0;
    }
    delLine[0]=d0;
    delLine[1]=d1;
    delLine[2]=d2;
    delLine[3]=d3;
} /* latrf_process4() */

#endif /* if HAVE_VFPU */
