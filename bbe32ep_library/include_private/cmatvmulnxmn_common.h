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
#ifndef CMATVMULMXNXN_COMMON_H__
#define CMATVMULMXNXN_COMMON_H__

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matop.h"
#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif
/* M<8 && N>=8 && !(L&7) */
void cmatvmulnxmn_Mlt8_Ngte8_L8x( 
                   void    *          pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q );

/* M<8 && N>=8 && !(L&1) */
void cmatvmulnxmn_Mlt8_Ngte8_L2x( 
                   void    *          pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q );

/* M>=8 && N<8 && !(L&7) */
void cmatvmulnxmn_Mgte8_Nlt8_L8x( 
                   void    *          pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q );

/* M>=8 && N<8 && !(L&1) */
void cmatvmulnxmn_Mgte8_Nlt8_L2x( 
                   void    *          pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q );

/* M>=8 && N>=8 && !(L&7) */
void cmatvmulnxmn_Mgte8_Ngte8_L8x( 
                   void    *          pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q );

/* M>=8 && N>=8 && !(L&1) */
void cmatvmulnxmn_Mgte8_Ngte8_L2x( 
                   void    *          pScr,
                   complex_fract16 * restrict z, 
             const complex_fract16 * restrict x, 
             const complex_fract16 * restrict y, 
             int N, int M, int L, int Q );

/* return scratch size */
size_t cmatvmulnxmn_Mlt8_Ngte8_L8x_getScratchSize(int N, int M);
size_t cmatvmulnxmn_Mlt8_Ngte8_L2x_getScratchSize(int N, int M);
size_t cmatvmulnxmn_Mgte8_Nlt8_L8x_getScratchSize(int N, int M);
size_t cmatvmulnxmn_Mgte8_Nlt8_L2x_getScratchSize(int N, int M);
size_t cmatvmulnxmn_Mgte8_Ngte8_L8x_getScratchSize(int N, int M);
size_t cmatvmulnxmn_Mgte8_Ngte8_L2x_getScratchSize(int N, int M);
#ifdef __cplusplus
}
#endif

#endif
