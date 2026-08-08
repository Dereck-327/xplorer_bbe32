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
#ifndef LTEPRS_COMMON_H__
#define LTEPRS_COMMON_H__
#include "NatureDSP_types.h"
#include "common.h"

typedef union { int16_t i[BBE_SIMD_WIDTH]; _vselN s; _vsaN v;} tlteprssel;


#define HAVE_LTEPRS (1 &(XCHAL_HAVE_BBEN_ALIGNING_LOAD_STORE && HAVE_LFSR && defined(BBE_ADDSAVSN) && defined(BBE_SHFLI_MMC2X4X4X2_M2_STEP_1)))

extern const uint32_t ALIGN(32) lteprs_poly[8][8];
extern const tlteprssel lteprs1_tbl[2]; 
extern const tlteprssel lteprs2_tbl[2]; 
extern const tlteprssel lteprs4_tbl[2]; 
extern const tlteprssel lteprs6_tbl[4]; 
extern const tlteprssel lteprs8_tbl[1]; 
extern const tlteprssel lteprs10_tbl[4];

// update LFSR state: takes c0,c1 and outputs combined value in 4 first elements of c0
#define UPDATE_LFSR(c0,c1,N)                                         \
{                                                                    \
    vsaN shr,shl;                                                    \
    vselN rot;                                                       \
    xb_vecNx16 x0,x1;                                                \
    /* combine 15:8 of c0,c1 */                                      \
    c0=BBE_SELNX16I(c1,c0,BBE_SELI_EXTRACT_HI_HALVES);               \
    /*c0=BBE_SELNX16I(c1,c0,BBE_SELI_INTERLEAVE_4_HI);*/             \
    /* first, rotate by 16-bit elements */                           \
    x0=BBE_SEQNX16();                                                \
    x1=BBE_MOVVA16(N>>4);                                            \
    x0=BBE_ADDNX16(x0,x1);                                           \
    rot=BBE_MOVVSELNX16(x0,0);                                       \
    c0=BBE_SELNX16(c0,c0,rot);                                       \
    /* next make bitwise rotation */                                 \
    N&=15;                                                           \
    shr=BBE_MOVVSA32(N);                                             \
    shl=BBE_MOVVSA32(16-N);                                          \
    x0=BBE_SRLNX16(c0,shr);                                          \
    x1=BBE_SLLNX16(c0,shl);                                          \
    x1=BBE_SELNX16I(x1,x1,BBE_SELI_ROTATE_RIGHT_1);                  \
    c0=BBE_ORNX16(x0,x1);                                            \
    /* finally, combine 0,1,8,9 */                                   \
    c0 = BBE_SHFLNX16I(c0,BBE_SHFLI_MMC1X4X4X4_M2_STEP_1_HIGH_HALF); \
    /*c0 = BBE_SHFLNX16I(c0,BBE_SHFLI_MMC2X4X4X2_M2_STEP_1);*/       \
}

// the same as above, but allows to rotate more than by 128 bits 
#define UPDATE_LFSRL(c0,c1,N)                                        \
{                                                                    \
    vsaN shr,shl;                                                    \
    vselN rot;                                                       \
    xb_vecNx16 x0,x1;                                                \
    /* first, rotate by 16-bit elements */                           \
    x0=BBE_SEQNX16();                                                \
    x1=BBE_MOVVA16(N>>4);                                            \
    x0=BBE_ADDNX16(x0,x1);                                           \
    rot=BBE_MOVVSELNX16(x0,0);                                       \
    c0=BBE_SELNX16(c0,c0,rot);                                       \
    c1=BBE_SELNX16(c1,c1,rot);                                       \
    /* combine c0_0,c0_1,...,c1_0,c1_1, */                           \
    c0=BBE_SELNX16I(c1,c0,BBE_SELI_EXTRACT_LO_HALVES);               \
    /*c0=BBE_SELNX16I(c1,c0,	BBE_SELI_INTERLEAVE_4_LO);*/         \
    /* next make bitwise rotation */                                 \
    N&=15;                                                           \
    shr=BBE_MOVVSA32(N);                                             \
    shl=BBE_MOVVSA32(16-N);                                          \
    x0=BBE_SRLNX16(c0,shr);                                          \
    x1=BBE_SLLNX16(c0,shl);                                          \
    x1=BBE_SELNX16I(x1,x1,BBE_SELI_ROTATE_RIGHT_1);                  \
    c0=BBE_ORNX16(x0,x1);                                            \
    /* finally, combine 0,1,16,17 */                                 \
    c0 = BBE_SHFLNX16I(c0,BBE_SHFLI_MMC1X4X4X4_M2_STEP_1_HIGH_HALF); \
    /*c0 = BBE_SHFLNX16I(c0,BBE_SHFLI_MMC2X4X4X2_M2_STEP_1);*/       \
}

#define UPDATE_LFSRL2(c0, c1, c0_prev, c1_prev, N, start_lane)                  \
{                                                                               \
    vsaN shr, shl;                                                              \
    vselN rot;                                                                  \
    xb_vecNx16 x0, x1;                                                          \
    /* first, rotate by 16-bit elements */                                      \
    x0 = BBE_SEQNX16();                                                         \
    /* start state in lanes start_lane and start_lane+1, plus (N >> 4) lanes */ \
    x1 = BBE_MOVVA16((N >> 4) + start_lane);                                    \
    x0 = BBE_ADDNX16(x0, x1);                                                   \
    rot = BBE_MOVVSELNX16(x0, 0);                                               \
    c0 = BBE_SELNX16(c0, c0_prev, rot);                                         \
    c1 = BBE_SELNX16(c1, c1_prev, rot);                                         \
    /* combine c0_0,c0_1,...,c1_0,c1_1, */                                      \
    c0=BBE_SELNX16I(c1,c0,BBE_SELI_EXTRACT_LO_HALVES);                          \
    /* next make bitwise rotation */                                            \
    N &= 15;                                                                    \
    shr = BBE_MOVVSA32(N);                                                      \
    shl = BBE_MOVVSA32(16 - N);                                                 \
    x0 = BBE_SRLNX16(c0, shr);                                                  \
    x1 = BBE_SLLNX16(c0, shl);                                                  \
    x1 = BBE_SELNX16I(x1, x1, BBE_SELI_ROTATE_RIGHT_1);                         \
    c0 = BBE_ORNX16(x0, x1);                                                    \
    /* finally, combine 0,1,16,17 */                                            \
    c0 = BBE_SHFLNX16I(c0,BBE_SHFLI_MMC1X4X4X4_M2_STEP_1_HIGH_HALF);            \
}

// apply BBE_BMUL32A instructions N times
#define BMUL32N(c,N)                                   \
{                                                      \
    xb_vecNx16 t=c;                                    \
    if(N>=2) { c=BBE_BMUL32A(c,0);t=BBE_BMUL32A(t,0);} \
    if(N>=4) { c=BBE_BMUL32A(c,0);t=BBE_BMUL32A(t,0);} \
    if(N>=6) { c=BBE_BMUL32A(c,0);t=BBE_BMUL32A(t,0);} \
    if(N>=8) { c=BBE_BMUL32A(c,0);t=BBE_BMUL32A(t,0);} \
    if(N>=10){ c=BBE_BMUL32A(c,0);t=BBE_BMUL32A(t,0);} \
    if(N>=12){ c=BBE_BMUL32A(c,0);t=BBE_BMUL32A(t,0);} \
    if(N>=14){ c=BBE_BMUL32A(c,0);t=BBE_BMUL32A(t,0);} \
    if(N>=16){ c=BBE_BMUL32A(c,0);t=BBE_BMUL32A(t,0);} \
    c=BBE_SELNX16I(t,c,BBE_SELI_INTERLEAVE_2_HI);      \
}

#endif
