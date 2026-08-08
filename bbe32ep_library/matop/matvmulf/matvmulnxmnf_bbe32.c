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
  NatureDSP_Baseband library. Matrix Operations
    Real Matrix-Matrix/Matrix-Vector Multiply
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"

#if !(HAVE_VFPU)
DISCARD_FUN(void, matvmulnxmnf,( void * pScr,
                            float32_t * restrict z,
                      const float32_t * restrict x,
                      const float32_t * restrict y,
                      int N, int M, int L ))
#else

#ifndef BBE_MOVN_2XF32_FROMN_4XCF32
#define BBE_MOVN_2XF32_FROMN_4XCF32(a) BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(a))
#endif
#ifndef BBE_SELN_2XF32
#define BBE_SELN_2XF32(a, b, c) (BBE_MOVN_2XF32_FROMNX16(BBE_SELNX16(BBE_MOVNX16_FROMN_2XF32(a), BBE_MOVNX16_FROMN_2XF32(b), (c) )))
#endif

/* N is multiple of 8 */
static void matvmulnxmnf_N8 ( float32_t * restrict z,
                        const float32_t * restrict x,
                        const float32_t * restrict y,
                        int N, int M, int L )
{
  int m, n, l;

  const xb_vecN_2xf32 * restrict px0;
  const xb_vecN_2xf32 * restrict px1;
  const xb_vecN_2xf32 * restrict px2;
  const xb_vecN_2xf32 * restrict px3;
  const xb_vecN_2xf32 * restrict px4;
  const xb_vecN_2xf32 * restrict px5;
  const xb_vecN_2xf32 * restrict px6;
  const xb_vecN_2xf32 * restrict px7;
  const xb_vecN_2xf32 * restrict py;
        xb_vecN_2xf32 * restrict pz;

  xb_vecN_2xf32 X0, X1, X2, X3, X4, X5, X6, X7;
  xb_vecN_2xf32 Y;
  xb_vecN_2xf32 Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7;
  valign al_z;
  int zstride;

  pz  = (      xb_vecN_2xf32 *)(z);
  px0 = (const xb_vecN_2xf32 *)(x+0*N);
  px1 = (const xb_vecN_2xf32 *)(x+1*N);
  px2 = (const xb_vecN_2xf32 *)(x+2*N);
  px3 = (const xb_vecN_2xf32 *)(x+3*N);
  px4 = (const xb_vecN_2xf32 *)(x+4*N);
  px5 = (const xb_vecN_2xf32 *)(x+5*N);
  px6 = (const xb_vecN_2xf32 *)(x+6*N);
  px7 = (const xb_vecN_2xf32 *)(x+7*N);

  al_z = BBE_ZALIGN();
  zstride = (M==4) ? 4*sizeof(float32_t) : 8*sizeof(float32_t);

  __Pragma("loop_count min=1");
  for (l=0; l<L; l++)
  {
      /* Compute vector by 8 values per iteration */
      for (m=0; m<(M>>3); m++)
      {
          py  = (const xb_vecN_2xf32 *)(y+l*N);
          Z0 = Z2 = Z4 = Z6 = BBE_ZERON_2XF32();

          __Pragma("loop_count min=1");
          for (n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-1)); n++)
          {
              /* Load matrix X */
              BBE_LVN_2XF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
              BBE_LVN_2XF32_IP(X1, px1, 2*BBE_SIMD_WIDTH);
              BBE_LVN_2XF32_IP(X2, px2, 2*BBE_SIMD_WIDTH);
              BBE_LVN_2XF32_IP(X3, px3, 2*BBE_SIMD_WIDTH);
              BBE_LVN_2XF32_IP(X4, px4, 2*BBE_SIMD_WIDTH);
              BBE_LVN_2XF32_IP(X5, px5, 2*BBE_SIMD_WIDTH);
              BBE_LVN_2XF32_IP(X6, px6, 2*BBE_SIMD_WIDTH);
              BBE_LVN_2XF32_IP(X7, px7, 2*BBE_SIMD_WIDTH);
              /* Load vector Y */
              BBE_LVN_2XF32_IP(Y, py, 2*BBE_SIMD_WIDTH);
              Y = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_MMC1X4X4X4_M2_STEP_1_LOW_HALF);

              BBE_DSELN_2XF32I(X1, X0, X1, X0, BBE_DSELI_INTERLEAVE_2);
              BBE_DSELN_2XF32I(X3, X2, X3, X2, BBE_DSELI_INTERLEAVE_2);
              BBE_DSELN_2XF32I(X5, X4, X5, X4, BBE_DSELI_INTERLEAVE_2);
              BBE_DSELN_2XF32I(X7, X6, X7, X6, BBE_DSELI_INTERLEAVE_2);

              BBE_MULMASN_2XF32(Z0, X0, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z0, X1, Y, 0, 0xD);
              BBE_MULMASN_2XF32(Z2, X2, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z2, X3, Y, 0, 0xD);
              BBE_MULMASN_2XF32(Z4, X4, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z4, X5, Y, 0, 0xD);
              BBE_MULMASN_2XF32(Z6, X6, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z6, X7, Y, 0, 0xD);
          }

          BBE_DSELN_2XF32I(Z2, Z0, Z2, Z0, BBE_DSELI_INTERLEAVE_4);
          BBE_DSELN_2XF32I(Z6, Z4, Z6, Z4, BBE_DSELI_INTERLEAVE_4);
          Z0 = BBE_ADDN_2XF32(Z0, Z2);
          Z4 = BBE_ADDN_2XF32(Z4, Z6);

          Z1 = BBE_SELN_2XF32I(Z4, Z0, BBE_SELI_EXTRACT_LO_HALVES);
          Z2 = BBE_SELN_2XF32I(Z4, Z0, BBE_SELI_EXTRACT_HI_HALVES);
          Z0 = BBE_ADDN_2XF32(Z1, Z2);

          /* Save outputs */
          BBE_SVN_2XF32_IP(Z0, pz, 2*BBE_SIMD_WIDTH);

          /* Jump to next 8 rows */
          px0 = (const xb_vecN_2xf32 *)((float32_t*)px7);
          px1 = (const xb_vecN_2xf32 *)((float32_t*)px0+N);
          px2 = (const xb_vecN_2xf32 *)((float32_t*)px1+N);
          px3 = (const xb_vecN_2xf32 *)((float32_t*)px2+N);
          px4 = (const xb_vecN_2xf32 *)((float32_t*)px3+N);
          px5 = (const xb_vecN_2xf32 *)((float32_t*)px4+N);
          px6 = (const xb_vecN_2xf32 *)((float32_t*)px5+N);
          px7 = (const xb_vecN_2xf32 *)((float32_t*)px6+N);
      }
      /* Process last M%8 values */
      if (M&4)
      {
          py  = (const xb_vecN_2xf32 *)(y+l*N);
          Z0 = Z1 = Z2 = Z3 = Z4 = Z5 = Z6 = Z7 = BBE_ZERON_2XF32();

          __Pragma("loop_count min=1");
          for (n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-1)); n++)
          {
              /* Load matrix X */
              BBE_LVN_2XF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
              BBE_LVN_2XF32_IP(X1, px1, 2*BBE_SIMD_WIDTH);
              BBE_LVN_2XF32_IP(X2, px2, 2*BBE_SIMD_WIDTH);
              BBE_LVN_2XF32_IP(X3, px3, 2*BBE_SIMD_WIDTH);
              /* Load vector Y */
              BBE_LVN_2XF32_IP(Y, py, 2*BBE_SIMD_WIDTH);
              Y = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_MMC1X4X4X4_M2_STEP_1_LOW_HALF);

              BBE_DSELN_2XF32I(X1, X0, X1, X0, BBE_DSELI_INTERLEAVE_2);
              BBE_DSELN_2XF32I(X3, X2, X3, X2, BBE_DSELI_INTERLEAVE_2);

              BBE_MULMASN_2XF32(Z0, X0, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z1, X1, Y, 0, 0xD);
              BBE_MULMASN_2XF32(Z2, X2, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z3, X3, Y, 0, 0xD);
          }

          Z0 = BBE_ADDN_2XF32(Z0, Z1);
          Z2 = BBE_ADDN_2XF32(Z2, Z3);

          BBE_DSELN_2XF32I(Z2, Z0, Z2, Z0, BBE_DSELI_INTERLEAVE_4);
          Z0 = BBE_ADDN_2XF32(Z0, Z2);

          Z1 = BBE_SHFLN_2XF32I(Z0, BBE_SHFLI_SWAP_8);
          Z0 = BBE_ADDN_2XF32(Z0, Z1);

          /* Save outputs */
          BBE_SAVN_2XF32_XP(Z0, al_z, pz, zstride);
          BBE_SAN_2XF32POS_FP(al_z, pz);

          /* Jump to next 4 rows */
          px0 = (const xb_vecN_2xf32 *)((float32_t*)px3);
          px1 = (const xb_vecN_2xf32 *)((float32_t*)px0+N);
          px2 = (const xb_vecN_2xf32 *)((float32_t*)px1+N);
          px3 = (const xb_vecN_2xf32 *)((float32_t*)px2+N);
          px4 = (const xb_vecN_2xf32 *)((float32_t*)px3+N);
          px5 = (const xb_vecN_2xf32 *)((float32_t*)px4+N);
          px6 = (const xb_vecN_2xf32 *)((float32_t*)px5+N);
          px7 = (const xb_vecN_2xf32 *)((float32_t*)px6+N);
      }
  }
} /* matvmulnxmnf_N8() */

/* N is multiple of 4 */
static void matvmulnxmnf_N4 ( float32_t * restrict z,
                        const float32_t * restrict x,
                        const float32_t * restrict y,
                        int N, int M, int L )
{
#if 1
  int m, n, l;

  const xb_vecN_2xf32 * restrict px0;
  const xb_vecN_2xf32 * restrict px1;
  const xb_vecN_2xf32 * restrict px2;
  const xb_vecN_2xf32 * restrict px3;
  const xb_vecN_2xf32 * restrict px4;
  const xb_vecN_2xf32 * restrict px5;
  const xb_vecN_2xf32 * restrict px6;
  const xb_vecN_2xf32 * restrict px7;
  const xb_vecN_2xf32 * restrict py;
        xb_vecN_2xf32 * restrict pz;

  xb_vecN_2xf32 X0, X1, X2, X3, X4, X5, X6, X7;
  xb_vecN_2xf32 Y;
  xb_vecN_2xf32 Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7;
  valign al_px1, al_px3, al_px5, al_px7, al_z;
  vboolN_2 mask_4n;
  int zstride, N_;

  pz  = (      xb_vecN_2xf32 *)(z);
  px0 = (const xb_vecN_2xf32 *)(x+0*N);
  px1 = (const xb_vecN_2xf32 *)(x+1*N);
  px2 = (const xb_vecN_2xf32 *)(x+2*N);
  px3 = (const xb_vecN_2xf32 *)(x+3*N);
  px4 = (const xb_vecN_2xf32 *)(x+4*N);
  px5 = (const xb_vecN_2xf32 *)(x+5*N);
  px6 = (const xb_vecN_2xf32 *)(x+6*N);
  px7 = (const xb_vecN_2xf32 *)(x+7*N);

  mask_4n = BBE_LTRN_2I(4);
  zstride = (M==4) ? 4*sizeof(float32_t) : 8*sizeof(float32_t);
  N_ = N+4;

  __Pragma("loop_count min=1");
  for (l=0; l<L; l++)
  {
      /* Compute vector by 8 values per iteration */
      for (m=0; m<(M>>3); m++)
      {
          py  = (const xb_vecN_2xf32 *)(y+l*N_);
          al_px1 = BBE_LAN_2XF32_PP(px1);
          al_px3 = BBE_LAN_2XF32_PP(px3);
          al_px5 = BBE_LAN_2XF32_PP(px5);
          al_px7 = BBE_LAN_2XF32_PP(px7);
          Z0 = Z2 = Z4 = Z6 = BBE_ZERON_2XF32();

          __Pragma("loop_count min=1");
          for (n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-1)); n++)
          {
              /* Load matrix X */
              BBE_LVN_2XF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
              BBE_LAN_2XF32_IP(X1, al_px1, px1);
              BBE_LVN_2XF32_IP(X2, px2, 2*BBE_SIMD_WIDTH);
              BBE_LAN_2XF32_IP(X3, al_px3, px3);
              BBE_LVN_2XF32_IP(X4, px4, 2*BBE_SIMD_WIDTH);
              BBE_LAN_2XF32_IP(X5, al_px5, px5);
              BBE_LVN_2XF32_IP(X6, px6, 2*BBE_SIMD_WIDTH);
              BBE_LAN_2XF32_IP(X7, al_px7, px7);
              /* Load vector Y */
              BBE_LVN_2XF32_IP(Y, py, 2*BBE_SIMD_WIDTH);
              Y = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_MMC1X4X4X4_M2_STEP_1_LOW_HALF);

              BBE_DSELN_2XF32I(X1, X0, X1, X0, BBE_DSELI_INTERLEAVE_2);
              BBE_DSELN_2XF32I(X3, X2, X3, X2, BBE_DSELI_INTERLEAVE_2);
              BBE_DSELN_2XF32I(X5, X4, X5, X4, BBE_DSELI_INTERLEAVE_2);
              BBE_DSELN_2XF32I(X7, X6, X7, X6, BBE_DSELI_INTERLEAVE_2);

              BBE_MULMASN_2XF32(Z0, X0, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z0, X1, Y, 0, 0xD);
              BBE_MULMASN_2XF32(Z2, X2, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z2, X3, Y, 0, 0xD);
              BBE_MULMASN_2XF32(Z4, X4, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z4, X5, Y, 0, 0xD);
              BBE_MULMASN_2XF32(Z6, X6, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z6, X7, Y, 0, 0xD);
          }
          /* Process last 4 values */
          {
              /* Load matrix X */
              BBE_LVN_2XF32T_IP(X0, px0, 2*BBE_SIMD_WIDTH, mask_4n);
              BBE_LAVN_2XF32_XP(X1, al_px1, px1, 4*sizeof(float32_t));
              BBE_LVN_2XF32T_IP(X2, px2, 2*BBE_SIMD_WIDTH, mask_4n);
              BBE_LAVN_2XF32_XP(X3, al_px3, px3, 4*sizeof(float32_t));
              BBE_LVN_2XF32T_IP(X4, px4, 2*BBE_SIMD_WIDTH, mask_4n);
              BBE_LAVN_2XF32_XP(X5, al_px5, px5, 4*sizeof(float32_t));
              BBE_LVN_2XF32T_IP(X6, px6, 2*BBE_SIMD_WIDTH, mask_4n);
              BBE_LAVN_2XF32_XP(X7, al_px7, px7, 4*sizeof(float32_t));
              /* Load vector Y */
              BBE_LVN_2XF32T_IP(Y, py, 2*BBE_SIMD_WIDTH, mask_4n);
              Y = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_MMC1X4X4X4_M2_STEP_1_LOW_HALF);

              BBE_DSELN_2XF32I(X1, X0, X1, X0, BBE_DSELI_INTERLEAVE_2);
              BBE_DSELN_2XF32I(X3, X2, X3, X2, BBE_DSELI_INTERLEAVE_2);
              BBE_DSELN_2XF32I(X5, X4, X5, X4, BBE_DSELI_INTERLEAVE_2);
              BBE_DSELN_2XF32I(X7, X6, X7, X6, BBE_DSELI_INTERLEAVE_2);

              BBE_MULMASN_2XF32(Z0, X0, Y, 0, 0x8);
              //BBE_MULMASN_2XF32(Z0, X1, Y, 0, 0xD);
              BBE_MULMASN_2XF32(Z2, X2, Y, 0, 0x8);
              //BBE_MULMASN_2XF32(Z2, X3, Y, 0, 0xD);
              BBE_MULMASN_2XF32(Z4, X4, Y, 0, 0x8);
              //BBE_MULMASN_2XF32(Z4, X5, Y, 0, 0xD);
              BBE_MULMASN_2XF32(Z6, X6, Y, 0, 0x8);
              //BBE_MULMASN_2XF32(Z6, X7, Y, 0, 0xD);
          }

          BBE_DSELN_2XF32I(Z2, Z0, Z2, Z0, BBE_DSELI_INTERLEAVE_4);
          BBE_DSELN_2XF32I(Z6, Z4, Z6, Z4, BBE_DSELI_INTERLEAVE_4);
          Z0 = BBE_ADDN_2XF32(Z0, Z2);
          Z4 = BBE_ADDN_2XF32(Z4, Z6);

          Z1 = BBE_SELN_2XF32I(Z4, Z0, BBE_SELI_EXTRACT_LO_HALVES);
          Z2 = BBE_SELN_2XF32I(Z4, Z0, BBE_SELI_EXTRACT_HI_HALVES);
          Z0 = BBE_ADDN_2XF32(Z1, Z2);

          /* Save outputs */
          BBE_SVN_2XF32_IP(Z0, pz, 2*BBE_SIMD_WIDTH);

          /* Jump to next 8 rows */
          px0 = (const xb_vecN_2xf32 *)((float32_t*)px7);
          px1 = (const xb_vecN_2xf32 *)((float32_t*)px0+N);
          px2 = (const xb_vecN_2xf32 *)((float32_t*)px1+N);
          px3 = (const xb_vecN_2xf32 *)((float32_t*)px2+N);
          px4 = (const xb_vecN_2xf32 *)((float32_t*)px3+N);
          px5 = (const xb_vecN_2xf32 *)((float32_t*)px4+N);
          px6 = (const xb_vecN_2xf32 *)((float32_t*)px5+N);
          px7 = (const xb_vecN_2xf32 *)((float32_t*)px6+N);
      }
      /* Process last M%8 values */
      if (M&4)
      {
          py  = (const xb_vecN_2xf32 *)(y+l*N_);
          Z0 = Z1 = Z2 = Z3 = Z4 = Z5 = Z6 = Z7 = BBE_ZERON_2XF32();
          al_px1 = BBE_LAN_2XF32_PP(px1);
          al_px3 = BBE_LAN_2XF32_PP(px3);

          __Pragma("loop_count min=1");
          for (n=0; n<(N>>(LOG2_BBE_SIMD_WIDTH-1)); n++)
          {
              /* Load matrix X */
              BBE_LVN_2XF32_IP(X0, px0, 2*BBE_SIMD_WIDTH);
              BBE_LAN_2XF32_IP(X1, al_px1, px1);
              BBE_LVN_2XF32_IP(X2, px2, 2*BBE_SIMD_WIDTH);
              BBE_LAN_2XF32_IP(X3, al_px3, px3);
              /* Load vector Y */
              BBE_LVN_2XF32_IP(Y, py, 2*BBE_SIMD_WIDTH);
              Y = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_MMC1X4X4X4_M2_STEP_1_LOW_HALF);

              BBE_DSELN_2XF32I(X1, X0, X1, X0, BBE_DSELI_INTERLEAVE_2);
              BBE_DSELN_2XF32I(X3, X2, X3, X2, BBE_DSELI_INTERLEAVE_2);

              BBE_MULMASN_2XF32(Z0, X0, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z1, X1, Y, 0, 0xD);
              BBE_MULMASN_2XF32(Z2, X2, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z3, X3, Y, 0, 0xD);
          }
          /* Process last 4 values */
          {
              /* Load matrix X */
              BBE_LVN_2XF32T_IP(X0, px0, 2*BBE_SIMD_WIDTH, mask_4n);
              BBE_LAVN_2XF32_XP(X1, al_px1, px1, 4*sizeof(float32_t));
              BBE_LVN_2XF32T_IP(X2, px2, 2*BBE_SIMD_WIDTH, mask_4n);
              BBE_LAVN_2XF32_XP(X3, al_px3, px3, 4*sizeof(float32_t));
              /* Load vector Y */
              BBE_LVN_2XF32T_IP(Y, py, 2*BBE_SIMD_WIDTH, mask_4n);
              Y = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_MMC1X4X4X4_M2_STEP_1_LOW_HALF);

              BBE_DSELN_2XF32I(X1, X0, X1, X0, BBE_DSELI_INTERLEAVE_2);
              BBE_DSELN_2XF32I(X3, X2, X3, X2, BBE_DSELI_INTERLEAVE_2);

              BBE_MULMASN_2XF32(Z0, X0, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z1, X1, Y, 0, 0xD);
              BBE_MULMASN_2XF32(Z2, X2, Y, 0, 0x8);
              BBE_MULMASN_2XF32(Z3, X3, Y, 0, 0xD);
          }

          Z0 = BBE_ADDN_2XF32(Z0, Z1);
          Z2 = BBE_ADDN_2XF32(Z2, Z3);

          BBE_DSELN_2XF32I(Z2, Z0, Z2, Z0, BBE_DSELI_INTERLEAVE_4);
          Z0 = BBE_ADDN_2XF32(Z0, Z2);

          Z1 = BBE_SHFLN_2XF32I(Z0, BBE_SHFLI_SWAP_8);
          Z0 = BBE_ADDN_2XF32(Z0, Z1);

          /* Save outputs */
          al_z = BBE_ZALIGN();
          BBE_SAVN_2XF32_XP(Z0, al_z, pz, zstride);
          BBE_SAN_2XF32POS_FP(al_z, pz);

          /* Jump to next 4 rows */
          px0 = (const xb_vecN_2xf32 *)((float32_t*)px3);
          px1 = (const xb_vecN_2xf32 *)((float32_t*)px0+N);
          px2 = (const xb_vecN_2xf32 *)((float32_t*)px1+N);
          px3 = (const xb_vecN_2xf32 *)((float32_t*)px2+N);
          px4 = (const xb_vecN_2xf32 *)((float32_t*)px3+N);
          px5 = (const xb_vecN_2xf32 *)((float32_t*)px4+N);
          px6 = (const xb_vecN_2xf32 *)((float32_t*)px5+N);
          px7 = (const xb_vecN_2xf32 *)((float32_t*)px6+N);
      }
  }
#else
  int m, n, l;
  int Sx = M*N;
  int Sy=(N>4)?((N+7)&~7):N;
  int Sz=(M>4)?((M+7)&~7):M;

  // block order, real data
  for (l=0; l<L; l++,x+=Sx,y+=Sy,z+=Sz)
  {
      for (m=0; m<M; m++)
      {
          float32_t A=0;

          for (n=0; n<N; n++) A+=x[m*N+n]*y[n];

          z[m]=A;
      }
  }
#endif
} /* matvmulnxmnf_N4() */

/* N=4 */
static void matvmul4xmnf ( float32_t * restrict z,
                     const float32_t * restrict x,
                     const float32_t * restrict y,
                           int M, int L )
{
#if 1
  const xb_vecN_2xf32 * restrict px;
  const xb_vecN_2xf32 * restrict py;
        xb_vecN_2xf32 * restrict pz;
  int l, m;
  int zstride;

  xb_vecN_2xf32 X00, X01, X10, X11, Y, Z;
  xb_vecN_2xf32 x0, x1, x2, x3, y01, y23, z0, z1;
  valign al_y, al_z;

  px = (const xb_vecN_2xf32 *)x;
  py = (const xb_vecN_2xf32 *)y;
  pz = (      xb_vecN_2xf32 *)z;
  al_y = BBE_LAN_2XF32_PP(py);
  al_z = BBE_ZALIGN();
  zstride = (M==4) ? 4*sizeof(float32_t) : 8*sizeof(float32_t);

  for (l=0; l<L; l++)
  {
      BBE_LAVN_2XF32_XP(Y, al_y, py, 4*sizeof(float32_t));
      y01 = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_REP_0X4);
      y23 = BBE_SHFLN_2XF32I(Y, BBE_SHFLI_REP_1X4);

      for (m=0; m<(M>>3); m++)
      {
          /* Load input matrices X and Y */
          BBE_LVN_2XF32_IP(X00, px, 2 * BBE_SIMD_WIDTH);
          BBE_LVN_2XF32_IP(X01, px, 2 * BBE_SIMD_WIDTH);
          BBE_LVN_2XF32_IP(X10, px, 2 * BBE_SIMD_WIDTH);
          BBE_LVN_2XF32_IP(X11, px, 2 * BBE_SIMD_WIDTH);

          BBE_DSELN_2XF32I(X01, X00, X01, X00, BBE_DSELI_DEINTERLEAVE_2);
          BBE_DSELN_2XF32I(X11, X10, X11, X10, BBE_DSELI_DEINTERLEAVE_2);
          BBE_DSELN_2XF32I(x2 , x0 , X10, X00, BBE_DSELI_DEINTERLEAVE_2);
          BBE_DSELN_2XF32I(x3 , x1 , X11, X01, BBE_DSELI_DEINTERLEAVE_2);

          /* Multiply input matrices X and Y */
          z0 = BBE_MULMN_2XF32( x0, y01, 0, 0x8);
          z1 = BBE_MULMN_2XF32( x1, y01, 0, 0xD);
          BBE_MULMASN_2XF32(z0, x2, y23, 0, 0x8);
          BBE_MULMASN_2XF32(z1, x3, y23, 0, 0xD);

          /* Save results */
          Z = BBE_ADDN_2XF32(z0, z1);
          BBE_SVN_2XF32_IP(Z, pz, 2 * BBE_SIMD_WIDTH);
      }
      if (M&4)
      {
          /* Load input matrix X */
          BBE_LVN_2XF32_IP(X00, px, 2 * BBE_SIMD_WIDTH);
          BBE_LVN_2XF32_IP(X01, px, 2 * BBE_SIMD_WIDTH);

          BBE_DSELN_2XF32I(X01, X00, X01, X00, BBE_DSELI_DEINTERLEAVE_2);
          BBE_DSELN_2XF32I(x2 , x0 , X00, X00, BBE_DSELI_DEINTERLEAVE_2);
          BBE_DSELN_2XF32I(x3 , x1 , X01, X01, BBE_DSELI_DEINTERLEAVE_2);

          /* Multiply input matrices X and Y */
          z0 = BBE_MULMN_2XF32( x0, y01, 0, 0x8);
          z1 = BBE_MULMN_2XF32( x1, y01, 0, 0xD);
          BBE_MULMASN_2XF32(z0, x2, y23, 0, 0x8);
          BBE_MULMASN_2XF32(z1, x3, y23, 0, 0xD);

          /* Save results */
          Z = BBE_ADDN_2XF32(z0, z1);
          /* Save outputs */
          BBE_SAVN_2XF32_XP(Z, al_z, pz, zstride);
          BBE_SAN_2XF32POS_FP(al_z, pz);
      }
  }
#else
  int m, n, l;
  int Sx=M*4;
  int Sy=4;
  int Sz=(M>4)?((M+7)&~7):M;

  // block order, real data
  for (l=0; l<L; l++,x+=Sx,y+=Sy,z+=Sz)
  {
      for (m=0; m<M; m++)
      {
          float32_t A=0;

          for (n=0; n<4; n++) A+=x[m*4+n]*y[n];

          z[m]=A;
      }
  }
#endif
} /* matvmul4xmnf() */

/*-------------------------------------------------------------------------
Real Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of real matrices or vectors. Both the block order and streaming 
order are allowed for input/output matrix sequences.

Data format and order options:
  Suffix   Data Order             Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Parameters:
Temporary: 
pScr        Scratch memory area. To determine the scratch area size required by
            a function <fun>, use the respective helper function 
            <fun>_getScratchSize()
Input:
x[L*Sx]     Sequence of left-hand input matrices
y[L*Sy]     Sequence of right-hand input matrices
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, Floating-Point, MxN*Nx1->Mx1, Sx=MxN, Sy=(N>4)?((N+7)&~7):N, Sz=(M>4)?((M+7)&~7):M
   Restrictions:
     L must be a multiple of 2
     N, M must be multiples of 4
*/
void matvmulnxmnf ( void * pScr,
                    float32_t * restrict z,
              const float32_t * restrict x,
              const float32_t * restrict y,
              int N, int M, int L )
{
  /* check restrictions */
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT(L % 2 == 0);
  NASSERT(M % 4 == 0);
  NASSERT(N % 4 == 0);

  if ( M<=0 || L<=0 ) return;
  if ( N<=0 )
  {
      int k;
      xb_vecN_2xf32 * pz;
      xb_vecN_2xf32 zero;

      pz = (xb_vecN_2xf32 *)z;
      zero = BBE_ZERON_2XF32();
      __Pragma("loop_count min=1");
      for (k=0; k<((M*L)>>(LOG2_BBE_SIMD_WIDTH-1)); k++)
      {
          BBE_SVN_2XF32_IP(zero, pz, 2*BBE_SIMD_WIDTH);
      }

      return;
  }

  if (4==N)
  {
      matvmul4xmnf(z, x, y, M, L);
  }
  else if (0!=(N&7))
  {
      matvmulnxmnf_N4(z, x, y, N, M, L);
  }
  else
  {
      matvmulnxmnf_N8(z, x, y, N, M, L);
  }

} /* matvmulnxmnf() */
#endif

/* Return the scratch area size, in bytes. */
size_t matvmulnxmnf_getScratchSize ( int N, int M, int L )
{
    return 0;
} /* matvmulnxmnf_getScratchSize() */
