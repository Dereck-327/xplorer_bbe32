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
DISCARD_FUN(void, matmulnxmnf,( void * pScr,
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

/* M is multiple of 8 */
static void matmulnxmnf_M8 ( float32_t * restrict z,
                       const float32_t * restrict x,
                       const float32_t * restrict y,
                       int N, int M, int L )
{
  int mrow, mcol, n, l;

  const xtcomplexfloat * restrict px0;
  const xtcomplexfloat * restrict px1;
  const xtcomplexfloat * restrict px2;
  const xtcomplexfloat * restrict px3;
  const xb_vecN_2xf32  * restrict py;
        xb_vecN_2xf32  * restrict pz0;
        xb_vecN_2xf32  * restrict pz1;
        xb_vecN_2xf32  * restrict pz2;
        xb_vecN_2xf32  * restrict pz3;
  xb_vecN_2xf32 x_sel001, x_sel023;
  xb_vecN_2xf32 x_sel101, x_sel123;
  xb_vecN_2xf32 x_sel201, x_sel223;
  xb_vecN_2xf32 x_sel301, x_sel323;
  xb_vecN_2xf32 y0, y1, y2, y3;
  xb_vecN_2xf32 z0, z1, z2, z3;
  xb_vecN_4xcf32 tmp0, tmp1, tmp2, tmp3;
  
  pz0 = (xb_vecN_2xf32 *)(z+0*M);
  pz1 = (xb_vecN_2xf32 *)(z+1*M);
  pz2 = (xb_vecN_2xf32 *)(z+2*M);
  pz3 = (xb_vecN_2xf32 *)(z+3*M);
  px0 = (const xtcomplexfloat *)(x+0*N);
  px1 = (const xtcomplexfloat *)(x+1*N);
  px2 = (const xtcomplexfloat *)(x+2*N);
  px3 = (const xtcomplexfloat *)(x+3*N);
  
  __Pragma("loop_count min=1");
  for (l=0; l<L; l++)
  {
      /* Compute matrices by 4x(BBE_SIMD_WIDTH/2) pieces */
      __Pragma("loop_count min=1");
      for (mrow=0; mrow<(M>>2); mrow++)
      {
          py  = (const xb_vecN_2xf32  *)(y+l*N*M);
          __Pragma("loop_count min=1");
          for (mcol=0; mcol<(M>>(LOG2_BBE_SIMD_WIDTH-1)); mcol++)
          {
              z0 = z1 = z2 = z3 = BBE_ZERON_2XF32();

              __Pragma("loop_count min=1");
              for (n=0; n<(N>>2); n++)
              {
                  /* Load 4x(BBE_SIMD_WIDTH/2) elements of matrix Y */
                  BBE_LVN_2XF32_XP(y0, py, sizeof(float32_t)*M);
                  BBE_LVN_2XF32_XP(y1, py, sizeof(float32_t)*M);
                  BBE_LVN_2XF32_XP(y2, py, sizeof(float32_t)*M);
                  BBE_LVN_2XF32_XP(y3, py, sizeof(float32_t)*M);

                  /* Load 4x4 elements of matrix X and interleave it */
                  BBE_LSN_4XCF32_IP(tmp0, px0, 2*sizeof(float32_t));
                  x_sel001 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp0, 0));
                  BBE_LSN_4XCF32_IP(tmp0, px0, 2*sizeof(float32_t));
                  x_sel023 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp0, 0));
                  BBE_LSN_4XCF32_IP(tmp1, px1, 2*sizeof(float32_t));
                  x_sel101 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp1, 0));
                  BBE_LSN_4XCF32_IP(tmp1, px1, 2*sizeof(float32_t));
                  x_sel123 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp1, 0));
                  BBE_LSN_4XCF32_IP(tmp2, px2, 2*sizeof(float32_t));
                  x_sel201 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp2, 0));
                  BBE_LSN_4XCF32_IP(tmp2, px2, 2*sizeof(float32_t));
                  x_sel223 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp2, 0));
                  BBE_LSN_4XCF32_IP(tmp3, px3, 2*sizeof(float32_t));
                  x_sel301 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp3, 0));
                  BBE_LSN_4XCF32_IP(tmp3, px3, 2*sizeof(float32_t));
                  x_sel323 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp3, 0));

                  /* Z[4x(BBE_SIMD_WIDTH/2)] = X[4x4]*Y[4x(BBE_SIMD_WIDTH/2)] */
                  BBE_MULMASN_2XF32(z0, x_sel001, y0, 0, 0x4);
                  BBE_MULMASN_2XF32(z0, x_sel001, y1, 0, 0xE);
                  BBE_MULMASN_2XF32(z0, x_sel023, y2, 0, 0x4);
                  BBE_MULMASN_2XF32(z0, x_sel023, y3, 0, 0xE);
                  BBE_MULMASN_2XF32(z1, x_sel101, y0, 0, 0x4);
                  BBE_MULMASN_2XF32(z1, x_sel101, y1, 0, 0xE);
                  BBE_MULMASN_2XF32(z1, x_sel123, y2, 0, 0x4);
                  BBE_MULMASN_2XF32(z1, x_sel123, y3, 0, 0xE);
                  BBE_MULMASN_2XF32(z2, x_sel201, y0, 0, 0x4);
                  BBE_MULMASN_2XF32(z2, x_sel201, y1, 0, 0xE);
                  BBE_MULMASN_2XF32(z2, x_sel223, y2, 0, 0x4);
                  BBE_MULMASN_2XF32(z2, x_sel223, y3, 0, 0xE);
                  BBE_MULMASN_2XF32(z3, x_sel301, y0, 0, 0x4);
                  BBE_MULMASN_2XF32(z3, x_sel301, y1, 0, 0xE);
                  BBE_MULMASN_2XF32(z3, x_sel323, y2, 0, 0x4);
                  BBE_MULMASN_2XF32(z3, x_sel323, y3, 0, 0xE);
              }
              /* Save outputs */
              BBE_SVN_2XF32_IP(z0, pz0, 2 * BBE_SIMD_WIDTH);
              BBE_SVN_2XF32_IP(z1, pz1, 2 * BBE_SIMD_WIDTH);
              BBE_SVN_2XF32_IP(z2, pz2, 2 * BBE_SIMD_WIDTH);
              BBE_SVN_2XF32_IP(z3, pz3, 2 * BBE_SIMD_WIDTH);

              /* Prepare pointers for next (BBE_SIMD_WIDTH/2) columns */
              px0 = (const xtcomplexfloat *)((float32_t*)px0-N);
              px1 = (const xtcomplexfloat *)((float32_t*)px1-N);
              px2 = (const xtcomplexfloat *)((float32_t*)px2-N);
              px3 = (const xtcomplexfloat *)((float32_t*)px3-N);
              py  = (const xb_vecN_2xf32  *)((float32_t*)py-N*M+(BBE_SIMD_WIDTH/2));
          }
          /* Jump to next 4 rows */
          pz0 = pz3;
          pz1 = (xb_vecN_2xf32 *)((float32_t*)pz0+M);
          pz2 = (xb_vecN_2xf32 *)((float32_t*)pz1+M);
          pz3 = (xb_vecN_2xf32 *)((float32_t*)pz2+M);
          
          px0 = (const xtcomplexfloat *)((float32_t*)px0+4*N);
          px1 = (const xtcomplexfloat *)((float32_t*)px1+4*N);
          px2 = (const xtcomplexfloat *)((float32_t*)px2+4*N);
          px3 = (const xtcomplexfloat *)((float32_t*)px3+4*N);

          py  = (const xb_vecN_2xf32  *)((float32_t*)py-M);
      }
  }
}

/* M is multiple of 4 */
static void matmulnxmnf_M4 ( float32_t * restrict z,
                       const float32_t * restrict x,
                       const float32_t * restrict y,
                       int N, int M, int L )
{
  static const uint16_t ALIGN(32) tblSel[BBE_SIMD_WIDTH]=
  {
    0x00, 0x01, 0x02, 0x03, 0x00, 0x01, 0x02, 0x03, 0x10, 0x11, 0x12, 0x13, 0x10, 0x11, 0x12, 0x13
  };
  int mrow, mcol, n, l;

  const xtcomplexfloat * restrict px0;
  const xtcomplexfloat * restrict px1;
  const xtcomplexfloat * restrict px2;
  const xtcomplexfloat * restrict px3;
  const xb_vecN_2xf32  * restrict py;
  const xb_vecN_2xf32  * restrict py_tmp;
        xb_vecN_2xf32  * restrict pz0;
        xb_vecN_2xf32  * restrict pz1;
        xb_vecN_2xf32  * restrict pz2;
        xb_vecN_2xf32  * restrict pz3;
  xb_vecN_2xf32 x_sel001, x_sel023;
  xb_vecN_2xf32 x_sel101, x_sel123;
  xb_vecN_2xf32 x_sel201, x_sel223;
  xb_vecN_2xf32 x_sel301, x_sel323;
  xb_vecN_2xf32 y0, y1, y2, y3;
  xb_vecN_2xf32 z0, z1, z2, z3;
  xb_vecN_4xcf32 tmp0, tmp1, tmp2, tmp3;
  valign aly1, aly3, alz1, alz3;
  xb_vecNx16 vTmp;
  vselN_2 sel0;
  vboolN_2 mask_4n;
  
  vTmp = BBE_LVNX16_I((const xb_vecNx16*)tblSel,0);
  sel0 = BBE_MOVVSELN_2NX16(vTmp,0);
  alz1 = BBE_ZALIGN();
  alz3 = BBE_ZALIGN();
  mask_4n = BBE_LTRN_2I(4);
  pz0 = (xb_vecN_2xf32 *)(z+0*M);
  pz1 = (xb_vecN_2xf32 *)(z+1*M);
  pz2 = (xb_vecN_2xf32 *)(z+2*M);
  pz3 = (xb_vecN_2xf32 *)(z+3*M);
  px0 = (const xtcomplexfloat *)(x+0*N);
  px1 = (const xtcomplexfloat *)(x+1*N);
  px2 = (const xtcomplexfloat *)(x+2*N);
  px3 = (const xtcomplexfloat *)(x+3*N);
  
  __Pragma("loop_count min=1");
  for (l=0; l<L; l++)
  {
      /* Compute matrices by 4x(BBE_SIMD_WIDTH/2) pieces */
      __Pragma("loop_count min=1");
      for (mrow=0; mrow<(M>>2); mrow++)
      {
          py  = (const xb_vecN_2xf32  *)(y+l*N*M);

          for (mcol=0; mcol<(M>>(LOG2_BBE_SIMD_WIDTH-1)); mcol++)
          {
              z0 = z1 = z2 = z3 = BBE_ZERON_2XF32();

              __Pragma("loop_count min=1");
              for (n=0; n<(N>>2); n++)
              {
                  /* Load 4x(BBE_SIMD_WIDTH/2) elements of matrix Y */
                  /* each odd load is unaligned */
                  BBE_LVN_2XF32_XP(y0, py, sizeof(float32_t)*M);
                  py_tmp = py; py = (const xb_vecN_2xf32 *)((intptr_t)py+sizeof(float32_t)*M);
                  aly1 = BBE_LAN_2XF32_PP(py_tmp);
                  BBE_LAN_2XF32_IP(y1, aly1, py_tmp);
                  BBE_LVN_2XF32_XP(y2, py, sizeof(float32_t)*M);
                  py_tmp = py; py = (const xb_vecN_2xf32 *)((intptr_t)py+sizeof(float32_t)*M);
                  aly3 = BBE_LAN_2XF32_PP(py_tmp);
                  BBE_LAN_2XF32_IP(y3, aly3, py_tmp);

                  /* Load 4x4 elements of matrix X and interleave it */
                  BBE_LSN_4XCF32_IP(tmp0, px0, 2*sizeof(float32_t));
                  x_sel001 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp0, 0));
                  BBE_LSN_4XCF32_IP(tmp0, px0, 2*sizeof(float32_t));
                  x_sel023 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp0, 0));
                  BBE_LSN_4XCF32_IP(tmp1, px1, 2*sizeof(float32_t));
                  x_sel101 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp1, 0));
                  BBE_LSN_4XCF32_IP(tmp1, px1, 2*sizeof(float32_t));
                  x_sel123 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp1, 0));
                  BBE_LSN_4XCF32_IP(tmp2, px2, 2*sizeof(float32_t));
                  x_sel201 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp2, 0));
                  BBE_LSN_4XCF32_IP(tmp2, px2, 2*sizeof(float32_t));
                  x_sel223 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp2, 0));
                  BBE_LSN_4XCF32_IP(tmp3, px3, 2*sizeof(float32_t));
                  x_sel301 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp3, 0));
                  BBE_LSN_4XCF32_IP(tmp3, px3, 2*sizeof(float32_t));
                  x_sel323 = BBE_MOVN_2XF32_FROMN_4XCF32(BBE_REPN_4XCF32(tmp3, 0));

                  /* Z[4x(BBE_SIMD_WIDTH/2)] = X[4x4]*Y[4x(BBE_SIMD_WIDTH/2)] */
                  BBE_MULMASN_2XF32(z0, x_sel001, y0, 0, 0x4);
                  BBE_MULMASN_2XF32(z0, x_sel001, y1, 0, 0xE);
                  BBE_MULMASN_2XF32(z0, x_sel023, y2, 0, 0x4);
                  BBE_MULMASN_2XF32(z0, x_sel023, y3, 0, 0xE);
                  BBE_MULMASN_2XF32(z1, x_sel101, y0, 0, 0x4);
                  BBE_MULMASN_2XF32(z1, x_sel101, y1, 0, 0xE);
                  BBE_MULMASN_2XF32(z1, x_sel123, y2, 0, 0x4);
                  BBE_MULMASN_2XF32(z1, x_sel123, y3, 0, 0xE);
                  BBE_MULMASN_2XF32(z2, x_sel201, y0, 0, 0x4);
                  BBE_MULMASN_2XF32(z2, x_sel201, y1, 0, 0xE);
                  BBE_MULMASN_2XF32(z2, x_sel223, y2, 0, 0x4);
                  BBE_MULMASN_2XF32(z2, x_sel223, y3, 0, 0xE);
                  BBE_MULMASN_2XF32(z3, x_sel301, y0, 0, 0x4);
                  BBE_MULMASN_2XF32(z3, x_sel301, y1, 0, 0xE);
                  BBE_MULMASN_2XF32(z3, x_sel323, y2, 0, 0x4);
                  BBE_MULMASN_2XF32(z3, x_sel323, y3, 0, 0xE);
              }
              /* Save outputs */
              BBE_SVN_2XF32_IP(z0, pz0, 2 * BBE_SIMD_WIDTH);
              BBE_SAN_2XF32_IP(z1, alz1, pz1);
              BBE_SVN_2XF32_IP(z2, pz2, 2 * BBE_SIMD_WIDTH);
              BBE_SAN_2XF32_IP(z3, alz3, pz3);

              /* Prepare pointers for next (BBE_SIMD_WIDTH/2) columns */
              px0 = (const xtcomplexfloat *)((float32_t*)px0-N);
              px1 = (const xtcomplexfloat *)((float32_t*)px1-N);
              px2 = (const xtcomplexfloat *)((float32_t*)px2-N);
              px3 = (const xtcomplexfloat *)((float32_t*)px3-N);
              py  = (const xb_vecN_2xf32  *)((float32_t*)py-N*M+(BBE_SIMD_WIDTH/2));
          }
          BBE_SAN_2XF32POS_FP(alz1, pz1);
          BBE_SAN_2XF32POS_FP(alz3, pz3);
          
          /* Compute last 4 columns of matrix */
          {
              z0 = z1 = z2 = z3 = BBE_ZERON_2XF32();

              __Pragma("loop_count min=1");
              for (n=0; n<(N>>2); n++)
              {
                  /* Load 4x4 elements of matrix Y */
                  /* each odd load is unaligned */
                  BBE_LVN_2XF32T_XP(y0, py, sizeof(float32_t)*(M-4), mask_4n);
                  BBE_LVN_2XF32F_XP(y1, py, sizeof(float32_t)*(M+4), mask_4n);
                  BBE_LVN_2XF32T_XP(y2, py, sizeof(float32_t)*(M-4), mask_4n);
                  BBE_LVN_2XF32F_XP(y3, py, sizeof(float32_t)*(M+4), mask_4n);
                  y0 = BBE_SHFLN_2XF32I(y0, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
                  y1 = BBE_SHFLN_2XF32I(y1, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);
                  y2 = BBE_SHFLN_2XF32I(y2, BBE_SHFLI_MMC4X4X4X4_M2_STEP_1);
                  y3 = BBE_SHFLN_2XF32I(y3, BBE_SHFLI_MMC4X4X4X4_M2_STEP_2);

                  /* Load 4x4 elements of matrix X and interleave it */
                  BBE_LSN_4XCF32_IP(tmp0, px0, 2*sizeof(float32_t));
                  BBE_LSN_4XCF32_IP(tmp1, px1, 2*sizeof(float32_t));
                  BBE_LSN_4XCF32_IP(tmp2, px2, 2*sizeof(float32_t));
                  BBE_LSN_4XCF32_IP(tmp3, px3, 2*sizeof(float32_t));
                  x_sel001 = BBE_MOVN_2XF32_FROMN_4XCF32(tmp0);
                  x_sel101 = BBE_MOVN_2XF32_FROMN_4XCF32(tmp1);
                  x_sel201 = BBE_MOVN_2XF32_FROMN_4XCF32(tmp2);
                  x_sel301 = BBE_MOVN_2XF32_FROMN_4XCF32(tmp3);
                  x_sel001 = BBE_SELN_2XF32(x_sel101, x_sel001, sel0);
                  x_sel201 = BBE_SELN_2XF32(x_sel301, x_sel201, sel0);

                  BBE_LSN_4XCF32_IP(tmp0, px0, 2*sizeof(float32_t));
                  BBE_LSN_4XCF32_IP(tmp1, px1, 2*sizeof(float32_t));
                  BBE_LSN_4XCF32_IP(tmp2, px2, 2*sizeof(float32_t));
                  BBE_LSN_4XCF32_IP(tmp3, px3, 2*sizeof(float32_t));
                  x_sel023 = BBE_MOVN_2XF32_FROMN_4XCF32(tmp0);
                  x_sel123 = BBE_MOVN_2XF32_FROMN_4XCF32(tmp1);
                  x_sel223 = BBE_MOVN_2XF32_FROMN_4XCF32(tmp2);
                  x_sel323 = BBE_MOVN_2XF32_FROMN_4XCF32(tmp3);
                  x_sel023 = BBE_SELN_2XF32(x_sel123, x_sel023, sel0);
                  x_sel223 = BBE_SELN_2XF32(x_sel323, x_sel223, sel0);

                  /* Z[4x(BBE_SIMD_WIDTH/2)] = X[4x4]*Y[4x(BBE_SIMD_WIDTH/2)] */
                  BBE_MULMASN_2XF32(z0, x_sel001, y0, 0, 0x4);
                  BBE_MULMASN_2XF32(z1, x_sel001, y1, 0, 0xE);
                  BBE_MULMASN_2XF32(z0, x_sel023, y2, 0, 0x4);
                  BBE_MULMASN_2XF32(z1, x_sel023, y3, 0, 0xE);
                  BBE_MULMASN_2XF32(z2, x_sel201, y0, 0, 0x4);
                  BBE_MULMASN_2XF32(z3, x_sel201, y1, 0, 0xE);
                  BBE_MULMASN_2XF32(z2, x_sel223, y2, 0, 0x4);
                  BBE_MULMASN_2XF32(z3, x_sel223, y3, 0, 0xE);
              }
              /* Save outputs */
              z0 = BBE_ADDN_2XF32(z0, z1);
              z2 = BBE_ADDN_2XF32(z2, z3);
              pz1 = (xb_vecN_2xf32*)((intptr_t)pz1-sizeof(float32_t)*4);
              pz3 = (xb_vecN_2xf32*)((intptr_t)pz3-sizeof(float32_t)*4);
              BBE_SVN_2XF32T_XP(z0, pz0, sizeof(float32_t)*4, mask_4n);
              BBE_SVN_2XF32F_IP(z0, pz1, sizeof(float32_t)*8, mask_4n);
              BBE_SVN_2XF32T_XP(z2, pz2, sizeof(float32_t)*4, mask_4n);
              BBE_SVN_2XF32F_IP(z2, pz3, sizeof(float32_t)*8, mask_4n);

              /* Prepare pointers for next (BBE_SIMD_WIDTH/2) columns */
              px0 = (const xtcomplexfloat *)((float32_t*)px0-N);
              px1 = (const xtcomplexfloat *)((float32_t*)px1-N);
              px2 = (const xtcomplexfloat *)((float32_t*)px2-N);
              px3 = (const xtcomplexfloat *)((float32_t*)px3-N);
              py  = (const xb_vecN_2xf32  *)((float32_t*)py-N*M+(BBE_SIMD_WIDTH/2));
          }
          /* Jump to next 4 rows */
          pz0 = pz3;
          pz1 = (xb_vecN_2xf32 *)((float32_t*)pz0+M);
          pz2 = (xb_vecN_2xf32 *)((float32_t*)pz1+M);
          pz3 = (xb_vecN_2xf32 *)((float32_t*)pz2+M);
          
          px0 = (const xtcomplexfloat *)((float32_t*)px0+4*N);
          px1 = (const xtcomplexfloat *)((float32_t*)px1+4*N);
          px2 = (const xtcomplexfloat *)((float32_t*)px2+4*N);
          px3 = (const xtcomplexfloat *)((float32_t*)px3+4*N);

          py  = (const xb_vecN_2xf32  *)((float32_t*)py-M);
      }
  }
}

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

/* Block Order, Floating-Point, MxN*NxM->MxM, Sx=MxN, Sy=NxM, Sz=MxM
   Restrictions:  
     N,M must be multiples of 4
*/
void matmulnxmnf ( void * pScr,
                   float32_t * restrict z,
             const float32_t * restrict x,
             const float32_t * restrict y,
             int N, int M, int L )
{
  /* check restrictions */
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
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
      __Pragma("loop_count min=2, factor=2");
      for (k=0; k<((M*M*L)>>(LOG2_BBE_SIMD_WIDTH-1)); k++)
      {
          BBE_SVN_2XF32_IP(zero, pz, 2*BBE_SIMD_WIDTH);
      }

      return;
  }

  if (0==(M&7))
  {
      matmulnxmnf_M8(z, x, y, N, M, L);
  }
  else
  {
      matmulnxmnf_M4(z, x, y, N, M, L);
  }

} /* matmulnxmnf() */
#endif

/* Return the scratch area size, in bytes. */
size_t matmulnxmnf_getScratchSize ( int N, int M, int L )
{
  return 0;
} /* matmulnxmnf_getScratchSize() */
