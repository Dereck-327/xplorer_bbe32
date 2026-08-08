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
    NatureDSP_Baseband library. Singular Value Decomposition
    Golub-Reinsch SVD algorithm for real upper bidiagonal matrices stored in
    stream order.
    Internal functions of real Golub-Kahan SVD step implementation for generic
    matrix sizes.
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#if 0
#include <math.h>
#endif
#include <float.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* SVD common declarations */
#include "svd_common.h"

#if HAVE_VFPU

#define MIN(a,b)     ((a)<(b)?(a):(b))
#define MAX(a,b)     ((a)>(b)?(a):(b))
#define SIDX(i,k)    ((i)*L+(k))

#define sz_f32       sizeof(float32_t)
#define sz_vbn2      sizeof(vboolN_2)

#if 0
/* Set p-th boolean to b
 * Input:
 *   p                         Boolean index, 0<=p<L
 *   b                         Boolean value, zero or non-zero
 * Input/Output:
 *   vb[L/(BBE_SIMD_WIDTH/2)]  Array of vector booleans */
static void vbn2Set( vboolN_2 * vb, int p, uint8_t b )
{
  vboolN vb0,vb1,vbm;
  int m,n;
  m = p/(BBE_SIMD_WIDTH/2); n = (p%(BBE_SIMD_WIDTH/2));
  vbm = BBE_EQNX16(BBE_SRAINX16(BBE_SEQNX16(),1), BBE_MOVVA16(n));
  vb0 = BBE_ANDNOTBN(BBE_MOVN_FROMN_2(vb[m]),vbm);
  vb1 = BBE_ANDBN(BBE_NEQNX16(BBE_MOVVA16(b),BBE_MOVVI16(0)),vbm);
  vb[m] = BBE_MOVN_2_FROMN(BBE_ORBN(vb0,vb1));
}

/* Get value of p-th boolean.
 * Input:
 *   p                         Boolean index, 0<=p<L
 *   vb[L/(BBE_SIMD_WIDTH/2)]  Array of vector booleans
 * Return value:
 *   Boolean value, 0 or 1. */
static uint8_t vbn2Get( const vboolN_2 * vb, int p )
{
  xb_vecNx16 v;
  int m,n;
  m = p/(BBE_SIMD_WIDTH/2); n = (p%(BBE_SIMD_WIDTH/2));
  v = BBE_MOVNX16T(BBE_MOVVI16(1), BBE_MOVVI16(0), BBE_MOVN_FROMN_2(vb[m]));
  return (uint8_t)BBE_EXTRANX16(v,2*n);
}
#endif

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * Summarize index data and fetch data for Wilkinson's shift computation.
 * Input:
 *   itsLim       Iterations count limit
 *   N            Number of columns in an input matrix
 *   L            Number of matrices
 *   D[N][L]      Main diagonal of input matrices
 *   F[N-1][L]    First superdiagonal of input matrices
 *   a_l[L]       Left index of working subblocks
 *   a_k[L]       Right index of working subblocks
 * Input/Output:
 *   a_its[L]     Iteration counters
 * Output:
 *   p_l_lo       Lowest left index over input matrices
 *   p_k_up       Uppermost right index over input matrices
 *   a_m[N][L/(BBE_SIMD_WIDTH/4)] (complex variant)
 *   a_m[N][L/(BBE_SIMD_WIDTH/2)] (real variant)
 *                Boolean labels of Working subblocks
 *   a_a[L]       Fetched data: D[a_k[0..L-1]-1] 
 *   a_b[L]       Fetched data: D[a_k[0..L-1]] 
 *   a_x[L]       Fetched data: D[a_l[0..L-1]] 
 *   a_y[L]       Fetched data: F[a_l[0..L-1]] 
 *   a_c[L]       Fetched data: F[a_k[0..L-1]-1] 
 *   a_s[L]       Fetched data: F[a_k[0..L-1]-2] 
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_its,rng_l,rng_k,a_m,a_a,a_b,a_x,a_y,a_c,a_s,a_l,a_k,D,F
 *      Must not overlap
 *   a_its,a_a,a_b,a_x,a_y,a_c,a_s,a_l,a_k,D,F
 *      Must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void rgrsvdsf_gks_fetch( 
                    int16_t   * restrict p_l_lo,
                    int16_t   * restrict p_k_up,
                    int16_t   * restrict a_its,
                    vboolN_2  * restrict a_m,
                    float32_t * restrict a_a,
                    float32_t * restrict a_b,
                    float32_t * restrict a_x,
                    float32_t * restrict a_y,
                    float32_t * restrict a_c,
                    float32_t * restrict a_s,
              const int16_t   *          a_l,
              const int16_t   *          a_k,
              const float32_t *          D,
              const float32_t *          F,
              int itsLim, int N, int L )
{
#if 1
  const xb_vecNx16    * restrict ITS_r;
        xb_vecNx16    * restrict ITS_w;
        vboolN_2      * restrict M_w;
        xb_vecN_2xf32 * restrict A_w;
        xb_vecN_2xf32 * restrict B_w;
        xb_vecN_2xf32 * restrict X_w;
        xb_vecN_2xf32 * restrict Y_w;
        xb_vecN_2xf32 * restrict C_w;
        xb_vecN_2xf32 * restrict S_w;
  const xb_vecNx16    * restrict L_r;
  const xb_vecNx16    * restrict K_r;
  const xb_vecN_2xf32 * restrict D0_r;
  const xb_vecN_2xf32 * restrict D1_r;
  const xb_vecN_2xf32 * restrict F0_r;
  const xb_vecN_2xf32 * restrict F1_r;

  int n,p;
  int l_lo,l_up,k_lo,k_up;
  int LW2 = L/(BBE_SIMD_WIDTH/2);

  NASSERT_ALIGN(a_its, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_y, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_l, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_k, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(D, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(F, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

#if 1
  {
    xb_vecNx16 vl_lo, vl_up, vk_lo, vk_up, vl, vk, its, vitsLim;
    vboolN blt;
    vl_lo = vk_lo = BBE_MOVVA16(N-1);
    vl_up = vk_up = BBE_ZERONX16();
    vitsLim = BBE_MOVVA16(itsLim);
    ITS_r = (xb_vecNx16*)a_its;
    L_r = (xb_vecNx16*)a_l;
    K_r = (xb_vecNx16*)a_k;
    for ( p=0; p<L/BBE_SIMD_WIDTH; p++ ) {
      BBE_LVNX16_IP(its, ITS_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
      blt = BBE_LTNX16(its, vitsLim);
      BBE_MINNX16T(vl_lo, vl_lo, vl, blt);
      BBE_MAXNX16T(vl_up, vl_up, vl, blt);
      BBE_MINNX16T(vk_lo, vk_lo, vk, blt);
      BBE_MAXNX16T(vk_up, vk_up, vk, blt);
    } /* p */
    if ((L&(BBE_SIMD_WIDTH/2))!=0) {
      BBE_LVNX16_IP(its, ITS_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
      its = BBE_SELNX16I(vitsLim, its, BBE_SELI_EXTRACT_LO_HALVES);
      blt = BBE_LTNX16(its, vitsLim);
      BBE_MINNX16T(vl_lo, vl_lo, vl, blt);
      BBE_MAXNX16T(vl_up, vl_up, vl, blt);
      BBE_MINNX16T(vk_lo, vk_lo, vk, blt);
      BBE_MAXNX16T(vk_up, vk_up, vk, blt);
    } /* L */
    vl_lo = BBE_MAXNX16(vl_lo, BBE_ZERONX16());
    vk_lo = BBE_MAXNX16(vk_lo, BBE_ZERONX16());
    l_lo = BBE_RMINA16VNX16(vl_lo); l_up = BBE_RMAXA16VNX16(vl_up);
    k_lo = BBE_RMINA16VNX16(vk_lo); k_up = BBE_RMAXA16VNX16(vk_up);
  }
#else
  l_lo = k_lo = N-1; l_up = k_up = 0;
  for ( p=0; p<L; p++ ) {
    if (a_its[p]<itsLim) {
      if (l_lo>a_l[p]) l_lo = a_l[p];
      if (l_up<a_l[p]) l_up = a_l[p];
      if (k_lo>a_k[p]) k_lo = a_k[p];
      if (k_up<a_k[p]) k_up = a_k[p];
    }
  } /* p */
  if (l_lo<0) l_lo = 0;
  if (k_lo<0) k_lo = 0;
#endif
  if (l_lo<k_up) {
#if 1
    xb_vecN_2xf32 a0,a1,b0,b1,c0,c1,s0,s1;
    xb_vecN_2xf32 x0,x1,y0,y1;
    xb_vecNx16 vk,vl,vits;
    vboolN bgt,blt,beq;
    vboolN_2 bm0,bm1,bs;
    D0_r = (xb_vecN_2xf32*)&D[SIDX(MAX(k_lo,1)-1,0)];
    D1_r = (xb_vecN_2xf32*)&D[SIDX(MAX(k_lo,1),0)];
    F0_r = (xb_vecN_2xf32*)&F[SIDX(MAX(k_lo,1)-1,0)];
    F1_r = (xb_vecN_2xf32*)&F[SIDX(MAX(k_lo,1)-2,0)];
    for ( n=MAX(k_lo,1); n<=k_up; n++ ) {
      bs = BBE_MOVN_2_FROMN(BBE_EXT0IB(BBE_MOVN_FROM1(BBE_MOVBA1(n>1)),BBE_SIMD_WIDTH));
      ITS_r = (xb_vecNx16*)a_its;
      L_r = (xb_vecNx16*)a_l;
      K_r = (xb_vecNx16*)a_k;
      A_w = (xb_vecN_2xf32*)a_a;
      B_w = (xb_vecN_2xf32*)a_b;
      C_w = (xb_vecN_2xf32*)a_c;
      S_w = (xb_vecN_2xf32*)a_s;
      for ( p=0; p<L/BBE_SIMD_WIDTH; p++ ) {
        BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vits, ITS_r, 2*BBE_SIMD_WIDTH);
        bgt = BBE_GTNX16(vk,vl);
        beq = BBE_EQNX16(vk, BBE_MOVVA16(n));
        blt = BBE_LTNX16(vits, BBE_MOVVA16(itsLim));
        BBE_EXTRACTBN(bm1, bm0, BBE_ANDBN(beq, BBE_ANDBN(bgt,blt)));
        /* D[SIDX(n-1,p)] */
        BBE_LVN_2XF32_IP(a0, D0_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(a1, D0_r, 2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32T_IP(a0, A_w, 2*BBE_SIMD_WIDTH, bm0);
        BBE_SVN_2XF32T_IP(a1, A_w, 2*BBE_SIMD_WIDTH, bm1);
        /* D[SIDX(n,p)] */
        BBE_LVN_2XF32_IP(b0, D1_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(b1, D1_r, 2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32T_IP(b0, B_w, 2*BBE_SIMD_WIDTH, bm0);
        BBE_SVN_2XF32T_IP(b1, B_w, 2*BBE_SIMD_WIDTH, bm1);
        /* F[SIDX(n-1,p)]; */
        BBE_LVN_2XF32_IP(c0, F0_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(c1, F0_r, 2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32T_IP(c0, C_w, 2*BBE_SIMD_WIDTH, bm0);
        BBE_SVN_2XF32T_IP(c1, C_w, 2*BBE_SIMD_WIDTH, bm1);
        /* ( n>1 ? F[SIDX(n-2,p)] : 0.f ) */
        BBE_LVN_2XF32T_IP(s0, F1_r, 2*BBE_SIMD_WIDTH, bs);
        BBE_LVN_2XF32T_IP(s1, F1_r, 2*BBE_SIMD_WIDTH, bs);
        BBE_SVN_2XF32T_IP(s0, S_w, 2*BBE_SIMD_WIDTH, bm0);
        BBE_SVN_2XF32T_IP(s1, S_w, 2*BBE_SIMD_WIDTH, bm1);
      } /* p */
      if ((L&(BBE_SIMD_WIDTH/2))!=0) {
        BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vits, ITS_r, 2*BBE_SIMD_WIDTH);
        bgt = BBE_GTNX16(vk,vl);
        beq = BBE_EQNX16(vk, BBE_MOVVA16(n));
        blt = BBE_LTNX16(vits, BBE_MOVVA16(itsLim));
        BBE_EXTRACTBN(bm1, bm0, BBE_ANDBN(beq, BBE_ANDBN(bgt,blt)));
        /* D[SIDX(n-1,p)] */
        BBE_LVN_2XF32_IP(a0, D0_r, 2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32T_IP(a0, A_w, 2*BBE_SIMD_WIDTH, bm0);
        /* D[SIDX(n,p)] */
        BBE_LVN_2XF32_IP(b0, D1_r, 2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32T_IP(b0, B_w, 2*BBE_SIMD_WIDTH, bm0);
        /* F[SIDX(n-1,p)]; */
        BBE_LVN_2XF32_IP(c0, F0_r, 2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32T_IP(c0, C_w, 2*BBE_SIMD_WIDTH, bm0);
        /* ( n>1 ? F[SIDX(n-2,p)] : 0.f ) */
        BBE_LVN_2XF32T_IP(s0, F1_r, 2*BBE_SIMD_WIDTH, bs);
        BBE_SVN_2XF32T_IP(s0, S_w, 2*BBE_SIMD_WIDTH, bm0);
      } /* L */
    } /* n */
    D0_r = (xb_vecN_2xf32*)&D[SIDX(l_lo,0)];
    F0_r = (xb_vecN_2xf32*)&F[SIDX(l_lo,0)];
    for ( n=l_lo; n<=MIN(l_up,N-2); n++ ) {
      ITS_r = (xb_vecNx16*)a_its;
      L_r = (xb_vecNx16*)a_l;
      K_r = (xb_vecNx16*)a_k;
      X_w = (xb_vecN_2xf32*)a_x;
      Y_w = (xb_vecN_2xf32*)a_y;
      for ( p=0; p<L/BBE_SIMD_WIDTH; p++ ) {
        BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vits, ITS_r, 2*BBE_SIMD_WIDTH);
        bgt = BBE_GTNX16(vk,vl);
        beq = BBE_EQNX16(vl, BBE_MOVVA16(n));
        blt = BBE_LTNX16(vits, BBE_MOVVA16(itsLim));
        BBE_EXTRACTBN(bm1, bm0, BBE_ANDBN(beq, BBE_ANDBN(bgt,blt)));
        /* D[SIDX(n,p)] */
        BBE_LVN_2XF32_IP(x0, D0_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(x1, D0_r, 2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32T_IP(x0, X_w, 2*BBE_SIMD_WIDTH, bm0);
        BBE_SVN_2XF32T_IP(x1, X_w, 2*BBE_SIMD_WIDTH, bm1);
        /* F[SIDX(n,p)] */
        BBE_LVN_2XF32_IP(y0, F0_r, 2*BBE_SIMD_WIDTH);
        BBE_LVN_2XF32_IP(y1, F0_r, 2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32T_IP(y0, Y_w, 2*BBE_SIMD_WIDTH, bm0);
        BBE_SVN_2XF32T_IP(y1, Y_w, 2*BBE_SIMD_WIDTH, bm1);
      } /* p */
      if ((L&(BBE_SIMD_WIDTH/2))!=0) {
        BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(vits, ITS_r, 2*BBE_SIMD_WIDTH);
        bgt = BBE_GTNX16(vk,vl);
        beq = BBE_EQNX16(vl, BBE_MOVVA16(n));
        blt = BBE_LTNX16(vits, BBE_MOVVA16(itsLim));
        BBE_EXTRACTBN(bm1, bm0, BBE_ANDBN(beq, BBE_ANDBN(bgt,blt)));
        /* D[SIDX(n,p)] */
        BBE_LVN_2XF32_IP(x0, D0_r, 2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32T_IP(x0, X_w, 2*BBE_SIMD_WIDTH, bm0);
        /* F[SIDX(n,p)] */
        BBE_LVN_2XF32_IP(y0, F0_r, 2*BBE_SIMD_WIDTH);
        BBE_SVN_2XF32T_IP(y0, Y_w, 2*BBE_SIMD_WIDTH, bm0);
      } /* L */
    } /* n */
#else
    for ( n=MAX(k_lo,1); n<=k_up; n++ ) {
      for ( p=0; p<L; p++ ) {
        if (a_k[p]>a_l[p] && a_its[p]<itsLim) {
          if (a_k[p]==n) {
            a_a[p] = D[SIDX(n-1,p)]; 
            a_b[p] = D[SIDX(n,p)]; 
            a_c[p] = F[SIDX(n-1,p)]; 
            a_s[p] = ( n>1 ? F[SIDX(n-2,p)] : 0.f );
          }
        } /* a_k[p], a_l[p], a_its[p] */
      } /* p */
    } /* n */
    for ( n=l_lo; n<=MIN(l_up,N-2); n++ ) {
      for ( p=0; p<L; p++ ) {
        if (a_k[p]>a_l[p] && a_its[p]<itsLim) {
          if (a_l[p]==n) {
            a_x[p] = D[SIDX(n,p)]; 
            a_y[p] = F[SIDX(n,p)]; 
          }
        } /* a_k[p], a_l[p], a_its[p] */
      } /* p */
    } /* n */
#endif
  } /* l_lo, k_up */
#if 1
  M_w = a_m + l_lo*LW2;
  for ( n=l_lo; n<k_up; n++ ) {
    xb_vecNx16 vl,vk,vn,vits,vitsLim;
    vboolN bl,bk,bi;
    vboolN_2 bm0,bm1;
    vn = BBE_MOVVA16(n);
    vitsLim = BBE_MOVVA16(itsLim);
    ITS_r = (xb_vecNx16*)a_its;
    L_r = (xb_vecNx16*)a_l;
    K_r = (xb_vecNx16*)a_k;
    for ( p=0; p<L/BBE_SIMD_WIDTH; p++ ) {
      BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vits, ITS_r, 2*BBE_SIMD_WIDTH);
      bl = BBE_LENX16(vl,vn);
      bk = BBE_LTNX16(vn,vk);
      bi = BBE_LTNX16(vits, vitsLim);
      BBE_EXTRACTBN(bm1, bm0, BBE_ANDBN(bl, BBE_ANDBN(bk,bi)));
      BBE_SBN_2_IP(bm0, M_w, sz_vbn2);
      BBE_SBN_2_IP(bm1, M_w, sz_vbn2);
    } /* p */
    if ((L&(BBE_SIMD_WIDTH/2))!=0) {
      BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vits, ITS_r, 2*BBE_SIMD_WIDTH);
      bl = BBE_LENX16(vl,vn);
      bk = BBE_LTNX16(vn,vk);
      bi = BBE_LTNX16(vits, vitsLim);
      BBE_EXTRACTBN(bm1, bm0, BBE_ANDBN(bl, BBE_ANDBN(bk,bi)));
      BBE_SBN_2_IP(bm0, M_w, sz_vbn2);
    } /* L */
  } /* n */
#else
  for ( n=l_lo; n<k_up; n++ ) {
    int LW2 = L/(BBE_SIMD_WIDTH/2);
    for ( p=0; p<L; p++ ) {
      vbn2Set(a_m+n*LW2, p, a_l[p]<=n && n<a_k[p] && a_its[p]<itsLim);
    } /* p */
  } /* n */
#endif
#if 1
  {
    xb_vecNx16 vl,vk,vits,vitsLim;
    vboolN bgt, blt;
    vitsLim = BBE_MOVVA16(itsLim);
    ITS_r = ITS_w = (xb_vecNx16*)a_its;
    L_r = (xb_vecNx16*)a_l;
    K_r = (xb_vecNx16*)a_k;
    for ( p=0; p<L/BBE_SIMD_WIDTH; p++ ) {
      BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vits, ITS_r, 2*BBE_SIMD_WIDTH);
      bgt = BBE_GTNX16(vk,vl);
      blt = BBE_LTNX16(vits,vitsLim);
      vits = BBE_ADDNX16(vits, BBE_MOVVINT16(1));
      BBE_SVNX16T_IP(vits, ITS_w, 2*BBE_SIMD_WIDTH, BBE_ANDBN(bgt,blt));
    } /* p */
    if ((L&(BBE_SIMD_WIDTH/2))!=0) {
      BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
      BBE_LVNX16_IP(vits, ITS_r, 2*BBE_SIMD_WIDTH);
      vits = BBE_SELNX16I(vitsLim, vits, BBE_SELI_EXTRACT_LO_HALVES);
      bgt = BBE_GTNX16(vk,vl);
      blt = BBE_LTNX16(vits,vitsLim);
      vits = BBE_ADDNX16(vits, BBE_MOVVINT16(1));
      BBE_SVNX16T_IP(vits, ITS_w, 2*BBE_SIMD_WIDTH, BBE_ANDBN(bgt,blt));
    } /* L */
  }
#else
  for ( p=0; p<L; p++ ) {
    if (a_k[p]>a_l[p] && a_its[p]<itsLim) {
      a_its[p]++;
    }
  } /* p */
#endif
  *p_l_lo = l_lo; *p_k_up = k_up;
#else
  int n,p;
  int l_lo,l_up,k_lo,k_up;
  int LW2 = L/(BBE_SIMD_WIDTH/2);

  NASSERT_ALIGN(a_its, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_y, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_l, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_k, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(D, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(F, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  l_lo = k_lo = N-1; l_up = k_up = 0;
  for ( p=0; p<L; p++ ) {
    if (a_its[p]<itsLim) {
      if (l_lo>a_l[p]) l_lo = a_l[p];
      if (l_up<a_l[p]) l_up = a_l[p];
      if (k_lo>a_k[p]) k_lo = a_k[p];
      if (k_up<a_k[p]) k_up = a_k[p];
    }
  } /* p */
  if (l_lo<0) l_lo = 0;
  if (k_lo<0) k_lo = 0;
  if (l_lo<k_up) {
    for ( n=MAX(k_lo,1); n<=k_up; n++ ) {
      for ( p=0; p<L; p++ ) {
        if (a_k[p]>a_l[p] && a_its[p]<itsLim) {
          if (a_k[p]==n) {
            a_a[p] = D[SIDX(n-1,p)]; 
            a_b[p] = D[SIDX(n,p)]; 
            a_c[p] = F[SIDX(n-1,p)]; 
            a_s[p] = ( n>1 ? F[SIDX(n-2,p)] : 0.f );
          }
        } /* a_k[p], a_l[p], a_its[p] */
      } /* p */
    } /* n */
    for ( n=l_lo; n<=MIN(l_up,N-2); n++ ) {
      for ( p=0; p<L; p++ ) {
        if (a_k[p]>a_l[p] && a_its[p]<itsLim) {
          if (a_l[p]==n) {
            a_x[p] = D[SIDX(n,p)]; 
            a_y[p] = F[SIDX(n,p)]; 
          }
        } /* a_k[p], a_l[p], a_its[p] */
      } /* p */
    } /* n */
  } /* l_lo, k_up */
  for ( n=l_lo; n<k_up; n++ ) {
    for ( p=0; p<L; p++ ) {
      vbn2Set(a_m+n*LW2, p, a_l[p]<=n && n<a_k[p] && a_its[p]<itsLim);
    } /* p */
  } /* n */
  for ( p=0; p<L; p++ ) {
    if (a_k[p]>a_l[p] && a_its[p]<itsLim) {
      a_its[p]++;
    }
  } /* p */
  *p_l_lo = l_lo; *p_k_up = k_up;
#endif
} /* rgrsvdsf_gks_fetch() */

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * Compute Wilkinson's shift.
 * Input:
 *   L            Number of matrices
 *   D[N][L]      Main diagonal of input matrices
 *   F[N-1][L]    First superdiagonal of input matrices
 *   a_c[L]       Fetched data: F[a_k[0..L-1]-1].  May be reused
 *                as a temporal storage of intermediate results.
 *   a_s[L]       Fetched data: F[a_k[0..L-1]-2]. May be reused
 *                as a temporal storage of intermediate results.
 *   a_x[L]       Fetched data: D[a_l[0..L-1]] 
 *   a_y[L]       Fetched data: F[a_l[0..L-1]] 
 * Input/Output:
 *   a_a[L]       In:  fetched data: D[a_k[0..L-1]-1] 
 *                Out: element (0,0) of B'*B-mu*I
 *   a_b[L]       In:  fetched data: D[a_k[0..L-1]] 
 *                Out: element (0,1) of B'*B-mu*I
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_a,a_b,a_x,a_y,a_c,a_s
 *      Must not overlap and must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void rgrsvdsf_gks_wilkShift(
                    float32_t * restrict a_a,
                    float32_t * restrict a_b,
                    float32_t * restrict a_c,
                    float32_t * restrict a_s,
              const float32_t *          a_x,
              const float32_t *          a_y,
              int L )
{
  /*
   * MATLAB outline:
   *   % Wilkinson's shift: for the bottom-right 2x2 subblock of B'*B, find an
   *   % eigenvalue which is closer to the last diagonal element. Derived from
   *   % Algol code provided in [2].
   *   x = d(k-1); y = d(k); z = d(l); v = f(k-1);
   *   if k>2, w = f(k-2); else w = 0; end;
   *   t = ((x-y)*(x+y)+(w-v)*(w+v))/(2*x*v);
   *   r = sqrt(t*t+1); if t<0, r = -r; end;
   *   % a <- t11-mu; b <- t12
   *   a = (z-y)*(z+y)+v*(x/(t+r)-v); b = z*f(l);
   */
#if 1
  const xb_vecN_2xf32 * restrict A_r;
        xb_vecN_2xf32 * restrict A_w;
  const xb_vecN_2xf32 * restrict B_r;
        xb_vecN_2xf32 * restrict B_w;
  const xb_vecN_2xf32 * restrict S_r;
        xb_vecN_2xf32 * restrict S_w;
  const xb_vecN_2xf32 * restrict C_r;
  const xb_vecN_2xf32 * restrict X_r;
  const xb_vecN_2xf32 * restrict Y_r;

  int p;

  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_y, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  /*
   * t <- ((x-y)*(x+y)+(w-v)*(w+v))/(2*x*v);
   */

#if 1
  {
    xb_vecN_2xf32 a,b,r,s,t,x,y,v,w;
    vboolN_2 blt;
    A_r = (xb_vecN_2xf32*)a_a; B_r = (xb_vecN_2xf32*)a_b;
    C_r = (xb_vecN_2xf32*)a_c; S_r = S_w = (xb_vecN_2xf32*)a_s;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
      BBE_LVN_2XF32_IP(x, A_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(y, B_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(v, C_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(w, S_r, 2*BBE_SIMD_WIDTH);
      BBE_ADDSUBN_2XF32(r,s,x,y);
      BBE_ADDSUBN_2XF32(a,b,w,v);
      /* r <- (x-y)*(x+y)+(w-v)*(w+v) */
      r = BBE_MULN_2XF32(r,s); BBE_MULAN_2XF32(r,a,b);
      /* s <- 2*x*v */
      s = BBE_MULN_2XF32(BBE_CONSTN_2XF32(2), BBE_MULN_2XF32(x,v));
      /* t <- r/s */
      blt = BBE_OLTN_2XF32(BBE_ABSN_2XF32(s), FLT_MIN);
      BBE_MULN_2XF32T(s,s,8388608.f,blt);
      a = BBE_RECIP0N_2XF32(s);
      b = BBE_CONSTN_2XF32(1); BBE_MULSN_2XF32(b,a,s); BBE_MULAN_2XF32(a,b,a);
      t = BBE_MULN_2XF32(r,a); BBE_MULSN_2XF32(r,t,s); BBE_MULAN_2XF32(t,r,a);
      BBE_MULN_2XF32T(t,t,8388608.f,blt);
      BBE_SVN_2XF32_IP(t, S_w, 2*BBE_SIMD_WIDTH);
    } /* p */
    __Pragma("no_reorder");
  }
#else
  {
    float32_t t,v,w,x,y;
    for ( p=0; p<L; p++ ) {
      x = a_a[p]; y = a_b[p]; v = a_c[p]; w = a_s[p];
      t = ((x-y)*(x+y) + (w-v)*(w+v))/(2*x*v);
      a_s[p] = t;
    } /* p */
  }
#endif

  /*
   * r <- sqrt(t*t+1);
   * s <- t<0 ? t-r : t+r;
   */

#if 1
  {
    xb_vecN_2xf32 d,t,r,s;
    S_r = S_w = (xb_vecN_2xf32*)a_s;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
      BBE_LVN_2XF32_IP(t, S_r, 2*BBE_SIMD_WIDTH);
      s = BBE_CONSTN_2XF32(1); BBE_MULAN_2XF32(s,t,t);
      r = BBE_FSQRTN_2XF32(s);
      BBE_ADDSUBN_2XF32(d,s,t,r);
      t = BBE_MOVN_2XF32T(d, s, BBE_OLTN_2XF32(t, BBE_CONSTN_2XF32(0)));
      BBE_SVN_2XF32_IP(t, S_w, 2*BBE_SIMD_WIDTH);
    } /* p */
  }
#else
  {
    float32_t r,t;
    for ( p=0; p<L; p++ ) {
      t = a_s[p];
      r = sqrtf(t*t+1.f); if (t<0) r = -r;
      a_s[p] = t+r;
    } /* p */
  }
#endif

  /*
   * a <- t11-mu == (z-y)*(z+y)+v*(x/s-v); b = z*f(l);
   * b <- t12 == z*f(l);
   */
#if 1
  {
    xb_vecN_2xf32 a,b,r,t,s,u,v,x,y,z;
    A_r = A_w = (xb_vecN_2xf32*)a_a;
    B_r = B_w = (xb_vecN_2xf32*)a_b;
    X_r = (xb_vecN_2xf32*)a_x;
    Y_r = (xb_vecN_2xf32*)a_y;
    C_r = (xb_vecN_2xf32*)a_c;
    S_r = (xb_vecN_2xf32*)a_s;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
      BBE_LVN_2XF32_IP(x, A_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(y, B_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(z, X_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(u, Y_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(v, C_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(s, S_r, 2*BBE_SIMD_WIDTH);
      /* r <- x/s; */
      a = BBE_RECIP0N_2XF32(s);
      b = BBE_CONSTN_2XF32(1); BBE_MULSN_2XF32(b,a,s); BBE_MULAN_2XF32(a,b,a);
      r = BBE_MULN_2XF32(x,a); BBE_MULSN_2XF32(x,r,s); BBE_MULAN_2XF32(r,x,a);
      /* a <- v*(x/s-v) */
      a = BBE_MULN_2XF32(v, BBE_SUBN_2XF32(r,v));
      /* r <- z-y; t <- z+y */
      BBE_ADDSUBN_2XF32(r,t,z,y);
      /* a <- (z-y)*(z+y) + v*(x/s-v); */
      BBE_MULAN_2XF32(a,r,t);
      /* b <- z*u; */
      b = BBE_MULN_2XF32(z,u);
      BBE_SVN_2XF32_IP(a, A_w, 2*BBE_SIMD_WIDTH);
      BBE_SVN_2XF32_IP(b, B_w, 2*BBE_SIMD_WIDTH);
    } /* p */
  }
#else
  {
    float32_t s,u,v,x,y,z;
    for ( p=0; p<L; p++ ) {
      x = a_a[p]; y = a_b[p]; z = a_x[p]; u = a_y[p]; v = a_c[p]; s = a_s[p];
      a_a[p] = (z-y)*(z+y) + v*(x/s-v); a_b[p] = z*u;
    } /* p */
  }
#endif
#else
  float32_t r,s,t,u,v,w,x,y,z;
  int p;

  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_x, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_y, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  for ( p=0; p<L; p++ ) {
    x = a_a[p]; y = a_b[p]; z = a_x[p]; v = a_c[p]; w = a_s[p];
    /* t <- ((x-y)*(x+y)+(w-v)*(w+v))/(2*x*v); */
    t = ((x-y)*(x+y) + (w-v)*(w+v))/(2*x*v);
    a_s[p] = t;
  } /* p */
  for ( p=0; p<L; p++ ) {
    t = a_s[p];
    /* r <- sqrt(t*t+1);
     * s <- t<0 ? t-r : t+r; */
    r = sqrtf(t*t+1.f); if (t<0) r = -r;
    a_s[p] = t+r;
  } /* p */
  for ( p=0; p<L; p++ ) {
    x = a_a[p]; y = a_b[p]; z = a_x[p]; u = a_y[p]; v = a_c[p]; s = a_s[p];
    /* a <- t11-mu == (z-y)*(z+y)+v*(x/s-v);
     * b <- t12 == z*f(l); */
    a_a[p] = (z-y)*(z+y) + v*(x/s-v); a_b[p] = z*u;
  } /* p */
#endif
} /* rgrsvdsf_gks_wilkShift() */

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * Compute Givens's rotating matrix.
 * For real data:
 *   G(a,b) <- [c,s;-s,c]: [a,b]*[c,s;-s,c] == [*,0], c^2+s^2 == 1
 * For complex data:
 *   G(a,b) <- [c,conj(s);-s,conj(c)]: [a,b]*G(a,b) == [*,0], c*conj(c)+s*conj(s) == 1
 * Temporary:
 *   a_x[L]       Scratch array of L entries
 * Input:
 *   L            Number of matrices
 *   a_a[L]       a-values
 *   a_b[L]       b-values
 * Output:
 *   a_c[L]       G(a,b) cosine values
 *   a_s[L]       G(a,b) sine values
 *  Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_a,a_b,a_c,a_s
 *      Must not overlap and must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void rgrsvdsf_gks_givens( 
                    float32_t * restrict a_c,
                    float32_t * restrict a_s,
              const float32_t *          a_a, 
              const float32_t *          a_b, 
              int L )
{
  /*
   * MATLAB reference code:
   *
   *   function [c,s] = givens(a,b)
   *   if b~=0
   *     if abs(b)>abs(a)
   *       t = a/b; s = 1/sqrt(1+t^2); c = -t*s;
   *     else
   *       t = b/a; c = 1/sqrt(1+t^2); s = -t*c;
   *     end
   *   else
   *     c = 1; s = 0;
   *   end;
   */
#if 1
  const xb_vecN_2xf32 * restrict C_r;
        xb_vecN_2xf32 * restrict C0_w;
        xb_vecN_2xf32 * restrict C1_w;
        xb_vecN_2xf32 * restrict S0_w;
        xb_vecN_2xf32 * restrict S1_w;
  const xb_vecN_2xf32 * restrict A_r;
  const xb_vecN_2xf32 * restrict B_r;

  xb_vecN_2xf32 a,b,f,g,h,t;
  vboolN_2 beq,blt;
  int p;

  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  C0_w = (xb_vecN_2xf32*)a_c;
  A_r = (xb_vecN_2xf32*)a_a;
  B_r = (xb_vecN_2xf32*)a_b;
  for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
    BBE_LVN_2XF32_IP(a, A_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(b, B_r, 2*BBE_SIMD_WIDTH);
    blt = BBE_ULTN_2XF32(BBE_ABSN_2XF32(a), BBE_ABSN_2XF32(b));
    f = BBE_MOVN_2XF32T(b,a,blt);
    g = BBE_MOVN_2XF32T(a,b,blt);
    /* if (b==0) f = 1.f; */
    beq = BBE_OEQN_2XF32(b, BBE_CONSTN_2XF32(0));
    BBE_CONSTN_2XF32T(f,1,beq);
    /* h <- ~1/f */
    h = BBE_RECIP0N_2XF32(f);
    /* Newton-Raphson refinement iteration for a reciprocal */
    t = BBE_CONSTN_2XF32(1); BBE_MULSN_2XF32(t,h,f); BBE_MULAN_2XF32(h,t,h);
    /* t <- ~g/f */
    t = BBE_MULN_2XF32(h,g);
    /* Modified Newton-Raphson iteration for a quotient */
    BBE_MULSN_2XF32(g,t,f); BBE_MULAN_2XF32(t,g,h);
    BBE_SVN_2XF32_IP(t, C0_w, 2*BBE_SIMD_WIDTH);
  } /* p */

  C_r = C0_w = C1_w = (xb_vecN_2xf32*)a_c;
  S0_w = S1_w = (xb_vecN_2xf32*)a_s;
  A_r = (xb_vecN_2xf32*)a_a;
  B_r = (xb_vecN_2xf32*)a_b;
  for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
    BBE_LVN_2XF32_IP(t, C_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(a, A_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(b, B_r, 2*BBE_SIMD_WIDTH);
    blt = BBE_ULTN_2XF32(BBE_ABSN_2XF32(a), BBE_ABSN_2XF32(b));
    /* f <- 1/sqrtf(1+t*t) */
    h = BBE_CONSTN_2XF32(1); BBE_MULAN_2XF32(h,t,t); f = BBE_RSQRTN_2XF32(h);
    /* g <- -f*t */
    g = BBE_MULN_2XF32(BBE_NEGN_2XF32(t), f);
    /* if (fabsf(a)<fabsf(b) || isnan(a) || isnan(b)) {
     *   *s = f; *c = g;
     * } else {
     *   *s = g; *c = f;
     * } */
    BBE_SVN_2XF32F_IP(f, C0_w, 2*BBE_SIMD_WIDTH, blt);
    BBE_SVN_2XF32F_IP(g, S0_w, 2*BBE_SIMD_WIDTH, blt);
    BBE_SVN_2XF32T_IP(g, C1_w, 2*BBE_SIMD_WIDTH, blt);
    BBE_SVN_2XF32T_IP(f, S1_w, 2*BBE_SIMD_WIDTH, blt);
  } /* p */
#else
  float32_t a,b,f,g,t;
  int p;

  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  for ( p=0; p<L; p++ ) {
    a = a_a[p]; b = a_b[p];
    if (fabsf(a)<fabsf(b) || isnan(a) || isnan(b)) {
      g = a; f = b;
    } else {
      g = b; f = a;
    }
    if (b==0) {
      f = 1.f;
    }
    a_c[p] = g/f;
  } /* p */

  for ( p=0; p<L; p++ ) {
    a = a_a[p]; b = a_b[p]; t = a_c[p];
    f = 1/sqrtf(1+t*t);
    g = -f*t;
    if (fabsf(a)<fabsf(b) || isnan(a) || isnan(b)) {
      a_s[p] = f; a_c[p] = g;
    } else {
      a_s[p] = g; a_c[p] = f;
    }
  } /* p */
#endif
} /* rgrsvdsf_gks_givens() */

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * 1st updating step of SVD step iteration: B <- B*G(a,b).
 * Input:
 *   n            Current position
 *   L            Number of matrices
 *   a_c[L]       G(a,b) cosine values
 *   a_s[L]       G(a,b) sine values
 *   a_l[L]       Left index of working subblocks
 *   a_m[N][L/(BBE_SIMD_WIDTH/4)] (complex variant)
 *   a_m[N][L/(BBE_SIMD_WIDTH/2)] (real variant)
 *                Boolean labels of Working subblocks
 * Input/Output:
 *   a_a[L]       In: B(n-1,n) Out: B(n,n)  
 *   a_b[L]       In: B(n-1,n+1) Out: B(n+1,n)
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_a,a_b,a_c,a_s,a_l,a_m
 *      Must not overlap
 *   a_a,a_b,a_c,a_s,a_l
 *      Must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void rgrsvdsf_gks_step1(
                    float32_t * restrict D,
                    float32_t * restrict F,
                    float32_t * restrict a_a,
                    float32_t * restrict a_b,
              const float32_t *          a_c,
              const float32_t *          a_s,
              const int16_t   *          a_l,
              const vboolN_2  *          a_m,
              int n, int L )
{
#if 1
        xb_vecN_2xf32 * restrict Fm_w;
  const xb_vecN_2xf32 * restrict Fn_r;
        xb_vecN_2xf32 * restrict Fn_w;
  const xb_vecN_2xf32 * restrict Dn_r;
  const xb_vecN_2xf32 * restrict Dq_r;
        xb_vecN_2xf32 * restrict Dq_w;
        xb_vecN_2xf32 * restrict A_rw;
        xb_vecN_2xf32 * restrict B_rw;
  const xb_vecN_2xf32 * restrict C_r;
  const xb_vecN_2xf32 * restrict S_r;
  const xb_vecNx16    * restrict L_r;
  const vboolN_2      * restrict M_r;

  xb_vecN_2xf32 a0,a1,b0,b1,c0,c1,f0,f1,g0,g1;
  xb_vecN_2xf32 h0,h1,s0,s1,x0,x1,y0,y1,z0,z1;
  xb_vecNx16 vl,vn;
  vboolN_2 bl0,bl1,bm0,bm1;
  int p, LW2 = L/(BBE_SIMD_WIDTH/2);

  NASSERT_ALIGN(D, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(F, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_l, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  vn = BBE_MOVVA16(n);
  Fm_w = (xb_vecN_2xf32*)&F[SIDX(n-1,0)];
  Fn_r = Fn_w = (xb_vecN_2xf32*)&F[SIDX(n,0)];
  Dn_r = (xb_vecN_2xf32*)&D[SIDX(n,0)];
  Dq_r = Dq_w = (xb_vecN_2xf32*)&D[SIDX(n+1,0)];
  A_rw = (xb_vecN_2xf32*)a_a;
  B_rw = (xb_vecN_2xf32*)a_b;
  C_r = (xb_vecN_2xf32*)a_c;
  S_r = (xb_vecN_2xf32*)a_s;
  L_r = (xb_vecNx16*)a_l;
  M_r = a_m + n*LW2;
  for ( p=0; p<L/BBE_SIMD_WIDTH; p++ ) {
    BBE_LBN_2_IP(bm0, M_r, sz_vbn2);
    BBE_LBN_2_IP(bm1, M_r, sz_vbn2);
    /* if (n>a_l[p]) F[SIDX(n-1,p)] = a_c[p]*a_a[p]-a_s[p]*a_b[p]; */
    BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
    BBE_EXTRACTBN(bl1, bl0, BBE_GTNX16(vn,vl));
    a0 = BBE_LVN_2XF32_I(A_rw, 0);
    a1 = BBE_LVN_2XF32_I(A_rw, 2*BBE_SIMD_WIDTH);
    b0 = BBE_LVN_2XF32_I(B_rw, 0);
    b1 = BBE_LVN_2XF32_I(B_rw, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(c0, C_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(c1, C_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(s0, S_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(s1, S_r, 2*BBE_SIMD_WIDTH);
    f0 = BBE_MULN_2XF32(c0,a0); BBE_MULSN_2XF32(f0,s0,b0);
    f1 = BBE_MULN_2XF32(c1,a1); BBE_MULSN_2XF32(f1,s1,b1);
    BBE_SVN_2XF32T_IP(f0, Fm_w, 2*BBE_SIMD_WIDTH, BBE_ANDBN_2(bm0,bl0));
    BBE_SVN_2XF32T_IP(f1, Fm_w, 2*BBE_SIMD_WIDTH, BBE_ANDBN_2(bm1,bl1));
    /* x = D[SIDX(n,p)]; y = F[SIDX(n,p)]; z = D[SIDX(q,p)]; */
    BBE_LVN_2XF32_IP(x0, Dn_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(x1, Dn_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, Fn_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y1, Fn_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(z0, Dq_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(z1, Dq_r, 2*BBE_SIMD_WIDTH);
    /* a_a[p] = a_c[p]*x-a_s[p]*y; */
    a0 = BBE_MULN_2XF32(c0,x0); BBE_MULSN_2XF32(a0,s0,y0);
    a1 = BBE_MULN_2XF32(c1,x1); BBE_MULSN_2XF32(a1,s1,y1);
    BBE_SVN_2XF32T_IP(a0, A_rw, 2*BBE_SIMD_WIDTH, bm0);
    BBE_SVN_2XF32T_IP(a1, A_rw, 2*BBE_SIMD_WIDTH, bm1);
    /* a_b[p] = -a_s[p]*z; */
    b0 = BBE_NEGN_2XF32(BBE_MULN_2XF32(s0,z0));
    b1 = BBE_NEGN_2XF32(BBE_MULN_2XF32(s1,z1));
    BBE_SVN_2XF32T_IP(b0, B_rw, 2*BBE_SIMD_WIDTH, bm0);
    BBE_SVN_2XF32T_IP(b1, B_rw, 2*BBE_SIMD_WIDTH, bm1);
    /* F[SIDX(n,p)] = a_s[p]*x+a_c[p]*y; */
    g0 = BBE_MULN_2XF32(s0,x0); BBE_MULAN_2XF32(g0,c0,y0);
    g1 = BBE_MULN_2XF32(s1,x1); BBE_MULAN_2XF32(g1,c1,y1);
    BBE_SVN_2XF32T_IP(g0, Fn_w, 2*BBE_SIMD_WIDTH, bm0);
    BBE_SVN_2XF32T_IP(g1, Fn_w, 2*BBE_SIMD_WIDTH, bm1);
    /* D[SIDX(q,p)] = a_c[p]*z; */
    h0 = BBE_MULN_2XF32(c0,z0);
    h1 = BBE_MULN_2XF32(c1,z1);
    BBE_SVN_2XF32T_IP(h0, Dq_w, 2*BBE_SIMD_WIDTH, bm0);
    BBE_SVN_2XF32T_IP(h1, Dq_w, 2*BBE_SIMD_WIDTH, bm1);
  } /* p */
  if ((L&(BBE_SIMD_WIDTH/2))!=0) {
    BBE_LBN_2_IP(bm0, M_r, sz_vbn2);
    /* if (n>a_l[p]) F[SIDX(n-1,p)] = a_c[p]*a_a[p]-a_s[p]*a_b[p]; */
    BBE_LVNX16_IP(vl, L_r, 2*BBE_SIMD_WIDTH);
    BBE_EXTRACTBN(bl1, bl0, BBE_GTNX16(vn,vl));
    a0 = BBE_LVN_2XF32_I(A_rw, 0);
    b0 = BBE_LVN_2XF32_I(B_rw, 0);
    BBE_LVN_2XF32_IP(c0, C_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(s0, S_r, 2*BBE_SIMD_WIDTH);
    f0 = BBE_MULN_2XF32(c0,a0); BBE_MULSN_2XF32(f0,s0,b0);
    BBE_SVN_2XF32T_IP(f0, Fm_w, 2*BBE_SIMD_WIDTH, BBE_ANDBN_2(bm0,bl0));
    /* x = D[SIDX(n,p)]; y = F[SIDX(n,p)]; z = D[SIDX(q,p)]; */
    BBE_LVN_2XF32_IP(x0, Dn_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y0, Fn_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(z0, Dq_r, 2*BBE_SIMD_WIDTH);
    /* a_a[p] = a_c[p]*x-a_s[p]*y; */
    a0 = BBE_MULN_2XF32(c0,x0); BBE_MULSN_2XF32(a0,s0,y0);
    BBE_SVN_2XF32T_IP(a0, A_rw, 2*BBE_SIMD_WIDTH, bm0);
    /* a_b[p] = -a_s[p]*z; */
    b0 = BBE_NEGN_2XF32(BBE_MULN_2XF32(s0,z0));
    BBE_SVN_2XF32T_IP(b0, B_rw, 2*BBE_SIMD_WIDTH, bm0);
    /* F[SIDX(n,p)] = a_s[p]*x+a_c[p]*y; */
    g0 = BBE_MULN_2XF32(s0,x0); BBE_MULAN_2XF32(g0,c0,y0);
    BBE_SVN_2XF32T_IP(g0, Fn_w, 2*BBE_SIMD_WIDTH, bm0);
    /* D[SIDX(q,p)] = a_c[p]*z; */
    h0 = BBE_MULN_2XF32(c0,z0);
    BBE_SVN_2XF32T_IP(h0, Dq_w, 2*BBE_SIMD_WIDTH, bm0);
  } /* L */
#else
  float32_t x,y,z;
  int p,q;
  int LW2 = L/(BBE_SIMD_WIDTH/2);

  NASSERT_ALIGN(D, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(F, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_l, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  q = n+1;
  for ( p=0; p<L; p++ ) {
    if (vbn2Get(a_m+n*LW2, p)) {
      /* B <- B*G(a,b) */
      if (n>a_l[p]) F[SIDX(n-1,p)] = a_c[p]*a_a[p]-a_s[p]*a_b[p];
      x = D[SIDX(n,p)]; y = F[SIDX(n,p)]; z = D[SIDX(q,p)];
      a_a[p] = a_c[p]*x-a_s[p]*y; a_b[p] = -a_s[p]*z;
      F[SIDX(n,p)] = a_s[p]*x+a_c[p]*y; D[SIDX(q,p)] = a_c[p]*z;
    } /* a_m */
  } /* p */
#endif
} /* rgrsvdsf_gks_step1() */

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * 2nd updating step of SVD step iteration: B <- G(a,b)'*B.
 * Input:
 *   n            Current position
 *   L            Number of matrices
 *   a_c[L]       G(a,b) cosine values
 *   a_s[L]       G(a,b) sine values
 *   a_k[L]       Right index of working subblocks
 *   a_m[N][L/(BBE_SIMD_WIDTH/4)] (complex variant)
 *   a_m[N][L/(BBE_SIMD_WIDTH/2)] (real variant)
 *                Boolean labels of Working subblocks
 * Input/Output:
 *   a_a[L]       In: B(n,n) Out: B(n,n+1)  
 *   a_b[L]       In: B(n+1,n) Out: B(n,n+2)
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   a_a,a_b,a_c,a_s,a_k,a_m
 *      Must not overlap
 *   a_a,a_b,a_c,a_s,a_k
 *      Must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void rgrsvdsf_gks_step2(
                    float32_t * restrict D,  
                    float32_t * restrict F,
                    float32_t * restrict a_a,
                    float32_t * restrict a_b,
              const float32_t *          a_c,
              const float32_t *          a_s,
              const int16_t   *          a_k,
              const vboolN_2  *          a_m,
              int n, int L )
{
#if 1
  const xb_vecN_2xf32 * restrict Fn_r;
        xb_vecN_2xf32 * restrict Fn_w;
  const xb_vecN_2xf32 * restrict Fq_r;
        xb_vecN_2xf32 * restrict Fq_w;
        xb_vecN_2xf32 * restrict Dn_w;
  const xb_vecN_2xf32 * restrict Dq_r;
        xb_vecN_2xf32 * restrict Dq_w;
        xb_vecN_2xf32 * restrict A_rw;
        xb_vecN_2xf32 * restrict B_rw;
  const xb_vecN_2xf32 * restrict C_r;
  const xb_vecN_2xf32 * restrict S_r;
  const xb_vecNx16    * restrict K_r;
  const vboolN_2      * restrict M_r;

  xb_vecN_2xf32 a0,a1,b0,b1,c0,c1,f0,f1,g0,g1;
  xb_vecN_2xf32 h0,h1,s0,s1,y0,y1,z0,z1,w0,w1;
  xb_vecNx16 vk,vn;
  vboolN_2 bk0,bk1,bm0,bm1;
  int p, LW2 = L/(BBE_SIMD_WIDTH/2);

  NASSERT_ALIGN(D, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(F, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_k, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  vn = BBE_MOVVA16(n+1);
  Fn_r = Fn_w = (xb_vecN_2xf32*)&F[SIDX(n,0)];
  Fq_r = Fq_w = (xb_vecN_2xf32*)&F[SIDX(n+1,0)];
  Dn_w = (xb_vecN_2xf32*)&D[SIDX(n,0)];
  Dq_r = Dq_w = (xb_vecN_2xf32*)&D[SIDX(n+1,0)];
  A_rw = (xb_vecN_2xf32*)a_a;
  B_rw = (xb_vecN_2xf32*)a_b;
  C_r = (xb_vecN_2xf32*)a_c;
  S_r = (xb_vecN_2xf32*)a_s;
  K_r = (xb_vecNx16*)a_k;
  M_r = a_m + n*LW2;
  for ( p=0; p<L/BBE_SIMD_WIDTH; p++ ) {
    BBE_LBN_2_IP(bm0, M_r, sz_vbn2);
    BBE_LBN_2_IP(bm1, M_r, sz_vbn2);
    /* D[SIDX(n,p)] = a_c[p]*a_a[p]-a_s[p]*a_b[p]; */
    a0 = BBE_LVN_2XF32_I(A_rw, 0);
    a1 = BBE_LVN_2XF32_I(A_rw, 2*BBE_SIMD_WIDTH);
    b0 = BBE_LVN_2XF32_I(B_rw, 0);
    b1 = BBE_LVN_2XF32_I(B_rw, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(c0, C_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(c1, C_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(s0, S_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(s1, S_r, 2*BBE_SIMD_WIDTH);
    f0 = BBE_MULN_2XF32(c0,a0); BBE_MULSN_2XF32(f0,s0,b0);
    f1 = BBE_MULN_2XF32(c1,a1); BBE_MULSN_2XF32(f1,s1,b1);
    BBE_SVN_2XF32T_IP(f0, Dn_w, 2*BBE_SIMD_WIDTH, bm0);
    BBE_SVN_2XF32T_IP(f1, Dn_w, 2*BBE_SIMD_WIDTH, bm1);
    /* y = F[SIDX(n,p)]; z = D[SIDX(q,p)]; */
    BBE_LVN_2XF32_IP(y0, Fn_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(y1, Fn_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(z0, Dq_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(z1, Dq_r, 2*BBE_SIMD_WIDTH);
    /* a_a[p] = F[SIDX(n,p)] = a_c[p]*y-a_s[p]*z;  */
    a0 = BBE_MULN_2XF32(c0,y0); BBE_MULSN_2XF32(a0,s0,z0);
    a1 = BBE_MULN_2XF32(c1,y1); BBE_MULSN_2XF32(a1,s1,z1);
    BBE_SVN_2XF32T_IP(a0, A_rw, 2*BBE_SIMD_WIDTH, bm0);
    BBE_SVN_2XF32T_IP(a1, A_rw, 2*BBE_SIMD_WIDTH, bm1);
    BBE_SVN_2XF32T_IP(a0, Fn_w, 2*BBE_SIMD_WIDTH, bm0);
    BBE_SVN_2XF32T_IP(a1, Fn_w, 2*BBE_SIMD_WIDTH, bm1);
    /* D[SIDX(q,p)] = a_s[p]*y+a_c[p]*z; */
    g0 = BBE_MULN_2XF32(s0,y0); BBE_MULAN_2XF32(g0,c0,z0);
    g1 = BBE_MULN_2XF32(s1,y1); BBE_MULAN_2XF32(g1,c1,z1);
    BBE_SVN_2XF32T_IP(g0, Dq_w, 2*BBE_SIMD_WIDTH, bm0);
    BBE_SVN_2XF32T_IP(g1, Dq_w, 2*BBE_SIMD_WIDTH, bm1);
    /* if (n<a_k[p]-1) {
     *   w = F[SIDX(q,p)]; a_b[p] = -a_s[p]*w; F[SIDX(q,p)] = a_c[p]*w;
     * } */
    BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
    BBE_EXTRACTBN(bk1, bk0, BBE_LTNX16(vn,vk));
    bk0 = BBE_ANDBN_2(bk0,bm0);
    bk1 = BBE_ANDBN_2(bk1,bm1);
    BBE_LVN_2XF32T_IP(w0, Fq_r, 2*BBE_SIMD_WIDTH, bk0);
    BBE_LVN_2XF32T_IP(w1, Fq_r, 2*BBE_SIMD_WIDTH, bk1);
    b0 = BBE_NEGN_2XF32(BBE_MULN_2XF32(s0,w0));
    b1 = BBE_NEGN_2XF32(BBE_MULN_2XF32(s1,w1));
    h0 = BBE_MULN_2XF32(c0,w0);
    h1 = BBE_MULN_2XF32(c1,w1);
    BBE_SVN_2XF32T_IP(b0, B_rw, 2*BBE_SIMD_WIDTH, bk0);
    BBE_SVN_2XF32T_IP(b1, B_rw, 2*BBE_SIMD_WIDTH, bk1);
    BBE_SVN_2XF32T_IP(h0, Fq_w, 2*BBE_SIMD_WIDTH, bk0);
    BBE_SVN_2XF32T_IP(h1, Fq_w, 2*BBE_SIMD_WIDTH, bk1);
  } /* p */
  if ((L&(BBE_SIMD_WIDTH/2))!=0) {
    BBE_LBN_2_IP(bm0, M_r, sz_vbn2);
    /* D[SIDX(n,p)] = a_c[p]*a_a[p]-a_s[p]*a_b[p]; */
    a0 = BBE_LVN_2XF32_I(A_rw, 0);
    b0 = BBE_LVN_2XF32_I(B_rw, 0);
    BBE_LVN_2XF32_IP(c0, C_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(s0, S_r, 2*BBE_SIMD_WIDTH);
    f0 = BBE_MULN_2XF32(c0,a0); BBE_MULSN_2XF32(f0,s0,b0);
    BBE_SVN_2XF32T_IP(f0, Dn_w, 2*BBE_SIMD_WIDTH, bm0);
    /* y = F[SIDX(n,p)]; z = D[SIDX(q,p)]; */
    BBE_LVN_2XF32_IP(y0, Fn_r, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(z0, Dq_r, 2*BBE_SIMD_WIDTH);
    /* a_a[p] = F[SIDX(n,p)] = a_c[p]*y-a_s[p]*z;  */
    a0 = BBE_MULN_2XF32(c0,y0); BBE_MULSN_2XF32(a0,s0,z0);
    BBE_SVN_2XF32T_IP(a0, A_rw, 2*BBE_SIMD_WIDTH, bm0);
    BBE_SVN_2XF32T_IP(a0, Fn_w, 2*BBE_SIMD_WIDTH, bm0);
    /* D[SIDX(q,p)] = a_s[p]*y+a_c[p]*z; */
    g0 = BBE_MULN_2XF32(s0,y0); BBE_MULAN_2XF32(g0,c0,z0);
    BBE_SVN_2XF32T_IP(g0, Dq_w, 2*BBE_SIMD_WIDTH, bm0);
    /* if (n<a_k[p]-1) {
     *   w = F[SIDX(q,p)]; a_b[p] = -a_s[p]*w; F[SIDX(q,p)] = a_c[p]*w;
     * } */
    BBE_LVNX16_IP(vk, K_r, 2*BBE_SIMD_WIDTH);
    BBE_EXTRACTBN(bk1, bk0, BBE_LTNX16(vn,vk));
    bk0 = BBE_ANDBN_2(bk0,bm0);
    BBE_LVN_2XF32T_IP(w0, Fq_r, 2*BBE_SIMD_WIDTH, bk0);
    b0 = BBE_NEGN_2XF32(BBE_MULN_2XF32(s0,w0));
    h0 = BBE_MULN_2XF32(c0,w0);
    BBE_SVN_2XF32T_IP(b0, B_rw, 2*BBE_SIMD_WIDTH, bk0);
    BBE_SVN_2XF32T_IP(h0, Fq_w, 2*BBE_SIMD_WIDTH, bk0);
  }/* L */
#else
  float32_t y,z,w;
  int p,q;
  int LW2 = L/(BBE_SIMD_WIDTH/2);
  
  NASSERT_ALIGN(D, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(F, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_a, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_b, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_k, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));
  
  q = n+1;
  for ( p=0; p<L; p++ ) {
    if (vbn2Get(a_m+n*LW2, p)) {
      /* B <- G(a,b)'*B */
      D[SIDX(n,p)] = a_c[p]*a_a[p]-a_s[p]*a_b[p];
      y = F[SIDX(n,p)]; z = D[SIDX(q,p)];
      a_a[p] = F[SIDX(n,p)] = a_c[p]*y-a_s[p]*z; 
      D[SIDX(q,p)] = a_s[p]*y+a_c[p]*z;
      if (n<a_k[p]-1) {
        w = F[SIDX(q,p)]; a_b[p] = -a_s[p]*w; F[SIDX(q,p)] = a_c[p]*w;
      }
    } /* a_m */
  } /* p */
#endif

} /* rgrsvdsf_gks_step2() */

/*
 * Internal function of Golub-Kahan SVD step implementation.
 * Accumulate transformations: U <- U*conj(G(a,b)) or V <- V*G(a,b),
 * where U (MxN) and V (NxN) are matrices comprised of left- and right
 * singular vectors. For real data, a single function is used for both
 * U and V updates. For complex data there are two separate functions
 * for matrices U and V.
 * Input:
 *   n            Current position
 *   M,N          Matrix dimensions
 *   L            Number of matrices
 *   a_c[L]       G(a,b) cosine values
 *   a_s[L]       G(a,b) sine values
 *   a_m[N][L/(BBE_SIMD_WIDTH/4)] (complex variant)
 *   a_m[N][L/(BBE_SIMD_WIDTH/2)] (real variant)
 *                Boolean labels of Working subblocks
 * Input/Output:
 *   W[M*N][L]    Matrix of left-singular (U) or right-singular (V) 
 *                orthonormal vectors
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2
 *   W,a_c,a_s,a_m
 *      Must not overlap
 *   W,a_c,a_s
 *      Must be 2*BBE_SIMD_WIDTH-byte aligned
 */

void rgrsvdsf_gks_accum(
                    float32_t * restrict W,
              const float32_t *          a_c,
              const float32_t *          a_s,
              const vboolN_2  *          a_m,
              int n, int M, int N, int L )
{
#if 1
  const xb_vecN_2xf32 * restrict Wn_r;
        xb_vecN_2xf32 * restrict Wn_w;
  const xb_vecN_2xf32 * restrict Wp_r;
        xb_vecN_2xf32 * restrict Wp_w;
  const xb_vecN_2xf32 * restrict C_r;
  const xb_vecN_2xf32 * restrict S_r;
  const vboolN_2      * restrict M_r;

  xb_vecN_2xf32 c,s,x,y,z,w;
  vboolN_2 bm;
  int i,p,LW2 = L/(BBE_SIMD_WIDTH/2);

  NASSERT_ALIGN(W, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  for ( i=0; i<M; i++ ) {
    Wn_r = Wn_w = (xb_vecN_2xf32*)&W[SIDX(i*N+n,0)];
    Wp_r = Wp_w = (xb_vecN_2xf32*)&W[SIDX(i*N+n+1,0)];
    C_r = (xb_vecN_2xf32*)a_c;
    S_r = (xb_vecN_2xf32*)a_s;
    M_r = a_m + n*LW2;
    for ( p=0; p<L/(BBE_SIMD_WIDTH/2); p++ ) {
      BBE_LBN_2_IP(bm, M_r, sz_vbn2);
      BBE_LVN_2XF32_IP(c, C_r, 2*BBE_SIMD_WIDTH);
      BBE_LVN_2XF32_IP(s, S_r, 2*BBE_SIMD_WIDTH);
      /* x = W[SIDX(i*N+n,p)]; y = W[SIDX(i*N+q,p)]; */
      BBE_LVN_2XF32_IP(x, Wn_r, 2*BBE_SIMD_WIDTH);  
      BBE_LVN_2XF32_IP(y, Wp_r, 2*BBE_SIMD_WIDTH);
      /* W[SIDX(i*N+n,p)] = a_c[p]*x-a_s[p]*y; */
      z = BBE_MULN_2XF32(c,x); BBE_MULSN_2XF32(z,s,y);
      BBE_SVN_2XF32T_IP(z, Wn_w, 2*BBE_SIMD_WIDTH, bm);
      /* W[SIDX(i*N+q,p)] = a_s[p]*x+a_c[p]*y; */
      w = BBE_MULN_2XF32(s,x); BBE_MULAN_2XF32(w,c,y);
      BBE_SVN_2XF32T_IP(w, Wp_w, 2*BBE_SIMD_WIDTH, bm);
    } /* p */
  } /* i */
#else
  float32_t x,y;
  int i,p,q;
  int LW2 = L/(BBE_SIMD_WIDTH/2);

  NASSERT_ALIGN(W, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_c, 2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(a_s, 2*BBE_SIMD_WIDTH);
  NASSERT(0==(L%(BBE_SIMD_WIDTH/2)));

  q = n+1;
  for ( i=0; i<M; i++ ) {
    for ( p=0; p<L; p++ ) {
      if (vbn2Get(a_m+n*LW2, p)) {
        /* W <- W*G(a,b) */
        x = W[SIDX(i*N+n,p)]; y = W[SIDX(i*N+q,p)];
        W[SIDX(i*N+n,p)] = a_c[p]*x-a_s[p]*y;
        W[SIDX(i*N+q,p)] = a_s[p]*x+a_c[p]*y;
      }
    } /* a_m */
  } /* p */
#endif
} /* rgrsvdsf_gks_accum() */

#endif /* HAVE_VFPU */
