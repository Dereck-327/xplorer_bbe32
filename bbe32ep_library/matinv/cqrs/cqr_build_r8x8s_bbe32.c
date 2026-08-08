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
/*          Copyright (C) 2009-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
  NatureDSP_Baseband library. QR-based matrix decomposition and inversion for streaming order
    cqr_build_rMxNs
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "cqr_common.h"

#if (HAVE_VSAMATH && HAVE_NSAENX40 && 1)

#define CURRENT_M 8

#define DEF_CHOUSEHOLDER_FN3_BODY                                   \
        int i = 0 ;                                                 \
        xb_vecNx16 k_norm,tmp0,tmp1;                                \
        xb_vecNx16 norm_x, Fi;                                      \
        xb_vecNx40 acc0;                                            \
        vsaN c_vec,v_exp;                                           \
        xb_vecNx16  x1, x2;                                         \
                                                                    \
        BBE_LVNX16_XP(x1, px, stride_bytes);                        \
        acc0 = BBE_MAGINX16C(x1, x1);                               \
        acc0=BBE_ADDNX40(acc0,acc0);                                \
        c_vec=BBE_NSAENX40(acc0);                                   \
        acc0=BBE_SLLNX40(acc0,c_vec);                               \
        BBE_RSQRTLUNX40_0(acc0,tmp0, tmp1, acc0);                   \
        BBE_MULUUSNX16( acc0, tmp1,  tmp0);                         \
        acc0=BBE_SRAINX40(acc0,24);                                 \
        tmp0=BBE_PACKLNX40(acc0);                                   \
        tmp0= BBE_SHFLNX16I(tmp0, BBE_SHFLI_DUPLICATE_1_EVEN);      \
        c_vec = BBE_SUBSR1SAVSN(18+1,c_vec);                        \
        acc0 = BBE_MULUSRNX16(tmp0,x1, c_vec);                      \
        norm_x = BBE_PACKVNX40(acc0, c_vec);                        \
        Fi = norm_x;                                                \
                                                                    \
        BBE_SVNX16_XP(norm_x, pfi, BBE_SIMD_WIDTH*sizeof(int16_t)); \
                                                                    \
        /* Build Householder's vector         */                    \
        /*v = x/sqrt(x'*x) + e1*x(1)/abs(x(1) */                    \
        acc0 = BBE_MAGINX16C( x1, x1);                              \
        if ((len_x&1)==0)                                           \
        {                                                           \
            BBE_LVNX16_XP(x2, px, stride_bytes);                    \
            BBE_MAGIANX16C(acc0, x2, x2);                           \
        }                                                           \
        for(i=2-(len_x&1); i<len_x; i+=2)                           \
        {                                                           \
          BBE_LVNX16_XP(x1, px, stride_bytes);                      \
          BBE_LVNX16_XP(x2, px, stride_bytes);                      \
          BBE_MAGIANX16C(acc0, x1, x1);                             \
          BBE_MAGIANX16C(acc0, x2, x2);                             \
        }                                                           \
        px -= (stride_bytes*(CURRENT_M-m))/sizeof(*px);             \
        BBE_LVNX16_XP(x1, px, stride_bytes);                        \
        acc0=BBE_ADDNX40(acc0,acc0);                                \
        c_vec=BBE_NSAENX40(acc0);                                   \
        acc0=BBE_SLLNX40(acc0,c_vec);                               \
        BBE_RSQRTLUNX40_0(acc0,tmp0, tmp1, acc0);                   \
        BBE_MULUUSNX16( acc0, tmp1,  tmp0);                         \
        acc0=BBE_SRAINX40(acc0,24);                                 \
        tmp0=BBE_PACKLNX40(acc0);                                   \
        mant= BBE_SHFLNX16I(tmp0, BBE_SHFLI_DUPLICATE_1_EVEN);      \
        v_exp = BBE_SUBSR1SAVSN(18+1,c_vec);                        \
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);                      \
        norm_x = BBE_PACKVNX40(acc0, v_exp);                        \
                                                                    \
        norm_x = BBE_ADDNX16(norm_x, Fi);                           \
                                                                    \
        acc0 = BBE_MULNX16J( norm_x, Fi);                           \
        acc0=BBE_ADDNX40(acc0,acc0);                                \
        c_vec=BBE_NSAENX40(acc0);                                   \
        acc0=BBE_SLLNX40(acc0,c_vec);                               \
        BBE_RSQRTLUNX40_0(acc0,tmp0, tmp1, acc0);                   \
        BBE_MULUUSNX16(acc0, tmp1,  tmp0);                          \
        c_vec = BBE_SUBSR1SAVSN(28,c_vec);                          \
        k_norm=BBE_PACKVNX40(acc0,c_vec);                           \
        k_norm = BBE_SHFLNX16I(k_norm, BBE_SHFLI_DUPLICATE_1_EVEN); \
                                                                    \
                                                                    \
        acc0 = BBE_MULNX16( norm_x, k_norm);                        \
        norm_x = BBE_PACKQNX40(acc0);                               \
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));              \

#if 1
inline_ void chouseholder_M8xL_m0(
                                  const int16_t *   A,        /* input matrix              */
                                  int16_t *            v ,        /* output streaming order */
                                  int16_t*            Fi_out, /* output phase rotator      */
                                  const int L
                                  )
{
    const int m = 0;
    int l;
    xb_vecNx16 mant;
    const int stride_bytes = L*2*CURRENT_M*sizeof(int16_t);
    const xb_vecNx16  * restrict px = (xb_vecNx16* )(A + L*2*(CURRENT_M+1)*m); 
    xb_vecNx16  * restrict pv  =  (xb_vecNx16*)(v);
    xb_vecNx16  * restrict pfi = (xb_vecNx16  *)Fi_out; 
    const int len_x = CURRENT_M-m; 

    NASSERT_ALIGN(px,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(v,2*BBE_SIMD_WIDTH);
    NASSERT(len_x>1 );
    NASSERT( CURRENT_M>1);
    NASSERT( m < CURRENT_M-1);

   __Pragma("loop_count min=1")
    for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2))
    {
        DEF_CHOUSEHOLDER_FN3_BODY;
        for(i=1; i<len_x; i++)
        {
            BBE_LVNX16_XP(x1, px, stride_bytes);
            acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
            norm_x = BBE_PACKVNX40(acc0, v_exp);
            norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
            BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        }
        px -= (stride_bytes*(CURRENT_M-m))/sizeof(*px)-1;
    } /* for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2)) */
}

inline_ void chouseholder_M8xL_m1(
                                  const int16_t *   A,        /* input matrix              */
                                  int16_t *            v ,        /* output streaming order */
                                  int16_t*            Fi_out, /* output phase rotator      */
                                  const int L
                                  )
{
  const int m = 1;
  int l;
  xb_vecNx16 mant;
  const int stride_bytes = L*2*CURRENT_M*sizeof(int16_t);
  const xb_vecNx16  * restrict px = (xb_vecNx16* )(A + L*2*(CURRENT_M+1)*m); 
  xb_vecNx16  * restrict pv  =  (xb_vecNx16*)(v);
  xb_vecNx16  * restrict pfi = (xb_vecNx16  *)Fi_out; 
  const int len_x = CURRENT_M-m; 

  NASSERT_ALIGN(px,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(v,2*BBE_SIMD_WIDTH);
  NASSERT(len_x>1 );
  NASSERT( CURRENT_M>1);
  NASSERT( m < CURRENT_M-1);

   __Pragma("loop_count min=1")
  for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2))
  {
        DEF_CHOUSEHOLDER_FN3_BODY;
        for(i=1; i<len_x; i++)
        {
            BBE_LVNX16_XP(x1, px, stride_bytes);
            acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
            norm_x = BBE_PACKVNX40(acc0, v_exp);
            norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
            BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        }
        px -= (stride_bytes*(CURRENT_M-m))/sizeof(*px)-1;
  } /* for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2)) */
}

inline_ void chouseholder_M8xL_m2(
                                  const int16_t *   A,        /* input matrix              */
                                  int16_t *            v ,        /* output streaming order */
                                  int16_t*            Fi_out, /* output phase rotator      */
                                  const int L
                                  )
{
  const int m = 2;
  int l;
  xb_vecNx16 mant;
  const int stride_bytes = L*2*CURRENT_M*sizeof(int16_t);
  const xb_vecNx16  * restrict px = (xb_vecNx16* )(A + L*2*(CURRENT_M+1)*m); 
  xb_vecNx16  * restrict pv  =  (xb_vecNx16*)(v);
  xb_vecNx16  * restrict pfi = (xb_vecNx16  *)Fi_out; 
  const int len_x = CURRENT_M-m; 

  NASSERT_ALIGN(px,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(v,2*BBE_SIMD_WIDTH);
  NASSERT(len_x>1 );
  NASSERT( CURRENT_M>1);
  NASSERT( m < CURRENT_M-1);

   __Pragma("loop_count min=1")
  for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2))
  {
        DEF_CHOUSEHOLDER_FN3_BODY;
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        px -= (stride_bytes*(CURRENT_M-m))/sizeof(*px)-1;
  } /* for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2)) */
}

inline_ void chouseholder_M8xL_m3(
                                  const int16_t *   A,        /* input matrix              */
                                  int16_t *            v ,        /* output streaming order */
                                  int16_t*            Fi_out, /* output phase rotator      */
                                  const int L
                                  )
{
  const int m = 3;
  int l;
  xb_vecNx16 mant;
  const int stride_bytes = L*2*CURRENT_M*sizeof(int16_t);
  const xb_vecNx16  * restrict px = (xb_vecNx16* )(A + L*2*(CURRENT_M+1)*m); 
  xb_vecNx16  * restrict pv  =  (xb_vecNx16*)(v);
  xb_vecNx16  * restrict pfi = (xb_vecNx16  *)Fi_out; 
  const int len_x = CURRENT_M-m; 

  NASSERT_ALIGN(px,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(v,2*BBE_SIMD_WIDTH);
  NASSERT(len_x>1 );
  NASSERT( CURRENT_M>1);
  NASSERT( m < CURRENT_M-1);

   __Pragma("loop_count min=1")
  for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2))
  {
        DEF_CHOUSEHOLDER_FN3_BODY;
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        px -= (stride_bytes*(CURRENT_M-m))/sizeof(*px)-1;
  } /* for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2)) */
}

inline_ void chouseholder_M8xL_m4(
                                  const int16_t *   A,        /* input matrix              */
                                  int16_t *            v ,        /* output streaming order */
                                  int16_t*            Fi_out, /* output phase rotator      */
                                  const int L
                                  )
{
  const int m = 4;
  int l;
  xb_vecNx16 mant;
  const int stride_bytes = L*2*CURRENT_M*sizeof(int16_t);
  const xb_vecNx16  * restrict px = (xb_vecNx16* )(A + L*2*(CURRENT_M+1)*m); 
  xb_vecNx16  * restrict pv  =  (xb_vecNx16*)(v);
  xb_vecNx16  * restrict pfi = (xb_vecNx16  *)Fi_out; 
  const int len_x = CURRENT_M-m; 

  NASSERT_ALIGN(px,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(v,2*BBE_SIMD_WIDTH);
  NASSERT(len_x>1 );
  NASSERT( CURRENT_M>1);
  NASSERT( m < CURRENT_M-1);

   __Pragma("loop_count min=1")
  for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2))
  {
        DEF_CHOUSEHOLDER_FN3_BODY;
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        px -= (stride_bytes*(CURRENT_M-m))/sizeof(*px)-1;
  } /* for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2)) */
}

inline_ void chouseholder_M8xL_m5(
                                  const int16_t *   A,        /* input matrix              */
                                  int16_t *            v ,        /* output streaming order */
                                  int16_t*            Fi_out, /* output phase rotator      */
                                  const int L
                                  )
{
  const int m = 5;
  int l;
  xb_vecNx16 mant;
  const int stride_bytes = L*2*CURRENT_M*sizeof(int16_t);
  const xb_vecNx16  * restrict px = (xb_vecNx16* )(A + L*2*(CURRENT_M+1)*m); 
  xb_vecNx16  * restrict pv  =  (xb_vecNx16*)(v);
  xb_vecNx16  * restrict pfi = (xb_vecNx16  *)Fi_out; 
  const int len_x = CURRENT_M-m; 

  NASSERT_ALIGN(px,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(v,2*BBE_SIMD_WIDTH);
  NASSERT(len_x>1 );
  NASSERT( CURRENT_M>1);
  NASSERT( m < CURRENT_M-1);

   __Pragma("loop_count min=1")
  for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2))
  {
        DEF_CHOUSEHOLDER_FN3_BODY;
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        px -= (stride_bytes*(CURRENT_M-m))/sizeof(*px)-1;
  } /* for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2)) */
}

inline_ void chouseholder_M8xL_m6(
                                  const int16_t *   A,        /* input matrix              */
                                  int16_t *            v ,        /* output streaming order */
                                  int16_t*            Fi_out, /* output phase rotator      */
                                  const int L
                                  )
{
  const int m = 6;
  int l;
  xb_vecNx16 mant;
  const int stride_bytes = L*2*CURRENT_M*sizeof(int16_t);
  const xb_vecNx16  * restrict px = (xb_vecNx16* )(A + L*2*(CURRENT_M+1)*m); 
  xb_vecNx16  * restrict pv  =  (xb_vecNx16*)(v);
  xb_vecNx16  * restrict pfi = (xb_vecNx16  *)Fi_out; 
  const int len_x = CURRENT_M-m; 

  NASSERT_ALIGN(px,2*BBE_SIMD_WIDTH);
  NASSERT_ALIGN(v,2*BBE_SIMD_WIDTH);
  NASSERT(len_x>1 );
  NASSERT( CURRENT_M>1);
  NASSERT( m < CURRENT_M-1);

   __Pragma("loop_count min=1")
  for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2))
  {
        DEF_CHOUSEHOLDER_FN3_BODY;
        BBE_LVNX16_XP(x1, px, stride_bytes);
        acc0 = BBE_MULUSRNX16(mant,x1, v_exp);
        norm_x = BBE_PACKVNX40(acc0, v_exp);
        norm_x = BBE_MULNX16PACKQ(norm_x, k_norm);
        BBE_SVNX16_IP(norm_x, pv, (2*BBE_SIMD_WIDTH));
        px -= (stride_bytes*(CURRENT_M-m))/sizeof(*px)-1;
  } /* for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2)) */
}


#else
DEF_CHOUSEHOLDER_FN2(8, 0)
DEF_CHOUSEHOLDER_FN2(8, 1)
DEF_CHOUSEHOLDER_FN2(8, 2)
DEF_CHOUSEHOLDER_FN2(8, 3)
DEF_CHOUSEHOLDER_FN2(8, 4)
DEF_CHOUSEHOLDER_FN2(8, 5)
DEF_CHOUSEHOLDER_FN2(8, 6)
#endif

#if 0
void cqr_reorder_V( void *scr,        /* CQR_SIZE_V bytes */
                    int16_t *V,       /* inOut */
                    int      M,
                    int      N,
                    int      L);

#define MAX_L 64

ALIGN(32) static int16_t tmp_buf[CQR_SIZE_V(8,8,MAX_L)/2];
#endif

/*
    Calculate product  vR = v' * R ,
    v - Householder's vectors, LxMx1,
    R - matrix MxPxL,
    vR -vector Px1xL,
    vR and R in a streaming format
    v - in a blocking format
*/

#define CALCVR8x8(vR,R,V,L,P,M,row_stride,q,pragma_text)   \
{                                                \
    vsaN _q=BBE_MOVVSA32(q);                     \
    int i;                                       \
    xb_vecNx16 r0, v0, vR0;                      \
    xb_vecNx40 acc0;                             \
    pR0 = (xb_vecNx16 *)R;                       \
    pV0 = (xb_vecNx16 *)V;                       \
    pVR = (xb_vecNx16 *)vR;                      \
    NASSERT_ALIGN(vR,2*BBE_SIMD_WIDTH);;         \
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);;          \
    NASSERT_ALIGN(vR,2*BBE_SIMD_WIDTH);;         \
    NASSERT(L%(BBE_SIMD_WIDTH/2) == 0 && L > 0); \
    NASSERT(M > 1 && P>=1);                      \
    WUR_CBEGIN( (unsigned) V);                   \
    WUR_CEND  ( (unsigned) (V + M*2*L));         \
    __Pragma(pragma_text)                        \
    for(i=0; i < P*(L/(BBE_SIMD_WIDTH/2)); i++)  \
    {                                            \
                  BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);acc0 =BBE_MULRNX16J(r0, v0,_q);  \
        if(M>2) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>3) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>4) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>5) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>6) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        if(M>7) { BBE_LVNX16_XP(r0, pR0, row_stride); BBE_LVNX16_IC(v0, pV0);BBE_MULANX16J(acc0, r0, v0); }   \
        BBE_LVNX16_XP(r0, pR0, -row_stride*(M-1) + 2*BBE_SIMD_WIDTH );                                        \
        BBE_LVNX16_IC(v0, pV0);                                                                               \
        BBE_MULANX16J(acc0, r0, v0);                                                                          \
        vR0 = BBE_PACKVNX40(acc0, _q);   /* Q13*/                                                             \
        BBE_SVNX16_IP(vR0, pVR, 2*BBE_SIMD_WIDTH);                                                            \
    }                                                                                                         \
}

inline_ void RotateRows8x8( int16_t * R, /*(io)*/
                         const int16_t * Fi,/* (i)*/
                         int L,
                         const int M,
                         const int P)
#if 1
{
    int i, j, l;
    xb_vecNx40 acc0;
    vsaN _14=BBE_MOVVSA32(14);
    xb_vecNx16* pFi0 = (xb_vecNx16 *)Fi;
    xb_vecNx16* pFi1 = (xb_vecNx16 *)(Fi + 1*2*L*M/4);
    xb_vecNx16* pFi2 = (xb_vecNx16 *)(Fi + 2*2*L*M/4);
    xb_vecNx16* pFi3 = (xb_vecNx16 *)(Fi + 3*2*L*M/4);

    xb_vecNx16* pR0 = (xb_vecNx16 *)R, *pR1, *pR2, *pR3;
    xb_vecNx16 r0, f0, f1, f2, f3, r1, r2, r3;

    NASSERT_ALIGN(Fi,2*BBE_SIMD_WIDTH);;
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);;

    ASSERT(M%4==0);

     pR0  =  (xb_vecNx16 *)(R +         0*2*L*P);
     pR1  =  (xb_vecNx16 *)(R + (0+1*M/4)*2*L*P);
     pR2  =  (xb_vecNx16 *)(R + (0+2*M/4)*2*L*P);
     pR3  =  (xb_vecNx16 *)(R + (0+3*M/4)*2*L*P);

#ifdef COMPILER_XTENSA
#pragma ymemory(pR0)
#pragma ymemory(pR1)
#pragma ymemory(pR2)
#pragma ymemory(pR3)
#endif
    for(i=0; i<M/4; i++)
    {
        for(l=0; l<L; l+=BBE_SIMD_WIDTH/2)
        {

            BBE_LVNX16_IP(f0, pFi0, 2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(f1, pFi1, 2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(f2, pFi2, 2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(f3, pFi3, 2*BBE_SIMD_WIDTH);

            for(j=0; j<P-1; j++)
            {
                BBE_LVNX16_XP(r0, pR0, 0);
                BBE_LVNX16_XP(r1, pR1, 0);
                BBE_LVNX16_XP(r2, pR2, 0);
                BBE_LVNX16_XP(r3, pR3, 0);

                acc0 = BBE_MULRNX16J(r0, f0,_14); r0=BBE_PACKVNX40(acc0,_14);
                acc0 = BBE_MULRNX16J(r1, f1,_14); r1=BBE_PACKVNX40(acc0,_14);
                acc0 = BBE_MULRNX16J(r2, f2,_14); r2=BBE_PACKVNX40(acc0,_14);
                acc0 = BBE_MULRNX16J(r3, f3,_14); r3=BBE_PACKVNX40(acc0,_14);

                BBE_SVNX16_XP(r0, pR0, 4*L);
                BBE_SVNX16_XP(r1, pR1, 4*L);
                BBE_SVNX16_XP(r2, pR2, 4*L);
                BBE_SVNX16_XP(r3, pR3, 4*L);
            }

            BBE_LVNX16_XP(r0, pR0, 0);
            BBE_LVNX16_XP(r1, pR1, 0);
            BBE_LVNX16_XP(r2, pR2, 0);
            BBE_LVNX16_XP(r3, pR3, 0);

            acc0 = BBE_MULRNX16J(r0, f0,_14); r0=BBE_PACKVNX40(acc0,_14);
            acc0 = BBE_MULRNX16J(r1, f1,_14); r1=BBE_PACKVNX40(acc0,_14);
            acc0 = BBE_MULRNX16J(r2, f2,_14); r2=BBE_PACKVNX40(acc0,_14);
            acc0 = BBE_MULRNX16J(r3, f3,_14); r3=BBE_PACKVNX40(acc0,_14);

            BBE_SVNX16_XP(r0, pR0, 2*BBE_SIMD_WIDTH - (P-1)*4*L);
            BBE_SVNX16_XP(r1, pR1, 2*BBE_SIMD_WIDTH - (P-1)*4*L);
            BBE_SVNX16_XP(r2, pR2, 2*BBE_SIMD_WIDTH - (P-1)*4*L);
            BBE_SVNX16_XP(r3, pR3, 2*BBE_SIMD_WIDTH - (P-1)*4*L);
        }
        // Adjust pointers to begin of a next rows
        pR0 = (xb_vecNx16*)( (P-1)*2*L*sizeof(int16_t) + (uintptr_t)pR0 );
        pR1 = (xb_vecNx16*)( (P-1)*2*L*sizeof(int16_t) + (uintptr_t)pR1 );
        pR2 = (xb_vecNx16*)( (P-1)*2*L*sizeof(int16_t) + (uintptr_t)pR2 );
        pR3 = (xb_vecNx16*)( (P-1)*2*L*sizeof(int16_t) + (uintptr_t)pR3 );
    }
}
#else
{
    int i, j, l;
    xb_vecNx40 acc0;
    vsaN _14=BBE_MOVVSA32(14);
    xb_vecNx16* pFi0 = (xb_vecNx16 *)Fi;
    xb_vecNx16* pR0 = (xb_vecNx16 *)R, *pR1, *pR2, *pR3;
    xb_vecNx16 r0, f0, f1, f2, f3, r1, r2, r3;

    NASSERT_ALIGN(Fi,2*BBE_SIMD_WIDTH);;
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);;

    ASSERT(M%4==0);

     pR0  =  (xb_vecNx16 *)(R +         0*2*L*P);
     pR1  =  (xb_vecNx16 *)(R + (0+1*M/4)*2*L*P);
     pR2  =  (xb_vecNx16 *)(R + (0+2*M/4)*2*L*P);
     pR3  =  (xb_vecNx16 *)(R + (0+3*M/4)*2*L*P);

#ifdef COMPILER_XTENSA
#pragma ymemory(pR0)
#pragma ymemory(pR1)
#pragma ymemory(pR2)
#pragma ymemory(pR3)
#endif
    for(i=0; i<M/4; i++)
    {
        int finc=0; 
        int rinc=4*L;
        #ifdef COMPILER_XTENSA
        #pragma loop_count min=8
        #endif
        for(j=1,l=0; l<L*P; l+=BBE_SIMD_WIDTH/2)
        {
                XT_MOVEQZ(finc,2*BBE_SIMD_WIDTH,j);
                XT_MOVEQZ(rinc,2*BBE_SIMD_WIDTH - (P-1)*4*L,j);
                f1=BBE_LVNX16_X ( pFi0, 1*4*L*M/4);
                f2=BBE_LVNX16_X ( pFi0, 2*4*L*M/4);
                f3=BBE_LVNX16_X ( pFi0, 3*4*L*M/4);
                BBE_LVNX16_XP(f0, pFi0, finc);
                finc=0;

                BBE_LVNX16_XP(r0, pR0, 0);
                BBE_LVNX16_XP(r1, pR1, 0);
                BBE_LVNX16_XP(r2, pR2, 0);
                BBE_LVNX16_XP(r3, pR3, 0);

                acc0 = BBE_MULRNX16J(r0, f0,_14); r0=BBE_PACKVNX40(acc0,_14);
                acc0 = BBE_MULRNX16J(r1, f1,_14); r1=BBE_PACKVNX40(acc0,_14);
                acc0 = BBE_MULRNX16J(r2, f2,_14); r2=BBE_PACKVNX40(acc0,_14);
                acc0 = BBE_MULRNX16J(r3, f3,_14); r3=BBE_PACKVNX40(acc0,_14);

                BBE_SVNX16_XP(r0, pR0, rinc);
                BBE_SVNX16_XP(r1, pR1, rinc);
                BBE_SVNX16_XP(r2, pR2, rinc);
                BBE_SVNX16_XP(r3, pR3, rinc);
                rinc=4*L;
                j=BBE_ADDMOD16U(j,(8<<16)|1);
        }
        // Adjust pointers to begin of a next rows
        pR0 = (xb_vecNx16*)( (P-1)*2*L*sizeof(int16_t) + (uintptr_t)pR0 );
        pR1 = (xb_vecNx16*)( (P-1)*2*L*sizeof(int16_t) + (uintptr_t)pR1 );
        pR2 = (xb_vecNx16*)( (P-1)*2*L*sizeof(int16_t) + (uintptr_t)pR2 );
        pR3 = (xb_vecNx16*)( (P-1)*2*L*sizeof(int16_t) + (uintptr_t)pR3 );
    }
}
#endif

#define Q 14
#undef MULCPACKQ_SUB
#define MULCPACKQ_SUB(inout, x, y)          \
{                                           \
   xb_vecNx16 _0x4000=BBE_MOVVA16(1<<Q);    \
   vsaN _14=BBE_MOVVSA32(Q);                \
   xb_vecNx40 acc0;                         \
   acc0=BBE_MULRNX16(inout,_0x4000,_14);    \
   BBE_MULSNX16C(acc0,x,y);                 \
   inout = BBE_PACKVNX40(acc0,_14);         \
}


#define UPDATE_X2(M,_t,row_stride)               \
if (M>_t)                                        \
{                                                \
            BBE_LVNX16_XP(r0, pX0,  row_stride); \
            BBE_LVNX16_XP(r1, pX1,  row_stride); \
            BBE_LVNX16_IC(v0, pV0);              \
            MULCPACKQ_SUB(r0, v0, vR0);          \
            MULCPACKQ_SUB(r1, v0, vR1);          \
            BBE_SVNX16_XP(r0, pY0,(row_stride)); \
            BBE_SVNX16_XP(r1, pY1,(row_stride)); \
}

#define UPDATE_X1(M,_t,row_stride)               \
if (M>_t)                                        \
{                                                \
            BBE_LVNX16_XP(r1, pX1,(row_stride)); \
            BBE_LVNX16_IC(v0, pV0);              \
            MULCPACKQ_SUB(r1, v0, vR1);          \
            BBE_SVNX16_XP(r1, pY1,(row_stride)); \
}


#define UPDATE8x8(X,Y,vR,V,L,P,M,row_stride,pragma_text)       \
{                                                              \
    int i;                                                     \
    xb_vecNx16 r0, v0, vR0, vR1, r1;                           \
    pX0 = (xb_vecNx16 *)(X);                                   \
    pY0 = (xb_vecNx16 *)(Y);                                   \
    pV0 = (xb_vecNx16 *)(V);                                   \
    pY1 = (xb_vecNx16 *)(Y+P/2*2*L);                           \
    pX1 = (xb_vecNx16 *)(X+P/2*2*L);                           \
    pVR  = (xb_vecNx16 *)vR;                                   \
    pVR1 = (xb_vecNx16 *)(vR+P/2*2*L);                         \
    NASSERT_ALIGN(vR,2*BBE_SIMD_WIDTH);;                       \
    NASSERT_ALIGN(X,2*BBE_SIMD_WIDTH);                         \
    NASSERT_ALIGN(vR,2*BBE_SIMD_WIDTH);;                       \
    NASSERT(L%(BBE_SIMD_WIDTH/2) == 0 && L > 0);               \
    NASSERT(M > 1 && P>=1);                                    \
                                                               \
    WUR_CBEGIN( (unsigned) V);                                 \
    WUR_CEND  ( (unsigned) (V + M*2*L));                       \
    __Pragma(pragma_text)                                      \
    for(i=0; i < (P/2)*(L/(BBE_SIMD_WIDTH/2)); i++)            \
    {                                                          \
        BBE_LVNX16_IP(vR0, pVR , 2*BBE_SIMD_WIDTH);            \
        BBE_LVNX16_IP(vR1, pVR1, 2*BBE_SIMD_WIDTH);            \
        UPDATE_X2(M,1,row_stride);                             \
        UPDATE_X2(M,2,row_stride);                             \
        UPDATE_X2(M,3,row_stride);                             \
        UPDATE_X2(M,4,row_stride);                             \
        UPDATE_X2(M,5,row_stride);                             \
        UPDATE_X2(M,6,row_stride);                             \
        UPDATE_X2(M,7,row_stride);                             \
        UPDATE_X2(M,0,-row_stride*(M-1) + 2*BBE_SIMD_WIDTH);   \
    }                                                          \
    if (P&1)                                                   \
    {                                                          \
        for(i=0; i < L/(BBE_SIMD_WIDTH/2); i++)                \
        {                                                      \
            BBE_LVNX16_IP(vR1, pVR1, 2*BBE_SIMD_WIDTH);        \
            UPDATE_X1(M,1,row_stride);                         \
            UPDATE_X1(M,2,row_stride);                         \
            UPDATE_X1(M,3,row_stride);                         \
            UPDATE_X1(M,4,row_stride);                         \
            UPDATE_X1(M,5,row_stride);                         \
            UPDATE_X1(M,6,row_stride);                         \
            UPDATE_X1(M,0,-row_stride*(M-1)+2*BBE_SIMD_WIDTH); \
        }                                                      \
    }                                                          \
}

/*-----------------------------------------------------------------------
[c]qr_build_rMxNs

QR decomposition of MxN complex matrices.
Instead of direct computation of Q factors, these functions produce a
set of N Householder vectors V for each of input matrices A. This
approach allow us to save CPU cycles and memory when solving a system
of linear equations: it is cheaper to perform N elementary reflections
for a right hand side vector if compared to explicit multiplication of
that vector by matrix Q.

Fixed point representation of output matrices R is the same as for
input matrices, but Householder vectors V are always Q14.

Data transform is performed in-place: upper triangular matrices R replace
input matrices A.

NOTE:
Data layout for matrices is selected as for other matrices written in a 
streaming order. 

Input:
R[M*N][L]                  matrices A (L matrices of size MxN)
Output:
V[((M*N+((N-1)*N)/2+M)*L]  L sets of Householder vectors
R[M*N][L]                  upper triangular matrices (L matrices
                           of size MxN)

Restrictions:
1. All matrices must not overlap and must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,N,L)
4. M must greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
-------------------------------------------------------------------------*/
void cqr_build_r8x8s  ( void* pScr, complex_fract16 * restrict _V, complex_fract16 * restrict _R, int L)
{
    int16_t * restrict V=(int16_t *)_V;
    int16_t * restrict R=(int16_t *)_R;
    xb_vecNx16 * restrict pR0 ;
    xb_vecNx16 * restrict pX0 ;
    xb_vecNx16 * restrict pY0 ;
    xb_vecNx16 * restrict pV0 ;
    xb_vecNx16 * restrict pY1 ;
    xb_vecNx16 * restrict pX1 ;
    xb_vecNx16 * restrict pVR ;
    xb_vecNx16 * restrict pVR1;

    int l;
    int m;
    int M = CURRENT_M;

    int16_t *pR = R;

    xb_vecNx16 * restrict pV;
    xb_vecNx16 * restrict pFi;

    int16_t *vR = (int16_t*) pScr;
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);;
    NASSERT_ALIGN(V,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);;
    NASSERT(L%(BBE_SIMD_WIDTH/2) == 0 && L > 0);

    pV =  (xb_vecNx16 *)(V);
    pFi = (xb_vecNx16 *)(V + 2*L*( SIZE_OF_V(CURRENT_M,CURRENT_M)));

    chouseholder_M8xL_m0(R, (int16_t*)pV, (int16_t*)pFi, L);
    CALCVR8x8(vR, pR, ((int16_t*)pV), L, 8, 8, (CURRENT_M*4*L),Q,"loop_count min=8");
    __Pragma("no_reorder");
    UPDATE8x8(pR, pR, vR, ((int16_t*)pV), L, 8, 8, (CURRENT_M*4*L),"loop_count min=4");
    __Pragma("no_reorder");
    pV  += 2*L*M*sizeof(int16_t)/sizeof(*pV);
    pFi += 2*L*sizeof(int16_t)/sizeof(*pFi);
    pR += 2*L*(CURRENT_M+1) ;
    M--;

    chouseholder_M8xL_m1(R, (int16_t*)pV, (int16_t*)pFi, L);
    CALCVR8x8(vR, pR, ((int16_t*)pV), L, 7, 7, (CURRENT_M*4*L),Q,"loop_count min=7");
    __Pragma("no_reorder");
    UPDATE8x8(pR, pR, vR, ((int16_t*)pV), L, 7, 7, (CURRENT_M*4*L),"loop_count min=3");
    __Pragma("no_reorder");
    pV  += 2*L*M*sizeof(int16_t)/sizeof(*pV);
    pFi += 2*L*sizeof(int16_t)/sizeof(*pFi);
    pR += 2*L*(CURRENT_M+1) ;
    M--;

    chouseholder_M8xL_m2(R, (int16_t*)pV, (int16_t*)pFi, L);
    CALCVR8x8(vR, pR, ((int16_t*)pV), L, 6, 6, (CURRENT_M*4*L),Q,"loop_count min=6");
    __Pragma("no_reorder");
    UPDATE8x8(pR, pR, vR, ((int16_t*)pV), L, 6, 6, (CURRENT_M*4*L),"loop_count min=3");
    __Pragma("no_reorder");
    pV  += 2*L*M*sizeof(int16_t)/sizeof(*pV);
    pFi += 2*L*sizeof(int16_t)/sizeof(*pFi);
    pR += 2*L*(CURRENT_M+1) ;
    M--;

    chouseholder_M8xL_m3(R,(int16_t*)pV, (int16_t*)pFi, L);
    CALCVR8x8(vR, pR, ((int16_t*)pV), L, 5, 5, (CURRENT_M*4*L),Q,"loop_count min=5");
    __Pragma("no_reorder");
    UPDATE8x8(pR, pR, vR, ((int16_t*)pV), L, 5, 5, (CURRENT_M*4*L),"loop_count min=2");
    __Pragma("no_reorder");
    pV  += 2*L*M*sizeof(int16_t)/sizeof(*pV);
    pFi += 2*L*sizeof(int16_t)/sizeof(*pFi);
    pR += 2*L*(CURRENT_M+1) ;
    M--;

    chouseholder_M8xL_m4(R, (int16_t*)pV, (int16_t*)pFi, L);
    CALCVR8x8(vR, pR, ((int16_t*)pV), L, 4, 4, (CURRENT_M*4*L),Q,"loop_count min=4");
    __Pragma("no_reorder");
    UPDATE8x8(pR, pR, vR, ((int16_t*)pV), L, 4, 4, (CURRENT_M*4*L),"loop_count min=2");
    __Pragma("no_reorder");
    pV  += 2*L*M*sizeof(int16_t)/sizeof(*pV);
    pFi += 2*L*sizeof(int16_t)/sizeof(*pFi);
    pR += 2*L*(CURRENT_M+1) ;
    M--;


    chouseholder_M8xL_m5(R, (int16_t*)pV, (int16_t*)pFi, L);
    CALCVR8x8(vR, pR, ((int16_t*)pV), L, 3, 3, (CURRENT_M*4*L),Q,"loop_count min=3");
    __Pragma("no_reorder");
    UPDATE8x8(pR, pR, vR, ((int16_t*)pV), L, 3, 3, (CURRENT_M*4*L),"loop_count min=1");
    __Pragma("no_reorder");
    pV  += 2*L*M*sizeof(int16_t)/sizeof(*pV);
    pFi += 2*L*sizeof(int16_t)/sizeof(*pFi);
    pR += 2*L*(CURRENT_M+1) ;
    M--;

    chouseholder_M8xL_m6(R, (int16_t*)pV, (int16_t*)pFi, L);
    CALCVR8x8(vR, pR, ((int16_t*)pV), L, 2, 2, (CURRENT_M*4*L),Q,"loop_count min=2");
    __Pragma("no_reorder");
    UPDATE8x8(pR, pR, vR, ((int16_t*)pV), L, 2, 2, (CURRENT_M*4*L),"loop_count min=1");
    __Pragma("no_reorder");
    pV  += 2*L*M*sizeof(int16_t)/sizeof(*pV);
    pFi += 2*L*sizeof(int16_t)/sizeof(*pFi);
    pR += 2*L*(CURRENT_M+1) ;
    M--;

    __Pragma("loop_count min=1")
    for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2) )
    {
        xb_vecNx16 tmp0,tmp1;
        vsaN c_vec,v_exp;
        xb_vecNx40 acc0;
        xb_vecNx16  mant;

        xb_vecNx16  x0, _fi;
        const xb_vecNx16 *  pa;
        m = CURRENT_M-1;
        pa = (xb_vecNx16 *)(R + m*L*2 + m*CURRENT_M*L*2+2*l);
        x0 = BBE_LVNX16_I( pa, 0 );
        acc0 = BBE_MAGINX16C( x0, x0);
        acc0=BBE_ADDNX40(acc0,acc0);
        c_vec=BBE_NSAENX40(acc0);
        acc0=BBE_SLLNX40(acc0,c_vec);
        BBE_RSQRTLUNX40_0(acc0,tmp0, tmp1, acc0);
        BBE_MULUUSNX16( acc0, tmp1,  tmp0);
        acc0=BBE_SRAINX40(acc0,24);
        tmp0=BBE_PACKLNX40(acc0);
        mant= BBE_SHFLNX16I(tmp0, BBE_SHFLI_DUPLICATE_1_EVEN);
        v_exp = BBE_SUBSR1SAVSN(18+1,c_vec);
        acc0 = BBE_MULUSRNX16(mant,x0, v_exp);
        _fi = BBE_PACKVNX40(acc0, v_exp);
        BBE_SVNX16_IP(_fi, pFi/*(Fi+8*2*m)*/, sizeof(*pFi) );
    }

    pFi -= CURRENT_M*L/(BBE_SIMD_WIDTH/2);

    RotateRows8x8(R, (const int16_t *)pFi, L, CURRENT_M, CURRENT_M);
} 

size_t cqr_build_r8x8s_getScratchSize (int M, int N, int L)
{
    (void)N;
    return M*2*L*sizeof(int16_t);
} /* cqr_build_r8x8s_getScratchSize() */

#else
DISCARD_FUN(void, cqr_build_r8x8s, ( void* pScr,complex_fract16 * restrict V, complex_fract16 * restrict R, int L))
size_t cqr_build_r8x8s_getScratchSize (int M, int N, int L) { (void)M;(void)N;(void)L; return 0; }
#endif
