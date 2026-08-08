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
    cqr_calc_qbMxNxPs/qr_calc_qbMxNxPs
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

#define CURRENT_M 16
#define CURRENT_P 1
#define UPDATE_COLUMN16x16(_v_step, _b_step12, _b_step13, _b_step14, _m,reload )    \
{                                                                       \
    xb_vecNx16 _0x4000=BBE_MOVPINT16(16);                               \
    xb_vecNx16 v0, vR0, ar[16];                                         \
    xb_vecNx40 acc0;                                                    \
    vsaN _14=BBE_MOVVSA32(14);                                          \
                v0= BBE_LVNX16_I(_pv,CLIP_TO_ZERO(15-_m)*2*BBE_SIMD_WIDTH);BBE_LVNX16_IP(ar[15],_pb15,     0);  acc0=BBE_MULRNX16J(  ar[15], v0,_14); \
                v0= BBE_LVNX16_I(_pv,CLIP_TO_ZERO(14-_m)*2*BBE_SIMD_WIDTH);BBE_LVNX16_IP(ar[14],_pb14,     0);  BBE_MULANX16J(acc0,  ar[14], v0); \
    if(_m<=13) {v0= BBE_LVNX16_I(_pv,CLIP_TO_ZERO(13-_m)*2*BBE_SIMD_WIDTH);BBE_LVNX16_IP(ar[13],_pb13,     0);  BBE_MULANX16J(acc0,  ar[13], v0);}\
    if(_m<=12) {v0= BBE_LVNX16_I(_pv,CLIP_TO_ZERO(12-_m)*2*BBE_SIMD_WIDTH);BBE_LVNX16_IP(ar[12],_pb12,     0);  BBE_MULANX16J(acc0,  ar[12], v0);}\
    if(_m<=0) {v0= BBE_LVNX16_I(_pv, CLIP_TO_ZERO(0-_m)*2*BBE_SIMD_WIDTH);    ar[ 0]=BBE_LVNX16_X( _pb12, _o_12); BBE_MULANX16J(acc0,  ar[ 0], v0);}\
    if(_m<=1) {v0= BBE_LVNX16_I(_pv, CLIP_TO_ZERO(1-_m)*2*BBE_SIMD_WIDTH);    ar[ 1]=BBE_LVNX16_X( _pb13, _o_12); BBE_MULANX16J(acc0,  ar[ 1], v0);}\
    if(_m<=2) {v0= BBE_LVNX16_I(_pv, CLIP_TO_ZERO(2-_m)*2*BBE_SIMD_WIDTH);    ar[ 2]=BBE_LVNX16_X( _pb14, _o_12); BBE_MULANX16J(acc0,  ar[ 2], v0);}\
    if(_m<=3) {v0= BBE_LVNX16_I(_pv, CLIP_TO_ZERO(3-_m)*2*BBE_SIMD_WIDTH);    ar[ 3]=BBE_LVNX16_X( _pb15, _o_12); BBE_MULANX16J(acc0,  ar[ 3], v0);}\
    if(_m<=4) {v0= BBE_LVNX16_I(_pv, CLIP_TO_ZERO(4-_m)*2*BBE_SIMD_WIDTH);    ar[ 4]=BBE_LVNX16_X( _pb12, _o_8);  BBE_MULANX16J(acc0,  ar[ 4], v0);}\
    if(_m<=5) {v0= BBE_LVNX16_I(_pv, CLIP_TO_ZERO(5-_m)*2*BBE_SIMD_WIDTH);    ar[ 5]=BBE_LVNX16_X( _pb13, _o_8);  BBE_MULANX16J(acc0,  ar[ 5], v0);}\
    if(_m<=6) {v0= BBE_LVNX16_I(_pv, CLIP_TO_ZERO(6-_m)*2*BBE_SIMD_WIDTH);    ar[ 6]=BBE_LVNX16_X( _pb14, _o_8);  BBE_MULANX16J(acc0,  ar[ 6], v0);}\
    if(_m<=7) {v0= BBE_LVNX16_I(_pv, CLIP_TO_ZERO(7-_m)*2*BBE_SIMD_WIDTH);    ar[ 7]=BBE_LVNX16_X( _pb15, _o_8);  BBE_MULANX16J(acc0,  ar[ 7], v0);}\
    if(_m<=8) {v0= BBE_LVNX16_I(_pv, CLIP_TO_ZERO(8-_m)*2*BBE_SIMD_WIDTH);    ar[ 8]=BBE_LVNX16_X( _pb12, _o_4);  BBE_MULANX16J(acc0,  ar[ 8], v0);}\
    if(_m<=9) {v0= BBE_LVNX16_I(_pv, CLIP_TO_ZERO(9-_m)*2*BBE_SIMD_WIDTH);    ar[ 9]=BBE_LVNX16_X( _pb13, _o_4);  BBE_MULANX16J(acc0,  ar[ 9], v0);}\
    if(_m<=10) {v0= BBE_LVNX16_I(_pv,CLIP_TO_ZERO(10-_m)*2*BBE_SIMD_WIDTH);ar[10]=BBE_LVNX16_X( _pb14,  _o_4);  BBE_MULANX16J(acc0,  ar[10], v0);}\
    if(_m<=11) {v0= BBE_LVNX16_I(_pv,CLIP_TO_ZERO(11-_m)*2*BBE_SIMD_WIDTH);ar[11]=BBE_LVNX16_X( _pb15,  _o_4);  BBE_MULANX16J(acc0,  ar[11], v0);}\
    vR0 = BBE_PACKVNX40(acc0, _14);                                                                                                               \
    if(_m<= 0) { BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 0],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[ 0]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_X (ar[ 0], _pb12, _o_12);}\
    if(_m<= 1) { BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 1],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[ 1]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_X (ar[ 1], _pb13, _o_12);}\
    if(_m<= 2) { BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 2],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[ 2]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_X (ar[ 2], _pb14, _o_12);}\
    if(_m<= 3) { BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 3],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[ 3]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_X (ar[ 3], _pb15, _o_12);}\
    if(_m<= 4) { BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 4],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[ 4]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_X (ar[ 4], _pb12, _o_8); }\
    if(_m<= 5) { BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 5],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[ 5]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_X (ar[ 5], _pb13, _o_8); }\
    if(_m<= 6) { BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 6],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[ 6]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_X (ar[ 6], _pb14, _o_8); }\
    if(_m<= 7) { BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 7],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[ 7]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_X (ar[ 7], _pb15, _o_8); }\
    if(_m<= 8) { BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 8],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[ 8]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_X (ar[ 8], _pb12, _o_4); }\
    if(_m<= 9) { BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[ 9],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[ 9]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_X (ar[ 9], _pb13, _o_4); }\
    if(_m<=10) { BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[10],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[10]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_X (ar[10], _pb14, _o_4); }\
    if(_m<=11) { BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[11],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[11]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_X (ar[11], _pb15, _o_4); }\
    if(_m<=12) { if(reload)ar[12]=BBE_LVNX16_I( _pb12,     0);BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[12],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[12]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[12], _pb12, _b_step12);}\
    if(_m<=13) { if(reload)ar[13]=BBE_LVNX16_I( _pb13,     0);BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[13],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[13]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[13], _pb13, _b_step13);}\
    if(reload)ar[14]=BBE_LVNX16_I( _pb14,     0);\
    if(reload)ar[15]=BBE_LVNX16_I( _pb15,     0);\
                 BBE_LVNX16_IP(v0,_pv,2*BBE_SIMD_WIDTH); acc0=BBE_MULRNX16(ar[14],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[14]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[14], _pb14, _b_step14); \
                 BBE_LVNX16_XP(v0,_pv,         _v_step); acc0=BBE_MULRNX16(ar[15],_0x4000,_14); BBE_MULSNX16C(acc0,vR0,v0); ar[15]=BBE_PACKVNX40(acc0,_14); BBE_SVNX16_XP(ar[15], _pb15, _b_step14); \
}

/*-------------------------------------------------------------------------
cqr_calc_qbMxNxPs/qr_calc_qbMxNxPs

These functions apply Householder reflections to L MxP matrices B in the
course of solving a set of complex-valued linear problems A*X=B through
the QR decomposition of matrices A: A*X=B, A=Q*R => Q*R*X=B => R*X=Q'*B.
Instead of direct multiplication of each matrix B by conjugate transpose
of the corresponding matrix Q, we use a set of vectors V to perform a
sequence of Householder  reflections (see QR decomposition routines
cqr_build_rMxN/qr_build_rMxN).

Fixed point representation of output matrices is the same as for input
matrices.

Data transform is performed in-place.

NOTE:
Data layout for matrices is selected as for other matrices written in a stream 
order. 

Input:
B[M*P]L]                   Matrices B (L matrices of size MxP)
V[((M*N+((N-1)*N)/2+M)*L]  L sets of Householder vectors
Output:
B[M*P][L]                  Matrices Q'*B (L matrices of size MxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(M,P,L)
4. M must be greater than or equal to N
5. Matrix sizes M,N,L must be greater than 1
---------------------------------------------------------------------------*/
void cqr_calc_qb16x16x1s (void *pScr,complex_fract16 *_B, const complex_fract16 *_V, int L)
{
          int16_t *B=(      int16_t *)_B;
    const int16_t *V=(const int16_t *)_V;
    int16_t *start_B = B; 
    const int16_t *start_V = V; 
    int l;

    xb_vecNx16 *  restrict   _pv = (xb_vecNx16 *)V; 
    xb_vecNx16 *  restrict  _pb12 = (xb_vecNx16 *)(B+(CURRENT_M-4)*2*L); 
    xb_vecNx16 *  restrict  _pb13 = (xb_vecNx16 *)(B+(CURRENT_M-3)*2*L); 
    xb_vecNx16 *  restrict  _pb14 = (xb_vecNx16 *)(B+(CURRENT_M-2)*2*L); 
    xb_vecNx16 *  restrict  _pb15 = (xb_vecNx16 *)(B+(CURRENT_M-1)*2*L); 

    const int     _o_4  =  -4*2*L*sizeof(int16_t);    // offset 4 rows
    const int     _o_8  =  -8*2*L*sizeof(int16_t);    // offset 8 rows
    const int     _o_12 =  -12*2*L*sizeof(int16_t);    // offset 12 rows
    const int     _v_step= sizeof(*_pv);
    const int     _v_last_step = BBE_SIMD_WIDTH*(1+SIZE_OF_FI(CURRENT_M))*sizeof(int16_t); 
    xb_vecNx16 *  restrict  pv ;
    xb_vecNx16 *  restrict  pb0; 
    xb_vecNx16 *  restrict  pb1; 
    xb_vecNx16 *  restrict  pb2; 
    xb_vecNx16 *  restrict  pb3; 

    xb_vecNx16 *  restrict  pb4; 
    xb_vecNx16 *  restrict  pb5; 
    xb_vecNx16 *  restrict  pb6; 
    xb_vecNx16 *  restrict  pb7; 

#ifdef COMPILER_XTENSA

#pragma ymemory( _pb12)
#pragma ymemory( _pb13)
#pragma ymemory( _pb14)
#pragma ymemory( _pb15)
#pragma ymemory(pb0)
#pragma ymemory(pb1)
#pragma ymemory(pb2)
#pragma ymemory(pb3)
#pragma ymemory(pb4)
#pragma ymemory(pb5)
#pragma ymemory(pb6)
#pragma ymemory(pb7)

#endif

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V,2*BBE_SIMD_WIDTH);
    NASSERT(L%(BBE_SIMD_WIDTH/2)==0);
    if (L<=0) return;
    __Pragma("loop_count min=1");
    for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2))
    {

      UPDATE_COLUMN16x16(_v_step,        0,    0,    0,    0, 1 );
      UPDATE_COLUMN16x16(_v_step,        0,    0,    0,    1, 1 );
      UPDATE_COLUMN16x16(_v_step,        0,    0,    0,    2, 1 );
      UPDATE_COLUMN16x16(_v_step,        0,    0,    0,    3, 1 );

      UPDATE_COLUMN16x16(_v_step,        0,    0,    0,    4, 1 );
      UPDATE_COLUMN16x16(_v_step,        0,    0,    0,    5, 1 );
      UPDATE_COLUMN16x16(_v_step,        0,    0,    0,    6, 1 );
      UPDATE_COLUMN16x16(_v_step,        0,    0,    0,    7, 1 );

      UPDATE_COLUMN16x16(_v_step,        0,    0,    0,    8, 0 );
      UPDATE_COLUMN16x16(_v_step,        0,    0,    0,    9, 0 );
      UPDATE_COLUMN16x16(_v_step,        0,    0,    0,    10, 0);
      UPDATE_COLUMN16x16(_v_step,        0,    0,    0,    11, 0);

      UPDATE_COLUMN16x16(_v_step,        2*BBE_SIMD_WIDTH,    0,    0,    12, 0);
      UPDATE_COLUMN16x16(_v_step,        0 , 2*BBE_SIMD_WIDTH,    0,    13, 0);
      UPDATE_COLUMN16x16(_v_last_step,    0 , 0,    2*BBE_SIMD_WIDTH,    14, 0);

      B += BBE_SIMD_WIDTH;
      V += BBE_SIMD_WIDTH*(SIZE_OF_V(CURRENT_M, CURRENT_M)+SIZE_OF_FI(CURRENT_M)); 
    }

    {
      vsaN _14 = BBE_MOVVSA32(14);
      pv = (xb_vecNx16 *)(start_V + BBE_SIMD_WIDTH*(SIZE_OF_V(CURRENT_M, CURRENT_M))+2*0*(SIZE_OF_V(CURRENT_M, CURRENT_M)+SIZE_OF_FI(CURRENT_M)));
      pb0 = (xb_vecNx16 *)(start_B + 2*L*0); 
      pb1 = (xb_vecNx16 *)(start_B + 2*L*1); 
      pb2 = (xb_vecNx16 *)(start_B + 2*L*2); 
      pb3 = (xb_vecNx16 *)(start_B + 2*L*3); 
      pb4 = (xb_vecNx16 *)(start_B + 2*L*4); 
      pb5 = (xb_vecNx16 *)(start_B + 2*L*5); 
      pb6 = (xb_vecNx16 *)(start_B + 2*L*6); 
      pb7 = (xb_vecNx16 *)(start_B + 2*L*7); 
      __Pragma("loop_count min=1");
      for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2))
      {
       
        xb_vecNx40 acc0;
        xb_vecNx16 b0,b1,b2,b3,b4,b5,b6,b7; 
        xb_vecNx16 v0,v1,v2,v3,v4,v5,v6,v7; 

        BBE_LVNX16_IP(b0, pb0, 0); BBE_LVNX16_IP(v0, pv, 2*BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b0, v0, _14); b0 = BBE_PACKVNX40(acc0, _14);    BBE_SVNX16_IP(b0, pb0, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(b1, pb1, 0); BBE_LVNX16_IP(v1, pv, 2 * BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b1, v1,_14); b1=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b1,pb1, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(b2, pb2, 0); BBE_LVNX16_IP(v2, pv, 2 * BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b2, v2,_14); b2=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b2,pb2, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(b3, pb3, 0); BBE_LVNX16_IP(v3, pv, 2 * BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b3, v3,_14); b3=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b3,pb3, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(b4, pb4, 0); BBE_LVNX16_IP(v4, pv, 2 * BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b4, v4,_14); b4=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b4,pb4, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(b5, pb5, 0); BBE_LVNX16_IP(v5, pv, 2 * BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b5, v5,_14); b5=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b5,pb5, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(b6, pb6, 0); BBE_LVNX16_IP(v6, pv, 2 * BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b6, v6,_14); b6=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b6,pb6, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(b7, pb7, 0);BBE_LVNX16_XP(v7, pv, 2*BBE_SIMD_WIDTH*(1+(8+SIZE_OF_V(CURRENT_M, CURRENT_M))));
        acc0 = BBE_MULRNX16J(b7, v7,_14); b7=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b7,pb7, 2*BBE_SIMD_WIDTH);
      }

      pv = (xb_vecNx16 *)(start_V + BBE_SIMD_WIDTH*(8+(SIZE_OF_V(CURRENT_M, CURRENT_M)))+2*0*(SIZE_OF_V(CURRENT_M, CURRENT_M)+SIZE_OF_FI(CURRENT_M))); 
      pb0 = (xb_vecNx16 *)(start_B + 2*L*8); 
      pb1 = (xb_vecNx16 *)(start_B + 2*L*9); 
      pb2 = (xb_vecNx16 *)(start_B + 2*L*10); 
      pb3 = (xb_vecNx16 *)(start_B + 2*L*11); 
      pb4 = (xb_vecNx16 *)(start_B + 2*L*12); 
      pb5 = (xb_vecNx16 *)(start_B + 2*L*13); 
      pb6 = (xb_vecNx16 *)(start_B + 2*L*14); 
      pb7 = (xb_vecNx16 *)(start_B + 2*L*15); 
      __Pragma("loop_count min=1");
      for(l=0; l<L; l+=(BBE_SIMD_WIDTH/2))
      {
        vsaN _14=BBE_MOVVSA32(14);
        xb_vecNx40 acc0;
        xb_vecNx16 b0,b1,b2,b3,b4,b5,b6,b7; 
        xb_vecNx16 v0,v1,v2,v3,v4,v5,v6,v7; 

        BBE_LVNX16_IP(b0, pb0, 0);
        BBE_LVNX16_IP(b1, pb1, 0);
        BBE_LVNX16_IP(b2, pb2, 0);
        BBE_LVNX16_IP(b3, pb3, 0);
        BBE_LVNX16_IP(b4, pb4, 0);
        BBE_LVNX16_IP(b5, pb5, 0);
        BBE_LVNX16_IP(b6, pb6, 0);
        BBE_LVNX16_IP(b7, pb7, 0);

        BBE_LVNX16_IP(v0, pv, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v1, pv, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v2, pv, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v3, pv, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v4, pv, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v5, pv, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(v6, pv, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(v7, pv, 2*BBE_SIMD_WIDTH*(1+(8+SIZE_OF_V(CURRENT_M, CURRENT_M))));

        acc0 = BBE_MULRNX16J(b0, v0,_14); b0=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b0,pb0, 2*BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b1, v1,_14); b1=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b1,pb1, 2*BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b2, v2,_14); b2=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b2,pb2, 2*BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b3, v3,_14); b3=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b3,pb3, 2*BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b4, v4,_14); b4=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b4,pb4, 2*BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b5, v5,_14); b5=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b5,pb5, 2*BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b6, v6,_14); b6=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b6,pb6, 2*BBE_SIMD_WIDTH);
        acc0 = BBE_MULRNX16J(b7, v7,_14); b7=BBE_PACKVNX40(acc0,_14);    BBE_SVNX16_IP(b7,pb7, 2*BBE_SIMD_WIDTH);
      }
    }

} /* cqr_calc_qb16x16x1s() */

size_t cqr_calc_qb16x16x1s_getScratchSize (int M, int P, int L)
{
    (void)M;(void)P;(void)L;
    return 0;
} /* cqr_calc_qb16x16x1s_getScratchSize() */
