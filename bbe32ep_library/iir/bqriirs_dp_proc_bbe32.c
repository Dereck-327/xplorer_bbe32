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
    Biquad real block IIR, streaming version
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_iir.h"
/* Filter processing functions. */
#include "bqriirs_common.h"

/* Biquad real block IIR processing function, Low Noise fixed-point implementation. */
void bqriirs_dp_proc( int16_t * restrict r,
                      int16_t * restrict sect, // 6*M
                const int16_t *          x,
                      int16_t * restrict scr,  // 4*N
                const int16_t *          coef,
                      int16_t            gain,
                      int N, int L, int M )
{
  unsigned int a12, b12;

  xb_vecNx16 G, B0, B1, B2, A1, A2;

        xb_vecNx16 *SECT, *SCR;
        xb_vecNx16 * restrict R;
  const xb_vecNx16 *X_;

  xb_vecNx16 D0_so, D0_bo, D0_s, D0_b, D1_s, D1_b, D2_s, D2_b;

  xb_vecNx16 X32_s, X32_b, X, Y, Gain, Y_s, Y_b;

  xb_vecNx40 Acc_fb, Acc_fb_s, Acc_fb_b;
  xb_vecNx40 Acc_ff, Acc_ff_s, Acc_ff_b;

  xb_vecNx40 z;

  vsaN rnd8, rnd12, rnd14, rnd16;

  int n, m, l;

  SECT = (xb_vecNx16 *)sect;
  SCR  = (xb_vecNx16 *)scr;
  R    = (xb_vecNx16 *)r;

  rnd8  = BBE_MOVVSA32(8);
  rnd12 = BBE_MOVVSA32(12);
  rnd14 = BBE_MOVVSA32(14);
  rnd16 = BBE_MOVVSA32(16);

  NASSERT(L%BBE_SIMD_WIDTH==0);

  Gain = BBE_LSNX16_I(&gain, 0);
  Gain = BBE_REPNX16(Gain, 0);

  if ( M == 1 )
  {
    BBE_LSNX16_IP(G,  coef, 2); 
    BBE_LSNX16_IP(B0, coef, 2); b12 = BBE_L32X((const int32_t *)coef, 0);
    BBE_LSNX16_IP(B1, coef, 2); 
    BBE_LSNX16_IP(B2, coef, 2); a12 = BBE_L32X((const int32_t *)coef, 0);
    BBE_LSNX16_IP(A1, coef, 2); 
    BBE_LSNX16_IP(A2, coef, 2);

    G  = BBE_REPNX16(G,  0); B0 = BBE_REPNX16(B0, 0); B1 = BBE_REPNX16(B1, 0);
    B2 = BBE_REPNX16(B2, 0); A1 = BBE_REPNX16(A1, 0); A2 = BBE_REPNX16(A2, 0);

    for ( l=0; l<(L>>LOG2_BBE_SIMD_WIDTH); l++ )
    {
      D0_so = BBE_LVNX16_I(SECT, 0*BBE_SIMD_WIDTH); D0_bo = BBE_LVNX16_I(SECT,  2*BBE_SIMD_WIDTH);
      D0_s  = BBE_LVNX16_I(SECT, 4*BBE_SIMD_WIDTH); D0_b  = BBE_LVNX16_I(SECT,  6*BBE_SIMD_WIDTH);
      D1_s  = BBE_LVNX16_I(SECT, 8*BBE_SIMD_WIDTH); D1_b  = BBE_LVNX16_I(SECT, 10*BBE_SIMD_WIDTH);
      
      X_ = (const xb_vecNx16 *)(x + BBE_SIMD_WIDTH*l);

      for ( n=0; n<N; n++ )
      {
        BBE_LVNX16_XP(X, X_, L*2);

        // Q29 <- Q15 + 14
        X32_s = BBE_SLLINX16(X, 14);
        X32_b = BBE_SRAINX16(X,  2);

        // Q43 <- Q29*Q14
        Acc_fb_s = BBE_MULUSNX16(D0_s, A2);
        Acc_fb_b = BBE_MULNX16(X32_b, G);
        Acc_ff_s = BBE_MULUSNX16(D0_so, B2);
        Acc_ff_b = BBE_MULRNX16(D1_b, B0, rnd12);

        // Q43 <- Q43 - Q29*Q14
        BBE_MULUSANX16(Acc_fb_s, D1_s, A1);
        BBE_MULUSANX16(Acc_fb_s, X32_s, G);
        BBE_MULUSANX16(Acc_ff_s, D0_s, B1);
        BBE_MULUSANX16(Acc_ff_s, D1_s, B0);
        BBE_MULANX16PR(Acc_fb_b, D0_b, D1_b, a12);
        BBE_MULANX16PR(Acc_ff_b, D0_bo, D0_b, b12);

        Acc_fb_b = BBE_SLLINX40(Acc_fb_b,  2);
        Acc_fb_s = BBE_SRAINX40(Acc_fb_s, 14);
        Acc_ff_b = BBE_SLLINX40(Acc_ff_b,  2);
        Acc_ff_s = BBE_SRAINX40(Acc_ff_s, 14);

        Acc_fb = BBE_ADDNX40(Acc_fb_s, Acc_fb_b);
        Acc_ff = BBE_ADDNX40(Acc_ff_s, Acc_ff_b);

        D2_s = BBE_PACKLNX40(Acc_fb); D2_b = BBE_PACKVNX40(Acc_fb, rnd16);

        D0_so = D0_s; D0_bo = D0_b;
        D0_s  = D1_s; D0_b  = D1_b;
        D1_s  = D2_s; D1_b  = D2_b;

        // Q15 <- ( Q43 - 14 ) - 14 w/ rounding
        Y = BBE_PACKVNX40(Acc_ff, rnd14);

        // Q15 <- Q15*Q8 - 8
        z = BBE_MULNX16(Y, Gain);
        Y = BBE_PACKVNX40(z, rnd8);

        BBE_SVNX16_X(Y, R, 2*n*L + 2*BBE_SIMD_WIDTH*l);
      }

      BBE_SVNX16_IP(D0_so, SECT, 2*BBE_SIMD_WIDTH); BBE_SVNX16_IP(D0_bo, SECT, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(D0_s,  SECT, 2*BBE_SIMD_WIDTH); BBE_SVNX16_IP(D0_b,  SECT, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(D1_s,  SECT, 2*BBE_SIMD_WIDTH); BBE_SVNX16_IP(D1_b,  SECT, 2*BBE_SIMD_WIDTH);
    }

    return;
  }

  BBE_LSNX16_IP(G,  coef, 2); 
  BBE_LSNX16_IP(B0, coef, 2); b12 = BBE_L32X((const int32_t *)coef, 0);
  BBE_LSNX16_IP(B1, coef, 2); 
  BBE_LSNX16_IP(B2, coef, 2); a12 = BBE_L32X((const int32_t *)coef, 0);
  BBE_LSNX16_IP(A1, coef, 2); 
  BBE_LSNX16_IP(A2, coef, 2);

  G  = BBE_REPNX16(G,  0); B0 = BBE_REPNX16(B0, 0); B1 = BBE_REPNX16(B1, 0);
  B2 = BBE_REPNX16(B2, 0); A1 = BBE_REPNX16(A1, 0); A2 = BBE_REPNX16(A2, 0);

  for ( l=0; l<L/BBE_SIMD_WIDTH; l++ )
  {
    D0_so = BBE_LVNX16_I(SECT, 0*BBE_SIMD_WIDTH); D0_bo = BBE_LVNX16_I(SECT,  2*BBE_SIMD_WIDTH);
    D0_s  = BBE_LVNX16_I(SECT, 4*BBE_SIMD_WIDTH); D0_b  = BBE_LVNX16_I(SECT,  6*BBE_SIMD_WIDTH);
    D1_s  = BBE_LVNX16_I(SECT, 8*BBE_SIMD_WIDTH); D1_b  = BBE_LVNX16_I(SECT, 10*BBE_SIMD_WIDTH);
        
    X_ = (const xb_vecNx16 *)(x + BBE_SIMD_WIDTH*l);

    for ( n=0; n<N; n++ )
    {
      BBE_LVNX16_XP(X, X_, L*2);

      // Q29 <- Q15 + 14
      X32_s = BBE_SLLINX16(X, 14);
      X32_b = BBE_SRAINX16(X,  2);

      // Q43 <- Q29*Q14
      // = BBE_MULAUSNX16(D0_s, A2);
      Acc_fb_s = BBE_MULUSNX16(D1_s, A1);
      BBE_MULUSANX16(Acc_fb_s,D0_s,A2);
      BBE_MULUSANX16(Acc_fb_s, X32_s, G);
      
      Acc_fb_b = BBE_MULNX16(X32_b, G);
      BBE_MULANX16PR(Acc_fb_b, D0_b, D1_b, a12);
      
      Acc_ff_s = BBE_MULUSNX16(D0_so, B2);
      BBE_MULUSANX16(Acc_ff_s, D0_s, B1);
      BBE_MULUSANX16(Acc_ff_s, D1_s, B0);
                
      Acc_ff_b = BBE_MULNX16(D1_b, B0);
      BBE_MULANX16PR(Acc_ff_b, D0_bo, D0_b, b12);

      // Q43 <- Q43 - Q29*Q14
      
      Acc_fb_b = BBE_SLLINX40(Acc_fb_b,  2);
      Acc_fb_s = BBE_SRAINX40(Acc_fb_s, 14);
      Acc_ff_b = BBE_SLLINX40(Acc_ff_b,  2);
      Acc_ff_s = BBE_SRAINX40(Acc_ff_s, 14);

      Acc_fb = BBE_ADDNX40(Acc_fb_s, Acc_fb_b);
      Acc_ff = BBE_ADDNX40(Acc_ff_s, Acc_ff_b);

      D2_s = BBE_PACKLNX40(Acc_fb); D2_b = BBE_PACKVNX40(Acc_fb, rnd16);
      Y_s  = BBE_PACKLNX40(Acc_ff); Y_b  = BBE_PACKVNX40(Acc_ff, rnd16);

      D0_so = D0_s; D0_bo = D0_b;
      D0_s  = D1_s; D0_b  = D1_b;
      D1_s  = D2_s; D1_b  = D2_b;

      BBE_SVNX16_X(Y_s, SCR, 2*n*L*2 + 2*BBE_SIMD_WIDTH*l*2 + 0*BBE_SIMD_WIDTH);
      BBE_SVNX16_X(Y_b, SCR, 2*n*L*2 + 2*BBE_SIMD_WIDTH*l*2 + 2*BBE_SIMD_WIDTH);
    }

    BBE_SVNX16_IP(D0_so, SECT, 2*BBE_SIMD_WIDTH); BBE_SVNX16_IP(D0_bo, SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D0_s,  SECT, 2*BBE_SIMD_WIDTH); BBE_SVNX16_IP(D0_b,  SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D1_s,  SECT, 2*BBE_SIMD_WIDTH); BBE_SVNX16_IP(D1_b,  SECT, 2*BBE_SIMD_WIDTH);
  }

  for ( m=1; m<M-1; m++ )
  {
    BBE_LSNX16_IP(G,  coef, 2); 
    BBE_LSNX16_IP(B0, coef, 2); b12 = BBE_L32X((const int32_t *)coef, 0);
    BBE_LSNX16_IP(B1, coef, 2); 
    BBE_LSNX16_IP(B2, coef, 2); a12 = BBE_L32X((const int32_t *)coef, 0);
    BBE_LSNX16_IP(A1, coef, 2); 
    BBE_LSNX16_IP(A2, coef, 2);

    G  = BBE_REPNX16(G,  0); B0 = BBE_REPNX16(B0, 0); B1 = BBE_REPNX16(B1, 0);
    B2 = BBE_REPNX16(B2, 0); A1 = BBE_REPNX16(A1, 0); A2 = BBE_REPNX16(A2, 0);

    for ( l=0; l<L/BBE_SIMD_WIDTH; l++ )
    {
      D0_so = BBE_LVNX16_I(SECT, 0*BBE_SIMD_WIDTH); D0_bo = BBE_LVNX16_I(SECT,  2*BBE_SIMD_WIDTH);
      D0_s  = BBE_LVNX16_I(SECT, 4*BBE_SIMD_WIDTH); D0_b  = BBE_LVNX16_I(SECT,  6*BBE_SIMD_WIDTH);
      D1_s  = BBE_LVNX16_I(SECT, 8*BBE_SIMD_WIDTH); D1_b  = BBE_LVNX16_I(SECT, 10*BBE_SIMD_WIDTH);
      
      X32_s = BBE_LVNX16_X(SCR, 2*BBE_SIMD_WIDTH*l*2 + 0*BBE_SIMD_WIDTH*2);
      X32_b = BBE_LVNX16_X(SCR, 2*BBE_SIMD_WIDTH*l*2 + 1*BBE_SIMD_WIDTH*2);
                  
      __Pragma("ymemory(SCR)")
      for ( n=0; n<N-1; n++ )
      {
        // Q43 <- Q29*Q14
        Acc_fb_s = BBE_MULUSNX16(D0_s, A2);
        Acc_fb_b = BBE_MULNX16(X32_b, G);
        Acc_ff_s = BBE_MULUSNX16(D0_so, B2);
        Acc_ff_b = BBE_MULNX16(D1_b, B0);

        // Q43 <- Q43 - Q29*Q14
        BBE_MULUSANX16(Acc_fb_s, D1_s, A1);
        BBE_MULUSANX16(Acc_fb_s, X32_s, G);
        BBE_MULUSANX16(Acc_ff_s, D0_s, B1);
        BBE_MULUSANX16(Acc_ff_s, D1_s, B0);
        BBE_MULANX16PR(Acc_fb_b, D0_b,  D1_b, a12);
        BBE_MULANX16PR(Acc_ff_b, D0_bo, D0_b, b12);

        Acc_fb_b = BBE_SLLINX40(Acc_fb_b,  2);
        Acc_fb_s = BBE_SRAINX40(Acc_fb_s, 14);
        Acc_ff_b = BBE_SLLINX40(Acc_ff_b,  2);
        Acc_ff_s = BBE_SRAINX40(Acc_ff_s, 14);

        Acc_fb = BBE_ADDNX40(Acc_fb_s, Acc_fb_b);
        Acc_ff = BBE_ADDNX40(Acc_ff_s, Acc_ff_b);

        D2_s = BBE_PACKLNX40(Acc_fb); D2_b = BBE_PACKVNX40(Acc_fb, rnd16);
        Y_s  = BBE_PACKLNX40(Acc_ff); Y_b  = BBE_PACKVNX40(Acc_ff, rnd16);

        D0_so = D0_s; D0_bo = D0_b;
        D0_s  = D1_s; D0_b  = D1_b;
        D1_s  = D2_s; D1_b  = D2_b;

        X32_s = BBE_LVNX16_X(SCR, 2*(n+1)*L*2 + 2*BBE_SIMD_WIDTH*l*2 + 0*BBE_SIMD_WIDTH*2);
        X32_b = BBE_LVNX16_X(SCR, 2*(n+1)*L*2 + 2*BBE_SIMD_WIDTH*l*2 + 1*BBE_SIMD_WIDTH*2);

        BBE_SVNX16_X(Y_s, SCR, 2*n*L*2 + 2*BBE_SIMD_WIDTH*l*2 + 0*BBE_SIMD_WIDTH);
        BBE_SVNX16_X(Y_b, SCR, 2*n*L*2 + 2*BBE_SIMD_WIDTH*l*2 + 2*BBE_SIMD_WIDTH);
      }

      // Q43 <- Q29*Q14
      Acc_fb_s = BBE_MULUSNX16(D0_s, A2);
      Acc_fb_b = BBE_MULNX16(X32_b, G);
      Acc_ff_s = BBE_MULUSNX16(D0_so, B2);
      Acc_ff_b = BBE_MULNX16(D1_b, B0);

      // Q43 <- Q43 - Q29*Q14
      BBE_MULUSANX16(Acc_fb_s, D1_s, A1);
      BBE_MULUSANX16(Acc_fb_s, X32_s, G);
      BBE_MULUSANX16(Acc_ff_s, D0_s, B1);
      BBE_MULUSANX16(Acc_ff_s, D1_s, B0);
      BBE_MULANX16PR(Acc_fb_b, D0_b,  D1_b, a12);
      BBE_MULANX16PR(Acc_ff_b, D0_bo, D0_b, b12);

      Acc_fb_b = BBE_SLLINX40(Acc_fb_b,  2);
      Acc_fb_s = BBE_SRAINX40(Acc_fb_s, 14);
      Acc_ff_b = BBE_SLLINX40(Acc_ff_b,  2);
      Acc_ff_s = BBE_SRAINX40(Acc_ff_s, 14);

      Acc_fb = BBE_ADDNX40(Acc_fb_s, Acc_fb_b);
      Acc_ff = BBE_ADDNX40(Acc_ff_s, Acc_ff_b);

      D2_s = BBE_PACKLNX40(Acc_fb); D2_b = BBE_PACKVNX40(Acc_fb, rnd16);
      Y_s  = BBE_PACKLNX40(Acc_ff); Y_b  = BBE_PACKVNX40(Acc_ff, rnd16);

      D0_so = D0_s; D0_bo = D0_b;
      D0_s  = D1_s; D0_b  = D1_b;
      D1_s  = D2_s; D1_b  = D2_b;

      BBE_SVNX16_X(Y_s, SCR, 2*(N-1)*L*2 + 2*BBE_SIMD_WIDTH*l*2 + 0*BBE_SIMD_WIDTH);
      BBE_SVNX16_X(Y_b, SCR, 2*(N-1)*L*2 + 2*BBE_SIMD_WIDTH*l*2 + 2*BBE_SIMD_WIDTH);

      BBE_SVNX16_IP(D0_so, SECT, 2*BBE_SIMD_WIDTH); BBE_SVNX16_IP(D0_bo, SECT, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(D0_s,  SECT, 2*BBE_SIMD_WIDTH); BBE_SVNX16_IP(D0_b,  SECT, 2*BBE_SIMD_WIDTH);
      BBE_SVNX16_IP(D1_s,  SECT, 2*BBE_SIMD_WIDTH); BBE_SVNX16_IP(D1_b,  SECT, 2*BBE_SIMD_WIDTH);
    }
  }

  BBE_LSNX16_IP(G,  coef, 2); 
  BBE_LSNX16_IP(B0, coef, 2); b12 = BBE_L32X((const int32_t *)coef, 0);
  BBE_LSNX16_IP(B1, coef, 2); 
  BBE_LSNX16_IP(B2, coef, 2); a12 = BBE_L32X((const int32_t *)coef, 0);
  BBE_LSNX16_IP(A1, coef, 2); 
  BBE_LSNX16_IP(A2, coef, 2);

  G  = BBE_REPNX16(G,  0); B0 = BBE_REPNX16(B0, 0); B1 = BBE_REPNX16(B1, 0);
  B2 = BBE_REPNX16(B2, 0); A1 = BBE_REPNX16(A1, 0); A2 = BBE_REPNX16(A2, 0);

  for ( l=0; l<L/BBE_SIMD_WIDTH; l++ )
  {
    D0_so = BBE_LVNX16_I(SECT, 0*BBE_SIMD_WIDTH); D0_bo = BBE_LVNX16_I(SECT,  2*BBE_SIMD_WIDTH);
    D0_s  = BBE_LVNX16_I(SECT, 4*BBE_SIMD_WIDTH); D0_b  = BBE_LVNX16_I(SECT,  6*BBE_SIMD_WIDTH);
    D1_s  = BBE_LVNX16_I(SECT, 8*BBE_SIMD_WIDTH); D1_b  = BBE_LVNX16_I(SECT, 10*BBE_SIMD_WIDTH);
        
    for ( n=0; n<N; n++ )
    {
      X32_s = BBE_LVNX16_X(SCR, 2*(n)*L*2 + 2*BBE_SIMD_WIDTH*l*2 + 0*BBE_SIMD_WIDTH*2);
      X32_b = BBE_LVNX16_X(SCR, 2*(n)*L*2 + 2*BBE_SIMD_WIDTH*l*2 + 1*BBE_SIMD_WIDTH*2);
      // Q43 <- Q29*Q14
      Acc_fb_s = BBE_MULUSNX16(D0_s, A2);
      Acc_fb_b = BBE_MULNX16(X32_b, G);
      Acc_ff_s = BBE_MULUSNX16(D0_so, B2);
      Acc_ff_b = BBE_MULRNX16(D1_b, B0, rnd12);

      // Q43 <- Q43 - Q29*Q14
      BBE_MULUSANX16(Acc_fb_s, D1_s, A1);
      BBE_MULUSANX16(Acc_fb_s, X32_s, G);
      BBE_MULUSANX16(Acc_ff_s, D0_s, B1);
      BBE_MULUSANX16(Acc_ff_s, D1_s, B0);
      BBE_MULANX16PR(Acc_fb_b, D0_b,  D1_b, a12);
      BBE_MULANX16PR(Acc_ff_b, D0_bo, D0_b, b12);

      Acc_fb_b = BBE_SLLINX40(Acc_fb_b,  2);
      Acc_fb_s = BBE_SRAINX40(Acc_fb_s, 14);
      Acc_ff_b = BBE_SLLINX40(Acc_ff_b,  2);
      Acc_ff_s = BBE_SRAINX40(Acc_ff_s, 14);

      Acc_fb = BBE_ADDNX40(Acc_fb_s, Acc_fb_b);
      Acc_ff = BBE_ADDNX40(Acc_ff_s, Acc_ff_b);

      D2_s = BBE_PACKLNX40(Acc_fb); D2_b = BBE_PACKVNX40(Acc_fb, rnd16);

      D0_so = D0_s; D0_bo = D0_b;
      D0_s  = D1_s; D0_b  = D1_b;
      D1_s  = D2_s; D1_b  = D2_b;

      // Q15 <- ( Q43 - 14 ) - 14 w/ rounding
      Y = BBE_PACKVNX40(Acc_ff, rnd14);

      // Q15 <- Q15*Q8 - 8
      z = BBE_MULNX16(Y, Gain);
      Y = BBE_PACKVNX40(z, rnd8);
      BBE_SVNX16_X(Y, R, 2*n*L + 2*BBE_SIMD_WIDTH*l);
    }
    BBE_SVNX16_IP(D0_so, SECT, 2*BBE_SIMD_WIDTH); BBE_SVNX16_IP(D0_bo, SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D0_s,  SECT, 2*BBE_SIMD_WIDTH); BBE_SVNX16_IP(D0_b,  SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D1_s,  SECT, 2*BBE_SIMD_WIDTH); BBE_SVNX16_IP(D1_b,  SECT, 2*BBE_SIMD_WIDTH);
  }
} // bqriirs_dp_proc()

