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
    Lattice real block IIR, streaming version
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Baseband_iir.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Processing functions declarations. */
#include "latrs_common.h"

#define sz_i16 sizeof(int16_t)

/* Lattice real block IIR processing function, Low Noise fixed-point implementation. */
void latrs_dp_proc4( int16_t * restrict r,
                     int16_t * restrict d,
               const int16_t *          x,
               const int16_t *          coef,
                     int16_t            gain,
                     int N, int L )
{
  xb_vecNx16 d0, d1;
  xb_vecNx16 X, Y, g, Coef0, Coef1, Coef2, Coef3;
  xb_vecNx16 D_l0, D_h0, D_l1, D_h1, D_l2, D_h2, D_l3, D_h3;
  int32_t Coef23;

  xb_vecNx40 T, Z;

  const xb_vecNx16 *_X;
        xb_vecNx16 *_D, *_R;

  vsaN rnd15, rnd16;

  int     j, l;

  ASSERT( r && d && x && coef );
    
  g = BBE_MOVVA16(gain);
  
  rnd15 = BBE_MOVVSA32(15);
  rnd16 = BBE_MOVVSA32(16);

  _D = (      xb_vecNx16 *)d;
  
  Coef23 = (((int32_t)(-coef[2]))<<16);
  Coef23 |= (uint16_t)(-coef[3]);
  Coef0 = BBE_LSNX16_I(coef, sz_i16*0);
  Coef1 = BBE_LSNX16_I(coef, sz_i16*1);
  Coef2 = BBE_LSNX16_I(coef, sz_i16*2);
  Coef3 = BBE_LSNX16_I(coef, sz_i16*3);
 
  Coef0 = BBE_REPNX16(Coef0, 0);
  Coef1 = BBE_REPNX16(Coef1, 0);
  Coef2 = BBE_REPNX16(Coef2, 0);
  Coef3 = BBE_REPNX16(Coef3, 0);
    
  for ( l=0; l<(L>>LOG2_BBE_SIMD_WIDTH); l++ )
  {
    _X = (const xb_vecNx16 *)(x + BBE_SIMD_WIDTH*l);
    _R = (      xb_vecNx16 *)(r + BBE_SIMD_WIDTH*l);
    
    D_l0 = BBE_LVNX16_I(_D, sz_i16*0*BBE_SIMD_WIDTH); 
    D_h0 = BBE_LVNX16_I(_D, sz_i16*1*BBE_SIMD_WIDTH);

    D_l1 = BBE_LVNX16_I(_D, sz_i16*2*BBE_SIMD_WIDTH); 
    D_h1 = BBE_LVNX16_I(_D, sz_i16*3*BBE_SIMD_WIDTH);

    D_l2 = BBE_LVNX16_I(_D, sz_i16*4*BBE_SIMD_WIDTH); 
    D_h2 = BBE_LVNX16_I(_D, sz_i16*5*BBE_SIMD_WIDTH);

    D_l3 = BBE_LVNX16_I(_D, sz_i16*6*BBE_SIMD_WIDTH); 
    D_h3 = BBE_LVNX16_I(_D, sz_i16*7*BBE_SIMD_WIDTH);

    for ( j=0; j<N; j++ )
    {
      BBE_LVNX16_XP(X, _X, sz_i16*L);
      // Q29 <- Q15*Q15 - 1
      T = BBE_MULNX16(X, g);
      T = BBE_SRAINX40(T, 1);
      
      // Q29 <- Q29 - ( Q29*Q15 - 15 w/ rounding )
      BBE_MULANX16PR(T, D_h2, D_h3, Coef23);
      BBE_MULANX16PR(T, D_h2, D_h3, Coef23);
      
      Z = BBE_MULUSRNX16(D_l3, Coef3, rnd15);
      BBE_MULUSANX16(Z,  D_l2, Coef2);

      Z = BBE_SRAINX40(Z, 15);
      T = BBE_SUBNX40(T, Z);
      Y = BBE_PACKQNX40(T);

      // Q29        
      d0 = BBE_SELNX16I(D_h2, D_l2, BBE_SELI_INTERLEAVE_1_LO);
      d1 = BBE_SELNX16I(D_h2, D_l2, BBE_SELI_INTERLEAVE_1_HI);

      Z = BBE_MOVSWV(d1, d0);

      BBE_MULANX16(Z, Y, Coef2);

      D_l3 = BBE_PACKLNX40(Z); D_h3 = BBE_PACKVNX40(Z, rnd16);
      
      // Q29 <- Q29 - ( Q29*Q15 - 15 w/ rounding )
      BBE_MULSNX16(T, D_h1, Coef1);
      BBE_MULSNX16(T, D_h1, Coef1);
      
      Z = BBE_MULUSRNX16(D_l1, Coef1, rnd15);
      Z = BBE_SRAINX40(Z, 15);

      T = BBE_SUBNX40(T, Z);

      Y = BBE_PACKQNX40(T);

      // Q29        
      d0 = BBE_SELNX16I(D_h1, D_l1, BBE_SELI_INTERLEAVE_1_LO);
      d1 = BBE_SELNX16I(D_h1, D_l1, BBE_SELI_INTERLEAVE_1_HI);

      Z = BBE_MOVSWV(d1, d0);

      BBE_MULANX16(Z, Y, Coef1);

      D_l2 = BBE_PACKLNX40(Z); D_h2 = BBE_PACKVNX40(Z, rnd16);
      
      // Q29 <- Q29 - ( Q29*Q15 - 15 w/ rounding )
      BBE_MULSNX16(T, D_h0, Coef0);
      BBE_MULSNX16(T, D_h0, Coef0);
      
      Z = BBE_MULUSRNX16(D_l0, Coef0, rnd15);
      Z = BBE_SRAINX40(Z, 15);

      T = BBE_SUBNX40(T, Z);
      
      // Q15 <- Q29 - 14 w/ rounding
      Z = BBE_ADDNX40(T, T);
      Y = BBE_PACKQNX40(Z);
      BBE_SVNX16_XP(Y, _R, sz_i16*L);
      
      // Q29
      Y = BBE_PACKQNX40(T);
      d0 = BBE_SELNX16I(D_h0, D_l0, BBE_SELI_INTERLEAVE_1_LO);
      d1 = BBE_SELNX16I(D_h0, D_l0, BBE_SELI_INTERLEAVE_1_HI);

      Z = BBE_MOVSWV(d1, d0);

      BBE_MULANX16(Z, Y, Coef0);

      D_l1 = BBE_PACKLNX40(Z); D_h1 = BBE_PACKVNX40(Z, rnd16);
      D_l0 = BBE_PACKLNX40(T); D_h0 = BBE_PACKVNX40(T, rnd16);
    }
    
    BBE_SVNX16_IP(D_l0, _D, sz_i16*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D_h0, _D, sz_i16*BBE_SIMD_WIDTH);
    
    BBE_SVNX16_IP(D_l1, _D, sz_i16*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D_h1, _D, sz_i16*BBE_SIMD_WIDTH);

    BBE_SVNX16_IP(D_l2, _D, sz_i16*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D_h2, _D, sz_i16*BBE_SIMD_WIDTH);

    BBE_SVNX16_IP(D_l3, _D, sz_i16*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D_h3, _D, sz_i16*BBE_SIMD_WIDTH);
  }
} // latr_dp_proc4()
