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

/* Lattice real block IIR processing function, Fast fixed-point implementation. */
void latrs_sp_proc7( int16_t * restrict r,
                     int16_t * restrict d,
               const int16_t *          x,
               const int16_t *          coef,
                     int16_t            gain,
                     int N, int L )
{
  xb_vecNx16 X, g, D, D0, D1, D2, Y, Y0;
  xb_vecNx16 Coef0,Coef1, Coef2, Coef3, Coef4, Coef5, Coef6;
  xb_vecNx16 c0, c1, c2, c3;
  xb_vecNx40 T, U;

  const xb_vecNx16 *_X;
        xb_vecNx16 *_D, *_R;
 
  int     j, l;

  ASSERT( r && d && x && coef );
  
  g = BBE_MOVVA16(gain);

  c1 = BBE_LPNX16_I( coef, 2*sz_i16);
  c2 = BBE_LPNX16_I( coef, 4*sz_i16);
  Coef2 = BBE_REPNX16(c1, 0);
  Coef3 = BBE_REPNX16(c1, 1);
  Coef4 = BBE_REPNX16(c2, 0);
  Coef5 = BBE_REPNX16(c2, 1);
  
  for ( l=0; l<(L>>LOG2_BBE_SIMD_WIDTH); l++)
  {
    _X = (const xb_vecNx16 *)(x +   BBE_SIMD_WIDTH*l);
    _D = (      xb_vecNx16 *)(d + 7*BBE_SIMD_WIDTH*l);
    _R = (      xb_vecNx16 *)(r +   BBE_SIMD_WIDTH*l);
    
    D0 = BBE_LVNX16_I(_D, sz_i16*0*BBE_SIMD_WIDTH);
    D1 = BBE_LVNX16_I(_D, sz_i16*1*BBE_SIMD_WIDTH);
    D2 = BBE_LVNX16_I(_D, sz_i16*2*BBE_SIMD_WIDTH);

    for ( j=0; j<N; j++ )
    {
      BBE_LVNX16_XP(X, _X, sz_i16*L);
      // Q29 <- Q15*Q15 - 1
      T = BBE_MULNX16(X, g);
      T = BBE_SRAINX40(T, 1);

      c0 = BBE_LPNX16_I( coef, 0);
      Coef0 = BBE_REPNX16(c0, 0);
      Coef1 = BBE_REPNX16(c0, 1);

      c3 = BBE_LSNX16_I( coef, 6*sz_i16);
      Coef6 = BBE_REPNX16(c3, 0);

      D = BBE_LVNX16_I(_D, sz_i16*6*BBE_SIMD_WIDTH);

      // Q29 <- Q29 - Q14*Q15
      BBE_MULSNX16(T, D, Coef6); 
    
      D = BBE_LVNX16_I(_D, sz_i16*5*BBE_SIMD_WIDTH);

      // Q29 <- Q29 - Q14*Q15
      BBE_MULSNX16(T, D, Coef5); 

      // Q29 <- (Q29 - 15 w/ rounding)*Q15
      Y = BBE_PACKQNX40(T);   

      // Q29 <- (Q29 - 15 w/ rounding)*Q15
      Y = BBE_MULNX16PACKQ(Y, Coef5);

      D = BBE_ADDNX16(D, Y);

      BBE_SVNX16_I(D, _D, sz_i16*6*BBE_SIMD_WIDTH);

      D = BBE_LVNX16_I(_D, sz_i16*4*BBE_SIMD_WIDTH);

      // Q29 <- Q29 - Q14*Q15
      BBE_MULSNX16(T, D, Coef4); 

      // Q29 <- (Q29 - 15 w/ rounding)*Q15
      Y = BBE_PACKQNX40(T);

      // Q29 <- (Q29 - 15 w/ rounding)*Q15
      Y = BBE_MULNX16PACKQ(Y, Coef4);

      D = BBE_ADDNX16(D, Y);

      BBE_SVNX16_I(D, _D, sz_i16*5*BBE_SIMD_WIDTH);

      D = BBE_LVNX16_I(_D, sz_i16*3*BBE_SIMD_WIDTH);

      // Q29 <- Q29 - Q14*Q15
      BBE_MULSNX16(T, D, Coef3); 

      // Q29 <- (Q29 - 15 w/ rounding)*Q15
      Y = BBE_PACKQNX40(T);

      // Q29 <- (Q29 - 15 w/ rounding)*Q15
      Y = BBE_MULNX16PACKQ(Y, Coef3);

      D = BBE_ADDNX16(D, Y);

      BBE_SVNX16_I(D, _D, sz_i16*4*BBE_SIMD_WIDTH);

      // Q29 <- Q29 - Q14*Q15
      BBE_MULSNX16(T, D2, Coef2); 

      // Q29 <- (Q29 - 15 w/ rounding)*Q15
      Y = BBE_PACKQNX40(T);

      // Q29 <- (Q29 - 15 w/ rounding)*Q15
      Y = BBE_MULNX16PACKQ(Y, Coef2);

      D = BBE_ADDNX16(D2, Y);

      BBE_SVNX16_I(D, _D, sz_i16*3*BBE_SIMD_WIDTH);

      // Q29 <- Q29 - Q14*Q15
      BBE_MULSNX16(T, D1, Coef1); 

      // Q29 <- (Q29 - 15 w/ rounding)*Q15
      Y = BBE_PACKQNX40(T);

      // Q29 <- (Q29 - 15 w/ rounding)*Q15
      Y = BBE_MULNX16PACKQ(Y, Coef1);

      D2 = BBE_ADDNX16(D1, Y);

      // Q29 <- Q29 - Q14*Q15
      BBE_MULSNX16(T, D0, Coef0); 

      // Q29 <- (Q29 - 15 w/ rounding)*Q15
      Y = BBE_PACKQNX40(T); 

      // Q29 <- (Q29 - 15 w/ rounding)*Q15
      D = BBE_MULNX16PACKQ(Y, Coef0);

      D1 = BBE_ADDNX16(D0, D); D0 = Y;
      
      // Q15 <- Q29 - 14 w/ rounding
      U = BBE_ADDNX40(T, T);
      Y0 = BBE_PACKQNX40(U);
      BBE_SVNX16_XP(Y0, _R, sz_i16*L);
    }
    
    BBE_SVNX16_I(D0, _D, 0*sz_i16*BBE_SIMD_WIDTH);    
    BBE_SVNX16_I(D1, _D, 1*sz_i16*BBE_SIMD_WIDTH);    
    BBE_SVNX16_I(D2, _D, 2*sz_i16*BBE_SIMD_WIDTH);
  }

} // latr_sp_proc7()
