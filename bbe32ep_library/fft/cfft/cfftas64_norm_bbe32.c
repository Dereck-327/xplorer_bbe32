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
  NatureDSP_Baseband library. FFT
    Radix-2 forward FFT on complex data, auto scaling
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
#include "fft_common.h"

#if !(HAVE_FFT && 1)
DISCARD_FUN(int, cfftas64_norm, (complex_fract16 * restrict y,complex_fract16 * restrict x))
#else
/*-------------------------------------------------------------------------
Radix-2 forward FFT on complex data, auto scaling

Description: These functions make forward FFT on complex data of power of 2
sizes: N=2^n, n=4..15. Functions with _norm suffix expect input data to be
normalized, i.e. the minimum number of redundant sign bits over x[]
(a.k.a the common block exponent) should be zero. Neglecting to normalize
data leads to significant loss in transform quality. On the contrary, regular
variants with no _norm suffix allow for non-zero common block exponent, but
they appear slightly slower due to internal data normalization.

Precision: 16-bit input, 16-bit output
Scaling  : Automatic data scaling at each stage

NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
     the call.

Parameters:
  Input:            
    x[N]        Complex input signal
    bexp        Common block exponent, that is the minimum number of redundant
                sign bits over input data x[]
  Output:            
    y[N]        Output spectrum samples
  Returned value:
                Total shift amount applied throughout the transform to scale
                the data, with positive numbers corresponding to the right
                shift. _norm-suffixed functions return strictly positive
                values, while for regular variants the total shift amount is
                bi-directional.
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
/* First stage radix 4, 
   suffix _mr mean : memory is input, registers is output  */

inline_ int R1_DFT4_L64_16_norm_mr(const int16_t *tw, int16_t *x, xb_vecNx16 *y,  int N)

{

  int count = N/4/(BBE_SIMD_WIDTH/2);                                        
  int stride = N/4*2*sizeof(int16_t); 
  VT * p_tw = (VT *)(tw);                                        
  xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                       
                                             
  VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
  xb_vecNx16 _t0, _t1, _t2, _t3;  
  const int bexp = 0; 
  int scaling = bexp;



  RANGE_BEGIN(4, -1, 1, 0, scaling);
    

  ASSERT(count==2); 
  if(count==2)
  {   
    NASSERT_ALIGN32(p_src);
    BBE_LVNX16_XP(_t0, p_src, stride);    
    NASSERT_ALIGN32(p_src);
    BBE_LVNX16_XP(_t1, p_src, stride);
    NASSERT_ALIGN32(p_src);
    BBE_LVNX16_XP(_t2, p_src, stride);
    NASSERT_ALIGN32(p_src);
    BBE_LVNX16_XP(_t3, p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
                                                                        
                                  
                                                                        
    t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
    t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

    BBE_MOVSAV(_t2);                                                  
    BBE_MOVSBV(_t3);                                                  
              
    BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
    BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
    BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

    __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           

    BBE_RANGENX16(_t0); 
    BBE_RANGENX16(_t1); 
    BBE_RANGENX16(_t2); 
    BBE_RANGENX16(_t3);


     
    BBE_DSELN_2XCQ15I(t1, t0, t1, t0, BBE_DSELI_INTERLEAVE_2);                                                 
    BBE_DSELN_2XCQ15I(t3, t2, t3, t2, BBE_DSELI_INTERLEAVE_2);          
    BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_4);

    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           


    y[0]=_t0; //BBE_SVRNX16_XP(_t0, p_dst,    2*BBE_SIMD_WIDTH);                   
    y[1]=_t2; //BBE_SVRNX16_XP(_t2, p_dst,    2*BBE_SIMD_WIDTH);            
    y[2]=_t1; //BBE_SVRNX16_XP(_t1, p_dst,    2*BBE_SIMD_WIDTH);                
    y[3]=_t3; //BBE_SVRNX16_XP(_t3, p_dst,    2*BBE_SIMD_WIDTH);
  }
  {   
    NASSERT_ALIGN32(p_src);
    BBE_LVNX16_XP(_t0, p_src, stride);    
    NASSERT_ALIGN32(p_src);
    BBE_LVNX16_XP(_t1, p_src, stride);
    NASSERT_ALIGN32(p_src);
    BBE_LVNX16_XP(_t2, p_src, stride);
    NASSERT_ALIGN32(p_src);
    BBE_LVNX16_XP(_t3, p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
                                                                        
                              
    t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
    t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

    BBE_MOVSAV(_t2);                                                  
    BBE_MOVSBV(_t3);                                                  
              
    BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
    BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
    BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

    __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
 
    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           

    BBE_RANGENX16(_t0); 
    BBE_RANGENX16(_t1); 
    BBE_RANGENX16(_t2); 
    BBE_RANGENX16(_t3);

    BBE_DSELN_2XCQ15I(t1, t0, t1, t0, BBE_DSELI_INTERLEAVE_2);                                                 
    BBE_DSELN_2XCQ15I(t3, t2, t3, t2, BBE_DSELI_INTERLEAVE_2);          
    BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_4);
    BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_4);

    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           

    y[4]=_t0; //BBE_SVRNX16_XP(_t0, p_dst,    2*BBE_SIMD_WIDTH);                   
    y[5]=_t2; //BBE_SVRNX16_XP(_t2, p_dst,    2*BBE_SIMD_WIDTH);            
    y[6]=_t1; //BBE_SVRNX16_XP(_t1, p_dst,    2*BBE_SIMD_WIDTH);                
    y[7]=_t3; //BBE_SVRNX16_XP(_t3, p_dst,    2*BBE_SIMD_WIDTH);
  }

  return scaling; 
}

int cfftas64_norm ( complex_fract16 * restrict y,complex_fract16 * restrict x           )
{
  int scaling;
  xb_vecNx16 a[8], b[8];

  NASSERT_ALIGN32(y);
  NASSERT_ALIGN32(x);

  scaling = R1_DFT4_L64_16_norm_mr(fft64_tw1, (int16_t*)x, /*y*/a, 64);
  scaling += R2_DFT4xI4_rr(fft64_tw2, a/* y*/, b, 64, BBE_RRANGE());
  scaling += R2_DFT4xIN_4_rm(b, (int16_t*)y, 64, BBE_RRANGE());

  return scaling;
} /* cfftas64_norm() */
#endif
