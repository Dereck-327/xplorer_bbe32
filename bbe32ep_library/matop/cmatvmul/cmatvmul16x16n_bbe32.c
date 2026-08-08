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
    Complex Matrix-Matrix/Matrix-Vector Multiply
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matop.h"

#if !(HAVE_MULPC && HAVE_PACKEDMUL && 1)
DISCARD_FUN(void, cmatvmul16x16n,(complex_fract16 * restrict z, 
            const complex_fract16 * restrict x, 
            const complex_fract16 * restrict y, 
            int L, int Q))
#else
/* multiply 8x8 by 8x1 */
#define CMUL8X8x8X1(r,X,Y,q,isUpdate)                          \
{                                                              \
    xb_vecNx16 XX,YY;                                          \
    /* X0=x00x01 x20x21 x02x03 x22x23 */                       \
    /* X2=x04x05 x24x25 x06x07 x26x27 */                       \
    /* X4=x40x41 x60x61 x42x43 x62x63 */                       \
    /* X6=x44x45 x64x65 x46x47 x66x67 */                       \
    /* X1=x10x11 x30x31 x12x13 x32x33 */                       \
    /* X3=x14x15 x34x35 x16x17 x36x37 */                       \
    /* X5=x50x51 x70x71 x52x53 x72x73 */                       \
    /* X7=x54x55 x74x75 x56x57 x76x77 */                       \
    BBE_DSELNX16I(X##2,X##0,X##2,X##0,BBE_DSELI_INTERLEAVE_4); \
    BBE_DSELNX16I(X##6,X##4,X##6,X##4,BBE_DSELI_INTERLEAVE_4); \
    BBE_DSELNX16I(X##3,X##1,X##3,X##1,BBE_DSELI_INTERLEAVE_4); \
    BBE_DSELNX16I(X##7,X##5,X##7,X##5,BBE_DSELI_INTERLEAVE_4); \
    /* multiply: */                                            \
    /* x00x01 x20x21 x40x41 x60x61  y0y1y0y1y0y1y0y1 */        \
    XX=BBE_SELNX16I(X##4,X##0,BBE_SELI_EXTRACT_LO_HALVES);     \
    YY=BBE_SHFLNX16I(Y,BBE_SHFLI_REP_0X4);                     \
    if (isUpdate) BBE_MULANX16PC_0(r,XX,YY);                   \
    else          r=BBE_MULRNX16PC_0(XX,YY,q);                 \
    XX=BBE_SELNX16I(X##5,X##1,BBE_SELI_EXTRACT_LO_HALVES);     \
    BBE_MULANX16PC_1(r,XX,YY);                                 \
    /* x02x03 x22x23 x42x43 x62x63  y2y3y2y3y2y3y2y3 */        \
    XX=BBE_SELNX16I(X##4,X##0,BBE_SELI_EXTRACT_HI_HALVES);     \
    YY=BBE_SHFLNX16I(Y,BBE_SHFLI_REP_1X4);                     \
    BBE_MULANX16PC_0(r,XX,YY);                                 \
    XX=BBE_SELNX16I(X##5,X##1,BBE_SELI_EXTRACT_HI_HALVES);     \
    BBE_MULANX16PC_1(r,XX,YY);                                 \
    /* x04x05 x24x25 x44x45 x64x65  y4y5y4y5y4y5y4y5 */        \
    XX=BBE_SELNX16I(X##6,X##2,BBE_SELI_EXTRACT_LO_HALVES);     \
    YY=BBE_SHFLNX16I(Y,BBE_SHFLI_REP_2X4);                     \
    BBE_MULANX16PC_0(r,XX,YY);                                 \
    XX=BBE_SELNX16I(X##7,X##3,BBE_SELI_EXTRACT_LO_HALVES);     \
    BBE_MULANX16PC_1(r,XX,YY);                                 \
    /* x06x07 x26x27 x46x47 x66x67  y6y7y6y7y6y7y6y7 */        \
    XX=BBE_SELNX16I(X##6,X##2,BBE_SELI_EXTRACT_HI_HALVES);     \
    YY=BBE_SHFLNX16I(Y,BBE_SHFLI_REP_3X4);                     \
    BBE_MULANX16PC_0(r,XX,YY);                                 \
    XX=BBE_SELNX16I(X##7,X##3,BBE_SELI_EXTRACT_HI_HALVES);     \
    BBE_MULANX16PC_1(r,XX,YY);                                 \
}

/*-------------------------------------------------------------------------
Complex Matrix-Matrix/Matrix-Vector Multiply

Description: These functions perform pairwise multiplication of two 
sequences of complex matrices or vectors. Both the block order and 
streaming order are allowed for input/output matrix sequences.

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
x[L*Sx]     Sequence of left-hand complex matrices
y[L*Sy]     Sequence of right-hand complex matrices
M           Matrix dimension 
N           Matrix dimension (columnar for MxN)
L           Number of matrices 
Q           Position of fractional point in matrix representation, 0..16
Output:
z[L*Sz]     Sequence of complex result matrices

Restrictions:
pScr,x,y,z  Aligned on 32-byte boundary
pScr,x,y,z  Must not overlap

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

/* Block Order, 16x16*16x1->16x1, Sx=256, Sy=16, Sz=16
   Restrictions:
     None
*/
void cmatvmul16x16n ( complex_fract16 * restrict z, 
                const complex_fract16 * restrict x, 
                const complex_fract16 * restrict y, 
                int L, int Q )
{
  int l, m, ystride;

  const xb_vecNx16 *px = (const xb_vecNx16 *)x;
  const xb_vecNx16 *py = (const xb_vecNx16 *)y;
  xb_vecNx16 X0, X1, Y0, Y1, Z;
  xb_vecNx40 r;
  vsaN  q = BBE_MOVVSA32(Q);

  /* check restrictions */
  NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
  NASSERT_ALIGN(z, 2 * BBE_SIMD_WIDTH);
  NASSERT(Q >= 0 && Q <= 16);
  if (L <= 0) return;

  px = (const xb_vecNx16 *)x;
  py = (const xb_vecNx16 *)y;
  __Pragma("loop_count min=8 factor=8");
  for (m = l = 0; l<L * 8; l++)
  {
    xb_c40 t;
    /* m=(m+2)&15; */
    m = BBE_ADDMOD16U(m, 0x100002);    
    ystride = -2 * BBE_SIMD_WIDTH;
    /* ystride= (m==0)? 2*BBE_SIMD_WIDTH:-2*BBE_SIMD_WIDTH; */
    XT_MOVEQZ(ystride, 2 * BBE_SIMD_WIDTH, m);

    BBE_LVNX16_IP(Y0, py, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_XP(Y1, py, ystride);

    BBE_LVNX16_IP(X0, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1, px, 2 * BBE_SIMD_WIDTH);
    r = BBE_MULNX16C(X0, Y0); BBE_MULANX16C(r, X1, Y1);
    t = BBE_RADDNX40C(r); r = BBE_MOVNX40_FROMC40(t);
    Z = BBE_PACKVNX40(r, q);
    BBE_SPNX16_IP(Z, z, 4);

    BBE_LVNX16_IP(X0, px, 2 * BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(X1, px, 2 * BBE_SIMD_WIDTH);
    r = BBE_MULNX16C(X0, Y0); BBE_MULANX16C(r, X1, Y1);
    t = BBE_RADDNX40C(r); r = BBE_MOVNX40_FROMC40(t);
    Z = BBE_PACKVNX40(r, q);
    BBE_SPNX16_IP(Z, z, 4);
  }
} /* cmatvmul16x16n() */
#endif
