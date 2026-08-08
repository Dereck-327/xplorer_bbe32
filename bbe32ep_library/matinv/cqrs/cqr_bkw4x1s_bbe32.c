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
    cqr_bkwNxPs/qr_bkwNxPs
    C code optimized for BBE32
  IntegrIT, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"

#if (HAVE_DIV && 1)

#define CURRENT_M 4
#define CURRENT_P 1

/*-------------------------------------------------------------------------
cqr_bkwNxPs/qr_bkwNxPs

Last stage of solving a set of L complex-valued linear problems A*X=B
through the QR decomposition by Householder reflections: back substitution
process for L systems of complex-valued linear equations R*X=QB, where R is
an MxM upper triangular matrix, X is an MxP matrix of unknowns, QB is an MxP
matrix resulting from Householder reflections being applied to the right
hand matrix B of the original linear problem: QB=Q'*B.

Fixed-point representation for output data is a function of fixed-point
format of input data: FPP(X) = FPP(QB)-FPP(R)+10, where FPP(x) stands for
the Fixed-Point Position of data item x.

Data transform is performed in-place.

NOTE:
1. Data layout for matrices is selected as for other matrices written 
   in a stream order. So, shorter dimension of output matrix B (NxP 
   instead of MxP as on input) does not require special management - 
   remaining (M-N)*P*L elements are kept unchanged

Input
B[M*P][L]  Matrices QB=Q'*B (L matrices of size MxP)
R[M*N][L]  upper triangular matrices R (L matrices of size MxN)
Output:
B[N*P][L]  Matrices X (L matrices of size NxP)

Restrictions:
1. All matrices must not overlap an must be aligned on 32-byte boundary 
2. Number of matrices L must be a multiple of 8 for complex data and 
   16 for real data
3. Scratch memory must be aligned on 32-byte boundary. Its size (in bytes)
   is defined by xxx_getScratchSize(N,P,L)
4. Matrix sizes N,L must be greater than 1
---------------------------------------------------------------------------*/
void cqr_bkw4x1s (void* pScr, complex_fract16* restrict _B, const complex_fract16* restrict _R, int L)
{
          int16_t* restrict B=(      int16_t*)_B;
    const int16_t* restrict R=(const int16_t*)_R;
    int i;
    int16_t *invDiagR = (int16_t *)pScr;
    int16_t *QB = B;

    const int16_t *start_R = R;
    xb_vecNx16 * restrict pf ;
    xb_vecNx40 z0;
    xb_vecNx16 z4;
    xb_vecNx16 c18 = BBE_MOVVINT16(19-1);
    xb_vecNx16 c1918= BBE_MOVVINT16(19);

    const xb_vecNx16 * restrict  pBr;
    xb_vecNx16 * restrict  pBw;

    xb_vecNx16 * restrict  pr12;
    xb_vecNx16 * restrict  pr13;
    xb_vecNx16 * restrict  pr14;
    xb_vecNx16 * restrict  pr23;
    xb_vecNx16 * restrict  pr24;
    xb_vecNx16 * restrict  pr34;
    xb_vecNx16 * restrict  pr11;
    xb_vecNx16 * restrict  pr22;
    xb_vecNx16 * restrict  pr33;
    xb_vecNx16 * restrict  pr44;
//cqr_bkwnxps(pScr,B,R,4,1,L);return 0;
    c1918=BBE_SELNX16I(c18,c1918,BBE_SELI_EXTRACT_LO_HALVES);

#ifdef COMPILER_XTENSA
#pragma ymemory( pr12 )
#pragma ymemory( pr13 )
#pragma ymemory( pr14 )
#pragma ymemory( pr23 )
#pragma ymemory( pr24 )
#pragma ymemory( pr34 )
#endif

    pf = (xb_vecNx16 *) invDiagR ;

    z0 = BBE_MOVWA32(0x2000000);
    z4 = BBE_MOVVINT16(4);
    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT(L>0 && L%(BBE_SIMD_WIDTH/2)==0);
//  return cqr_bkwnxps(pScr, B, R, 4, 1, L);
    pr11 = (xb_vecNx16 *)R;
    pr22 = (xb_vecNx16 *)(R+1*(CURRENT_M+1)*L*2);
    pr33 = (xb_vecNx16 *)(R+2*(CURRENT_M+1)*L*2);
    pr44 = (xb_vecNx16 *)(R+3*(CURRENT_M+1)*L*2);
    
    __Pragma("loop_count min=1");
    for(i=0; i<(L>>(LOG2_BBE_SIMD_WIDTH-1)); i++)
    {
        xb_vecNx16 f0,nsa,ex;
        vsaN nsa_vsa;
        xb_vecNx16 x0, x1, x2, x3;

        BBE_LVNX16_IP(x0, pr44, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_IP(x1, pr33, (2*BBE_SIMD_WIDTH));
        x0 = BBE_SELNX16I ( x1, x0,  BBE_SELI_EXTRACT_1_OF_2_OFF_0);
        nsa_vsa = BBE_NSANX16(x0);
        nsa = BBE_MOVVVS(nsa_vsa);
        ex = BBE_SUBNX16(c1918, nsa);
        nsa = BBE_SUBNX16(nsa,z4);
        nsa_vsa = BBE_MOVVSV(nsa, 0);
        x0 =  BBE_SLANX16(x0,nsa_vsa);
        f0 = BBE_QUONX32(z0,x0);
        // note: here we have ex lower by 1 in half word to minimize operations in the next loop!
        BBE_SVNX16_IP(f0, pf,(2*BBE_SIMD_WIDTH));
        BBE_SVNX16_IP(ex, pf,(2*BBE_SIMD_WIDTH));

        BBE_LVNX16_IP(x2, pr22, (2*BBE_SIMD_WIDTH));
        BBE_LVNX16_IP(x3, pr11, (2*BBE_SIMD_WIDTH));
        x2 = BBE_SELNX16I ( x3, x2,  BBE_SELI_EXTRACT_1_OF_2_OFF_0);
        nsa_vsa = BBE_NSANX16(x2);
        nsa = BBE_MOVVVS(nsa_vsa);
        ex = BBE_SUBNX16(c18, nsa);
        nsa = BBE_SUBNX16(nsa,z4);
        nsa_vsa = BBE_MOVVSV(nsa, 0);
        x2 =  BBE_SLANX16(x2,nsa_vsa);
        f0 = BBE_QUONX32(z0,x2);
        BBE_SVNX16_IP(f0, pf,(2*BBE_SIMD_WIDTH));
        BBE_SVNX16_IP(ex, pf,(2*BBE_SIMD_WIDTH));
    }

    R = start_R;
    {

        const int start_m = CURRENT_M-1;

        pr34 = (xb_vecNx16 *) (R + (start_m-1)*CURRENT_M*L*2/*row = start_m-1*/ +start_m*L*2 /*column = start_m*/);

        pr14 = pr34 -  2*CURRENT_M*L*2*sizeof(int16_t)/sizeof(*pr24);
        pr13 = pr34 - (2*CURRENT_M+1)*L*2*sizeof(int16_t)/sizeof(*pr24);
        pr12 = pr34 - (2*CURRENT_M+2)*L*2*sizeof(int16_t)/sizeof(*pr24);

        pr24 = pr34 -  1*CURRENT_M   *L*2*sizeof(int16_t)/sizeof(*pr24);
        pr23 = pr34 - (1*CURRENT_M+1)*L*2*sizeof(int16_t)/sizeof(*pr24);

        pf = (xb_vecNx16 *) invDiagR ;

        pBr=(const xb_vecNx16 *) (QB+3*CURRENT_P*L*2);
        pBw=(      xb_vecNx16 *) (QB+3*CURRENT_P*L*2);
        __Pragma("loop_count min=1");
        for(i=0; i<(L>>(LOG2_BBE_SIMD_WIDTH-1)); i++)
        {
            vsaN _11=BBE_MOVVSA32(11);
            vsaN v_ex;
            xb_vecNx16 dp;
            xb_vecNx16 r0,  ex, f0;
            xb_vecNx40 acc0;
            xb_vecNx40 acc02;
            xb_vecNx40 acc01;
            xb_vecNx40 acc00;
            xb_vecNx16 b0, b1, b2, b3;

            BBE_LVNX16_XP(b3,pBr, -CURRENT_P*L*4);
            BBE_LVNX16_XP(b2,pBr, -CURRENT_P*L*4);
            BBE_LVNX16_XP(b1,pBr, -CURRENT_P*L*4);
            BBE_LVNX16_XP(b0,pBr,  (2*BBE_SIMD_WIDTH)+3*CURRENT_P*L*4);

            /************** row 3 *********************/
            f0=BBE_LVNX16_X( pf,  0);
            ex=BBE_LVNX16_X( pf, (2*BBE_SIMD_WIDTH));

            f0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_LO);
            ex = BBE_SHFLNX16I(ex, BBE_SHFLI_DOUBLE_1_LO);
            v_ex=BBE_MOVVSV(ex,0);
            acc0 = BBE_MULRNX16( f0, b3, v_ex);
            b3 = BBE_PACKVNX40(acc0, v_ex);
            /************** row 2 *********************/
            BBE_LVNX16_XP(f0, pf, (2*BBE_SIMD_WIDTH));
            BBE_LVNX16_XP(ex, pf, (2*BBE_SIMD_WIDTH));
            f0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_HI);
            ex = BBE_SHFLNX16I(ex, BBE_SHFLI_DOUBLE_1_HI);
            BBE_LVNX16_XP( r0, pr34, sizeof(*pr34));
            acc02 = BBE_MULRNX16C(b3, r0,_11);
            BBE_LVNX16_XP( r0, pr24, sizeof(*pr34));
            acc01 = BBE_MULRNX16C(b3, r0,_11);
            BBE_LVNX16_XP( r0, pr14, sizeof(*pr14));
            acc00 = BBE_MULRNX16C(b3, r0,_11);
            BBE_SVNX16_XP(b3,pBw, -CURRENT_P*L*4);
            b2=BBE_SRAINX16(b2,1);
            dp = BBE_PACKVNX40(acc02,_11);
            acc0 = BBE_MULRNX16(f0, b2,ex);
            BBE_MULSNX16(acc0,f0,dp);
            b2 = BBE_PACKVNX40(acc0, ex);
            /************** row 1 *********************/
            f0=BBE_LVNX16_I( pf,  0);
            ex=BBE_LVNX16_I( pf, (2*BBE_SIMD_WIDTH));
            f0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_LO);
            ex = BBE_SHFLNX16I(ex, BBE_SHFLI_DOUBLE_1_LO);
            BBE_LVNX16_XP( r0, pr23, sizeof(*pr34));
            BBE_MULANX16C( acc01, b2, r0);
            BBE_LVNX16_XP( r0, pr13, sizeof(*pr14));
            BBE_MULANX16C( acc00, b2, r0);
            BBE_SVNX16_XP(b2,pBw, -CURRENT_P*L*4);
            dp = BBE_PACKPNX40(acc01);
            b1=BBE_SRAINX16(b1,1);
            dp = BBE_PACKVNX40( acc01,_11);
            acc0 = BBE_MULRNX16(f0, b1,ex);
            BBE_MULSNX16(acc0,f0,dp);
            b1 = BBE_PACKVNX40(acc0, ex);

            /************** row 0 *********************/
            BBE_LVNX16_XP(f0, pf, (2*BBE_SIMD_WIDTH));
            BBE_LVNX16_XP(ex, pf, (2*BBE_SIMD_WIDTH));
            f0 = BBE_SHFLNX16I(f0, BBE_SHFLI_DOUBLE_1_HI);
            ex = BBE_SHFLNX16I(ex, BBE_SHFLI_DOUBLE_1_HI);
            BBE_LVNX16_XP( r0, pr12, sizeof(*pr14));
            BBE_MULANX16C(acc00, b1, r0);
            BBE_SVNX16_XP(b1,pBw, -CURRENT_P*L*4);
            b0=BBE_SRAINX16(b0,1);
            dp = BBE_PACKVNX40( acc00,_11);
            acc0 = BBE_MULRNX16( f0, b0,ex);
            BBE_MULSNX16(acc0,f0,dp);
            b0 = BBE_PACKVNX40(acc0, ex);

            BBE_SVNX16_XP(b0,pBw, (2*BBE_SIMD_WIDTH)+3*CURRENT_P*L*4);
        }
    }
} /* cqr_bkw4x1s() */

size_t cqr_bkw4x1s_getScratchSize (int N, int P, int L)
{
    (void)P;
    return L*4*N*sizeof(int16_t);
} /* cqr_bkw4x1s_getScratchSize() */

#else
DISCARD_FUN(void, cqr_bkw4x1s,  (void* pScr, complex_fract16* restrict B, const complex_fract16* restrict R, int L))
size_t cqr_bkw4x1s_getScratchSize (int N, int P, int L) { return (0); }
#endif
