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
/* Common utility declarations. */
#include "qr_common.h"

#if HAVE_VSAMATH && 1

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

void  qr_bkw2x2x2s (void* pScr, int16_t* restrict B, const int16_t* restrict Q, const int16_t* restrict  R, int L)
{
    int l;
    xb_vecNx16 f0; 

    xb_vecNx16 * restrict pBr0; 
    xb_vecNx16 * restrict pBr1; 
    xb_vecNx16 * restrict pBw; 

    xb_vecNx16 * restrict pQ;

    xb_vecNx16 * restrict pR00; 

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Q,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT(L>0 && L%BBE_SIMD_WIDTH==0);
    (void)pScr;

    /*
    Calculate QB = Q' * B; 

    Q = [q00, q01
         q10, q11 ]

    Q' = [ q00', q10'
           q01', q11' ]

    Q' * B =  [ q00', q10'      *  [b00, b01
                q01', q11' ]         b10, b11 ]
    */


    pR00 = (xb_vecNx16*) (R); 
    pBw = (      xb_vecNx16*) (B); 
    pQ = (xb_vecNx16*) (Q); 
    pBr0 = (xb_vecNx16*) (B+2*L); 
    pBr1 = (xb_vecNx16*) (B+0*L);

    {
//        const xb_vecNx40 z0 = BBE_MOVWA32(0x2000000);
        xb_vecNx40 z0;
        int _0x2000000=0x2000000;
        xb_vecNx16 t,r00, r01, r11, zero=0; 
        xb_vecNx16 b00, b01, b11, b10; 
        xb_vecNx16 q00, q01, q11, q10; 
        xb_vecNx40 acc0,w0,w1; 
        vsaN shft;
        for(l = 0; l<L; l+=BBE_SIMD_WIDTH)
        {
            XT_MOVEQZ(_0x2000000,_0x2000000,_0x2000000);    // prevent z0 from caching
            z0= BBE_MOVWA32(_0x2000000);
            b11=BBE_LVNX16_X ( pBr0, 1*2*L); 
            BBE_LVNX16_IP(b10,pBr0, (2*BBE_SIMD_WIDTH)); 
            b01=BBE_LVNX16_X ( pBr1, 1*2*L);
            BBE_LVNX16_IP(b00,pBr1, (2*BBE_SIMD_WIDTH));
            q11 = BBE_LVNX16_X(pQ, 3*2*L);
            q10 = BBE_LVNX16_X(pQ, 2*2*L); 
            q01 = BBE_LVNX16_X(pQ, 1*2*L); 
            BBE_LVNX16_IP(q00, pQ, (2*BBE_SIMD_WIDTH)); 

            w0 = BBE_MULNX16(b00,q00);
            BBE_MULANX16(w0, b10,q10);
            w1 = BBE_MULNX16(b00,q01);
            BBE_MULANX16(w1, b10,q11);
            b00=BBE_PACKQNX40(w0);
            b10=BBE_PACKQNX40(w1);

            w0 =BBE_MULNX16(b01, q00); 
            BBE_MULANX16(w0,b11, q10); 
            w1 = BBE_MULNX16(b01, q01);
            BBE_MULANX16(w1, b11, q11);
            b01=BBE_PACKQNX40(w0);
            b11=BBE_PACKQNX40(w1);

            r01=BBE_LVNX16_X ( pR00, 1*2*L ); 
            r11=BBE_LVNX16_X ( pR00, 3*2*L ); 
            BBE_LVNX16_IP(r00, pR00, (2*BBE_SIMD_WIDTH) ); 

            t=r11;
            shft = BBE_NSANX16(t);
            shft = BBE_ADDSAVSN(-4,shft);
            r11 =  BBE_SLANX16(r11,shft);
            QUO32X32(f0,z0,r11);


            shft=BBE_SUBSAVSN(15,shft);
            acc0 = BBE_MULNX16( f0, b10); 
            b10 = BBE_PACKVNX40(acc0, shft); 
            acc0 = BBE_MULNX16( f0, b11); 
            b11 = BBE_PACKVNX40(acc0, shft); 
            shft=BBE_MOVVSA32(11);
            acc0=BBE_MULRNX16(r01, b10,shft); 
            t = BBE_PACKVNX40(acc0, shft); 
            b00 = BBE_ADDSR1RNX16(b00, zero);
            b00 = BBE_SUBNX16(b00, t);
            acc0=BBE_MULRNX16(r01, b11,shft); 
            t = BBE_PACKVNX40(acc0, shft); 
            b01 = BBE_ADDSR1RNX16(b01, zero);
            b01 = BBE_SUBNX16(b01, t);

            t=r00;
            shft = BBE_NSANX16(t);
            shft = BBE_ADDSAVSN(-4,shft);
            r11 =  BBE_SLANX16(r11,shft);
            r00 =  BBE_SLANX16(r00,shft);
            QUO32X32(f0,z0,r00);

            shft=BBE_SUBSAVSN(14,shft);

            acc0 = BBE_MULNX16( f0, b00); 
            b00 = BBE_PACKVNX40(acc0, shft);
            acc0 = BBE_MULNX16( f0, b01); 
            b01 = BBE_PACKVNX40(acc0, shft);

            BBE_SVNX16_X (b11, pBw, 3*2*L);
            BBE_SVNX16_X (b10, pBw, 2*2*L);
            BBE_SVNX16_X (b01, pBw, 1*2*L);
            BBE_SVNX16_IP(b00, pBw, (2*BBE_SIMD_WIDTH));
        }
    }
} /* qr_bkw2x2x2s() */
#else
DISCARD_FUN(void,qr_bkw2x2x2s,(void* pScr, int16_t* restrict B, const int16_t* restrict  Q, const int16_t* restrict R, int L))
#endif

size_t  qr_bkw2x2x2s_getScratchSize (int N, int P, int L)
{
    (void)N;(void)P;(void)L;
    return 0;
} /* qr_bkw2x2x2s_getScratchSize() */
