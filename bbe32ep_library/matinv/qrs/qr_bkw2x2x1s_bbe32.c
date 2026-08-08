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

void  qr_bkw2x2x1s (void* pScr, int16_t* restrict B, const int16_t* restrict Q, const int16_t* restrict  R, int L)
{
    int l;

    xb_vecNx16 f0; 
    xb_vecNx16 * restrict pB00; 
    xb_vecNx16 * restrict pB10; 

    xb_vecNx16 * restrict  pQ00; 
    xb_vecNx16 * restrict  pQ01; 
    xb_vecNx16 * restrict  pQ10; 
    xb_vecNx16 * restrict  pQ11; 

    xb_vecNx16 * restrict  pR00; 
    xb_vecNx16 * restrict  pR01; 
    xb_vecNx16 * restrict  pR11; 

    NASSERT_ALIGN(pScr,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(B,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(Q,2*BBE_SIMD_WIDTH);
    NASSERT_ALIGN(R,2*BBE_SIMD_WIDTH);
    NASSERT(L%BBE_SIMD_WIDTH==0);
    (void)pScr;

    pQ00 = (xb_vecNx16*) (Q); 
    pQ01 = (xb_vecNx16*) (Q+L);
    pQ10 = (xb_vecNx16*) (Q+2*L); 
    pQ11 = (xb_vecNx16*) (Q+3*L); 

    /*
    Calculate QB = Q' * B; 

    Q = [q00, q01
        q10, q11 ]

    Q' = [ q00', q10'
           q01', q11' ]

    Q' * B =  [ q00', q10'      *  [b00, b01
                q01', q11' ]        b10, b11 ]
    */

    pB00 = (xb_vecNx16*) (B); 
    pB10 = (xb_vecNx16*) (B); 

    pR00 = (xb_vecNx16*) (R); 
    pR01 = (xb_vecNx16*) (R+L);
    pR11 = (xb_vecNx16*) (R+3*L); 

    {
        const xb_vecNx40 z0 = BBE_MOVWA32(0x2000000);
        xb_vecNx16 t,r00, r01, r11; 
        xb_vecNx16 b00, b10; 
        xb_vecNx16 q00, q01, q11, q10; 
        xb_vecNx40 acc0,w0,w1; 
        vsaN shft;
        for(l = 0; l<L; l+=BBE_SIMD_WIDTH)
        {
            BBE_LVNX16_IP(b00, pB00, (2*BBE_SIMD_WIDTH)); 
            b10=BBE_LVNX16_X ( pB00, 1*2*L-(2*BBE_SIMD_WIDTH)); 

            BBE_LVNX16_IP(q00, pQ00, sizeof(*pQ00)); 
            BBE_LVNX16_IP(q01, pQ01, sizeof(*pQ00)); 
            BBE_LVNX16_IP(q10, pQ10, sizeof(*pQ00)); 
            BBE_LVNX16_IP(q11, pQ11, sizeof(*pQ00)); 

            w0=BBE_MULNX16 (b00,q00);
            BBE_MULANX16(w0,b10,q10);
            w1=BBE_MULNX16 (b00,q01);
            BBE_MULANX16(w1,b10,q11);
            b00=BBE_PACKQNX40(w0);
            b10=BBE_PACKQNX40(w1);

            BBE_LVNX16_IP(r00, pR00, sizeof(*pR00) ); 
            BBE_LVNX16_IP(r01, pR01, sizeof(*pR01) ); 
            BBE_LVNX16_IP(r11, pR11, sizeof(*pR11) ); 


            t = r11;
            shft = BBE_NSANX16(t);
            shft = BBE_ADDSAVSN(-4,shft);
            r11 =  BBE_SLANX16(r11,shft);
            QUO32X32(f0,z0,r11);


            shft=BBE_SUBSAVSN(15,shft);
            acc0 = BBE_MULRNX16( f0, b10, shft); 
            b10 = BBE_PACKVNX40(acc0, shft); 

            w0=BBE_MULNX16(r01, b10); 
            w0=BBE_SRAINX40(w0,1);
            t=BBE_PACKPNX40(w0);
            b00=BBE_SRAINX16(b00,1);
            b00 = BBE_SUBNX16(b00, t);


            t=r00;
            shft = BBE_NSANX16(t);
            shft = BBE_ADDSAVSN(-4,shft);
            r00 =  BBE_SLANX16(r00,shft);
            QUO32X32(f0,z0,r00);


            shft=BBE_SUBSAVSN(14,shft);
            acc0 = BBE_MULRNX16( f0, b00, shft); 
            b00 = BBE_PACKVNX40(acc0, shft);

            BBE_SVNX16_IP(b00, pB10, (2*BBE_SIMD_WIDTH)); 
            BBE_SVNX16_X (b10, pB10, 1*2*L-(2*BBE_SIMD_WIDTH)); 
        }
    }
} /* qr_bkw2x2x1s() */
#else
DISCARD_FUN(void,qr_bkw2x2x1s,(void* pScr, int16_t* restrict B, const int16_t* restrict  Q, const int16_t* restrict R, int L))
#endif

size_t  qr_bkw2x2x1s_getScratchSize (int N, int P, int L)
{
    (void)N; (void)P; (void)L;
    return 0;
} /* qr_bkw2x2x1s_getScratchSize() */
