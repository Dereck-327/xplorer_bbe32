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
    Biquad complex block IIR
    C code optimized for BBE32
    Integrit, 2006-2016
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_iir.h"
/* Filter processing functions. */
#include "bqciir_common.h"

/* Biquad complex block IIR processing function, Fast fixed-point implementation. */
void bqciir_sp_proc( int16_t * restrict r,
                     int16_t * restrict sect, // 4*M
               const int16_t *          x,
               const int16_t *          coef,
                     int16_t            gain,
                     int N, int M )
{
    xb_vecNx16 G,  A1, A2;
    xb_vecNx16 B0, B1, B2;

    xb_vecNx16 D0o, D0, D1, D2, Delay, X, Gain, Y;

    const xb_vecNx16 *COEF = (const xb_vecNx16 *)coef;
          xb_vecNx16 *SECT = (      xb_vecNx16 *)sect;

    xb_vecNx40 Acc_fb, Acc_ff, z;

    vsaN rnd14 = BBE_MOVVSA32(14);
    vsaN rnd8  = BBE_MOVVSA32(8);

    void * restrict R;

    int n, m;

    Gain = BBE_MOVVA16(gain); 

    for ( m=0; m<(M>>3)-1; m++ )
    {
        BBE_LVNX16_IP(G,  COEF, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(B0, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(B1, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(B2, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(A1, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(A2, COEF, 2*BBE_SIMD_WIDTH);

        D0o   = BBE_LVNX16_I(SECT, 0*2*BBE_SIMD_WIDTH);
        D0    = BBE_LVNX16_I(SECT, 1*2*BBE_SIMD_WIDTH);
        D1    = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
        Delay = BBE_LVNX16_I(SECT, 3*2*BBE_SIMD_WIDTH);

        R = (void *)r;

        for ( n=0; n<N; n++ )
        {
            BBE_LPNX16_IP(X, x, 4);

            Delay = BBE_SELNX16I(Delay, X, BBE_SELI_PACK_2);

            // Q29 <- Q14*Q15
            Acc_fb = BBE_MULRNX16(A1, D1,  rnd14);
            Acc_ff = BBE_MULRNX16(B2, D0o, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16(Acc_fb, A2, D0); 
            BBE_MULANX16(Acc_ff, B1, D0);
            BBE_MULANX16(Acc_fb, Delay, G);
            BBE_MULANX16(Acc_ff, B0, D1);

            // Q15 <- Q29 - 14 w/ rounding
            D2    = BBE_PACKVNX40(Acc_fb, rnd14);
            Delay = BBE_PACKVNX40(Acc_ff, rnd14);

            D0o = D0; D0 = D1; D1 = D2; 

            Y = BBE_REPNX16C(Delay, 7);
            BBE_SPNX16_IP(Y, R, 4);
        }

        x = r; // Use output of the last stage further

        BBE_SVNX16_IP(D0o,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D0,    SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1,    SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(Delay, SECT, 2*BBE_SIMD_WIDTH);
    }

    BBE_LVNX16_IP(G,  COEF, 2*BBE_SIMD_WIDTH); 
    BBE_LVNX16_IP(B0, COEF, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(B1, COEF, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(B2, COEF, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(A1, COEF, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(A2, COEF, 2*BBE_SIMD_WIDTH);

    D0o   = BBE_LVNX16_I(SECT, 0*2*BBE_SIMD_WIDTH);
    D0    = BBE_LVNX16_I(SECT, 1*2*BBE_SIMD_WIDTH);
    D1    = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
    Delay = BBE_LVNX16_I(SECT, 3*2*BBE_SIMD_WIDTH);

    R = (void *)r;

    for ( n=0; n<N; n++ )
    {
        BBE_LPNX16_IP(X, x, 4);

        Delay = BBE_SELNX16I(Delay, X, BBE_SELI_PACK_2);

        // Q29 <- Q14*Q15
        Acc_fb = BBE_MULRNX16(A1, D1,  rnd14);
        Acc_ff = BBE_MULRNX16(B2, D0o, rnd14);

        // Q29 <- Q29 - Q14*Q15
        BBE_MULANX16(Acc_ff, B1, D0);
        BBE_MULANX16(Acc_fb, A2, D0); 
        BBE_MULANX16(Acc_fb, Delay, G);
        BBE_MULANX16(Acc_ff, B0, D1);

        // Q15 <- Q29 - 14 w/ rounding
        D2    = BBE_PACKVNX40(Acc_fb, rnd14);
        Delay = BBE_PACKVNX40(Acc_ff, rnd14);

        D0o = D0; D0 = D1; D1 = D2; 
        
        Y = BBE_REPNX16C(Delay, 7);

        // Q15 <- Q15*Q8 - 8
        z = BBE_MULNX16(Y, Gain);
        Y = BBE_PACKVNX40(z, rnd8);

        BBE_SPNX16_IP(Y, R, 4);
    }

    x = r; // Use output of the last stage further

    BBE_SVNX16_IP(D0o,   SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D0,    SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D1,    SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Delay, SECT, 2*BBE_SIMD_WIDTH);
} // bqciir_sp_proc()
