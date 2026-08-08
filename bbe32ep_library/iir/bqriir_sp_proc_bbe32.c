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
    Biquad real block IIR
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
#include "bqriir_common.h"

/* Biquad real block IIR processing function, Fast fixed-point implementation. */
void bqriir_sp_proc( int16_t * restrict r,
                     int16_t * restrict sect, // 2*M
               const int16_t *          x,
               const int16_t *          coef,
                     int16_t            gain,
                     int N, int M )
{
    xb_vecNx16 GB0, AB1, AB2;
    xb_vecNx16 D0,  D1,  D2;

    xb_vecNx16 X, Y, Gain;

    xb_vecNx40 Acc_fb_ff, z;

    static const union {int16_t i[16]; _vselN s;} ALIGN(32) Sel3 = {{16, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7}};

    const xb_vecNx16 *COEF = (const xb_vecNx16 *)coef;
          xb_vecNx16 *SECT = (      xb_vecNx16 *)sect;

    vsaN rnd14 = BBE_MOVVSA32(14);
    vsaN rnd8  = BBE_MOVVSA32(8);

    void * restrict R;

    int n;

    Gain = BBE_MOVVA16(gain); 

    if ( M == 8 )
    {
        BBE_LVNX16_IP(GB0, COEF, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(AB1, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(AB2, COEF, 2*BBE_SIMD_WIDTH);

        D0 = BBE_LVNX16_I(SECT, 0*2*BBE_SIMD_WIDTH);
        D1 = BBE_LVNX16_I(SECT, 1*2*BBE_SIMD_WIDTH);
        D2 = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);

        R = (void *)r;
        
        for ( n=0; n<N; n++ )
        {
            BBE_LSNX16_IP(X, x, 2);
            D2 = BBE_SELNX16(X, D2, _S(Sel3.s));

            // Q29 <- Q14*Q15
            Acc_fb_ff = BBE_MULRNX16(AB2, D0, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16(Acc_fb_ff, AB1, D1);
            BBE_MULANX16(Acc_fb_ff, GB0, D2);

            // Q15 <- Q29 - 14 w/ rounding
            D2 = BBE_PACKVNX40(Acc_fb_ff, rnd14);

            D0 = D1; D1 = BBE_SELNX16I(D1, D2, BBE_SELI_PACK_8); Y = BBE_REPNX16(D2, 15);

            // Q15 <- Q15*Q8 - 8
            z = BBE_MULNX16(Y, Gain);
            Y = BBE_PACKVNX40(z, rnd8);

            BBE_SSNX16_IP(Y, R, 2);
        }

        BBE_SVNX16_IP(D0, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D2, SECT, 2*BBE_SIMD_WIDTH);
    }
    else if ( M == 16 )
    {
        xb_vecNx16 GB0_1, AB1_1, AB2_1;
        xb_vecNx16 D0_1,  D1_1,  D2_1;

        BBE_LVNX16_IP(GB0,   COEF, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(AB1,   COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(AB2,   COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(GB0_1, COEF, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(AB1_1, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(AB2_1, COEF, 2*BBE_SIMD_WIDTH);

        D0   = BBE_LVNX16_I(SECT, 0*2*BBE_SIMD_WIDTH);
        D1   = BBE_LVNX16_I(SECT, 1*2*BBE_SIMD_WIDTH);
        D2   = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
        D0_1 = BBE_LVNX16_I(SECT, 3*2*BBE_SIMD_WIDTH);
        D1_1 = BBE_LVNX16_I(SECT, 4*2*BBE_SIMD_WIDTH);
        D2_1 = BBE_LVNX16_I(SECT, 5*2*BBE_SIMD_WIDTH);

        R = (void *)r;
        
        for ( n=0; n<N; n++ )
        {
            BBE_LSNX16_IP(X, x, 2);
            D2 = BBE_SELNX16(X, D2, _S(Sel3.s));

            // Q29 <- Q14*Q15
            Acc_fb_ff = BBE_MULRNX16(AB2, D0, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16(Acc_fb_ff, AB1, D1);
            BBE_MULANX16(Acc_fb_ff, GB0, D2);

            // Q15 <- Q29 - 14 w/ rounding
            D2 = BBE_PACKVNX40(Acc_fb_ff, rnd14);

            D0 = D1; D1 = BBE_SELNX16I(D1, D2, BBE_SELI_PACK_8); Y = BBE_REPNX16(D2, 15);

            D2_1 = BBE_SELNX16(Y, D2_1, _S(Sel3.s));

            // Q29 <- Q14*Q15
            Acc_fb_ff = BBE_MULRNX16(AB2_1, D0_1, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16(Acc_fb_ff, AB1_1, D1_1);
            BBE_MULANX16(Acc_fb_ff, GB0_1, D2_1);

            // Q15 <- Q29 - 14 w/ rounding
            D2_1 = BBE_PACKVNX40(Acc_fb_ff, rnd14);

            D0_1 = D1_1; D1_1 = BBE_SELNX16I(D1_1, D2_1, BBE_SELI_PACK_8); Y = BBE_REPNX16(D2_1, 15);

            // Q15 <- Q15*Q8 - 8
            z = BBE_MULNX16(Y, Gain);
            Y = BBE_PACKVNX40(z, rnd8);

            BBE_SSNX16_IP(Y, R, 2);
        }

        BBE_SVNX16_IP(D0,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D2,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D0_1, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1_1, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D2_1, SECT, 2*BBE_SIMD_WIDTH);
    }
    else if ( M == 24 )
    {
        xb_vecNx16 GB0_1, AB1_1, AB2_1;
        xb_vecNx16 D0_1,  D1_1,  D2_1;

        BBE_LVNX16_IP(GB0,   COEF, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(AB1,   COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(AB2,   COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(GB0_1, COEF, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(AB1_1, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(AB2_1, COEF, 2*BBE_SIMD_WIDTH);

        D0   = BBE_LVNX16_I(SECT, 0*2*BBE_SIMD_WIDTH);
        D1   = BBE_LVNX16_I(SECT, 1*2*BBE_SIMD_WIDTH);
        D2   = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
        D0_1 = BBE_LVNX16_I(SECT, 3*2*BBE_SIMD_WIDTH);
        D1_1 = BBE_LVNX16_I(SECT, 4*2*BBE_SIMD_WIDTH);
        D2_1 = BBE_LVNX16_I(SECT, 5*2*BBE_SIMD_WIDTH);

        R = (void *)r;
        
        for ( n=0; n<N; n++ )
        {
            BBE_LSNX16_IP(X, x, 2);
            D2 = BBE_SELNX16(X, D2, _S(Sel3.s));

            // Q29 <- Q14*Q15
            Acc_fb_ff = BBE_MULRNX16(AB2, D0, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16(Acc_fb_ff, AB1, D1);
            BBE_MULANX16(Acc_fb_ff, GB0, D2);

            // Q15 <- Q29 - 14 w/ rounding
            D2 = BBE_PACKVNX40(Acc_fb_ff, rnd14);

            D0 = D1; D1 = BBE_SELNX16I(D1, D2, BBE_SELI_PACK_8); Y = BBE_REPNX16(D2, 15);

            D2_1 = BBE_SELNX16(Y, D2_1, _S(Sel3.s));

            // Q29 <- Q14*Q15
            Acc_fb_ff = BBE_MULRNX16(AB2_1, D0_1, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16(Acc_fb_ff, AB1_1, D1_1);
            BBE_MULANX16(Acc_fb_ff, GB0_1, D2_1);

            // Q15 <- Q29 - 14 w/ rounding
            D2_1 = BBE_PACKVNX40(Acc_fb_ff, rnd14);

            D0_1 = D1_1; D1_1 = BBE_SELNX16I(D1_1, D2_1, BBE_SELI_PACK_8); Y = BBE_REPNX16(D2_1, 15);

            BBE_SSNX16_IP(Y, R, 2);
        }

        x = r; // Use output of the last stage further

        BBE_SVNX16_IP(D0,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D2,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D0_1, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1_1, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D2_1, SECT, 2*BBE_SIMD_WIDTH);

        BBE_LVNX16_IP(GB0, COEF, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(AB1, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(AB2, COEF, 2*BBE_SIMD_WIDTH);

        D0 = BBE_LVNX16_I(SECT, 0*2*BBE_SIMD_WIDTH);
        D1 = BBE_LVNX16_I(SECT, 1*2*BBE_SIMD_WIDTH);
        D2 = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);

        R = (void *)r;
        
        for ( n=0; n<N; n++ )
        {
            BBE_LSNX16_IP(X, x, 2);
            D2 = BBE_SELNX16(X, D2, _S(Sel3.s));

            // Q29 <- Q14*Q15
            Acc_fb_ff = BBE_MULRNX16(AB2, D0, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16(Acc_fb_ff, AB1, D1);
            BBE_MULANX16(Acc_fb_ff, GB0, D2);

            // Q15 <- Q29 - 14 w/ rounding
            D2 = BBE_PACKVNX40(Acc_fb_ff, rnd14);

            D0 = D1; D1 = BBE_SELNX16I(D1, D2, BBE_SELI_PACK_8); Y = BBE_REPNX16(D2, 15);

            // Q15 <- Q15*Q8 - 8
            z = BBE_MULNX16(Y, Gain);
            Y = BBE_PACKVNX40(z, rnd8);

            BBE_SSNX16_IP(Y, R, 2);
        }

        BBE_SVNX16_IP(D0, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D2, SECT, 2*BBE_SIMD_WIDTH);
    }
    else 
    {
        xb_vecNx16 GB0_1, AB1_1, AB2_1;
        xb_vecNx16 D0_1,  D1_1,  D2_1;

        BBE_LVNX16_IP(GB0,   COEF, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(AB1,   COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(AB2,   COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(GB0_1, COEF, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(AB1_1, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(AB2_1, COEF, 2*BBE_SIMD_WIDTH);

        D0   = BBE_LVNX16_I(SECT, 0*2*BBE_SIMD_WIDTH);
        D1   = BBE_LVNX16_I(SECT, 1*2*BBE_SIMD_WIDTH);
        D2   = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
        D0_1 = BBE_LVNX16_I(SECT, 3*2*BBE_SIMD_WIDTH);
        D1_1 = BBE_LVNX16_I(SECT, 4*2*BBE_SIMD_WIDTH);
        D2_1 = BBE_LVNX16_I(SECT, 5*2*BBE_SIMD_WIDTH);

        R = (void *)r;
        
        for ( n=0; n<N; n++ )
        {
            BBE_LSNX16_IP(X, x, 2);
            D2 = BBE_SELNX16(X, D2, _S(Sel3.s));

            // Q29 <- Q14*Q15
            Acc_fb_ff = BBE_MULRNX16(AB2, D0, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16(Acc_fb_ff, AB1, D1);
            BBE_MULANX16(Acc_fb_ff, GB0, D2);

            // Q15 <- Q29 - 14 w/ rounding
            D2 = BBE_PACKVNX40(Acc_fb_ff, rnd14);

            D0 = D1; D1 = BBE_SELNX16I(D1, D2, BBE_SELI_PACK_8); Y = BBE_REPNX16(D2, 15);

            D2_1 = BBE_SELNX16(Y, D2_1, _S(Sel3.s));

            // Q29 <- Q14*Q15
            Acc_fb_ff = BBE_MULRNX16(AB2_1, D0_1, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16(Acc_fb_ff, AB1_1, D1_1);
            BBE_MULANX16(Acc_fb_ff, GB0_1, D2_1);

            // Q15 <- Q29 - 14 w/ rounding
            D2_1 = BBE_PACKVNX40(Acc_fb_ff, rnd14);

            D0_1 = D1_1; D1_1 = BBE_SELNX16I(D1_1, D2_1, BBE_SELI_PACK_8); Y = BBE_REPNX16(D2_1, 15);

            BBE_SSNX16_IP(Y, R, 2);
        }

        x = r; // Use output of the last stage further

        BBE_SVNX16_IP(D0,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D2,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D0_1, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1_1, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D2_1, SECT, 2*BBE_SIMD_WIDTH);

        BBE_LVNX16_IP(GB0,   COEF, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(AB1,   COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(AB2,   COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(GB0_1, COEF, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(AB1_1, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(AB2_1, COEF, 2*BBE_SIMD_WIDTH);

        D0   = BBE_LVNX16_I(SECT, 0*2*BBE_SIMD_WIDTH);
        D1   = BBE_LVNX16_I(SECT, 1*2*BBE_SIMD_WIDTH);
        D2   = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
        D0_1 = BBE_LVNX16_I(SECT, 3*2*BBE_SIMD_WIDTH);
        D1_1 = BBE_LVNX16_I(SECT, 4*2*BBE_SIMD_WIDTH);
        D2_1 = BBE_LVNX16_I(SECT, 5*2*BBE_SIMD_WIDTH);

        R = (void *)r;
        
        for ( n=0; n<N; n++ )
        {
            BBE_LSNX16_IP(X, x, 2);
            D2 = BBE_SELNX16(X, D2, _S(Sel3.s));

            // Q29 <- Q14*Q15
            Acc_fb_ff = BBE_MULRNX16(AB2, D0, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16(Acc_fb_ff, AB1, D1);
            BBE_MULANX16(Acc_fb_ff, GB0, D2);

            // Q15 <- Q29 - 14 w/ rounding
            D2 = BBE_PACKVNX40(Acc_fb_ff, rnd14);

            D0 = D1; D1 = BBE_SELNX16I(D1, D2, BBE_SELI_PACK_8); Y = BBE_REPNX16(D2, 15);

            D2_1 = BBE_SELNX16(Y, D2_1, _S(Sel3.s));

            // Q29 <- Q14*Q15
            Acc_fb_ff = BBE_MULRNX16(AB2_1, D0_1, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16(Acc_fb_ff, AB1_1, D1_1);
            BBE_MULANX16(Acc_fb_ff, GB0_1, D2_1);

            // Q15 <- Q29 - 14 w/ rounding
            D2_1 = BBE_PACKVNX40(Acc_fb_ff, rnd14);

            D0_1 = D1_1; D1_1 = BBE_SELNX16I(D1_1, D2_1, BBE_SELI_PACK_8); Y = BBE_REPNX16(D2_1, 15);

            // Q15 <- Q15*Q8 - 8
            z = BBE_MULNX16(Y, Gain);
            Y = BBE_PACKVNX40(z, rnd8);

            BBE_SSNX16_IP(Y, R, 2);
        }

        BBE_SVNX16_IP(D0,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D2,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D0_1, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1_1, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D2_1, SECT, 2*BBE_SIMD_WIDTH);
    }
} // bqriir_sp_proc()
