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


/* Biquad complex block IIR processing function, Low Noise fixed-point implementation. */
void bqciir_dp_proc( int16_t * restrict r,
                     int16_t * restrict sect, // 10*M
               const int16_t *          x,
                     int16_t * restrict scr,  // 4*N
               const int16_t *          coef,
                     int16_t            gain,
                     int N, int M )
{
    xb_vecNx16 G,  A1, A2;
    xb_vecNx16 B0, B1, B2;

    xb_vecNx16 D0_lo, D0_ho, D0_l, D0_h, D1_l, D1_h, D2_l, D2_h, Delay_l, Delay_h;

    xb_vecNx16 X, X_l, X_h;
    xb_vecNx16 Y, Y_l, Y_h;

    xb_vecNx16 Gain;

    xb_vecNx40 Acc_fb, Acc_fb_l, Acc_fb_h;
    xb_vecNx40 Acc_ff, Acc_ff_l, Acc_ff_h;

    xb_vecNx40 z;

    const xb_vecNx16 *COEF = (const xb_vecNx16 *)coef;
          xb_vecNx16 *SECT = (      xb_vecNx16 *)sect;

    void * restrict SCR;

    vsaN rnd16 = BBE_MOVVSA32(16);
    vsaN rnd8  = BBE_MOVVSA32(8);

    int n, m;
    NASSERT_ALIGN32(scr);

    Gain = BBE_MOVVA16(gain); 

    if ( M == 8 )
    {
        BBE_LVNX16_IP(G,  COEF, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(B0, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(B1, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(B2, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(A1, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(A2, COEF, 2*BBE_SIMD_WIDTH);

        D0_l    = BBE_LVNX16_I(SECT, 0*2*BBE_SIMD_WIDTH);
        D0_h    = BBE_LVNX16_I(SECT, 1*2*BBE_SIMD_WIDTH);
        D1_l    = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
        D1_h    = BBE_LVNX16_I(SECT, 3*2*BBE_SIMD_WIDTH);
        Delay_l = BBE_LVNX16_I(SECT, 4*2*BBE_SIMD_WIDTH);
        Delay_h = BBE_LVNX16_I(SECT, 5*2*BBE_SIMD_WIDTH);

        for ( n=0; n<N; n++ )
        {
            BBE_LPNX16_IP(X, x, 4);

            // Q29 <- Q15 + 14
            X_l = BBE_SLLINX16(X, 14);
            X_h = BBE_SRAINX16(X,  2);

            Delay_l = BBE_SELNX16I(Delay_l, X_l, BBE_SELI_PACK_2);
            Delay_h = BBE_SELNX16I(Delay_h, X_h, BBE_SELI_PACK_2);
            
            // Q43 <- Q29*Q14
            Acc_fb_l = BBE_MULUSNX16(D0_l, A2);
            Acc_fb_h = BBE_MULNX16(D0_h, A2);

            // Q43 <- Q43 - Q29*Q14
            BBE_MULUSANX16(Acc_fb_l, D1_l, A1);
            BBE_MULUSANX16(Acc_fb_l, Delay_l, G);
            BBE_MULANX16(Acc_fb_h, D1_h, A1);
            BBE_MULANX16(Acc_fb_h, Delay_h, G);

            Acc_fb_l = BBE_SRAINX40(Acc_fb_l, 14);
            Acc_fb_h = BBE_SLLINX40(Acc_fb_h,  2);

            Acc_fb = BBE_ADDNX40(Acc_fb_h, Acc_fb_l);

            D2_l = BBE_PACKLNX40(Acc_fb); D2_h = BBE_PACKVNX40(Acc_fb, rnd16);

            // Q43 <- Q29*Q14
            Acc_ff_l = BBE_MULUSNX16(D0_l, B2);
            Acc_ff_h = BBE_MULNX16(D0_h, B2);

            // Q43 <- Q43 - Q29*Q14
            BBE_MULUSANX16(Acc_ff_l, D1_l, B1);
            BBE_MULUSANX16(Acc_ff_l, D2_l, B0);
            BBE_MULANX16(Acc_ff_h, D1_h, B1);
            BBE_MULANX16(Acc_ff_h, D2_h, B0);

            Acc_ff_l = BBE_SRAINX40(Acc_ff_l, 14);
            Acc_ff_h = BBE_SLLINX40(Acc_ff_h,  2);

            Acc_ff = BBE_ADDNX40(Acc_ff_h, Acc_ff_l);

            D0_l = D1_l; D1_l = D2_l; 
            D0_h = D1_h; D1_h = D2_h;
            
            Delay_l = BBE_PACKLNX40(Acc_ff); Delay_h = BBE_PACKVNX40(Acc_ff, rnd16);

            // Q15 <- ( Q43 - 14 ) - 14 w/ rounding
            Acc_ff = BBE_ADDNX40(Acc_ff, Acc_ff);
            Y = BBE_PACKQNX40(Acc_ff);

            // Q15 <- Q15*Q8 - 8
            z = BBE_MULNX16(Y, Gain);
            Y = BBE_PACKVNX40(z, rnd8);
            Y = BBE_REPNX16C(Y, 7);

            BBE_SPNX16_IP(Y, r, 4);
        }

        BBE_SVNX16_IP(D0_l,    SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D0_h,    SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1_l,    SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1_h,    SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(Delay_l, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(Delay_h, SECT, 2*BBE_SIMD_WIDTH);

        return;
    }

    SCR = (void *)scr; 

    BBE_LVNX16_IP(G,  COEF, 2*BBE_SIMD_WIDTH); 
    BBE_LVNX16_IP(B0, COEF, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(B1, COEF, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(B2, COEF, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(A1, COEF, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(A2, COEF, 2*BBE_SIMD_WIDTH);

    D0_lo   = BBE_LVNX16_I(SECT, 0*2*BBE_SIMD_WIDTH);
    D0_ho   = BBE_LVNX16_I(SECT, 1*2*BBE_SIMD_WIDTH);
    D0_l    = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
    D0_h    = BBE_LVNX16_I(SECT, 3*2*BBE_SIMD_WIDTH);
    D1_l    = BBE_LVNX16_I(SECT, 4*2*BBE_SIMD_WIDTH);
    D1_h    = BBE_LVNX16_I(SECT, 5*2*BBE_SIMD_WIDTH);
    Delay_l = BBE_LVNX16_I(SECT, 6*2*BBE_SIMD_WIDTH);
    Delay_h = BBE_LVNX16_I(SECT, 7*2*BBE_SIMD_WIDTH);

    for ( n=0; n<N; n++ )
    {
        BBE_LPNX16_IP(X, x, 4);

        // Q29 <- Q15 + 14
        X_l = BBE_SLLINX16(X, 14);
        X_h = BBE_SRAINX16(X,  2);

        Delay_l = BBE_SELNX16I(Delay_l, X_l, BBE_SELI_PACK_2);
        Delay_h = BBE_SELNX16I(Delay_h, X_h, BBE_SELI_PACK_2);
        
        // Q43 <- Q29*Q14
        Acc_fb_l = BBE_MULUSNX16(D0_l, A2);
        Acc_fb_h = BBE_MULNX16(D0_h, A2);
        Acc_ff_l = BBE_MULUSNX16(D0_lo, B2);
        Acc_ff_h = BBE_MULNX16(D0_ho, B2);

        // Q43 <- Q43 - Q29*Q14
        BBE_MULUSANX16(Acc_fb_l, D1_l, A1);
        BBE_MULUSANX16(Acc_fb_l, Delay_l, G);
        BBE_MULANX16(Acc_fb_h, D1_h, A1);
        BBE_MULANX16(Acc_fb_h, Delay_h, G);
        BBE_MULUSANX16(Acc_ff_l, D0_l, B1);
        BBE_MULANX16(Acc_ff_h, D0_h, B1);
        
        Acc_fb_l = BBE_SRAINX40(Acc_fb_l, 14);
        Acc_fb_h = BBE_SLLINX40(Acc_fb_h,  2);

        Acc_fb = BBE_ADDNX40(Acc_fb_h, Acc_fb_l);

        D2_l = BBE_PACKLNX40(Acc_fb); D2_h = BBE_PACKVNX40(Acc_fb, rnd16);

        // Q43 <- Q43 - Q29*Q14
        BBE_MULUSANX16(Acc_ff_l, D1_l, B0);
        BBE_MULANX16(Acc_ff_h, D1_h, B0);

        Acc_ff_l = BBE_SRAINX40(Acc_ff_l, 14);
        Acc_ff_h = BBE_SLLINX40(Acc_ff_h,  2);

        Acc_ff = BBE_ADDNX40(Acc_ff_h, Acc_ff_l);

        D0_lo = D0_l; D0_l = D1_l; D1_l = D2_l; 
        D0_ho = D0_h; D0_h = D1_h; D1_h = D2_h;
        
        Delay_l = BBE_PACKLNX40(Acc_ff); Delay_h = BBE_PACKVNX40(Acc_ff, rnd16);

        Y_l = BBE_REPNX16C(Delay_l, 7);
        Y_h = BBE_REPNX16C(Delay_h, 7);

        BBE_SPNX16_IP(Y_l, SCR, 4);
        BBE_SPNX16_IP(Y_h, SCR, 4);
    }

    BBE_SVNX16_IP(D0_lo,   SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D0_ho,   SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D0_l,    SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D0_h,    SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D1_l,    SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D1_h,    SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Delay_l, SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Delay_h, SECT, 2*BBE_SIMD_WIDTH);

    for ( m=1; m<(M>>3) - 1; m++ )
    {
        SCR = (void *)scr; 

        BBE_LVNX16_IP(G,  COEF, 2*BBE_SIMD_WIDTH); 
        BBE_LVNX16_IP(B0, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(B1, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(B2, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(A1, COEF, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(A2, COEF, 2*BBE_SIMD_WIDTH);

        D0_lo   = BBE_LVNX16_I(SECT, 0*2*BBE_SIMD_WIDTH);
        D0_ho   = BBE_LVNX16_I(SECT, 1*2*BBE_SIMD_WIDTH);
        D0_l    = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
        D0_h    = BBE_LVNX16_I(SECT, 3*2*BBE_SIMD_WIDTH);
        D1_l    = BBE_LVNX16_I(SECT, 4*2*BBE_SIMD_WIDTH);
        D1_h    = BBE_LVNX16_I(SECT, 5*2*BBE_SIMD_WIDTH);
        Delay_l = BBE_LVNX16_I(SECT, 6*2*BBE_SIMD_WIDTH);
        Delay_h = BBE_LVNX16_I(SECT, 7*2*BBE_SIMD_WIDTH);

        for ( n=0; n<N; n++ )
        {
            X_l = BBE_LPNX16_I(SCR, 0);
            X_h = BBE_LPNX16_I(SCR, 4);

            Delay_l = BBE_SELNX16I(Delay_l, X_l, BBE_SELI_PACK_2);
            Delay_h = BBE_SELNX16I(Delay_h, X_h, BBE_SELI_PACK_2);
            
            // Q43 <- Q29*Q14
            Acc_fb_l = BBE_MULUSNX16(D0_l, A2);
            Acc_fb_h = BBE_MULNX16(D0_h, A2);
            Acc_ff_l = BBE_MULUSNX16(D0_lo, B2);
            Acc_ff_h = BBE_MULNX16(D0_ho, B2);

            // Q43 <- Q43 - Q29*Q14
            BBE_MULUSANX16(Acc_fb_l, D1_l, A1);
            BBE_MULUSANX16(Acc_fb_l, Delay_l, G);
            BBE_MULANX16(Acc_fb_h, D1_h, A1);
            BBE_MULANX16(Acc_fb_h, Delay_h, G);
            BBE_MULUSANX16(Acc_ff_l, D0_l, B1);
            BBE_MULANX16(Acc_ff_h, D0_h, B1);
            
            Acc_fb_l = BBE_SRAINX40(Acc_fb_l, 14);
            Acc_fb_h = BBE_SLLINX40(Acc_fb_h,  2);

            Acc_fb = BBE_ADDNX40(Acc_fb_h, Acc_fb_l);

            D2_l = BBE_PACKLNX40(Acc_fb); D2_h = BBE_PACKVNX40(Acc_fb, rnd16);

            // Q43 <- Q43 - Q29*Q14
            BBE_MULUSANX16(Acc_ff_l, D1_l, B0);
            BBE_MULANX16(Acc_ff_h, D1_h, B0);

            Acc_ff_l = BBE_SRAINX40(Acc_ff_l, 14);
            Acc_ff_h = BBE_SLLINX40(Acc_ff_h,  2);

            Acc_ff = BBE_ADDNX40(Acc_ff_h, Acc_ff_l);

            D0_lo = D0_l; D0_l = D1_l; D1_l = D2_l; 
            D0_ho = D0_h; D0_h = D1_h; D1_h = D2_h;
            
            Delay_l = BBE_PACKLNX40(Acc_ff); Delay_h = BBE_PACKVNX40(Acc_ff, rnd16);

            Y_l = BBE_REPNX16C(Delay_l, 7);
            Y_h = BBE_REPNX16C(Delay_h, 7);

            BBE_SPNX16_IP(Y_l, SCR, 4);
            BBE_SPNX16_IP(Y_h, SCR, 4);
        }

        BBE_SVNX16_IP(D0_lo,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D0_ho,   SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D0_l,    SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D0_h,    SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1_l,    SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1_h,    SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(Delay_l, SECT, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(Delay_h, SECT, 2*BBE_SIMD_WIDTH);
    } 

    SCR = (void *)scr; 

    BBE_LVNX16_IP(G,  COEF, 2*BBE_SIMD_WIDTH); 
    BBE_LVNX16_IP(B0, COEF, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(B1, COEF, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(B2, COEF, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(A1, COEF, 2*BBE_SIMD_WIDTH);
    BBE_LVNX16_IP(A2, COEF, 2*BBE_SIMD_WIDTH);

    D0_l    = BBE_LVNX16_I(SECT, 0*2*BBE_SIMD_WIDTH);
    D0_h    = BBE_LVNX16_I(SECT, 1*2*BBE_SIMD_WIDTH);
    D1_l    = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
    D1_h    = BBE_LVNX16_I(SECT, 3*2*BBE_SIMD_WIDTH);
    Delay_l = BBE_LVNX16_I(SECT, 4*2*BBE_SIMD_WIDTH);
    Delay_h = BBE_LVNX16_I(SECT, 5*2*BBE_SIMD_WIDTH);

    for ( n=0; n<N; n++ )
    {
        BBE_LPNX16_IP(X_l, SCR, 4);
        BBE_LPNX16_IP(X_h, SCR, 4);

        Delay_l = BBE_SELNX16I(Delay_l, X_l, BBE_SELI_PACK_2);
        Delay_h = BBE_SELNX16I(Delay_h, X_h, BBE_SELI_PACK_2);
        
        // Q43 <- Q29*Q14
        Acc_fb_l = BBE_MULUSNX16(D0_l, A2);
        Acc_fb_h = BBE_MULNX16(D0_h, A2);

        // Q43 <- Q43 - Q29*Q14
        BBE_MULUSANX16(Acc_fb_l, D1_l, A1);
        BBE_MULUSANX16(Acc_fb_l, Delay_l, G);
        BBE_MULANX16(Acc_fb_h, D1_h, A1);
        BBE_MULANX16(Acc_fb_h, Delay_h, G);

        Acc_fb_l = BBE_SRAINX40(Acc_fb_l, 14);
        Acc_fb_h = BBE_SLLINX40(Acc_fb_h,  2);

        Acc_fb = BBE_ADDNX40(Acc_fb_h, Acc_fb_l);

        D2_l = BBE_PACKLNX40(Acc_fb); D2_h = BBE_PACKVNX40(Acc_fb, rnd16);

        // Q43 <- Q29*Q14
        Acc_ff_l = BBE_MULUSNX16(D0_l, B2);
        Acc_ff_h = BBE_MULNX16(D0_h, B2);

        // Q43 <- Q43 - Q29*Q14
        BBE_MULUSANX16(Acc_ff_l, D1_l, B1);
        BBE_MULUSANX16(Acc_ff_l, D2_l, B0);
        BBE_MULANX16(Acc_ff_h, D1_h, B1);
        BBE_MULANX16(Acc_ff_h, D2_h, B0);

        Acc_ff_l = BBE_SRAINX40(Acc_ff_l, 14);
        Acc_ff_h = BBE_SLLINX40(Acc_ff_h,  2);

        Acc_ff = BBE_ADDNX40(Acc_ff_h, Acc_ff_l);

        D0_l = D1_l; D1_l = D2_l; 
        D0_h = D1_h; D1_h = D2_h;
        
        Delay_l = BBE_PACKLNX40(Acc_ff); Delay_h = BBE_PACKVNX40(Acc_ff, rnd16);

        // Q15 <- ( Q43 - 14 ) - 14 w/ rounding
        Acc_ff = BBE_ADDNX40(Acc_ff, Acc_ff);
        Y = BBE_PACKQNX40(Acc_ff);

        // Q15 <- Q15*Q8 - 8
        z = BBE_MULNX16(Y, Gain);
        Y = BBE_PACKVNX40(z, rnd8);
        Y = BBE_REPNX16C(Y, 7);

        BBE_SPNX16_IP(Y, r, 4);
    }

    BBE_SVNX16_IP(D0_l,    SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D0_h,    SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D1_l,    SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(D1_h,    SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Delay_l, SECT, 2*BBE_SIMD_WIDTH);
    BBE_SVNX16_IP(Delay_h, SECT, 2*BBE_SIMD_WIDTH);
} // bqciir_dp_proc()
