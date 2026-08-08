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
    Biquad real block IIR, streaming version
    C code optimized for BBE32
    IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_iir.h"
/* Filter processing functions. */
#include "bqriirs_common.h"

/* Biquad real block IIR processing function, Fast fixed-point implementation. */
void bqriirs_sp_proc( int16_t * restrict r,
                      int16_t * restrict sect, // 2*M
                const int16_t *          x,
                const int16_t *          coef,
                      int16_t            gain,
                      int N, int L, int M )
{
    unsigned int a12, b12;

    xb_vecNx16 G, B0, B1, B2, A1, A2;

    xb_vecNx16 D0, D1, D2;
    xb_vecNx16 D0o;
    xb_vecNx16 D0oo;

    xb_vecNx16 X0, X1, Y0, Y1, Gain;

    xb_vecNx40 Acc_fb0, Acc_fb1, Acc_ff0, Acc_ff1, z0, z1;

          xb_vecNx16 *SECT, *R;
    const xb_vecNx16 *_X, *SECT_ld;

    vsaN rnd14 = BBE_MOVVSA32(14);
    vsaN rnd8  = BBE_MOVVSA32(8);

    int n, m, l1;

    NASSERT( L%BBE_SIMD_WIDTH==0 );

    SECT    = (      xb_vecNx16 *)sect;
    SECT_ld = (const xb_vecNx16 *)sect;

    Gain = BBE_LSNX16_I(&gain, 0);
    Gain = BBE_REPNX16(Gain, 0);

    if ( N == 1 )
    {
        for ( m=0; m<M-1; m++ )
        {
            BBE_LSNX16_IP(G,  coef, 2); BBE_LSNX16_IP(B0, coef, 2); BBE_LSNX16_IP(B1, coef, 2);
            BBE_LSNX16_IP(B2, coef, 2); BBE_LSNX16_IP(A1, coef, 2); BBE_LSNX16_IP(A2, coef, 2);

            G  = BBE_REPNX16(G,  0); B0 = BBE_REPNX16(B0, 0); B1 = BBE_REPNX16(B1, 0);
            B2 = BBE_REPNX16(B2, 0); A1 = BBE_REPNX16(A1, 0); A2 = BBE_REPNX16(A2, 0);
                
            _X = (const xb_vecNx16 *)x;
            R  = (      xb_vecNx16 *)r;

#ifdef COMPILER_XTENSA 
    #pragma loop_count min=1
#endif
            for ( l1=0; l1<(L>>LOG2_BBE_SIMD_WIDTH); l1++)
            {
                BBE_LVNX16_IP(D0  , SECT_ld, 2*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(D1  , SECT_ld, 2*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(D0o , SECT_ld, 2*BBE_SIMD_WIDTH);
                BBE_LVNX16_IP(D0oo, SECT_ld, 2*BBE_SIMD_WIDTH);

                BBE_LVNX16_IP(X0, _X, 2*BBE_SIMD_WIDTH);

                // Q29 <- Q14*Q15
                Acc_fb0 = BBE_MULRNX16(A2, D0,   rnd14); 
                Acc_ff0 = BBE_MULRNX16(B2, D0oo, rnd14);

                // Q29 <- Q29 - Q14*Q15
                BBE_MULANX16(Acc_fb0, A1, D1);
                BBE_MULANX16(Acc_ff0, B1, D0o);
                BBE_MULANX16(Acc_fb0, G, X0);
                BBE_MULANX16(Acc_ff0, B0, D0);

                D0oo = D0o; D0o = D0; D0 = D1;  
                
                // Q15 <- Q29 - 14 w/ rounding
                D1 = BBE_PACKVNX40(Acc_fb0, rnd14);
                Y0 = BBE_PACKVNX40(Acc_ff0, rnd14);

                BBE_SVNX16_IP(Y0, R, 2*BBE_SIMD_WIDTH);

                BBE_SVNX16_IP(D0,   SECT, 2*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(D1,   SECT, 2*BBE_SIMD_WIDTH); 
                BBE_SVNX16_IP(D0o,  SECT, 2*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(D0oo, SECT, 2*BBE_SIMD_WIDTH); 
            }

            x = r; // Use output of the last stage further
        }

        BBE_LSNX16_IP(G,  coef, 2); BBE_LSNX16_IP(B0, coef, 2); BBE_LSNX16_IP(B1, coef, 2);
        BBE_LSNX16_IP(B2, coef, 2); BBE_LSNX16_IP(A1, coef, 2); BBE_LSNX16_IP(A2, coef, 2);

        G  = BBE_REPNX16(G,  0); B0 = BBE_REPNX16(B0, 0); B1 = BBE_REPNX16(B1, 0);
        B2 = BBE_REPNX16(B2, 0); A1 = BBE_REPNX16(A1, 0); A2 = BBE_REPNX16(A2, 0);

        _X = (const xb_vecNx16 *)x;
        R  = (      xb_vecNx16 *)r;
     
#ifdef COMPILER_XTENSA 
    #pragma loop_count min=1
#endif
        for ( l1=0; l1<(L>>LOG2_BBE_SIMD_WIDTH); l1++)
        {
            BBE_LVNX16_IP(D0  , SECT_ld, 2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(D1  , SECT_ld, 2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(D0o , SECT_ld, 2*BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(D0oo, SECT_ld, 2*BBE_SIMD_WIDTH);

            BBE_LVNX16_IP(X0, _X, 2*BBE_SIMD_WIDTH);

            // Q29 <- Q14*Q15
            Acc_fb0 = BBE_MULRNX16(A2, D0,   rnd14);
            Acc_ff0 = BBE_MULRNX16(B2, D0oo, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16(Acc_fb0, A1, D1);
            BBE_MULANX16(Acc_fb0, G, X0);
            BBE_MULANX16(Acc_ff0, B1, D0o);
            BBE_MULANX16(Acc_ff0, B0, D0);
            
            D0oo = D0o; D0o = D0; D0 = D1;  
            
            // Q15 <- Q29 - 14 w/ rounding
            D1 = BBE_PACKVNX40(Acc_fb0, rnd14);
            Y0 = BBE_PACKVNX40(Acc_ff0, rnd14);

            // Q15 <- Q15*Q8 - 8
            z0 = BBE_MULNX16(Y0, Gain);
            Y0 = BBE_PACKVNX40(z0, rnd8);

            BBE_SVNX16_IP(Y0, R, 2*BBE_SIMD_WIDTH);

            BBE_SVNX16_IP(D0,   SECT, 2*BBE_SIMD_WIDTH); 
            BBE_SVNX16_IP(D1,   SECT, 2*BBE_SIMD_WIDTH); 
            BBE_SVNX16_IP(D0o,  SECT, 2*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(D0oo, SECT, 2*BBE_SIMD_WIDTH);
        }

    }
    else if ( N&1 )
    {
        for ( m=0; m<M-1; m++ )
        {
            BBE_LSNX16_IP(G,  coef, 2); 
            BBE_LSNX16_IP(B0, coef, 2); b12 = BBE_L32X((const int32_t *)coef, 0); 
            BBE_LSNX16_IP(B1, coef, 2); 
            BBE_LSNX16_IP(B2, coef, 2); a12 = BBE_L32X((const int32_t *)coef, 0); 
            BBE_LSNX16_IP(A1, coef, 2); 
            BBE_LSNX16_IP(A2, coef, 2);

            G  = BBE_REPNX16(G,  0); 
            B0 = BBE_REPNX16(B0, 0); 

#ifdef COMPILER_XTENSA 
    #pragma loop_count min=1
#endif
            for ( l1=0; l1<(L>>LOG2_BBE_SIMD_WIDTH); l1++)
            {
                _X = (const xb_vecNx16 *)(x+BBE_SIMD_WIDTH*l1);
                R  = (      xb_vecNx16 *)(r+BBE_SIMD_WIDTH*l1);

                D0   = BBE_LVNX16_I(SECT, 2*0*BBE_SIMD_WIDTH); 
                D1   = BBE_LVNX16_I(SECT, 2*1*BBE_SIMD_WIDTH);
                D0o  = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
                D0oo = BBE_LVNX16_I(SECT, 2*3*BBE_SIMD_WIDTH);
             
#ifdef COMPILER_XTENSA 
    #pragma ymemory(_X)
    #pragma loop_count min=1
#endif
                for ( n=0; n<N-1; n+=2 )
                {                
                    BBE_LVNX16_XP(X0, _X, L*2);
                    BBE_LVNX16_XP(X1, _X, L*2);

                    // Q29 <- Q14*Q15
                    Acc_ff0 = BBE_MULRNX16(B0, D0, rnd14);
                    Acc_ff1 = BBE_MULRNX16(B0, D1, rnd14);
                    Acc_fb0 = BBE_MULRNX16(G,  X0, rnd14); 
                    Acc_fb1 = BBE_MULRNX16(G,  X1, rnd14); 

                    // Q29 <- Q29 - Q14*Q15
                    BBE_MULANX16PR(Acc_fb0, D0,   D1,  a12);
                    BBE_MULANX16PR(Acc_ff0, D0oo, D0o, b12);
                    BBE_MULANX16PR(Acc_ff1, D0o,  D0,  b12);

                    // Q15 <- Q29 - 14 w/ rounding
                    D2 = BBE_PACKVNX40(Acc_fb0, rnd14);
                    Y0 = BBE_PACKVNX40(Acc_ff0, rnd14);

                    // Q29 <- Q29 - Q14*Q15
                    BBE_MULANX16PR(Acc_fb1, D1, D2, a12);

                    D0oo = D0; D0o = D1; D0 = D2;
                    
                    // Q15 <- Q29 - 14 w/ rounding
                    D1 = BBE_PACKVNX40(Acc_fb1, rnd14);
                    Y1 = BBE_PACKVNX40(Acc_ff1, rnd14);

                    BBE_SVNX16_XP(Y0, R, 2*L);
                    BBE_SVNX16_XP(Y1, R, 2*L);
                }

                BBE_LVNX16_XP(X0, _X, L*2);

                // Q29 <- Q14*Q15
                Acc_fb0 = BBE_MULRNX16(G,  X0, rnd14); 
                Acc_ff0 = BBE_MULRNX16(B0, D0, rnd14);

                // Q29 <- Q29 - Q14*Q15
                BBE_MULANX16PR(Acc_fb0, D0,   D1,  a12);
                BBE_MULANX16PR(Acc_ff0, D0oo, D0o, b12);

                D0oo = D0o; D0o = D0; D0 = D1;  
                
                // Q15 <- Q29 - 14 w/ rounding
                D1 = BBE_PACKVNX40(Acc_fb0, rnd14);
                Y0 = BBE_PACKVNX40(Acc_ff0, rnd14);

                BBE_SVNX16_XP(Y0, R, 2*L);

                BBE_SVNX16_IP(D0,   SECT, 2*BBE_SIMD_WIDTH); 
                BBE_SVNX16_IP(D1,   SECT, 2*BBE_SIMD_WIDTH); 
                BBE_SVNX16_IP(D0o,  SECT, 2*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(D0oo, SECT, 2*BBE_SIMD_WIDTH);
            }

            x = r; // Use output of the last stage further
        }

        BBE_LSNX16_IP(G,  coef, 2); 
        BBE_LSNX16_IP(B0, coef, 2); b12 = BBE_L32X((const int32_t *)coef, 0); 
        BBE_LSNX16_IP(B1, coef, 2); 
        BBE_LSNX16_IP(B2, coef, 2); a12 = BBE_L32X((const int32_t *)coef, 0);
        BBE_LSNX16_IP(A1, coef, 2); 
        BBE_LSNX16_IP(A2, coef, 2);

        G  = BBE_REPNX16(G,  0); 
        B0 = BBE_REPNX16(B0, 0); 

#ifdef COMPILER_XTENSA 
    #pragma loop_count min=1
#endif
        for ( l1=0; l1<(L>>LOG2_BBE_SIMD_WIDTH); l1++)
        {
            _X = (const xb_vecNx16 *)(x+BBE_SIMD_WIDTH*l1);
            R  = (      xb_vecNx16 *)(r+BBE_SIMD_WIDTH*l1);
         
            D0   = BBE_LVNX16_I(SECT, 2*0*BBE_SIMD_WIDTH); 
            D1   = BBE_LVNX16_I(SECT, 2*1*BBE_SIMD_WIDTH);
            D0o  = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
            D0oo = BBE_LVNX16_I(SECT, 2*3*BBE_SIMD_WIDTH);

            BBE_LVNX16_XP(X0, _X, L*2);
            BBE_LVNX16_XP(X1, _X, L*2);

#ifdef COMPILER_XTENSA 
    #pragma ymemory(_X)
#endif
            for ( n=0; n<N-3; n+=2 )
            {                
                // Q29 <- Q14*Q15
                Acc_ff0 = BBE_MULRNX16(B0, D0, rnd14);
                Acc_ff1 = BBE_MULRNX16(B0, D1, rnd14);
                Acc_fb0 = BBE_MULRNX16(G,  X0, rnd14); 
                Acc_fb1 = BBE_MULRNX16(G, X1, rnd14); 

                // Q29 <- Q29 - Q14*Q15
                BBE_MULANX16PR(Acc_fb0, D0,   D1,  a12);
                BBE_MULANX16PR(Acc_ff0, D0oo, D0o, b12);
                BBE_MULANX16PR(Acc_ff1, D0o,  D0,  b12);

                // Q15 <- Q29 - 14 w/ rounding
                D2 = BBE_PACKVNX40(Acc_fb0, rnd14);
                Y0 = BBE_PACKVNX40(Acc_ff0, rnd14);

                // Q29 <- Q29 - Q14*Q15
                BBE_MULANX16PR(Acc_fb1, D1, D2, a12);
                
                D0oo = D0; D0o = D1; D0 = D2; 
                
                // Q15 <- Q29 - 14 w/ rounding
                D1 = BBE_PACKVNX40(Acc_fb1, rnd14);
                Y1 = BBE_PACKVNX40(Acc_ff1, rnd14);

                // Q15 <- Q15*Q8 - 8
                z0 = BBE_MULNX16(Y0, Gain);
                z1 = BBE_MULNX16(Y1, Gain);
                Y0 = BBE_PACKVNX40(z0, rnd8);
                Y1 = BBE_PACKVNX40(z1, rnd8);

                BBE_LVNX16_XP(X0, _X, L*2);
                BBE_LVNX16_XP(X1, _X, L*2);

                BBE_SVNX16_XP(Y0, R, 2*L);
                BBE_SVNX16_XP(Y1, R, 2*L);
            }

            // Q29 <- Q14*Q15
            Acc_ff0 = BBE_MULRNX16(B0, D0, rnd14);
            Acc_ff1 = BBE_MULRNX16(B0, D1, rnd14);
            Acc_fb0 = BBE_MULRNX16(G,  X0, rnd14); 
            Acc_fb1 = BBE_MULRNX16(G, X1, rnd14); 

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16PR(Acc_fb0, D0,   D1,  a12);
            BBE_MULANX16PR(Acc_ff0, D0oo, D0o, b12);
            BBE_MULANX16PR(Acc_ff1, D0o,  D0,  b12);

            // Q15 <- Q29 - 14 w/ rounding
            D2 = BBE_PACKVNX40(Acc_fb0, rnd14);
            Y0 = BBE_PACKVNX40(Acc_ff0, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16PR(Acc_fb1, D1, D2, a12);
            
            D0oo = D0; D0o = D1; D0 = D2; 
            
            // Q15 <- Q29 - 14 w/ rounding
            D1 = BBE_PACKVNX40(Acc_fb1, rnd14);
            Y1 = BBE_PACKVNX40(Acc_ff1, rnd14);

            // Q15 <- Q15*Q8 - 8
            z0 = BBE_MULNX16(Y0, Gain);
            z1 = BBE_MULNX16(Y1, Gain);
            Y0 = BBE_PACKVNX40(z0, rnd8);
            Y1 = BBE_PACKVNX40(z1, rnd8);

            BBE_SVNX16_XP(Y0, R, 2*L);
            BBE_SVNX16_XP(Y1, R, 2*L);

            BBE_LVNX16_XP(X0, _X, L*2);

            // Q29 <- Q14*Q15
            Acc_fb0 = BBE_MULRNX16(G,  X0, rnd14);
            Acc_ff0 = BBE_MULRNX16(B0, D0, rnd14);
                                         
            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16PR(Acc_fb0, D0,   D1,  a12);
            BBE_MULANX16PR(Acc_ff0, D0oo, D0o, b12);
            
            D0oo = D0o; D0o = D0; D0 = D1;  
            
            // Q15 <- Q29 - 14 w/ rounding
            D1 = BBE_PACKVNX40(Acc_fb0, rnd14);
            Y0 = BBE_PACKVNX40(Acc_ff0, rnd14);

            // Q15 <- Q15*Q8 - 8
            z0 = BBE_MULNX16(Y0, Gain);
            Y0 = BBE_PACKVNX40(z0, rnd8);

            BBE_SVNX16_XP(Y0, R, 2*L);

            BBE_SVNX16_IP(D0,   SECT, 2*BBE_SIMD_WIDTH); 
            BBE_SVNX16_IP(D1,   SECT, 2*BBE_SIMD_WIDTH); 
            BBE_SVNX16_IP(D0o,  SECT, 2*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(D0oo, SECT, 2*BBE_SIMD_WIDTH);
        }
    }
    else
    {
        for ( m=0; m<M-1; m++ )
        {
            BBE_LSNX16_IP(G,  coef, 2); 
            BBE_LSNX16_IP(B0, coef, 2); b12 = BBE_L32X((const int32_t *)coef, 0);
            BBE_LSNX16_IP(B1, coef, 2); 
            BBE_LSNX16_IP(B2, coef, 2); a12 = BBE_L32X((const int32_t *)coef, 0);
            BBE_LSNX16_IP(A1, coef, 2); 
            BBE_LSNX16_IP(A2, coef, 2);

            G  = BBE_REPNX16(G,  0); 
            B0 = BBE_REPNX16(B0, 0); 

#ifdef COMPILER_XTENSA 
    #pragma loop_count min=1
#endif
            for ( l1=0; l1<(L>>LOG2_BBE_SIMD_WIDTH); l1++)
            {
                _X = (const xb_vecNx16 *)(x+BBE_SIMD_WIDTH*l1);
                R  = (      xb_vecNx16 *)(r+BBE_SIMD_WIDTH*l1);

                D0   = BBE_LVNX16_I(SECT, 2*0*BBE_SIMD_WIDTH); 
                D1   = BBE_LVNX16_I(SECT, 2*1*BBE_SIMD_WIDTH);
                D0o  = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
                D0oo = BBE_LVNX16_I(SECT, 2*3*BBE_SIMD_WIDTH);

#ifdef COMPILER_XTENSA 
    #pragma ymemory(_X)
    #pragma loop_count min=1
#endif
                for ( n=0; n<N; n+=2 )
                {                
                    BBE_LVNX16_XP(X0, _X, L*2);
                    BBE_LVNX16_XP(X1, _X, L*2);

                    // Q29 <- Q14*Q15
                    Acc_fb0 = BBE_MULRNX16(G,  X0, rnd14);
                    Acc_ff0 = BBE_MULRNX16(B0, D0, rnd14);
                    Acc_fb1 = BBE_MULRNX16(G,  X1, rnd14);
                    Acc_ff1 = BBE_MULRNX16(B0, D1, rnd14);

                    // Q29 <- Q29 - Q14*Q15
                    BBE_MULANX16PR(Acc_fb0, D0,   D1,  a12);
                    BBE_MULANX16PR(Acc_ff0, D0oo, D0o, b12);
                    BBE_MULANX16PR(Acc_ff1, D0o,  D0,  b12);
                    
                    // Q15 <- Q29 - 14 w/ rounding
                    D2 = BBE_PACKVNX40(Acc_fb0, rnd14);
                    Y0 = BBE_PACKVNX40(Acc_ff0, rnd14);

                    // Q29 <- Q29 - Q14*Q15
                    BBE_MULANX16PR(Acc_fb1, D1, D2, a12);

                    D0oo = D0; D0o = D1; D0 = D2;
                    
                    // Q15 <- Q29 - 14 w/ rounding
                    D1 = BBE_PACKVNX40(Acc_fb1, rnd14);
                    Y1 = BBE_PACKVNX40(Acc_ff1, rnd14);

                    BBE_SVNX16_XP(Y0, R, 2*L);
                    BBE_SVNX16_XP(Y1, R, 2*L);
                }

                BBE_SVNX16_IP(D0,   SECT, 2*BBE_SIMD_WIDTH); 
                BBE_SVNX16_IP(D1,   SECT, 2*BBE_SIMD_WIDTH); 
                BBE_SVNX16_IP(D0o,  SECT, 2*BBE_SIMD_WIDTH);
                BBE_SVNX16_IP(D0oo, SECT, 2*BBE_SIMD_WIDTH);
            }

            x = r; // Use output of the last stage further
        }

        BBE_LSNX16_IP(G,  coef, 2); 
        BBE_LSNX16_IP(B0, coef, 2); b12 = BBE_L32X((const int32_t *)coef, 0);
        BBE_LSNX16_IP(B1, coef, 2); 
        BBE_LSNX16_IP(B2, coef, 2); a12 = BBE_L32X((const int32_t *)coef, 0);
        BBE_LSNX16_IP(A1, coef, 2); 
        BBE_LSNX16_IP(A2, coef, 2);

        G  = BBE_REPNX16(G,  0); 
        B0 = BBE_REPNX16(B0, 0); 

#ifdef COMPILER_XTENSA 
    #pragma loop_count min=1
#endif
        for ( l1=0; l1<(L>>LOG2_BBE_SIMD_WIDTH); l1++)
        {
            _X = (const xb_vecNx16 *)(x+BBE_SIMD_WIDTH*l1);
            R  = (      xb_vecNx16 *)(r+BBE_SIMD_WIDTH*l1);
         
            D0   = BBE_LVNX16_I(SECT, 2*0*BBE_SIMD_WIDTH); 
            D1   = BBE_LVNX16_I(SECT, 2*1*BBE_SIMD_WIDTH);
            D0o  = BBE_LVNX16_I(SECT, 2*2*BBE_SIMD_WIDTH);
            D0oo = BBE_LVNX16_I(SECT, 2*3*BBE_SIMD_WIDTH);

            BBE_LVNX16_XP(X0, _X, L*2);
            BBE_LVNX16_XP(X1, _X, L*2);
     
            // Q29 <- Q14*Q15
            Acc_fb0 = BBE_MULRNX16(G,  X0, rnd14);
            Acc_ff0 = BBE_MULRNX16(B0, D0, rnd14);
            Acc_fb1 = BBE_MULRNX16(G,  X1, rnd14);
            Acc_ff1 = BBE_MULRNX16(B0, D1, rnd14);
            
#ifdef COMPILER_XTENSA 
    #pragma ymemory(_X)
#endif
            for ( n=0; n<N-2; n+=2 )
            {                
                // Q29 <- Q29 - Q14*Q15
                BBE_MULANX16PR(Acc_fb0, D0,   D1,  a12);
                BBE_MULANX16PR(Acc_ff0, D0oo, D0o, b12);
                BBE_MULANX16PR(Acc_ff1, D0o,  D0,  b12);
                
                // Q15 <- Q29 - 14 w/ rounding
                D2 = BBE_PACKVNX40(Acc_fb0, rnd14);
                Y0 = BBE_PACKVNX40(Acc_ff0, rnd14);

                // Q29 <- Q29 - Q14*Q15
                BBE_MULANX16PR(Acc_fb1, D1, D2, a12);
                
                D0oo = D0; D0o = D1; D0 = D2; 
                
                // Q15 <- Q29 - 14 w/ rounding
                D1 = BBE_PACKVNX40(Acc_fb1, rnd14);
                Y1 = BBE_PACKVNX40(Acc_ff1, rnd14);

                BBE_LVNX16_XP(X0, _X, L*2);
                BBE_LVNX16_XP(X1, _X, L*2);

                // Q29 <- Q14*Q15
                Acc_fb0 = BBE_MULRNX16(G,  X0, rnd14);
                Acc_ff0 = BBE_MULRNX16(B0, D0, rnd14);
                Acc_fb1 = BBE_MULRNX16(G,  X1, rnd14);
                Acc_ff1 = BBE_MULRNX16(B0, D1, rnd14);

                // Q15 <- Q15*Q8 - 8
                z0 = BBE_MULNX16(Y0, Gain);
                z1 = BBE_MULNX16(Y1, Gain);
                Y0 = BBE_PACKVNX40(z0, rnd8);
                Y1 = BBE_PACKVNX40(z1, rnd8);

                BBE_SVNX16_XP(Y0, R, 2*L);
                BBE_SVNX16_XP(Y1, R, 2*L);
            }

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16PR(Acc_fb0, D0,   D1,  a12);
            BBE_MULANX16PR(Acc_ff0, D0oo, D0o, b12);
            BBE_MULANX16PR(Acc_ff1, D0o,  D0,  b12);
            
            // Q15 <- Q29 - 14 w/ rounding
            D2 = BBE_PACKVNX40(Acc_fb0, rnd14);
            Y0 = BBE_PACKVNX40(Acc_ff0, rnd14);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULANX16PR(Acc_fb1, D1, D2, a12);
            
            D0oo = D0; D0o = D1; D0 = D2; 
            
            // Q15 <- Q29 - 14 w/ rounding
            D1 = BBE_PACKVNX40(Acc_fb1, rnd14);
            Y1 = BBE_PACKVNX40(Acc_ff1, rnd14);

            // Q15 <- Q15*Q8 - 8
            z0 = BBE_MULNX16(Y0, Gain);
            z1 = BBE_MULNX16(Y1, Gain);
            Y0 = BBE_PACKVNX40(z0, rnd8);
            Y1 = BBE_PACKVNX40(z1, rnd8);

            BBE_SVNX16_XP(Y0, R, 2*L);
            BBE_SVNX16_XP(Y1, R, 2*L);

            BBE_SVNX16_IP(D0,   SECT, 2*BBE_SIMD_WIDTH); 
            BBE_SVNX16_IP(D1,   SECT, 2*BBE_SIMD_WIDTH); 
            BBE_SVNX16_IP(D0o,  SECT, 2*BBE_SIMD_WIDTH);
            BBE_SVNX16_IP(D0oo, SECT, 2*BBE_SIMD_WIDTH);
        }
    }
} // bqriirs_sp_proc()
