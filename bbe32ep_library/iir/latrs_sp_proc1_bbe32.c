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
void latrs_sp_proc1( int16_t * restrict r,
                     int16_t * restrict d,
               const int16_t *          x,
               const int16_t *          coef,
                     int16_t            gain,
                     int N, int L )
{
    xb_vecNx16 X0, X1, D0, D1, Y0, Y1, Coef, g;

    xb_vecNx40 T0, T1;

    const xb_vecNx16 *_X;
          xb_vecNx16 *_D, *_R;

    int     j, l;

    ASSERT( r && d && x && coef );
    
    g = BBE_MOVVA16(gain>>1);
    Coef = BBE_LSNX16_I(coef, 0);
    Coef = BBE_REPNX16(Coef, 0);

    _D = (xb_vecNx16 *)d;

    for ( l=0; l<(L>>(LOG2_BBE_SIMD_WIDTH+1)); l++)
    {
        _X = (const xb_vecNx16 *)(x + sz_i16*BBE_SIMD_WIDTH*l);
        _R = (      xb_vecNx16 *)(r + sz_i16*BBE_SIMD_WIDTH*l);

        D0 = BBE_LVNX16_I(_D, 0*sz_i16*BBE_SIMD_WIDTH);
        D1 = BBE_LVNX16_I(_D, 1*sz_i16*BBE_SIMD_WIDTH);

        X1 = BBE_LVNX16_I(_X, sz_i16*BBE_SIMD_WIDTH);
        BBE_LVNX16_XP(X0, _X, sz_i16*L);

        for ( j=0; j<N-1; j++ )
        {
            // Q29 <- Q15*Q14
            T0 = BBE_MULNX16(X0, g);
            T1 = BBE_MULNX16(X1, g);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULSNX16(T0, D0, Coef); 
            BBE_MULSNX16(T1, D1, Coef); 

            D0 = BBE_PACKQNX40(T0);
            D1 = BBE_PACKQNX40(T1);

            // Q15 <- Q29 - 14 w/ rounding
            T0 = BBE_ADDNX40(T0, T0);
            T1 = BBE_ADDNX40(T1, T1);

            Y0 = BBE_PACKQNX40(T0);
            Y1 = BBE_PACKQNX40(T1);

            X1 = BBE_LVNX16_I(_X, sz_i16*BBE_SIMD_WIDTH);
            BBE_LVNX16_XP(X0, _X, sz_i16*L);

            BBE_SVNX16_X (Y1, _R, sz_i16*BBE_SIMD_WIDTH);
            BBE_SVNX16_XP(Y0, _R, sz_i16*L);
        }

        // Q29 <- Q15*Q14
        T0 = BBE_MULNX16(X0, g);
        T1 = BBE_MULNX16(X1, g);

        // Q29 <- Q29 - Q14*Q15
        BBE_MULSNX16(T0, D0, Coef); 
        BBE_MULSNX16(T1, D1, Coef); 

        D0 = BBE_PACKQNX40(T0);
        D1 = BBE_PACKQNX40(T1);

        // Q15 <- Q29 - 14 w/ rounding
        T0 = BBE_ADDNX40(T0, T0);
        T1 = BBE_ADDNX40(T1, T1);

        Y0 = BBE_PACKQNX40(T0);
        Y1 = BBE_PACKQNX40(T1);

        BBE_SVNX16_X (Y1, _R, sz_i16*BBE_SIMD_WIDTH);
        BBE_SVNX16_XP(Y0, _R, sz_i16*L);

        BBE_SVNX16_IP(D0, _D, sz_i16*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(D1, _D, sz_i16*BBE_SIMD_WIDTH);
    }

    if ( L&(2*BBE_SIMD_WIDTH-1) )
    {
        _X = (const xb_vecNx16 *)(x + L - BBE_SIMD_WIDTH);
        _R = (      xb_vecNx16 *)(r + L - BBE_SIMD_WIDTH);

        D0 = BBE_LVNX16_I(_D, 0);

        if ( N&1 )
        {

            BBE_LVNX16_XP(X0, _X, sz_i16*L);
            // Q29 <- Q15*Q14
            T0 = BBE_MULNX16(X0, g);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULSNX16(T0, D0, Coef); 

            D0 = BBE_PACKQNX40(T0);

            // Q15 <- Q29 - 14 w/ rounding
            T0 = BBE_ADDNX40(T0, T0);

            Y0 = BBE_PACKQNX40(T0);

            BBE_SVNX16_XP(Y0, _R, sz_i16*L);
        }

        for ( j=0; j<(N>>1); j++ )
        {

            BBE_LVNX16_XP(X0, _X, sz_i16*L);
            BBE_LVNX16_XP(X1, _X, sz_i16*L);
            // Q29 <- Q15*Q14
            T0 = BBE_MULNX16(X0, g);
            T1 = BBE_MULNX16(X1, g);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULSNX16(T0, D0, Coef); 

            D1 = BBE_PACKQNX40(T0);

            // Q29 <- Q29 - Q14*Q15
            BBE_MULSNX16(T1, D1, Coef); 

            D0 = BBE_PACKQNX40(T1);

            // Q15 <- Q29 - 14 w/ rounding
            T0 = BBE_ADDNX40(T0, T0);
            T1 = BBE_ADDNX40(T1, T1);

            Y0 = BBE_PACKQNX40(T0);
            Y1 = BBE_PACKQNX40(T1);

            BBE_SVNX16_XP(Y0, _R, sz_i16*L);
            BBE_SVNX16_XP(Y1, _R, sz_i16*L);
        }

        BBE_SVNX16_IP(D0, _D, sz_i16*BBE_SIMD_WIDTH);
    }
} // latr_sp_proc1()
