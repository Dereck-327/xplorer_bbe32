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
  NatureDSP_Baseband library. FIR filters and Related Functions
    Despreading
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fir.h"

/*-------------------------------------------------------------------------
Despreading

Perform a special kind of correlation utilized in Direct-Sequence Spread Spectrum
(DSSS) modulation technique, with Spreading Factor (SF) 4 or 8. Despreading
functions multiply the input sequence of 16-bit complex samples x[] by a coded
pseudo-noise sequence of 2-bit QPSK symbols pn[], then reduce the product sequence
by adding together 4 or 8 contiguous products. Resulting sums are shifted to the
right by rsh bit positions, saturated to -32768 or 32767 and stored to the output
array y[].

For the coded sequence multiplicand pn[], there are two variants of 2-bit QPSK
codeset identified with the qpsk_type argument:

         | QPKS symbol for | QPKS symbol for 
   Dibit | qpsk_type == 0  | qpsk_type == 1
  -------+-----------------+-----------------
     00  |      1 + j      |       1
     01  |     -1 + j      |      -j
     10  |      1 - j      |      -1
     11  |     -1 - j      |       j

First 8 consecutive dibits of the coded sequence are concatenated in a 16-bit
word pn[0], with the very first dibit going into 2 LSBs of pn[0] and the last
dibit - into 2 MSBs. The next 8 dibits are stored in pn[1], and so on.

Parameters:
  Input:
  x[N]          Input signal, 16-bit complex samples
  pn[N/8]       Pseudo-noise coded sequence of 2-bit QPSK symbols
  qpsk_type     QPSK codeset selection, 0 or 1
  rsh           Right shift amount for reduced product sums, 0..31
  Output:
  y[N/SF]       Output signal, 16-bit complex samples. Fixed point position
                for the output signal is Qx-rsh, where Qx is the fixed point
                position for the input signal.
Restrictions:
  x[],y[],pn[]  Must not overlap
  x[],y[],pn[]  Must be aligned on 32-byte boundary
  N             Multiple of 8
 -------------------------------------------------------------------------*/

#if !HAVE_DSPR
DISCARD_FUN(void, despread8, (complex_fract16 * restrict y,
         const complex_fract16 * restrict x,
         const int16_t * restrict pn,
         int N, int qpsk_type, int rsh))
#else

#define NSTEP   (2*BBE_SIMD_WIDTH/2)

#define DESPREAD8_ITERN( qpsk_type )                         \
{                                                            \
    xb_vecNx16 p0;                                           \
    xb_vecNx16 x0, x1, x2, x3, x4, x5, x6, x7;               \
    xb_vecNx16 x8, x9, x10, x11, x12, x13, x14, x15;         \
    xb_vecNx16 y0, y1, y2, y3, y4, y5, y6, y7;               \
    xb_vecNx16 y8, y9, y10, y11, y12, y13, y14, y15;         \
    xb_vecNx40 w0, w1, w2, w3, w4, w5, w6, w7;               \
    xb_vecNx40 w8, w9, w10, w11, w12, w13, w14, w15;         \
                                                             \
    BBE_LVNX16_IP( p0, PN, +4*BBE_SIMD_WIDTH/2 );            \
                                                             \
    BBE_LVNX16_IP( x0 , X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x1 , X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x2 , X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x3 , X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x4 , X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x5 , X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x6 , X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x7 , X, +4*BBE_SIMD_WIDTH/2 );            \
                                                             \
    BBE_LVNX16_IP( x8 , X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x9 , X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x10, X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x11, X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x12, X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x13, X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x14, X, +4*BBE_SIMD_WIDTH/2 );            \
    BBE_LVNX16_IP( x15, X, +4*BBE_SIMD_WIDTH/2 );            \
                                                             \
    w0  = BBE_DSPR1DNX16CSF8( x0 , p0, 0 , qpsk_type );      \
    w1  = BBE_DSPR1DNX16CSF8( x1 , p0, 1 , qpsk_type );      \
    w2  = BBE_DSPR1DNX16CSF8( x2 , p0, 2 , qpsk_type );      \
    w3  = BBE_DSPR1DNX16CSF8( x3 , p0, 3 , qpsk_type );      \
    w4  = BBE_DSPR1DNX16CSF8( x4 , p0, 4 , qpsk_type );      \
    w5  = BBE_DSPR1DNX16CSF8( x5 , p0, 5 , qpsk_type );      \
    w6  = BBE_DSPR1DNX16CSF8( x6 , p0, 6 , qpsk_type );      \
    w7  = BBE_DSPR1DNX16CSF8( x7 , p0, 7 , qpsk_type );      \
                                                            \
    w8  = BBE_DSPR1DNX16CSF8( x8 , p0, 8 , qpsk_type );      \
    w9  = BBE_DSPR1DNX16CSF8( x9 , p0, 9 , qpsk_type );      \
    w10 = BBE_DSPR1DNX16CSF8( x10, p0, 10, qpsk_type );      \
    w11 = BBE_DSPR1DNX16CSF8( x11, p0, 11, qpsk_type );      \
    w12 = BBE_DSPR1DNX16CSF8( x12, p0, 12, qpsk_type );      \
    w13 = BBE_DSPR1DNX16CSF8( x13, p0, 13, qpsk_type );      \
    w14 = BBE_DSPR1DNX16CSF8( x14, p0, 14, qpsk_type );      \
    w15 = BBE_DSPR1DNX16CSF8( x15, p0, 15, qpsk_type );      \
                                                             \
    y0  = BBE_PACKVNX40( w0 , vsa0 );                        \
    y1  = BBE_PACKVNX40( w1 , vsa0 );                        \
    y2  = BBE_PACKVNX40( w2 , vsa0 );                        \
    y3  = BBE_PACKVNX40( w3 , vsa0 );                        \
    y4  = BBE_PACKVNX40( w4 , vsa0 );                        \
    y5  = BBE_PACKVNX40( w5 , vsa0 );                        \
    y6  = BBE_PACKVNX40( w6 , vsa0 );                        \
    y7  = BBE_PACKVNX40( w7 , vsa0 );                        \
                                                             \
    y8  = BBE_PACKVNX40( w8 , vsa0 );                        \
    y9  = BBE_PACKVNX40( w9 , vsa0 );                        \
    y10 = BBE_PACKVNX40( w10, vsa0 );                        \
    y11 = BBE_PACKVNX40( w11, vsa0 );                        \
    y12 = BBE_PACKVNX40( w12, vsa0 );                        \
    y13 = BBE_PACKVNX40( w13, vsa0 );                        \
    y14 = BBE_PACKVNX40( w14, vsa0 );                        \
    y15 = BBE_PACKVNX40( w15, vsa0 );                        \
                                                             \
    y0 = BBE_SELNX16I( y1 , y0 , BBE_SELI_INTERLEAVE_2_LO ); \
    y1 = BBE_SELNX16I( y3 , y2 , BBE_SELI_INTERLEAVE_2_LO ); \
    y2 = BBE_SELNX16I( y5 , y4 , BBE_SELI_INTERLEAVE_2_LO ); \
    y3 = BBE_SELNX16I( y7 , y6 , BBE_SELI_INTERLEAVE_2_LO ); \
    y4 = BBE_SELNX16I( y9 , y8 , BBE_SELI_INTERLEAVE_2_LO ); \
    y5 = BBE_SELNX16I( y11, y10, BBE_SELI_INTERLEAVE_2_LO ); \
    y6 = BBE_SELNX16I( y13, y12, BBE_SELI_INTERLEAVE_2_LO ); \
    y7 = BBE_SELNX16I( y15, y14, BBE_SELI_INTERLEAVE_2_LO ); \
                                                             \
    y0 = BBE_SELNX16( y1, y0, sel0 );                        \
    y1 = BBE_SELNX16( y3, y2, sel0 );                        \
    y2 = BBE_SELNX16( y5, y4, sel0 );                        \
    y3 = BBE_SELNX16( y7, y6, sel0 );                        \
                                                             \
    y0 = BBE_SELNX16I( y1, y0, BBE_SELI_EXTRACT_LO_HALVES ); \
    y1 = BBE_SELNX16I( y3, y2, BBE_SELI_EXTRACT_LO_HALVES ); \
                                                             \
    BBE_SVNX16_IP( y0, Y, +4*BBE_SIMD_WIDTH/2 );             \
    BBE_SVNX16_IP( y1, Y, +4*BBE_SIMD_WIDTH/2 );             \
}

#define DESPREAD8_ITER1( qpsk_type )                         \
{                                                            \
  xb_vecNx16 p0, x0, y0;                                     \
  xb_vecNx40 w0;                                             \
                                                             \
  BBE_LSNX16_IP( p0, PN_, +2 );                              \
                                                             \
  BBE_LVNX16_IP( x0, X, +4*BBE_SIMD_WIDTH/2 );               \
                                                             \
  w0 = BBE_DSPR1DNX16CSF8( x0, p0, 0, qpsk_type );           \
                                                             \
  y0 = BBE_PACKVNX40( w0, vsa0 );                            \
                                                             \
  BBE_SAVNX16_XP( y0, Y_va, Y, +4*(BBE_SIMD_WIDTH/2)/8 );    \
}

void despread8(complex_fract16 * restrict y,
         const complex_fract16 * restrict x,
         const int16_t * restrict pn,
         int N, int qpsk_type, int rsh)
{
    xb_vecNx16 * restrict Y;
    const xb_vecNx16 *          X;
    const xb_vecNx16 *          PN;

    valign Y_va;

    vsaN  vsa0;
    vselN sel0;

    static const union
    {
        int16_t i[16];
        _vselN s;
        _vsaN v;
    } ALIGN(32) cst = { { 0, 1, 2, 3, 16, 17, 18, 19 } };

    int n;

    const short * PN_;

    vsa0 = BBE_MOVVSA32(rsh);

    sel0 = _S(cst.s);

    Y = (xb_vecNx16*)y;
    X = (const xb_vecNx16*)x;
    PN = (const xb_vecNx16*)pn;

    if (qpsk_type == 0)
    {
        __Pragma("ymemory(PN)");
        __Pragma("ymemory(X)");
        for (n = 0; n<(N / 8) / NSTEP; n++)
        {
            DESPREAD8_ITERN(0);
        }

        Y_va = BBE_ZALIGN();
        PN_ = (const short *)PN;

        __Pragma("concurrent");
        __Pragma("ymemory(PN)");
        __Pragma("ymemory(X)");
        for (n = 0; n<((N / 8) & (NSTEP - 1)); n++)
        {
            DESPREAD8_ITER1(0);
        }

        BBE_SAVNX16POS_FP(Y_va, Y);
    }
    else
    {
        __Pragma("ymemory(PN)");
        __Pragma("ymemory(X)");
        for (n = 0; n<(N / 8) / NSTEP; n++)
        {
            DESPREAD8_ITERN(1);
        }

        Y_va = BBE_ZALIGN();
        PN_ = (const short *)PN;

        __Pragma("concurrent");
        __Pragma("ymemory(PN)");
        __Pragma("ymemory(X)");
        for (n = 0; n<((N / 8) & (NSTEP - 1)); n++)
        {
            DESPREAD8_ITER1(1);
        }

        BBE_SAVNX16POS_FP(Y_va, Y);
    }
} /* despread8() */

#endif
