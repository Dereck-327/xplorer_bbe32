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

/* Portable data types. */
#include "NatureDSP_types.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_fft.h"
/* Common utility declarations. */
#include "common.h"
/* Twiddles tables for float point FFT */
#include "fft_fp_tw.h"
/* Internal components for the floating point FFT */
#include "fft_fp_common.h"
/* Internal components for the floating point FFT */
#include "fft_fp_common.h"

#if XCHAL_HAVE_BBEN_VECTORFFT && HAVE_VFPU
/* The in-place real-to-complex spectrum conversion, N = 16*/
inline_ void blkrfft_spec_conv_packed_N16_fp(complex_float *_y,   /* in/out  complex array, size is N/2+1 */
    const complex_float *twiddle_table,
    int N, 
    int L)
{
    valign u, u1;
    int l;
    xb_vecN_4xcf32 *px0;
    xb_vecN_4xcf32 *px1;
    xb_vecN_4xcf32 *py1;
    xb_vecN_4xcf32 *py0;

    const int v_sz = (BBE_SIMD_WIDTH * 2 / sizeof(*_y)); //Vector size (number of complex elements)
    NASSERT_ALIGN32(_y);
    NASSERT(N==16); 
    vboolN_4 mask = BBE_NOTBN_4(BBE_MOVN_4_FROMN(BBE_LTRNI(4)));

    {
        xb_vecN_4xcf32 _a0, _a1, _b0, _b1, tw;
        xb_vecN_4xcf32 *ptw = (xb_vecN_4xcf32 *)(twiddle_table + N / 4 - v_sz);
     
        px1 = (xb_vecN_4xcf32 *)(_y + N / 2 - v_sz + 1);
        py1 = (xb_vecN_4xcf32 *)(_y + N / 2 - v_sz + 1);
        px0 = (xb_vecN_4xcf32 *)(_y);
        py0 = (xb_vecN_4xcf32 *)(_y );
             
        BBE_LVN_4XCF32_XP(tw, ptw, -2 * BBE_SIMD_WIDTH);

        /* 11 cycles per pipeline stage in steady state with unroll=2 */
        for (l = 0; l<L; l++)
        {

            u1 = BBE_ZALIGN();
            u = BBE_LAN_4XCF32_PP(px1);
                    

            BBE_LVN_4XCF32_XP(_a0, px0, N / 2 * sizeof(*_y));
            BBE_LAN_4XCF32_IP(_a1, u, px1);

            px1 = (xb_vecN_4xcf32*)(N / 2 * sizeof(*_y) - 2 * BBE_SIMD_WIDTH + (uintptr_t)px1);

            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            _a1 = BBE_CONJN_4XCF32(_a1);

            BBE_ADDSUBN_4XCF32(_b1, _b0, _a0, _a1);

            _b1 = BBE_MULN_4XCF32(_b1, tw);
            xb_vecN_2xf32 tmp = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(_b0));
            tmp = BBE_MULN_2XF32(tmp, 0.5f);
            _b0 = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(tmp));

            BBE_ADDSUBN_4XCF32(_a1, _a0, _b0, _b1);

            _a1 = BBE_SHFLN_4XCF32I(_a1, BBE_SHFLI_REVERSE_4);
            _a1 = BBE_CONJN_4XCF32(_a1);

            BBE_SAN_4XCF32_IP(_a1, u1, py1);

            py1 = (xb_vecN_4xcf32*)(N / 2 * sizeof(*_y) - 2 * BBE_SIMD_WIDTH + (uintptr_t)py1);
            BBE_SVN_4XCF32T_XP(_a0, py0, N / 2 * sizeof(*_y), mask);
                 
        } // for (l=0; l<L; l++)
    }   

    __Pragma("no_reorder");

    px0 = (xb_vecN_4xcf32 *)(_y);
    py0 = (xb_vecN_4xcf32 *)(_y);
     /* 16 cycles per pipeline stage in steady state with unroll = 8 */
    for (l = 0; l<L; l++)
    {
        xb_vecN_4xcf32 a0, aN_4;

        BBE_LVN_4XCF32_XP(a0, px0, N / 4 * sizeof(*_y));
        BBE_LVN_4XCF32_XP(aN_4, px0, N / 4 * sizeof(*_y));
        aN_4 = BBE_CONJN_4XCF32(aN_4);
        a0 = BBE_MULJN_4XCF32(BBE_CONSTN_4XCF32(0x1), a0);
        BBE_SVN_4XCF32F_XP(a0, py0, N / 4 * sizeof(*_y), mask);
        BBE_SVN_4XCF32F_XP(aN_4, py0, N / 4 * sizeof(*_y), mask);
    }
} /* blkrfft_spec_conv_packed_N16_fp */



/*
Table of the twiddles for cifftf8, generated by matlab:

tw_cifftf8 = [exp(1j*2*pi/8*(0:3)'); exp(1j*2*pi/4*fix((0:3)'/2))]
print_twd({tw_cifftf8}, 'cifftf8_twd', 1)
%=================================================================
function print_twd(all_twd, name, fid)
for i=1:length(all_twd)
t = all_twd{i};
L = length(t);
tmp =  reshape( [real(t), imag(t)].', 2*L, 1 ) ;
s = num2hex(single(tmp));
fprintf(2, 'extern const union ufloat32uint32 ');
fprintf(2, '%s%d[%d];\n', name, i, 2*L);

fprintf(fid, '\nALIGN(32) const union ufloat32uint32 ');
fprintf(fid, '%s%d[%d] = \n{\n', name, i, 2*L);
for k=1:2*L
fprintf(fid, '  { 0x%s },', s(k, :));
if(mod(k, 4)==0)
fprintf(fid, '\n');
end
end
fprintf(fid, '};\n');
end
*/
static ALIGN(32) const union ufloat32uint32 cifftf8_twd1[16] =
{
    { 0x3f800000 }, { 0x00000000 }, { 0x3f3504f3 }, { 0x3f3504f3 },
    { 0x248d3132 }, { 0x3f800000 }, { 0xbf3504f3 }, { 0x3f3504f3 },
    { 0x3f800000 }, { 0x00000000 }, { 0x3f800000 }, { 0x00000000 },
    { 0x248d3132 }, { 0x3f800000 }, { 0x248d3132 }, { 0x3f800000 },
};



static void blk_cfftf8_fp(complex_float *x,
    complex_float *y,
    int L)
{

    int i;


    xb_vecN_4xcf32 *px = (xb_vecN_4xcf32 *)x;
    xb_vecN_4xcf32 *py = (xb_vecN_4xcf32 *)y;

    xb_vecN_4xcf32 _a0, _a1, _t1, _t2, _b0, _b1;




    _t1 = BBE_LVN_4XCF32_I((xb_vecN_4xcf32*)cifftf8_twd1, 0);
    _t2 = BBE_LVN_4XCF32_I((xb_vecN_4xcf32*)cifftf8_twd1, 2*BBE_SIMD_WIDTH);

    vboolN_2 m = BBE_NOTBN_2(BBE_MOVN_2_FROMN(BBE_LTRNI(8)));
    for (i = 0; i < L; i++)
    {
        xb_vecN_2xf32 z0;
       // xb_vecN_4xcf32 test; 
        /* 13 cycles per pipeline stage in steady state with unroll=2 */

        /*
        complex_float a0, a1;
        for (i = 0; i < N/2; i++)
        {
        a0.s.re = x[i].s.re + x[i + N / 2].s.re;
        a1.s.re = x[i].s.re - x[i + N / 2].s.re;
        a0.s.im = x[i].s.im + x[i + N / 2].s.im;
        a1.s.im = x[i].s.im - x[i + N / 2].s.im;
        a1 = mul_fl32c(tw1[i], a1);

        y[2*i  ] = a0;
        y[2*i+1] = a1;
        }

        */
        BBE_LVN_4XCF32_IP(_a0, px, sizeof(*px));
        BBE_LVN_4XCF32_IP(_a1, px, sizeof(*px));
        BBE_ADDSUBN_4XCF32(_a1, _a0, _a0, _a1);
        _a1 = BBE_MULJN_4XCF32(_a1, _t1);

        BBE_DSELN_4XCF32I(_b1, _b0, _a1, _a0, BBE_DSELI_INTERLEAVE_4);
        BBE_ADDSUBN_4XCF32(_b1, _b0, _b0, _b1);
              
        z0 = BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(_b1) );
        /*  muliply by {1, 1, -j, -j} */
        BBE_MULMN_2XF32T(z0, z0, BBE_CONSTN_2XF32(1), 2, 6, m); 
        _b1 = BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(z0)); 

        _a0 = BBE_SELN_4XCF32I(_b1, _b0, BBE_SELI_EXTRACT_LO_HALVES);
        _a1 = BBE_SELN_4XCF32I(_b1, _b0, BBE_SELI_EXTRACT_HI_HALVES);
        /*


        for (i = 0; i < N/4; i++)
        {
        int j;
        for (j = 0; j < 2; j++)
        {
        a0.s.re = y[(2 * i + j)].s.re + y[(2 * i + j) + N / 2].s.re;
        a1.s.re = y[(2 * i + j)].s.re - y[(2 * i + j) + N / 2].s.re;
        a0.s.im = y[(2 * i + j)].s.im + y[(2 * i + j) + N / 2].s.im;
        a1.s.im = y[(2 * i + j)].s.im - y[(2 * i + j) + N / 2].s.im;
        a1 = mul_fl32c(tw2[i*2+j], a1);

        x[4 * i + j] = a0;
        x[4 * i + j + 2] = a1;
        }

        }
        */

        BBE_ADDSUBN_4XCF32(_a1, _a0, _a0, _a1);

        /*
        for (i = 0; i < N / 2; i++)
        {
        a0.s.re = (x[i].s.re + x[i + N / 2].s.re);
        a1.s.re = (x[i].s.re - x[i + N / 2].s.re);
        a0.s.im = (x[i].s.im + x[i + N / 2].s.im);
        a1.s.im = (x[i].s.im - x[i + N / 2].s.im);
        y[ i]      = a0;
        y[ i + N/2] = a1;
        }
        */
        BBE_SVN_4XCF32_IP(_a0, py, sizeof(*py));
        BBE_SVN_4XCF32_IP(_a1, py, sizeof(*py));
    }//for(i=0; i<L; i++)
} //blk_cfftf8_fp

#endif //#if XCHAL_HAVE_BBEN_VECTORFFT && HAVE_VFPU

/*-------------------------------------------------------------------------
Blockwise radix-2 floating point forward FFT on real data

Description: These functions make forward real FFT on L blocks, each of N=2^n
complex samples, where n=4..7. 
FFT implementation for real signal exploits the symmetry properties of the 
Fourier Transform: first, a blockwise complex FFT of half the original size 
is applied to input data, and then the resulting spectrum undergoes a 
conversion procedure which results in complex spectrum of real input data.
NOTES:
1. Bit-reversing permutation is done here. 
2. FFT runs an in-place algorithm, so INPUT DATA WILL APPEAR DAMAGED after 
   the call. 
3. As the complex spectrum of a real FFT is conjugate symmetric about the 
   midpoint, the brfftf functions generate only the first (N/2)+1 points 
   of the FFT, with the first point being the DC component, and the last 
   point ÿ the Nyquist frequency component. For real signal these two 
   spectral components have zero imaginary part, thus they can be packed 
   in a single 'complex' value. Finally, the the m-th row of the output 
   matrix will contain the following values:
    - DC component of the signal in y[m][0].re
    - Nyquist frequency component in y[m][0].im
    - first half of the complex spectrum in y[m][1]...y[m][(N/2)-1]
   These data are used to reconstruct N values that are stored into the 
   m-th row of real output matrix y[m][N]. The reconstruction is 
   accomplished through a special conversion of input spectrum, such 
   that subsequent invocation of blockwise complex inverse FFT of size 
   N/2 actually produces the desired real signal.

Representation: floating point


Parameters:
  Input:            
    x[L][N]   Real input signal
  Output:          
    y[L][N/2] Half of output spectrum, see note above
  Returned value:
                None
Restrictions:
  x,y           Must not overlap and must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
#if !XCHAL_HAVE_BBEN_VECTORFFT || !HAVE_VFPU
DISCARD_FUN(void, brfft16f, (complex_float * restrict y, float32_t * restrict x, int L) )
#else
void brfft16f(complex_float * restrict y, float32_t * restrict x, int L)
{   
    complex_float *twiddle_table = (complex_float*)rfftf16_twd1;

    NASSERT_ALIGN(x, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(y, 2 * BBE_SIMD_WIDTH);
    NASSERT((uintptr_t)x != (uintptr_t)y);
    NASSERT(L>0);

    blk_cfftf8_fp((complex_float*)x, y, L); 
    blkrfft_spec_conv_packed_N16_fp(y, twiddle_table, 16, L); 
 }
#endif //#if !XCHAL_HAVE_BBEN_VECTORFFT || !HAVE_VFPU






