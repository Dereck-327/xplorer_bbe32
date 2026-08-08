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
    NatureDSP_Baseband library. Eigenvalues and eigenvectors
    Real 2x2 stream ordered matrices
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#include <float.h>
/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* Baseband Library API */
#include "NatureDSP_Baseband_id.h"
#include "NatureDSP_Baseband_matinv.h"
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"
#include "vfpu_math.h"

#if HAVE_VFPU

#if 1

#define sz_f32   sizeof(float32_t)
#define sz_f32c  sizeof(complex_float)

#define MAX(a,b)    ( (a)>(b) ? (a) : (b) )
#define EPS       FLT_EPSILON

#if 0

#if !defined(COMPILER_MSVC)
#include <complex.h>
#endif

#include <math.h>

static complex_float _makecomplexf(float32_t re, float32_t im)
{
    union { float32_t r[2]; complex_float c; } u = { { re, im } };
    return (u.c);
}

/* Complex floating-point division, single precision. */
static complex_float _cdivf(complex_float x, complex_float y)
{
    /*
    * Based on cdiv code from "Similarity Reduction of a General Matrix to
    * Hessenberg Form" by R.S. Martin and J.H. Wilkinson, Handbook for
    * Automatie Computation, Vol.II Linear Algebra, Contribution II/13.
    * MATLAB reference code:
    *   function z = cdiv(x,y)
    *   xr = real(x); xi = imag(x);
    *   yr = real(y); yi = imag(y);
    *   if abs(yr)>abs(yi)
    *     h = yi/yr; yr = h*yi+yr;
    *     zr = (xr+h*xi)/yr;
    *     zi = (xi-h*xr)/yr;
    *   else
    *     h = yr/yi; yi = h*yr+yi;
    *     zr = (h*xr+xi)/yi;
    *     zi = (h*xi-xr)/yi;
    *   end
    *   z = complex(zr,zi);
    */

    float32_t xr, xi, yr, yi, zr, zi, h;
    xr = crealf(x); xi = cimagf(x);
    yr = crealf(y); yi = cimagf(y);
    if (fabsf(yr)>fabsf(yi)) {
        h = yi / yr; yr = h*yi + yr;
        zr = (xr + h*xi) / yr;
        zi = (xi - h*xr) / yr;
    }
    else {
        h = yr / yi; yi = h*yr + yi;
        zr = (h*xr + xi) / yi;
        zi = (h*xi - xr) / yi;
    }
    return (_makecomplexf(zr, zi));

} /* _cdivf() */

#endif

static void reigen_hqr_2x2_sf(complex_float * restrict e,
                              float32_t     * restrict H,
                              float32_t     * restrict P,
                              int L, int calcEigenVec)
#if 0
{
    const int N = 2;
    float32_t p, q, r;
    float32_t x, y, z, w;

    int k;

    for (k = 0; k < L; k++){
        x = H[(2 * 1 + 0)*L];
        y = H[(2 * 0 + 0)*L];
        z = H[(2 * 1 + 1)*L];

        if (fabsf(x) <= EPS*(fabsf(y) + fabsf(z))){
            /* Revealed a standalone eigenvalue. */
            e[0] = _makecomplexf(H[(2 * 0 + 0)*L], 0.f);
            e[L] = _makecomplexf(H[(2 * 1 + 1)*L], 0.f);
        } else {
            /* Revealed a pair of eigenvalues. */
            x = H[(2 * 1 + 1)*L];
            y = H[(2 * 0 + 0)*L];
            w = H[(2 * 0 + 1)*L] * H[(2 * 1 + 0)*L];

            p = (y - x) / 2; q = p*p + w; z = sqrtf(fabsf(q));
            if (q >= 0.f) {
                /* Pair of real eigenvalues */
                z = (p < 0.f ? p - z : p + z);
                e[0] = _makecomplexf(w != 0.f ? x + z : y, 0.f);
                e[L] = _makecomplexf(w != 0.f ? x - w / z : x, 0.f);
                if (P) {
                    /* Annihilate T(en,en-1) entry by similarity transformation involving
                    * a rotation in the en,en-1 plane. */
                    x = H[(2 * 1 + 0)*L];
                    r = sqrtf(x*x + z*z); p = x / r; q = z / r;
                    /* Row modification */
                    z = H[(2 * 0 + 0)*L];
                    H[(2 * 0 + 0)*L] = q*z + p*H[(2 * 1 + 0)*L];
                    H[(2 * 1 + 0)*L] = q*H[(2 * 1 + 0)*L] - p*z;
                    z = H[(2 * 0 + 1)*L];
                    H[(2 * 0 + 1)*L] = q*z + p*H[(2 * 1 + 1)*L];
                    H[(2 * 1 + 1)*L] = q*H[(2 * 1 + 1)*L] - p*z;
                    /* Column modification */
                    z = H[(2 * 0 + 0)*L];
                    H[(2 * 0 + 0)*L] = q*z + p*H[(2 * 0 + 1)*L];
                    H[(2 * 0 + 1)*L] = q*H[(2 * 0 + 1)*L] - p*z;
                    z = H[(2 * 1 + 0)*L];
                    H[(2 * 1 + 0)*L] = q*z + p*H[(2 * 1 + 1)*L];
                    H[(2 * 1 + 1)*L] = q*H[(2 * 1 + 1)*L] - p*z;
                    /* Accumulate transformations: right-multiply P by the rotator. */
                    P[0] = q;
                    P[L] = -p;
                    P[2*L] = p;
                    P[3*L] = q;
                }
            } else {
                /* Complex eigenvalue and its conjugate */
                e[0] = _makecomplexf(x + p, z);
                e[L] = _makecomplexf(x + p, -z);
            }
        }
        e++; H++; if (P) P++;
    } /* k */
} /* reigen_hqr_2x2_sf() */
#else
{
    const xb_vecN_2xf32  * restrict H_r;
    const xb_vecN_2xf32  * restrict Pd_r;
    const xb_vecN_2xf32  * restrict W_r;
    const xb_vecN_2xf32  * restrict Q_r;
    const xb_vecN_2xf32  * restrict Z_r;
    const vboolN_2       * B_r;
          xb_vecN_2xf32  * restrict H_w;
          xb_vecN_2xf32  * restrict Pd_w;
          xb_vecN_2xf32  * restrict W_w;
          xb_vecN_2xf32  * restrict Q_w;
          xb_vecN_2xf32  * restrict Z_w;
          vboolN_2       * B_w;
          xb_vecN_2xf32  * P_w;
          xb_vecN_4xcf32 * E_w;

    const xb_vecN_2xf32 _EPS = EPS;
    const xb_vecN_2xf32 c0f = BBE_CONSTN_2XF32(0);
    const xb_vecN_2xf32 c1f = BBE_CONSTN_2XF32(1);
    const xb_vecN_2xf32 c05f = BBE_CONSTN_2XF32(3);

    float32_t * Pd = &P[0 * L];
    float32_t * W  = &P[1 * L];
    float32_t * Q  = &P[2 * L];
    float32_t * Z  = &P[3 * L];
    vboolN_2  * B = (vboolN_2 *)&P[1 * L] - L / (BBE_SIMD_WIDTH / 2);

    xb_vecN_2xf32 h00, h01, h10, h11;
    xb_vecN_2xf32 p00, p01, p10, p11;
    xb_vecN_2xf32 _x, _y, _z, _q, _p, _r, _w;
    vboolN_2 ble, bqmez, bnez;   

    int k;

    NASSERT(0 == (L % (BBE_SIMD_WIDTH / 2)));

    H_r = (xb_vecN_2xf32  *)H;
    E_w = (xb_vecN_4xcf32 *)e;

    Pd_w = (xb_vecN_2xf32 *)Pd;
    Q_w = (xb_vecN_2xf32 *)Q;
    W_w = (xb_vecN_2xf32 *)W;
    Z_w = (xb_vecN_2xf32 *)Z;

    for (k = 0; k < L / (BBE_SIMD_WIDTH / 2); k++){
        xb_vecN_4xcf32 e00, e01, e10, e11;
        BBE_LVN_2XF32_XP(h00, H_r, 1 * L*sz_f32);
        BBE_LVN_2XF32_XP(h01, H_r, 1 * L*sz_f32);
        BBE_LVN_2XF32_XP(h10, H_r, 1 * L*sz_f32);
        BBE_LVN_2XF32_XP(h11, H_r, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);

        //x = H[(2 * 1 + 1)*L];
        //y = H[(2 * 0 + 0)*L];
        //w = H[(2 * 0 + 1)*L] * H[(2 * 1 + 0)*L];
        //p = (y - x) / 2; q = p*p + w; z = sqrtf(fabsf(q));
        // e[0] = _makecomplexf(x + p, z);
        // e[L] = _makecomplexf(x + p, -z);
        _x = h11;
        _y = h00;
        _w = _q = BBE_MULN_2XF32(h01, h10);
        _p = BBE_MULN_2XF32(BBE_SUBN_2XF32(_y, _x), c05f);
        BBE_MULAN_2XF32(_q, _p, _p);

        BBE_SVN_2XF32_IP(_w,  W_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_SVN_2XF32_IP(_p, Pd_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_SVN_2XF32_IP(_q,  Q_w, (BBE_SIMD_WIDTH / 2)*sz_f32);
        _z = BBE_SQRTN_2XF32(BBE_ABSN_2XF32(_q));
        BBE_SVN_2XF32_IP(_z,  Z_w, (BBE_SIMD_WIDTH / 2)*sz_f32);

        /* Complex eigenvalue and its conjugate */
        BBE_DSELN_2XF32I(p01, p00, _z, BBE_ADDN_2XF32(_x, _p), BBE_DSELI_INTERLEAVE_2);
        e00 = BBE_MOVN_4XCF32_FROMN_2XF32(p00);
        e01 = BBE_MOVN_4XCF32_FROMN_2XF32(p01);
        e10 = BBE_CONJN_4XCF32(e00);
        e11 = BBE_CONJN_4XCF32(e01);

        BBE_SVN_4XCF32_IP(e00, E_w, (BBE_SIMD_WIDTH / 4)*sz_f32c);
        BBE_SVN_4XCF32_XP(e01, E_w, L*sz_f32c);
        BBE_SVN_4XCF32_IP(e11, E_w, -(int)((BBE_SIMD_WIDTH / 4)*sz_f32c));
        BBE_SVN_4XCF32_XP(e10, E_w, (2 * (BBE_SIMD_WIDTH / 4) - L)*sz_f32c);
    }

    H_r = (xb_vecN_2xf32  *)H;
    E_w = (xb_vecN_4xcf32 *)e;

    Pd_r = (xb_vecN_2xf32 *)Pd;
    Q_r = (xb_vecN_2xf32 *)Q;
    W_r = (xb_vecN_2xf32 *)W;
    Z_r = Z_w = (xb_vecN_2xf32 *)Z;
    for (k = 0; k < L / (BBE_SIMD_WIDTH / 2); k++){
        xb_vecN_4xcf32 e00, e01, e10, e11;
        vboolN_2 ble_2_0, ble_2_1;
        vboolN_4 ble_4_0, ble_4_1;

        BBE_LVN_2XF32_XP(h00, H_r, 3 * L*sz_f32);
        BBE_LVN_2XF32_XP(h11, H_r, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);

        _x = h11;
        _y = h00;

        BBE_LVN_2XF32_IP(_p, Pd_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_q, Q_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_z, Z_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        BBE_LVN_2XF32_IP(_w, W_r, (BBE_SIMD_WIDTH / 2)*sz_f32);

        //e[0] = _makecomplexf(w != 0.f ? x + z : y, 0.f);
        //e[L] = _makecomplexf(w != 0.f ? x - w / z : x, 0.f);
        bqmez = BBE_OLEN_2XF32(c0f, _q); // q >= 0.f
        BBE_EXTRACTBN(ble_2_1, ble_2_0, BBE_MOVN_FROMN_2(bqmez));
        ble_4_0 = BBE_MOVN_4_FROMN_2(ble_2_0); ble_4_1 = BBE_MOVN_4_FROMN_2(ble_2_1);
        //z = (p < 0.f ? p - z : p + z);
        _z = BBE_MOVN_2XF32T(BBE_SUBN_2XF32(_p, _z), BBE_ADDN_2XF32(_p, _z), BBE_OLTN_2XF32(_p, c0f));
        BBE_SVN_2XF32_IP(_z, Z_w, (BBE_SIMD_WIDTH / 2)*sz_f32);

        bnez = BBE_NOTBN_2(BBE_UEQN_2XF32(c0f, _w));
        _y = BBE_MOVN_2XF32T(BBE_ADDN_2XF32(_x, _z), _y, bnez);
        _x = BBE_MOVN_2XF32T(BBE_SUBN_2XF32(_x, IT_FDIVN_2XF32(_w, _z, 1)), _x, bnez);

        BBE_DSELN_2XF32I(p01, p00, c0f, _y, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(p11, p10, c0f, _x, BBE_DSELI_INTERLEAVE_2);
        e00 = BBE_MOVN_4XCF32_FROMN_2XF32(p00);
        e01 = BBE_MOVN_4XCF32_FROMN_2XF32(p01);
        e10 = BBE_MOVN_4XCF32_FROMN_2XF32(p10);
        e11 = BBE_MOVN_4XCF32_FROMN_2XF32(p11);

        BBE_SVN_4XCF32T_IP(e00, E_w, (BBE_SIMD_WIDTH / 4)*sz_f32c, ble_4_0);
        BBE_SVN_4XCF32T_XP(e01, E_w, L*sz_f32c, ble_4_1);
        BBE_SVN_4XCF32T_IP(e11, E_w, -(int)((BBE_SIMD_WIDTH / 4)*sz_f32c), ble_4_1);
        BBE_SVN_4XCF32T_XP(e10, E_w, (2 * (BBE_SIMD_WIDTH / 4) - L)*sz_f32c, ble_4_0);
    }

    B_w = B;

    H_r = (xb_vecN_2xf32  *)H;
    E_w = (xb_vecN_4xcf32 *)e;
    for (k = 0; k < L / (BBE_SIMD_WIDTH / 2); k++){
        xb_vecN_4xcf32 e00, e01, e10, e11;
        vboolN_2 ble_2_0, ble_2_1;
        vboolN_4 ble_4_0, ble_4_1;
#if 0
        BBE_LVN_2XF32_XP(h00, H_r, 2 * L*sz_f32);
        BBE_LVN_2XF32_XP(h10, H_r, 1 * L*sz_f32);
        BBE_LVN_2XF32_XP(h11, H_r, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);
#elif 0
        h10 = BBE_LVN_2XF32_X(H_r, 2 * L*sz_f32);
        h11 = BBE_LVN_2XF32_X(H_r, 3 * L*sz_f32);
        //BBE_LVN_2XF32_IP(h00, H_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        h00 = BBE_LVN_2XF32_I(H_r, 0);
#else
        BBE_LVN_2XF32_XP(h00, H_r, 3 * L*sz_f32);
        BBE_LVN_2XF32_XP(h11, H_r, -1 * L*sz_f32);
        //BBE_LVN_2XF32_XP(h10, H_r, ((BBE_SIMD_WIDTH / 2) - 2 * L)*sz_f32);
        BBE_LVN_2XF32_XP(h10, H_r, -2 * L*sz_f32);
#endif

        //x = H[(2 * 1 + 0)*L];
        //y = H[(2 * 0 + 0)*L];
        //z = H[(2 * 1 + 1)*L];
        //fabsf(x) <= EPS*(fabsf(y) + fabsf(z))
        _x = BBE_ABSN_2XF32(h10);
        _y = BBE_MULN_2XF32(_EPS, BBE_ADDN_2XF32(BBE_ABSN_2XF32(h00), BBE_ABSN_2XF32(h11)));

#if 1
        ble = BBE_OLEN_2XF32(_x, _y);
        BBE_SBN_2_IP(ble, B_w, sizeof(vboolN_2));
        BBE_EXTRACTBN(ble_2_1, ble_2_0, BBE_MOVN_FROMN_2(ble));
        ble_4_0 = BBE_MOVN_4_FROMN_2(ble_2_0); ble_4_1 = BBE_MOVN_4_FROMN_2(ble_2_1);
#else
        ble_2_0 = BBE_OLEN_2XF32(BBE_SHFLN_2XF32I(_x, BBE_SHFLI_DOUBLE_2_LO), BBE_SHFLN_2XF32I(_y, BBE_SHFLI_DOUBLE_2_LO));
        ble_4_0 = BBE_MOVN_4_FROMN_2(ble_2_0);
        ble_2_1 = BBE_OLEN_2XF32(BBE_SHFLN_2XF32I(_x, BBE_SHFLI_DOUBLE_2_HI), BBE_SHFLN_2XF32I(_y, BBE_SHFLI_DOUBLE_2_HI));
        ble_4_1 = BBE_MOVN_4_FROMN_2(ble_2_1);
#endif

#if 1
        h11 = BBE_LVN_2XF32_X(H_r, 3 * L*sz_f32);
        BBE_LVN_2XF32_IP(h00, H_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
#else
        BBE_LVN_2XF32_XP(h11, H_r, -3 * L*sz_f32);
        BBE_LVN_2XF32_IP(h00, H_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
#endif

        //e[0] = _makecomplexf(H[(2 * 0 + 0)*L], 0.f);
        //e[L] = _makecomplexf(H[(2 * 1 + 1)*L], 0.f);
#if 0
        e00 = BBE_MOVN_4XCF32_FROMN_2XF32(BBE_SELN_2XF32I(c0f, h00, BBE_SELI_INTERLEAVE_2_LO));
        e01 = BBE_MOVN_4XCF32_FROMN_2XF32(BBE_SELN_2XF32I(c0f, h00, BBE_SELI_INTERLEAVE_2_HI));
        e10 = BBE_MOVN_4XCF32_FROMN_2XF32(BBE_SELN_2XF32I(c0f, h11, BBE_SELI_INTERLEAVE_2_LO));
        e11 = BBE_MOVN_4XCF32_FROMN_2XF32(BBE_SELN_2XF32I(c0f, h11, BBE_SELI_INTERLEAVE_2_HI));
#else
        BBE_DSELN_2XF32I(p01, p00, c0f, h00, BBE_DSELI_INTERLEAVE_2);
        BBE_DSELN_2XF32I(p11, p10, c0f, h11, BBE_DSELI_INTERLEAVE_2);
        e00 = BBE_MOVN_4XCF32_FROMN_2XF32(p00);
        e01 = BBE_MOVN_4XCF32_FROMN_2XF32(p01);
        e10 = BBE_MOVN_4XCF32_FROMN_2XF32(p10);
        e11 = BBE_MOVN_4XCF32_FROMN_2XF32(p11);
#endif

#if 0
        BBE_SVN_4XCF32T_IP(e00, E_w, (BBE_SIMD_WIDTH / 4)*sz_f32c, ble_4_0);
        BBE_SVN_4XCF32T_XP(e01, E_w, L*sz_f32c, ble_4_1);
        BBE_SVN_4XCF32T_IP(e11, E_w, -(int)((BBE_SIMD_WIDTH / 4)*sz_f32c), ble_4_1);
        BBE_SVN_4XCF32T_XP(e10, E_w, (2 * (BBE_SIMD_WIDTH / 4) - L)*sz_f32c, ble_4_0);
#else
        BBE_SVN_4XCF32T_X(e10, E_w, L                      *sz_f32c, ble_4_0);
        BBE_SVN_4XCF32T_X(e11, E_w, (L + (BBE_SIMD_WIDTH / 4))*sz_f32c, ble_4_1);
        BBE_SVN_4XCF32T_I(e01, E_w, (BBE_SIMD_WIDTH / 4)*sz_f32c, ble_4_1);
        BBE_SVN_4XCF32T_IP(e00, E_w, 2 * (BBE_SIMD_WIDTH / 4)*sz_f32c, ble_4_0);
#endif
    }

    if (calcEigenVec){
        H_r = H_w = (xb_vecN_2xf32  *)H;
        E_w = (xb_vecN_4xcf32 *)e;
        P_w = (xb_vecN_2xf32 *)P;

        Q_r  = (xb_vecN_2xf32 *)Q;
        Z_r  = (xb_vecN_2xf32 *)Z;
        B_r  = B;

        for (k = 0; k < L / (BBE_SIMD_WIDTH / 2); k++){
            BBE_LVN_2XF32_XP(h00, H_r, 1 * L*sz_f32);
            BBE_LVN_2XF32_XP(h01, H_r, 1 * L*sz_f32);
            BBE_LVN_2XF32_XP(h10, H_r, 1 * L*sz_f32);
            BBE_LVN_2XF32_XP(h11, H_r, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);

            // Load q, z
            BBE_LVN_2XF32_IP(_q, Q_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
            BBE_LVN_2XF32_IP(_z, Z_r, (BBE_SIMD_WIDTH / 2)*sz_f32);

            BBE_LBN_2_IP(ble, B_r, sizeof(vboolN_2));

            bqmez = BBE_OLEN_2XF32(c0f, _q); // q >= 0.f
            bqmez = BBE_ANDBN_2(bqmez, BBE_NOTBN_2(ble));

            //x = H[(2 * 1 + 0)*L];
            //r = sqrtf(x*x + z*z); p = x / r; q = z / r;
            _x = h10;
            _r = BBE_MULN_2XF32(_x,_x); BBE_MULAN_2XF32(_r, _z, _z);
            _r = BBE_RSQRTN_2XF32(_r);
            _p = BBE_MULN_2XF32(_x, _r);
            _q = BBE_MULN_2XF32(_z, _r);

            /* Row modification */
            _z = h00;
            h00 = BBE_MULN_2XF32(_q, _z ); BBE_MULAN_2XF32(h00, _p, h10);
            h10 = BBE_MULN_2XF32(_q, h10); BBE_MULSN_2XF32(h10, _p, _z);
            _z = h01;
            h01 = BBE_MULN_2XF32(_q, _z ); BBE_MULAN_2XF32(h01, _p, h11);
            h11 = BBE_MULN_2XF32(_q, h11); BBE_MULSN_2XF32(h11, _p, _z);

            /* Column modification */
            _z = h00;
            h00 = BBE_MULN_2XF32(_q, _z ); BBE_MULAN_2XF32(h00, _p, h01);
            h01 = BBE_MULN_2XF32(_q, h01); BBE_MULSN_2XF32(h01, _p, _z);
            _z = h10;
            h10 = BBE_MULN_2XF32(_q, _z); BBE_MULAN_2XF32(h10, _p, h11);
            h11 = BBE_MULN_2XF32(_q, h11); BBE_MULSN_2XF32(h11, _p, _z);

            BBE_SVN_2XF32T_XP(h00, H_w, 1 * L*sz_f32, bqmez);
            BBE_SVN_2XF32T_XP(h01, H_w, 1 * L*sz_f32, bqmez);
            BBE_SVN_2XF32T_XP(h10, H_w, 1 * L*sz_f32, bqmez);
            BBE_SVN_2XF32T_XP(h11, H_w, ((BBE_SIMD_WIDTH / 2) - 3*L)*sz_f32, bqmez);

            /* Accumulate transformations: right-multiply P by the rotator. */
            p00 = BBE_MOVN_2XF32T(_q, c1f , bqmez);
            p01 = BBE_MOVN_2XF32T(BBE_NEGN_2XF32(_p), c0f, bqmez);
            p10 = BBE_MOVN_2XF32T(_p, c0f, bqmez);
            p11 = BBE_MOVN_2XF32T(_q, c1f, bqmez);

            BBE_SVN_2XF32_XP(p00, P_w, 1 * L*sz_f32);
            BBE_SVN_2XF32_XP(p01, P_w, 1 * L*sz_f32);
            BBE_SVN_2XF32_XP(p10, P_w, 1 * L*sz_f32);
            BBE_SVN_2XF32_XP(p11, P_w, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);
        } /* k */
    }
} /* reigen_hqr_2x2_sf() */
#endif


static void reigen_bksubst_2x2_sf(float32_t     * restrict T,
                            const complex_float * restrict e,
                            int L)
#if 0
{
    const float32_t small = FLT_MIN;
    const float32_t big = 1e3f;
    const int N = 2;
    float32_t p, q, r, t;
    float32_t w, wr, wi;
    complex_float c;
    int k, exp;

    for (k = 0; k < L; k++) {
        p = crealf(e[L]); q = cimagf(e[L]);
        if (q == 0.f) {
            /* Real eigenvector */
            T[(2 * 1 + 1)*L] = 1.f;
            w = T[(2 * 0 + 0)*L] - p;
            r = T[(2 * 0 + 1)*L] * T[(2 * 1 + 1)*L];
            wr = crealf(e[0]); wi = cimagf(e[0]);
            T[(2 * 0 + 1)*L] = (w != 0.f ? -r / w : -r / MAX(EPS*fabsf(wr), small));
            t = fabsf(T[(2 * 0 + 1)*L]);
            /* Conditionally rescale the vector to prevent overflow. */
            if (t >= big) {
                exp = ilogbf(t);
                T[(2 * 0 + 1)*L] = ldexpf(T[(2 * 0 + 1)*L], -exp);
                T[(2 * 1 + 1)*L] = ldexpf(T[(2 * 1 + 1)*L], -exp);
            }

            if (wi == 0.f) {
                /* Real eigenvector */
                T[(2 * 0 + 0)*L] = 1.f; T[(2 * 1 + 0)*L] = 0.f;
            }
        } else if (q<0.f) {
            /* Complex vector associated with p-1j*q */
            if (fabsf(T[(2 * 1 + 0)*L]) > fabsf(T[(2 * 0 + 1)*L])) {
                T[(2 * 0 + 0)*L] = -(T[(2 * 1 + 1)*L] - p) / T[(2 * 1 + 0)*L];
                T[(2 * 0 + 1)*L] = -q / T[(2 * 1 + 0)*L];
            } else {
                c = _cdivf(_makecomplexf(-T[(2 * 0 + 1)*L], 0.f),
                    _makecomplexf(T[(2 * 0 + 0)*L] - p, q));
                T[(2 * 0 + 0)*L] = crealf(c); T[(2 * 0 + 1)*L] = cimagf(c);
            }
            T[(2 * 1 + 0)*L] = 1.f; T[(2 * 1 + 1)*L] = 0.f;
        }
        e++; T++;
    }
} /* reigen_bksubst_2x2_sf() */
#else
{
    int k;

#if 0
    const float32_t small = FLT_MIN;
    const float32_t big = 1e3f;
    const int N = 2;
    const complex_float * _e = e;
    float32_t * _T = T;

    float32_t p, q, r, t;
    float32_t w, wr, wi;
    complex_float c;
    int exp;
    for (k = 0; k < L; k++) {
        p = crealf(e[L]); q = cimagf(e[L]);
        if (q<0.f){
            /* Complex vector associated with p-1j*q */
            if (fabsf(T[(2 * 1 + 0)*L]) > fabsf(T[(2 * 0 + 1)*L])) {
                T[(2 * 0 + 0)*L] = -(T[(2 * 1 + 1)*L] - p) / T[(2 * 1 + 0)*L];
                T[(2 * 0 + 1)*L] = -q / T[(2 * 1 + 0)*L];
            } else {
                c = _cdivf(_makecomplexf(-T[(2 * 0 + 1)*L], 0.f),
                    _makecomplexf(T[(2 * 0 + 0)*L] - p, q));
                T[(2 * 0 + 0)*L] = crealf(c); T[(2 * 0 + 1)*L] = cimagf(c);
            }
            T[(2 * 1 + 0)*L] = 1.f; T[(2 * 1 + 1)*L] = 0.f;
        }
        e++; T++;
    }
    e = _e;
    T = _T;

    for (k = 0; k < L; k++) {
        p = crealf(e[L]); q = cimagf(e[L]);
        if (q == 0.f) {
            /* Real eigenvector */
            T[(2 * 1 + 1)*L] = 1.f;
            w = T[(2 * 0 + 0)*L] - p;
            r = T[(2 * 0 + 1)*L];
            wr = crealf(e[0]); wi = cimagf(e[0]);
            T[(2 * 0 + 1)*L] = (w != 0.f ? -r / w : -r / MAX(EPS*fabsf(wr), small));
            t = fabsf(T[(2 * 0 + 1)*L]);
            /* Conditionally rescale the vector to prevent overflow. */
            if (t >= big) {
                exp = ilogbf(t);
                T[(2 * 0 + 1)*L] = ldexpf(T[(2 * 0 + 1)*L], -exp);
                T[(2 * 1 + 1)*L] = ldexpf(1.f, -exp);
            }

            if (wi == 0.f) {
                /* Real eigenvector */
                T[(2 * 0 + 0)*L] = 1.f; T[(2 * 1 + 0)*L] = 0.f;
            }
        }
        e++; T++;
    }
#elif 1
    const xb_vecN_2xf32 * E_r;
    const xb_vecN_2xf32 * restrict T_r;
          xb_vecN_2xf32 * restrict T_w;

    const xb_vecN_2xf32 c0f = BBE_CONSTN_2XF32(0);
    const xb_vecN_2xf32 c1f = BBE_CONSTN_2XF32(1);
    const xb_vecNx16 exp_mask = BBE_MOVVINX16(BBE_FLOAT_MASK_EXP);
    const xb_vecN_2xf32 _small = FLT_MIN;
    const xb_vecN_2xf32 _big = 1e3f;
    const xb_vecN_2xf32 _EPS = EPS;

    xb_vecN_2xf32 _p, _q;
    xb_vecN_2xf32 e0, e1;
    xb_vecN_2xf32 _x, _y, _c0, _c1;
    xb_vecN_2xf32 t00, t01, t10, t11;
    xb_vecN_2xf32 d0, d1;
    xb_vecN_2xf32 n0, n1;
    xb_vecN_2xf32 _wr, _wi;
    xb_vecN_2xf32 _w, _r, _t;
    xb_vecN_2xf32 exp;

    vboolN_2 blz, blt;
    vboolN_2 beqz, ble, bwieqz;

    E_r = (xb_vecN_2xf32 *)&e[L];
    T_r = T_w = (xb_vecN_2xf32 *)T;
#if 0
    for (k = 0; k < L / (BBE_SIMD_WIDTH / 2); k++) {
        xb_vecN_2xf32 _t00, _t01;
        //p = crealf(e[L]); q = cimagf(e[L]);
        e1 = BBE_LVN_2XF32_I(E_r, (BBE_SIMD_WIDTH / 4)*sz_f32c);
        BBE_LVN_2XF32_XP(e0, E_r, 2 * (BBE_SIMD_WIDTH / 4)*sz_f32c);

        _p = BBE_SELN_2XF32I(e1, e0, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
        _q = BBE_SELN_2XF32I(e1, e0, BBE_SELI_EXTRACT_2_OF_4_OFF_2);

        blz = BBE_OLTN_2XF32(_q, c0f); // q<0.f
                
        t01 = BBE_LVN_2XF32_X(T_r, 1 * L*sz_f32);
        t10 = BBE_LVN_2XF32_X(T_r, 2 * L*sz_f32);
        t11 = BBE_LVN_2XF32_X(T_r, 3 * L*sz_f32);
        BBE_LVN_2XF32_IP(t00, T_r, (BBE_SIMD_WIDTH / 2)*sz_f32);

        //fabsf(T[(2 * 1 + 0)*L]) > fabsf(T[(2 * 0 + 1)*L])
        _x = BBE_ABSN_2XF32(t10); _y = BBE_ABSN_2XF32(t01);
        blt = BBE_OLTN_2XF32(_y, _x);

        //c = _cdivf(_makecomplexf(-T[(2 * 0 + 1)*L], 0.f),
        //    _makecomplexf(T[(2 * 0 + 0)*L] - p, q));
        //T[(2 * 0 + 0)*L] = crealf(c); T[(2 * 0 + 1)*L] = cimagf(c);
        _x = BBE_NEGN_2XF32(t01);
        n0 = BBE_SELN_2XF32I(c0f, _x, BBE_SELI_INTERLEAVE_2_LO);
        n1 = BBE_SELN_2XF32I(c0f, _x, BBE_SELI_INTERLEAVE_2_HI);
        _y = BBE_SUBN_2XF32(t00, _p);
        d0 = BBE_SELN_2XF32I(_q, _y, BBE_SELI_INTERLEAVE_2_LO);
        d1 = BBE_SELN_2XF32I(_q, _y, BBE_SELI_INTERLEAVE_2_HI);
        _c0 = IT_CDIVN_2XF32(n0, d0, 1);
        _c1 = IT_CDIVN_2XF32(n1, d1, 1);
        _t00 = BBE_SELN_2XF32I(_c1, _c0, BBE_SELI_EXTRACT_2_OF_4_OFF_0); // real 
        _t01 = BBE_SELN_2XF32I(_c1, _c0, BBE_SELI_EXTRACT_2_OF_4_OFF_2); // img

        // T[(2 * 0 + 0)*L] = -(T[(2 * 1 + 1)*L] - p) / T[(2 * 1 + 0)*L];
        // T[(2 * 0 + 1)*L] = -q / T[(2 * 1 + 0)*L];
        _x = BBE_SUBN_2XF32(t11, _p);
        _x = BBE_NEGN_2XF32(IT_FDIVN_2XF32(_x, t10, 1));
        _y = BBE_NEGN_2XF32(IT_FDIVN_2XF32(_q, t10, 1));
        t00 = BBE_MOVN_2XF32T(_x, _t00, blt);
        t01 = BBE_MOVN_2XF32T(_y, _t01, blt);

        //T[(2 * 1 + 0)*L] = 1.f; T[(2 * 1 + 1)*L] = 0.f;
        t10 = BBE_CONSTN_2XF32(1); t11 = c0f;        
        BBE_SVN_2XF32T_X(t01, T_w, 1 * L*sz_f32, blz);
        BBE_SVN_2XF32T_X(t10, T_w, 2 * L*sz_f32, blz);
        BBE_SVN_2XF32T_X(t11, T_w, 3 * L*sz_f32, blz);
        BBE_SVN_2XF32T_IP(t00, T_w, (BBE_SIMD_WIDTH / 2)*sz_f32, blz);
    }
#else
    for (k = 0; k < L / (BBE_SIMD_WIDTH / 2); k++) {
        xb_vecN_2xf32 nr, ni;
        xb_vecN_2xf32 dr, di;

        //p = crealf(e[L]); q = cimagf(e[L]);
        e1 = BBE_LVN_2XF32_I(E_r, (BBE_SIMD_WIDTH / 4)*sz_f32c);
        BBE_LVN_2XF32_XP(e0, E_r, 2 * (BBE_SIMD_WIDTH / 4)*sz_f32c);

        _p = BBE_SELN_2XF32I(e1, e0, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
        _q = BBE_SELN_2XF32I(e1, e0, BBE_SELI_EXTRACT_2_OF_4_OFF_2);
        blz = BBE_OLTN_2XF32(_q, c0f); // q<0.f

        t01 = BBE_LVN_2XF32_X(T_r, 1 * L*sz_f32);
        t10 = BBE_LVN_2XF32_X(T_r, 2 * L*sz_f32);
        t11 = BBE_LVN_2XF32_X(T_r, 3 * L*sz_f32);
        BBE_LVN_2XF32_IP(t00, T_r, (BBE_SIMD_WIDTH / 2)*sz_f32);

        //fabsf(T[(2 * 1 + 0)*L]) > fabsf(T[(2 * 0 + 1)*L])
        _x = BBE_ABSN_2XF32(t10); _y = BBE_ABSN_2XF32(t01);
        blt = BBE_OLTN_2XF32(_y, _x);

        nr = BBE_MOVN_2XF32T(BBE_SUBN_2XF32(_p, t11), BBE_NEGN_2XF32(t01), blt);
        ni = BBE_MOVN_2XF32T(BBE_NEGN_2XF32(_q)     , c0f                , blt);
        dr = BBE_MOVN_2XF32T(t10, BBE_SUBN_2XF32(t00, _p), blt);
        di = BBE_MOVN_2XF32T(c0f, _q                     , blt);

        n0 = BBE_SELN_2XF32I(ni, nr, BBE_SELI_INTERLEAVE_2_LO);
        d0 = BBE_SELN_2XF32I(di, dr, BBE_SELI_INTERLEAVE_2_LO);
        n1 = BBE_SELN_2XF32I(ni, nr, BBE_SELI_INTERLEAVE_2_HI);
        d1 = BBE_SELN_2XF32I(di, dr, BBE_SELI_INTERLEAVE_2_HI);
        _c0 = IT_CDIVN_2XF32(n0, d0, 1);
        _c1 = IT_CDIVN_2XF32(n1, d1, 1);

        t00 = BBE_SELN_2XF32I(_c1, _c0, BBE_SELI_EXTRACT_2_OF_4_OFF_0); // real 
        t01 = BBE_SELN_2XF32I(_c1, _c0, BBE_SELI_EXTRACT_2_OF_4_OFF_2); // img

        //T[(2 * 1 + 0)*L] = 1.f; T[(2 * 1 + 1)*L] = 0.f;
        t10 = BBE_CONSTN_2XF32(1); t11 = c0f;
        BBE_SVN_2XF32T_X(t01, T_w, 1 * L*sz_f32, blz);
        BBE_SVN_2XF32T_X(t10, T_w, 2 * L*sz_f32, blz);
        BBE_SVN_2XF32T_X(t11, T_w, 3 * L*sz_f32, blz);
        BBE_SVN_2XF32T_IP(t00, T_w, (BBE_SIMD_WIDTH / 2)*sz_f32, blz);
    }
#endif

    E_r = (xb_vecN_2xf32 *)e;
    T_r = T_w = (xb_vecN_2xf32 *)T;

    for (k = 0; k < L / (BBE_SIMD_WIDTH / 2); k++) {
        // p = crealf(e[L]); q = cimagf(e[L]);
        e0 = BBE_LVN_2XF32_X(E_r, L*sz_f32c);
        e1 = BBE_LVN_2XF32_X(E_r, L*sz_f32c + (BBE_SIMD_WIDTH / 4)*sz_f32c);
        _p = BBE_SELN_2XF32I(e1, e0, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
        _q = BBE_SELN_2XF32I(e1, e0, BBE_SELI_EXTRACT_2_OF_4_OFF_2);
        // wr = crealf(e[0]); wi = cimagf(e[0]);
        e1 = BBE_LVN_2XF32_I(E_r, (BBE_SIMD_WIDTH / 4)*sz_f32c);
        BBE_LVN_2XF32_XP(e0, E_r, 2 * (BBE_SIMD_WIDTH / 4)*sz_f32c);
        _wr = BBE_SELN_2XF32I(e1, e0, BBE_SELI_EXTRACT_2_OF_4_OFF_0);
        _wi = BBE_SELN_2XF32I(e1, e0, BBE_SELI_EXTRACT_2_OF_4_OFF_2);

        beqz = BBE_OEQN_2XF32(_q, c0f); // q==0.f

        t01 = BBE_LVN_2XF32_X(T_r, 1 * L*sz_f32);
        t10 = BBE_LVN_2XF32_X(T_r, 2 * L*sz_f32);
        t11 = BBE_LVN_2XF32_X(T_r, 3 * L*sz_f32);
        BBE_LVN_2XF32_IP(t00, T_r, (BBE_SIMD_WIDTH / 2)*sz_f32);
        //T[(2 * 1 + 1)*L] = 1.f;
        //w = T[(2 * 0 + 0)*L] - p;
        //r = T[(2 * 0 + 1)*L];
        t11 = c1f;
        _w = BBE_SUBN_2XF32(t00, _p);
        _r = t01;
        //T[(2 * 0 + 1)*L] = (w != 0.f ? -r / w : -r / MAX(EPS*fabsf(wr), small));
        _x = BBE_MAXN_2XF32(BBE_MULN_2XF32(_EPS, BBE_ABSN_2XF32(_wr)), _small);
        d0 = BBE_MOVN_2XF32T(_w, _x, BBE_NOTBN_2(BBE_UEQN_2XF32(c0f, _w)));
        t01 = BBE_NEGN_2XF32(IT_FDIVN_2XF32(_r, d0, 1));
        //t = fabsf(T[(2 * 0 + 1)*L]);
        _t = BBE_ABSN_2XF32(t01);
        // t >= big
        ble = BBE_OLEN_2XF32(_big, _t);
        // exp = ilogbf(t);
        //T[(2 * 0 + 1)*L] = ldexpf(T[(2 * 0 + 1)*L], -exp);
        //T[(2 * 1 + 1)*L] = ldexpf(1.f, -exp);
        exp = BBE_RECIP0N_2XF32(BBE_MOVN_2XF32_FROMNX16(BBE_ANDNX16(BBE_MOVNX16_FROMN_2XF32(_t), exp_mask)));
        BBE_MULN_2XF32T(t01, t01, exp, ble);
        //BBE_MULN_2XF32T(t11, t11, exp, ble); 
        t11 = BBE_MOVN_2XF32T(exp, t11, ble);
        //wi == 0.f
        bwieqz = BBE_OEQN_2XF32(_wi, c0f);
        //T[(2 * 0 + 0)*L] = 1.f; T[(2 * 1 + 0)*L] = 0.f;
        t00 = BBE_MOVN_2XF32T(c1f, t00, bwieqz);
        t10 = BBE_MOVN_2XF32T(c0f, t10, bwieqz);

        BBE_SVN_2XF32T_X(t01, T_w, 1 * L*sz_f32, beqz);
        BBE_SVN_2XF32T_X(t10, T_w, 2 * L*sz_f32, beqz);
        BBE_SVN_2XF32T_X(t11, T_w, 3 * L*sz_f32, beqz);
        BBE_SVN_2XF32T_IP(t00, T_w, (BBE_SIMD_WIDTH / 2)*sz_f32, beqz);
    }
#endif
} /* reigen_bksubst_2x2_sf() */
#endif

static void reigen_hmulp_2x2sf(float32_t * restrict U, /* U[N*N][L]         */
                         const float32_t * restrict P, /* P[N*N][L]         */
                         const float32_t * restrict H, /* H[N*(N+3)/2-1][L] */
                         int L){
    int k;

    const xb_vecN_2xf32 * P_r;
    const xb_vecN_2xf32 * H_r;
          xb_vecN_2xf32 * U_w;

    NASSERT_ALIGN(U, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(H, 2 * BBE_SIMD_WIDTH);
    NASSERT(0 == (L % (BBE_SIMD_WIDTH / 2)));

    P_r = (xb_vecN_2xf32*)P;
    H_r = (xb_vecN_2xf32*)H;
    U_w = (xb_vecN_2xf32*)U;

    for (k = 0; k<L / (BBE_SIMD_WIDTH / 2); k++) {
        xb_vecN_2xf32 s;
        xb_vecN_2xf32 p00, p01, p10, p11;
        xb_vecN_2xf32 h00, h01, h10, h11;

        BBE_LVN_2XF32_XP(p00, P_r, L*sz_f32);
        BBE_LVN_2XF32_XP(p01, P_r, L*sz_f32);
        BBE_LVN_2XF32_XP(p10, P_r, L*sz_f32);
        BBE_LVN_2XF32_XP(p11, P_r, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);
        BBE_LVN_2XF32_XP(h00, H_r, L*sz_f32);
        BBE_LVN_2XF32_XP(h01, H_r, L*sz_f32);
        BBE_LVN_2XF32_XP(h10, H_r, L*sz_f32);
        BBE_LVN_2XF32_XP(h11, H_r, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);
        
        s = BBE_MULN_2XF32(p00, h00);
        BBE_MULAN_2XF32(s, p01, h10);
        BBE_SVN_2XF32_XP(s, U_w, L*sz_f32);
        s = BBE_MULN_2XF32(p00, h01);
        BBE_MULAN_2XF32(s, p01, h11);
        BBE_SVN_2XF32_XP(s, U_w, L*sz_f32);
        s = BBE_MULN_2XF32(p10, h00);
        BBE_MULAN_2XF32(s, p11, h10);
        BBE_SVN_2XF32_XP(s, U_w, L*sz_f32);
        s = BBE_MULN_2XF32(p10, h01);
        BBE_MULAN_2XF32(s, p11, h11);
        BBE_SVN_2XF32_XP(s, U_w, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);
    }
} /* reigen_hmulp_2x2sf() */

static void reigen_evcomb_2x2sf(complex_float * restrict V, /* V[N*N][L] */
                          const float32_t     * restrict U, /* U[N*N][L] */
                          const complex_float * restrict e, /* e[N][L]   */
                          int L)
#if 0
{
    const int N = 2;
    float32_t f, g, s;
    int i, k;
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(U, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(e, 2 * BBE_SIMD_WIDTH);
    NASSERT(0 == (L % (BBE_SIMD_WIDTH / 2)));

    for (k = 0; k<L; k++) {
        if (cimagf(e[0])>0.f){
            // First iteration cimagf(e[0])>0.f
            for (s = 0.f, i = 0; i < N; i++) {
                f = U[(i*N + 0)*L];
                g = U[(i*N + 1)*L];
                s += f*f + g*g;
            }
            s = sqrtf(s);
            for (i = 0; i < N; i++) {
                f = U[(i*N + 0)*L];
                g = U[(i*N + 1)*L];
                V[i*N*L] = _makecomplexf(f / s, g / s);
            }
            // Second iteration cimagf(e[0])<0.f
            for (i = 0; i < N; i++) {
                f = U[(i*N + 0)*L];
                g = U[(i*N + 1)*L];
                V[(i*N + 1)*L] = _makecomplexf(f / s, -g / s);
            }
        } else { // cimagf(e[0]) == 0
            // First iteration
            for (s = 0.f, i = 0; i < N; i++) {
                f = U[(i*N)*L];
                s += f*f;
            }
            s = sqrtf(s);
            for (i = 0; i < N; i++) {
                f = U[i*N*L];
                V[i*N*L] = _makecomplexf(f / s, 0.f);
            }
            // Second iteration cimagf(e[0]) == 0
            for (s = 0.f, i = 0; i < N; i++) {
                f = U[(i*N + 1)*L];
                s += f*f;
            }
            s = sqrtf(s);
            for (i = 0; i < N; i++) {
                f = U[(i*N + 1)*L];
                V[(i*N + 1)*L] = _makecomplexf(f / s, 0.f);
            }
        }
        V++; U++; e++;
    }
} /* reigen_evcomb_2x2sf() */
#elif 1
{
    int k;

    const xb_vecN_2xf32  * U_r;
    const xb_vecN_4xcf32 * E_r;
          xb_vecN_4xcf32 * V_w;

    const xb_vecN_2xf32 c0f = BBE_CONSTN_2XF32(0);

    xb_vecN_2xf32  h0, h1, f, _f, g, s, _s;
    xb_vecN_4xcf32 e0, e1, v0, v1;
    vboolN_2 bmoz;

    s = BBE_ZERON_2XF32();

    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(U, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(e, 2 * BBE_SIMD_WIDTH);
    NASSERT(0 == (L % (BBE_SIMD_WIDTH / 2)));
   
    E_r = (xb_vecN_4xcf32 *)e;
    U_r = (xb_vecN_2xf32  *)U;
    V_w = (xb_vecN_4xcf32 *)V;

    for (k = 0; k<L / (BBE_SIMD_WIDTH / 2); k++) {
        e1 = BBE_LVN_4XCF32_I(E_r, (BBE_SIMD_WIDTH / 4)*sz_f32c);
        BBE_LVN_4XCF32_XP(e0, E_r, (BBE_SIMD_WIDTH / 2)*sz_f32c);
        h0 = BBE_MOVN_2XF32_FROMN_4XCF32(e0); h1 = BBE_MOVN_2XF32_FROMN_4XCF32(e1);
        h0 = BBE_SELN_2XF32I(h1, h0, BBE_SELI_EXTRACT_2_OF_4_OFF_2);
        // cimagf(e[j*stride]) > 0.f
        bmoz = BBE_OLTN_2XF32(c0f, h0);

        f = BBE_LVN_2XF32_I(U_r, 0);
        g = BBE_LVN_2XF32T_X(U_r, L*sz_f32, bmoz);
        s = BBE_MULN_2XF32(f, f);
        BBE_MULAN_2XF32(s, g, g);
        f = BBE_LVN_2XF32_X (U_r, 2*L*sz_f32);
        g = BBE_LVN_2XF32T_X(U_r, 3*L*sz_f32, bmoz);
        BBE_MULAN_2XF32(s, f, f);
        BBE_MULAN_2XF32(s, g, g);
        s = BBE_RSQRTN_2XF32(s);

        // Update V00
        f = BBE_LVN_2XF32_I(U_r, 0);
        g = BBE_LVN_2XF32T_X(U_r, L*sz_f32, bmoz);
        f = BBE_MULN_2XF32(f, s); g = BBE_MULN_2XF32(g, s);
        h1 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_HI);
        h0 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_LO);
        v0 = BBE_MOVN_4XCF32_FROMN_2XF32(h0);
        v1 = BBE_MOVN_4XCF32_FROMN_2XF32(h1);
        BBE_SVN_4XCF32_IP(v0, V_w, (BBE_SIMD_WIDTH / 4) * sz_f32c);
        BBE_SVN_4XCF32_XP(v1, V_w, 2 * L*sz_f32c);
        // Update V01
        f = BBE_LVN_2XF32_X(U_r, 2 * L*sz_f32);
        g = BBE_LVN_2XF32T_X(U_r, 3 * L*sz_f32, bmoz);
        f = BBE_MULN_2XF32(f, s); g = BBE_MULN_2XF32(g, s);
        h1 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_HI);
        h0 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_LO);
        v0 = BBE_MOVN_4XCF32_FROMN_2XF32(h0);
        v1 = BBE_MOVN_4XCF32_FROMN_2XF32(h1);
        BBE_SVN_4XCF32_IP(v1, V_w, -(int)((BBE_SIMD_WIDTH / 4) * sz_f32c));
        BBE_SVN_4XCF32_XP(v0, V_w, -L*sz_f32c);

        // Update second column
        f = BBE_LVN_2XF32_X(U_r, 1 * L*sz_f32);
        _s = BBE_MULN_2XF32(f, f);
        f = BBE_LVN_2XF32_X(U_r, 3 * L*sz_f32);
        BBE_MULAN_2XF32(_s, f, f);
        _s = BBE_RSQRTN_2XF32(_s);

        s = BBE_MOVN_2XF32T(s, _s, bmoz);
        
        // Update V10
        BBE_LVN_2XF32_XP(f, U_r, L*sz_f32);
        BBE_LVN_2XF32_XP(_f, U_r, L*sz_f32);
        g = BBE_MOVN_2XF32T(_f, c0f, bmoz); f = BBE_MOVN_2XF32T(f, _f, bmoz);
        f = BBE_MULN_2XF32(f, s); g = BBE_MULN_2XF32(g, s);
        BBE_NEGN_2XF32T(g, g, bmoz);
        h1 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_HI);
        h0 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_LO);
        v0 = BBE_MOVN_4XCF32_FROMN_2XF32(h0);
        v1 = BBE_MOVN_4XCF32_FROMN_2XF32(h1);
        BBE_SVN_4XCF32_IP(v0, V_w, (BBE_SIMD_WIDTH / 4) * sz_f32c);
        BBE_SVN_4XCF32_XP(v1, V_w, 2 * L*sz_f32c);
        // Update V11
        BBE_LVN_2XF32_XP(f, U_r, L*sz_f32);
        BBE_LVN_2XF32_XP(_f, U_r, ((BBE_SIMD_WIDTH / 2) - 3 * L)*sz_f32);
        g = BBE_MOVN_2XF32T(_f, c0f, bmoz); f = BBE_MOVN_2XF32T(f, _f, bmoz);
        f = BBE_MULN_2XF32(f, s); g = BBE_MULN_2XF32(g, s);
        BBE_NEGN_2XF32T(g, g, bmoz);
        h1 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_HI);
        h0 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_LO);
        v0 = BBE_MOVN_4XCF32_FROMN_2XF32(h0);
        v1 = BBE_MOVN_4XCF32_FROMN_2XF32(h1);
        BBE_SVN_4XCF32_IP(v1, V_w, -(int)((BBE_SIMD_WIDTH / 4) * sz_f32c));
        BBE_SVN_4XCF32_XP(v0, V_w, (2*(BBE_SIMD_WIDTH / 4) - 3 * L)*sz_f32c);
    }
} /* reigen_evcomb_2x2sf() */
#else
{}
#endif
#endif

/*-------------------------------------------------------------------------
Eigenvalues And Eigenvectors Of Real/Complex Stream Ordered Matrices

Description: for each complex/real input matrix A of size NxN, compute N
(possibly repeated) eigenvalues s[N], and (optonally) N right eigenvectors
of size Nx1 V[N]. Input and output data are stored in stream order.

Data format: IEEE-754 Std single precision floating-point

Notes:
1. Functions may perform in-place transformations of input matrices, so that
   INPUT DATA MAY APPEAR DAMAGED after the call.
2. Once the eigenvectors are not required, set the corresponding output pointer
   V to zero, so that a lower complexity algorithm will be used.
3. Floating-point functions assume that input data are reasonably scaled. That
   is, the base-2 exponent e of the maximum absolute value over an input matrix
   belongs to the range -E<e<E, where E = 63-log2(N)/2.
4. In order to reduce the computational complexity, a preprocessing step known
   as "matrix balancing" is omitted from the implementation.

Temporary:
  pScr        Scratch area. Required size (in bytes) is defined by 
              functions [r]eigen<size>sf_getScratchSize(N,L)
Input:
  N           Matrix size
  L           Number of matrices
  A[N*N][L]   NxN input matrices
Output:
  e[N][L]     Nx1 vectors of eigenvalues. In an exceptional case when the
              iterative algorithm fails to converge for a particular matrix,
              all elements of the respective vector are set to NaN.
  V[N*N][L]   NxN matrices comprised of N column eigenvectors (optional)
Restrictions:
  pScr,e,V,A  Must not overlap and must be aligned on 32-byte boundary 
  N           N>1
  L           Must be a multiple of 8 for real data, or a multiple of 4 for
              complex data
---------------------------------------------------------------------------*/
#if 0
void reigen2x2sf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            float32_t     * restrict A,
            int L )
{
  NASSERT_ALIGN( pScr, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( e   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A   , 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );

  reigennxnsf(pScr,e,V,A,2,L);

} /* reigen2x2sf() */

size_t reigen2x2sf_getScratchSize ( int N, int L )
{
  NASSERT( 2==N );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
  return ( reigennxnsf_getScratchSize(N,L) );
}
#else
void reigen2x2sf(
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            float32_t     * restrict A,
            int L )
{
    const int N = 2;
    float32_t *Ps, *Us;
    //int p;

    NASSERT_ALIGN( pScr, 2*BBE_SIMD_WIDTH );
    NASSERT_ALIGN( e   , 2*BBE_SIMD_WIDTH );
    NASSERT_ALIGN( V   , 2*BBE_SIMD_WIDTH );
    NASSERT_ALIGN( A   , 2*BBE_SIMD_WIDTH );
    NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );

    {
        void * p = pScr;
        /* Partition the scratch area. */
        Ps = (float32_t    *)p; p = Ps + L*N*N;
        Us = (float32_t    *)p; p = Us + L*N*N;
        /* Make sure that scratch arrays fit into the reserved space. */
        NASSERT((uint8_t*)p - (uint8_t*)pScr <= (int)reigen2x2sf_getScratchSize(N, L));
    }

    if (V) {
#if 0
        /* Initialize transformation matrices with 2x2 identity matrix. */
        for (p = 0; p < L; p++) {
            Ps[0 * L + p] = Ps[3 * L + p] = 1.f;
            Ps[1 * L + p] = Ps[2 * L + p] = 0.f;
        }
#endif
        reigen_hqr_2x2_sf(e, A, Ps, L, 1);
        reigen_bksubst_2x2_sf(A, e, L);
        /* Left-multiply eigenvectors by transformation matrices to obtain re/im
        * components of eigenvectors for original input matrices. */
        reigen_hmulp_2x2sf(Us, Ps, A, L);
        /* Combine re/im components and rescale eigenvectors. */
        reigen_evcomb_2x2sf(V, Us, e, L);
    } else {
        reigen_hqr_2x2_sf(e, A, Ps, L, 0);
    }
} /* reigen2x2sf() */

size_t reigen2x2sf_getScratchSize(int N, int L)
{
    NASSERT(2 == N);
    NASSERT(0 == (L % (BBE_SIMD_WIDTH / 2)));
    return ( L*N*N*sz_f32 +  /* Ps: transformation matrix in stream order */
             L*N*N*sz_f32);  /* Us: Real matrices of re/im components of eigenvectors in stream order */
}

#endif

#else /* HAVE_VFPU */

DISCARD_FUN( void, reigen2x2sf, ( void * pScr,
                         complex_float * restrict e,
                         complex_float * restrict V,
                         float32_t     * restrict A,
                         int L ) )

size_t reigen2x2sf_getScratchSize ( int N, int L )
{
  NASSERT( 2==N );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
  return (0);
}

#endif /* HAVE_VFPU */
