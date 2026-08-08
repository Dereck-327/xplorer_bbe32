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
    Complex 2x2 stream ordered matrices
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

#define sz_f32c  sizeof(complex_float)
#define EPS       FLT_EPSILON
#define MAX(a,b)    ( (a)>(b) ? (a) : (b) )

#if 0

#include <math.h>
#if !defined(COMPILER_MSVC)
#include <complex.h>
#endif

static complex_float _makecomplexf(float32_t re, float32_t im)
{
    union { float32_t r[2]; complex_float c; } u = { { re, im } };
    return (u.c);
}

/* Complex floating-point negation, single precision. */
static complex_float cnegf(complex_float x)
{
    return (_makecomplexf(-crealf(x), -cimagf(x)));
}

/* Complex floating-point conjugate, single precision. */
static complex_float _conjf(complex_float x)
{
    return (_makecomplexf(crealf(x), -cimagf(x)));
}

/* Cheap estimation of complex number's magnitude, single precision
* floating-point. */
static float32_t cmagf(complex_float x) 
{
    return (fabsf(crealf(x)) + fabsf(cimagf(x)));
}

/* Squared absolute of a complex floating-point number, single precision. */
static float32_t cabs2f(complex_float x)
{
    return (crealf(x)*crealf(x) + cimagf(x)*cimagf(x));
}

/* Complex floating-point multiplication, single precision */
static complex_float cmulf(complex_float x, complex_float y)
{
    return (_makecomplexf(crealf(x)*crealf(y) - cimagf(x)*cimagf(y),
        cimagf(x)*crealf(y) + crealf(x)*cimagf(y)));
}

/* Complex floating-point multiplication with conjugation of the second
* argument, single precision */
static complex_float cmuljf(complex_float x, complex_float y)
{
    return (_makecomplexf(crealf(x)*crealf(y) + cimagf(x)*cimagf(y),
        cimagf(x)*crealf(y) - crealf(x)*cimagf(y)));
}

/* Real by complex floating-point multiplication, single precision. */
static complex_float rcmulf(float32_t x, complex_float y)
{
    return (_makecomplexf(x*crealf(y), x*cimagf(y)));
}

/* Complex floating-point addition, single precision */
static complex_float caddf(complex_float x, complex_float y)
{
    return (_makecomplexf(crealf(x) + crealf(y), cimagf(x) + cimagf(y)));
}

/* Complex floating-point subtraction, single precision */
static complex_float csubf(complex_float x, complex_float y)
{
    return (_makecomplexf(crealf(x) - crealf(y), cimagf(x) - cimagf(y)));
}

/* Multiply a complex floating-point number by an integral power
* of two, single precision. */
static complex_float cldexpf(complex_float x, int y)
{
    return (_makecomplexf(ldexpf(crealf(x), y), ldexpf(cimagf(x), y)));
}

/* Complex floating-point absolute value, single precision. */
static float32_t _cabsf(complex_float x)
{
    /*
    * Based on cabs code from "Similarity Reduction of a General Matrix to
    * Hessenberg Form" by R.S. Martin and J.H. Wilkinson, Handbook for
    * Automatie Computation, Vol.II Linear Algebra, Contribution II/13.
    * MATLAB reference code:
    *   function r = cabs(x)
    *   xr = abs(real(x));
    *   xi = abs(imag(x));
    *   if xi>0
    *     if xi>xr, t = xr; xr = xi; xi = t; end;
    *     if xi==xr, t = 1; else t = xi/xr; end;
    *     r = xr*sqrt(1+t^2);
    *   else
    *     r = xr;
    *   end
    */

    float32_t t, xr, xi;
    xr = fabsf(crealf(x));
    xi = fabsf(cimagf(x));
    if (xi>0) {
        if (xi>xr) { t = xr; xr = xi; xi = t; }
        t = (xi == xr ? 1 : xi / xr);  /* Avoid Inf/Inf */
        return (xr*sqrtf(1 + t*t));
    }
    else {
        return (xi == xi ? xr : xi);
    }

} /* _cabsf() */

/* Complex by real floating-point division, single precision. */
static complex_float crdivf(complex_float x, float32_t y)
{
    return (_makecomplexf(crealf(x) / y, cimagf(x) / y));
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

/* Complex floating-point square root, single precision. Returns the
* square root that lies in the right half of the complex plane. */
static complex_float _csqrtf(complex_float x)
{
    /*
    * Based on csqrt code from "Similarity Reduction of a General Matrix to
    * Hessenberg Form" by R.S. Martin and J.H. Wilkinson, Handbook for
    * Automatie Computation, Vol.II Linear Algebra, Contribution II/13.
    * MATLAB reference code:
    *   function z = csqrt(x)
    *   xr = real(x); xi = imag(x);
    *   h = sqrt((abs(xr)+cabs(x))/2);
    *   if xi~=0, xi = xi/(2*h); end;
    *   if xr<0
    *     if xi<0
    *       xr = -xi; xi = -h;
    *     else
    *       xr = xi; xi = h;
    *     end
    *   else
    *     xr = h;
    *   end
    *   z = complex(xr,xi);
    */

    float32_t xr, xi, h;
    xr = crealf(x); xi = cimagf(x);
    h = sqrtf((fabsf(xr) + _cabsf(x)) / 2);
    if (xi != 0) xi /= 2 * h;
    if (xr<0) {
        xr = (xi<0 ? -xi : xi);
        xi = (xi<0 ? -h : h);
    }
    else {
        xr = h;
    }
    return (_makecomplexf(xr, xi));

} /* _csqrtf() */
#endif

static void eigen_hqr_2x2sf(complex_float * restrict e,
    complex_float * restrict H,
    complex_float * restrict P,
    int L, int calcEigenVec)
#if 0
{
    /*
    * This implementation is based on the following articles:
    * [1] "The QR Algorithm for Real Hessenberg Matrices" by R.S. Martin,
    *     G. Petern and J.H. Wilkinson, Handbook for Automatie Computation,
    *     Vol.II Linear Algebra, Contribution II/14.
    * [2] "Eigenvectors of Real and Complex Matrices by LR and QR
    *     triangularizations" by G. Petern and J.H. Wilkinson, Handbook for
    *     Automatie Computation, Vol.II Linear Algebra, Contribution II/15.
    */
    int k;
    complex_float p, q, x, y, z, w;
    float32_t s;

    for (k = 0; k < L; k++){
        x = H[(2 * 1 + 0)*L];
        y = H[(2 * 0 + 0)*L];
        z = H[(2 * 1 + 1)*L];

        if (cmagf(x) <= EPS*(cmagf(y) + cmagf(z))){
            /* Revealed a standalone eigenvalue. */
            e[(0)*L] = H[(2 * 0 + 0)*L];
            e[(1)*L] = H[(2 * 1 + 1)*L];
        } else {
            /* Revealed a pair of eigenvalues. */
            x = H[(2 * 1 + 1)*L];
            y = H[(2 * 0 + 0)*L];
            w = cmulf(H[(2 * 0 + 1)*L], H[(2 * 1 + 0)*L]);
            p = rcmulf(.5f, csubf(y, x));
            z = _csqrtf(caddf(cmulf(p, p), w));
            if (crealf(p)*crealf(z) + cimagf(p)*cimagf(z) < 0) {
                z = csubf(p, z);
            } else {
                z = caddf(p, z);
            }
            if (cmagf(w)>0) {
                /* Quadratic divisor */
                e[(0)*L] = caddf(x, z);
                e[(1)*L] = csubf(x, _cdivf(w, z));
            } else {
                /* Two linear divisors */
                e[(0)*L] = y; e[(1)*L] = x;
            }
            if (P) {
                /* Annihilate T(en,en-1) entry by similarity transformation involving
                * a rotation in the en,en-1 plane. */
                x = H[(2 * 1 + 0)*L]; s = sqrtf(cabs2f(x) + cabs2f(z));
                p = crdivf(x, s); q = crdivf(z, s);
                /* Row modification */
                z = H[(2 * 0 + 0)*L];
                H[(2 * 0 + 0)*L] = caddf(cmuljf(z, q), cmuljf(H[(2 * 1 + 0)*L], p));
                H[(2 * 1 + 0)*L] = csubf(cmulf(q, H[(2 * 1 + 0)*L]), cmulf(p, z));
                z = H[(2 * 0 + 1)*L];
                H[(2 * 0 + 1)*L] = caddf(cmuljf(z, q), cmuljf(H[(2 * 1 + 1)*L], p));
                H[(2 * 1 + 1)*L] = csubf(cmulf(q, H[(2 * 1 + 1)*L]), cmulf(p, z));
                /* Column modification */
                z = H[(2 * 0 + 0)*L];
                H[(2 * 0 + 0)*L] = caddf(cmulf(q, z), cmulf(p, H[(2 * 0 + 1)*L]));
                H[(2 * 0 + 1)*L] = csubf(cmuljf(H[(2 * 0 + 1)*L], q), cmuljf(z, p));
                z = H[(2 * 1 + 0)*L];
                H[(2 * 1 + 0)*L] = caddf(cmulf(q, z), cmulf(p, H[(2 * 1 + 1)*L]));
                H[(2 * 1 + 1)*L] = csubf(cmuljf(H[(2 * 1 + 1)*L], q), cmuljf(z, p));
                /* Accumulate transformations: right-multiply P by the rotator. */
                P[(2 * 0 + 0)*L] = q;
                P[(2 * 0 + 1)*L] = cnegf(_conjf(p));
                P[(2 * 1 + 0)*L] = p;
                P[(2 * 1 + 1)*L] = _conjf(q);
            }
        }
        e++; H++; if (P) P++;
    }
} /* eigen_hqr_2x2sf() */
#else
{
    /*
    * This implementation is based on the following articles:
    * [1] "The QR Algorithm for Real Hessenberg Matrices" by R.S. Martin,
    *     G. Petern and J.H. Wilkinson, Handbook for Automatie Computation,
    *     Vol.II Linear Algebra, Contribution II/14.
    * [2] "Eigenvectors of Real and Complex Matrices by LR and QR
    *     triangularizations" by G. Petern and J.H. Wilkinson, Handbook for
    *     Automatie Computation, Vol.II Linear Algebra, Contribution II/15.
    */
    int k;

    const xb_vecN_2xf32  _EPS = EPS;
    const xb_vecN_2xf32  _c0f = BBE_CONSTN_2XF32(0);
    const xb_vecN_4xcf32 _cplx0f = BBE_CONSTN_4XCF32(0);
    const xb_vecN_4xcf32 _c1f = BBE_SELN_4XCF32I(BBE_CONSTN_4XCF32(0), BBE_CONSTN_4XCF32(1), BBE_SELI_INTERLEAVE_2_LO);

    const xb_vecN_4xcf32 * restrict H_r;
    const xb_vecN_4xcf32 * restrict P1_r;
    const xb_vecN_4xcf32 * restrict Q_r;
    const xb_vecN_4xcf32 * restrict Z_r;
    const vboolN_2       * B_r;
          xb_vecN_4xcf32 * restrict H_w;
          xb_vecN_4xcf32 * restrict E_w;
          xb_vecN_4xcf32 * restrict P_w;
          xb_vecN_4xcf32 * restrict P1_w;
          xb_vecN_4xcf32 * restrict Q_w;
          xb_vecN_4xcf32 * restrict Z_w;
          vboolN_2       * B_w;

    complex_float * P1 = &P[0 * L];
    complex_float * Q  = &P[1 * L];
    complex_float * Z  = &P[2 * L];
    vboolN_2      * B = (vboolN_2 *)&P[4 * L] - L / (BBE_SIMD_WIDTH / 4);

    xb_vecN_4xcf32 h00, h01, h10, h11;
    xb_vecN_4xcf32 p00, p01, p10, p11;
    xb_vecN_4xcf32 e0, e1;
    xb_vecN_4xcf32 _x, _y, _z, _p, _q, _w;
    xb_vecN_2xf32 _a, _b, _s;
    vboolN_2 ble;
    vboolN_2 blz;    

    Z_w = (xb_vecN_4xcf32 *)Z;
    H_r = H_w = (xb_vecN_4xcf32 *)H;
    E_w = (xb_vecN_4xcf32 *)e;

    for (k = 0; k < L / (BBE_SIMD_WIDTH / 4); k++){
        BBE_LVN_4XCF32_XP(h00, H_r, 1 * L*sz_f32c);
        BBE_LVN_4XCF32_XP(h01, H_r, 1 * L*sz_f32c);
        BBE_LVN_4XCF32_XP(h10, H_r, 1 * L*sz_f32c);
        BBE_LVN_4XCF32_XP(h11, H_r, ((BBE_SIMD_WIDTH / 4) - 3 * L)*sz_f32c);

        /* Revealed a pair of eigenvalues. */
        //x = H[(2 * 1 + 1)*L];
        //y = H[(2 * 0 + 0)*L];
        //w = cmulf(H[(2 * 0 + 1)*L], H[(2 * 1 + 0)*L]);
        //p = rcmulf(.5f, csubf(y, x));
        //z = _csqrtf(caddf(cmulf(p, p), w));
        _x = h11; _y = h00;
        _w = BBE_MULN_4XCF32(h01, h10);
        _p = IT_RCMULN_4XCF32(BBE_CONSTN_2XF32(3), BBE_SUBN_4XCF32(_y, _x));
        _z = _w; BBE_MULAN_4XCF32(_z, _p, _p); _z = IT_SQRTN_4XCF32(_z);
        //if (crealf(p)*crealf(z) + cimagf(p)*cimagf(z) < 0) {
        //    z = csubf(p, z);
        //} else {
        //    z = caddf(p, z);
        //}
        _a = BBE_MULN_2XF32(BBE_MOVN_2XF32_FROMN_4XCF32(_p), BBE_MOVN_2XF32_FROMN_4XCF32(_z));
        _a = BBE_ADDN_2XF32(_a, BBE_SHFLN_2XF32I(_a, BBE_SHFLI_SWAP_2));
        blz = BBE_OLTN_2XF32(_a, _c0f);
        _z = BBE_MOVN_4XCF32T(BBE_SUBN_4XCF32(_p, _z), BBE_ADDN_4XCF32(_p, _z), BBE_MOVN_4_FROMN_2(blz));
        //if (cmagf(w)>0) {
        //    /* Quadratic divisor */
        //    e[(0)*L] = caddf(x, z);
        //    e[(1)*L] = csubf(x, _cdivf(w, z));
        //} else {
        //    /* Two linear divisors */
        //    e[(0)*L] = y; e[(1)*L] = x;
        //}
        blz = BBE_OLTN_2XF32(_c0f, IT_CMAGN_4XCF32(_w));
        e0 = BBE_MOVN_4XCF32T(BBE_ADDN_4XCF32(_x, _z), _y, BBE_MOVN_4_FROMN_2(blz));
        e1 = BBE_MOVN_4XCF32T(BBE_SUBN_4XCF32(_x, IT_DIVN_4XCF32(_w, _z, 1)), _x, BBE_MOVN_4_FROMN_2(blz));

        BBE_SVN_4XCF32_XP(e0, E_w, 1 * L*sz_f32c);
        BBE_SVN_4XCF32_XP(e1, E_w, ((BBE_SIMD_WIDTH / 4) - L)*sz_f32c);

        BBE_SVN_4XCF32_IP(_z, Z_w, (BBE_SIMD_WIDTH / 4)*sz_f32c);
    }

    H_r = (xb_vecN_4xcf32 *)H;
    E_w = (xb_vecN_4xcf32 *)e;
    B_w = B;

    for (k = 0; k < L / (BBE_SIMD_WIDTH / 4); k++){
        BBE_LVN_4XCF32_XP(h00, H_r, 3 * L*sz_f32c);
        BBE_LVN_4XCF32_XP(h11, H_r,-1 * L*sz_f32c);
        BBE_LVN_4XCF32_XP(h10, H_r,-2 * L*sz_f32c);
        
        // x = H[(2 * 1 + 0)*L];
        // y = H[(2 * 0 + 0)*L];
        // z = H[(2 * 1 + 1)*L];
        // !(cmagf(x) <= EPS*(cmagf(y) + cmagf(z)))
        _x = h10; _y = h00; _z = h11;
        _a = IT_CMAGN_4XCF32(_x);
        _b = BBE_ADDN_2XF32(BBE_ABSN_2XF32(BBE_MOVN_2XF32_FROMN_4XCF32(_y)), BBE_ABSN_2XF32(BBE_MOVN_2XF32_FROMN_4XCF32(_z)));
        _b = BBE_ADDN_2XF32(_b, BBE_SHFLN_2XF32I(_b, BBE_SHFLI_SWAP_2));
        _b = BBE_MULN_2XF32(_EPS, _b);
        ble = BBE_OLEN_2XF32(_a, _b);

        BBE_SBN_2_IP( ble, B_w, sizeof(vboolN_2));

        BBE_LVN_4XCF32_XP(h00, H_r, 3 * L*sz_f32c);
        BBE_LVN_4XCF32_XP(h11, H_r, ((BBE_SIMD_WIDTH / 4) - 3*L)*sz_f32c);

        /* Revealed a standalone eigenvalue. */
        //e[(0)*L] = H[(2 * 0 + 0)*L];
        //e[(1)*L] = H[(2 * 1 + 1)*L];
        e0 = h00; e1 = h11;
#if 1
        BBE_SVN_4XCF32T_X(e1, E_w, 1 * L*sz_f32c, BBE_MOVN_4_FROMN_2(ble));
        BBE_SVN_4XCF32T_IP(e0, E_w, (BBE_SIMD_WIDTH / 4)*sz_f32c, BBE_MOVN_4_FROMN_2(ble));
#else
        BBE_SVN_4XCF32T_XP(e0, E_w, 1 * L*sz_f32c, BBE_MOVN_4_FROMN_2(ble));
        BBE_SVN_4XCF32T_XP(e1, E_w, ((BBE_SIMD_WIDTH / 4) - 1 * L)*sz_f32c, BBE_MOVN_4_FROMN_2(ble));
#endif
    }

    if (calcEigenVec){
#if 0
        // Update H and P
        for (k = 0; k < L / (BBE_SIMD_WIDTH / 4); k++){
            H_r = H_w = (xb_vecN_4xcf32 *)H;
            P_w = (xb_vecN_4xcf32 *)P;
            h00 = BBE_LVN_4XCF32_I(H_r, 0);
            h01 = BBE_LVN_4XCF32_X(H_r, 1 * L*sz_f32c);
            h10 = BBE_LVN_4XCF32_X(H_r, 2 * L*sz_f32c);
            h11 = BBE_LVN_4XCF32_X(H_r, 3 * L*sz_f32c);

            // x = H[(2 * 1 + 0)*L];
            // y = H[(2 * 0 + 0)*L];
            // z = H[(2 * 1 + 1)*L];
            // !(cmagf(x) <= EPS*(cmagf(y) + cmagf(z)))
            _x = h10; _y = h00; _z = h11;
            _a = IT_CMAGN_4XCF32(_x);
            _b = BBE_ADDN_2XF32(BBE_ABSN_2XF32(BBE_MOVN_2XF32_FROMN_4XCF32(_y)), BBE_ABSN_2XF32(BBE_MOVN_2XF32_FROMN_4XCF32(_z)));
            _b = BBE_ADDN_2XF32(_b, BBE_SHFLN_2XF32I(_b, BBE_SHFLI_SWAP_2));
            _b = BBE_MULN_2XF32(_EPS, _b);
            bnle = BBE_OLEN_2XF32(_a, _b);
            // load temporary stored _z
            BBE_LVN_4XCF32_IP(_z, P_w, 0);

            /* Annihilate T(en,en-1) entry by similarity transformation involving
            * a rotation in the en,en-1 plane. */
            //x = H[(2 * 1 + 0)*L]; s = sqrtf(cabs2f(x) + cabs2f(z));
            //p = crdivf(x, s); q = crdivf(z, s);
            _x = h10;
            _s = BBE_MOVN_2XF32_FROMN_4XCF32(_x); _s = BBE_MULN_2XF32(_s, _s);
            _a = BBE_MOVN_2XF32_FROMN_4XCF32(_z); BBE_MULAN_2XF32(_s, _a, _a);
            _s = BBE_ADDN_2XF32(_s, BBE_SHFLN_2XF32I(_s, BBE_SHFLI_SWAP_2));
            _s = BBE_RSQRTN_2XF32(_s);
            _p = IT_RCMULN_4XCF32(_s, _x); _q = IT_RCMULN_4XCF32(_s, _z);
            /* Row modification */
            _z = h00;
            h00 = BBE_MULJN_4XCF32(_z, _q); BBE_MULJAN_4XCF32(h00, h10, _p);
            h10 = BBE_MULN_4XCF32(_q, h10); BBE_MULSN_4XCF32(h10, _p, _z);
            _z = h01;
            h01 = BBE_MULJN_4XCF32(_z, _q); BBE_MULJAN_4XCF32(h01, h11, _p);
            h11 = BBE_MULN_4XCF32(_q, h11); BBE_MULSN_4XCF32(h11, _p, _z);
            /* Column modification */
            _z = h00;
            h00 = BBE_MULN_4XCF32(_q, _z); BBE_MULAN_4XCF32(h00, _p, h01);
            h01 = BBE_MULJN_4XCF32(h01, _q); BBE_MULJSN_4XCF32(h01, _z, _p);
            _z = h10;
            h10 = BBE_MULN_4XCF32(_q, _z); BBE_MULAN_4XCF32(h10, _p, h11);
            h11 = BBE_MULJN_4XCF32(h11, _q); BBE_MULJSN_4XCF32(h11, _z, _p);
            /* Accumulate transformations: right-multiply P by the rotator. */
            p00 = _q;
            p01 = BBE_NEGN_4XCF32(BBE_CONJN_4XCF32(_p));
            p10 = _p;
            p11 = BBE_CONJN_4XCF32(_q);

            BBE_SVN_4XCF32F_I(h00, H_w,             0, BBE_MOVN_4_FROMN_2(bnle));
            BBE_SVN_4XCF32F_X(h01, H_w, 1 * L*sz_f32c, BBE_MOVN_4_FROMN_2(bnle));
            BBE_SVN_4XCF32F_X(h10, H_w, 2 * L*sz_f32c, BBE_MOVN_4_FROMN_2(bnle));
            BBE_SVN_4XCF32F_X(h11, H_w, 3 * L*sz_f32c, BBE_MOVN_4_FROMN_2(bnle));

            BBE_SVN_4XCF32F_I(p00, P_w,            0 , BBE_MOVN_4_FROMN_2(bnle));
            BBE_SVN_4XCF32F_X(p01, P_w, 1 * L*sz_f32c, BBE_MOVN_4_FROMN_2(bnle));
            BBE_SVN_4XCF32F_X(p10, P_w, 2 * L*sz_f32c, BBE_MOVN_4_FROMN_2(bnle));
            BBE_SVN_4XCF32F_X(p11, P_w, 3 * L*sz_f32c, BBE_MOVN_4_FROMN_2(bnle));

            BBE_SVN_4XCF32T_I(_c1f, P_w, 0, BBE_MOVN_4_FROMN_2(bnle));
            BBE_SVN_4XCF32T_X(BBE_CONSTN_4XCF32(0), P_w, 1 * L*sz_f32c, BBE_MOVN_4_FROMN_2(bnle));
            BBE_SVN_4XCF32T_X(BBE_CONSTN_4XCF32(0), P_w, 2 * L*sz_f32c, BBE_MOVN_4_FROMN_2(bnle));
            BBE_SVN_4XCF32T_X(_c1f, P_w, 3 * L*sz_f32c, BBE_MOVN_4_FROMN_2(bnle));

            e += (BBE_SIMD_WIDTH / 4); H += (BBE_SIMD_WIDTH / 4); P += (BBE_SIMD_WIDTH / 4);
        }
#else
        H_r = (xb_vecN_4xcf32 *)&H[2 * L];
        Z_r = (xb_vecN_4xcf32 *)Z;
        P1_w = (xb_vecN_4xcf32 *)P1;
        Q_w  = (xb_vecN_4xcf32 *)Q;

        for (k = 0; k < L / (BBE_SIMD_WIDTH / 4); k++){
            BBE_LVN_4XCF32_IP(h10, H_r, (BBE_SIMD_WIDTH / 4)*sz_f32c);
            BBE_LVN_4XCF32_IP(_z, Z_r, (BBE_SIMD_WIDTH / 4)*sz_f32c);

            /* Annihilate T(en,en-1) entry by similarity transformation involving
            * a rotation in the en,en-1 plane. */
            //x = H[(2 * 1 + 0)*L]; s = sqrtf(cabs2f(x) + cabs2f(z));
            //p = crdivf(x, s); q = crdivf(z, s);
            _x = h10;
            _s = BBE_MOVN_2XF32_FROMN_4XCF32(_x); _s = BBE_MULN_2XF32(_s, _s);
            _a = BBE_MOVN_2XF32_FROMN_4XCF32(_z); BBE_MULAN_2XF32(_s, _a, _a);
            _s = BBE_ADDN_2XF32(_s, BBE_SHFLN_2XF32I(_s, BBE_SHFLI_SWAP_2));
            _s = BBE_RSQRTN_2XF32(_s);
            _p = IT_RCMULN_4XCF32(_s, _x); _q = IT_RCMULN_4XCF32(_s, _z);

            BBE_SVN_4XCF32_IP(_p, P1_w, (BBE_SIMD_WIDTH / 4)*sz_f32c);
            BBE_SVN_4XCF32_IP(_q,  Q_w, (BBE_SIMD_WIDTH / 4)*sz_f32c);
        }

        H_r = H_w = (xb_vecN_4xcf32 *)H;
        P_w = (xb_vecN_4xcf32 *)P;
        P1_r = (xb_vecN_4xcf32 *)P1;
        Q_r  = (xb_vecN_4xcf32 *)Q;
        B_r = B;

        for (k = 0; k < L / (BBE_SIMD_WIDTH / 4); k++){
            BBE_LVN_4XCF32_XP(h00, H_r, 1 * L*sz_f32c);
            BBE_LVN_4XCF32_XP(h01, H_r, 1 * L*sz_f32c);
            BBE_LVN_4XCF32_XP(h10, H_r, 1 * L*sz_f32c);
            BBE_LVN_4XCF32_XP(h11, H_r, ((BBE_SIMD_WIDTH / 4) - 3 * L)*sz_f32c);

            BBE_LVN_4XCF32_IP(_p, P1_r, (BBE_SIMD_WIDTH / 4)*sz_f32c);
            BBE_LVN_4XCF32_IP(_q,  Q_r, (BBE_SIMD_WIDTH / 4)*sz_f32c);

            BBE_LBN_2_IP(ble, B_r, sizeof(vboolN_2));

            /* Row modification */
            _z = h00;
            h00 = BBE_MULJN_4XCF32(_z, _q); BBE_MULJAN_4XCF32(h00, h10, _p);
            h10 = BBE_MULN_4XCF32(_q, h10); BBE_MULSN_4XCF32(h10, _p, _z);
            _z = h01;
            h01 = BBE_MULJN_4XCF32(_z, _q); BBE_MULJAN_4XCF32(h01, h11, _p);
            h11 = BBE_MULN_4XCF32(_q, h11); BBE_MULSN_4XCF32(h11, _p, _z);
            /* Column modification */
            _z = h00;
            h00 = BBE_MULN_4XCF32(_q, _z); BBE_MULAN_4XCF32(h00, _p, h01);
            h01 = BBE_MULJN_4XCF32(h01, _q); BBE_MULJSN_4XCF32(h01, _z, _p);
            _z = h10;
            h10 = BBE_MULN_4XCF32(_q, _z); BBE_MULAN_4XCF32(h10, _p, h11);
            h11 = BBE_MULJN_4XCF32(h11, _q); BBE_MULJSN_4XCF32(h11, _z, _p);
            /* Accumulate transformations: right-multiply P by the rotator. */
            p00 = BBE_MOVN_4XCF32T(_c1f, _q, BBE_MOVN_4_FROMN_2(ble));
            p01 = BBE_MOVN_4XCF32T(_cplx0f, BBE_NEGN_4XCF32(BBE_CONJN_4XCF32(_p)), BBE_MOVN_4_FROMN_2(ble));
            p10 = BBE_MOVN_4XCF32T(_cplx0f, _p, BBE_MOVN_4_FROMN_2(ble));
            p11 = BBE_MOVN_4XCF32T(_c1f, BBE_CONJN_4XCF32(_q), BBE_MOVN_4_FROMN_2(ble));

            BBE_SVN_4XCF32F_XP(h00, H_w, 1 * L*sz_f32c, BBE_MOVN_4_FROMN_2(ble));
            BBE_SVN_4XCF32F_XP(h01, H_w, 1 * L*sz_f32c, BBE_MOVN_4_FROMN_2(ble));
            BBE_SVN_4XCF32F_XP(h10, H_w, 1 * L*sz_f32c, BBE_MOVN_4_FROMN_2(ble));
            BBE_SVN_4XCF32F_XP(h11, H_w, ((BBE_SIMD_WIDTH / 4) - 3 * L)*sz_f32c, BBE_MOVN_4_FROMN_2(ble));

            BBE_SVN_4XCF32_XP(p00, P_w, 1 * L*sz_f32c);
            BBE_SVN_4XCF32_XP(p01, P_w, 1 * L*sz_f32c);
            BBE_SVN_4XCF32_XP(p10, P_w, 1 * L*sz_f32c);
            BBE_SVN_4XCF32_XP(p11, P_w, ((BBE_SIMD_WIDTH / 4) - 3 * L)*sz_f32c);
        }
    }
#endif
} /* eigen_hqr_2x2sf() */
#endif

static void eigen_bksubst_2x2sf(complex_float * restrict T, int L)
#if 0
{
    float32_t w;
    complex_float y, z;
    int k, exp;

    const complex_float c1f = _makecomplexf(1.f, 0.f);
    const float32_t small = FLT_MIN;
    const float32_t big = 1e3f;

    for (k = 0; k < L; k++){
        y = csubf(T[(2 * 1 + 1)*L], T[(2 * 0 + 0)*L]);
        T[(2 * 1 + 1)*L] = c1f;
        if (cabs2f(y) == 0) {
            y = _makecomplexf(MAX(EPS*cmagf(T[(2 * 0 + 0)*L]), small), 0.f);
        }
        z = T[(2 * 0 + 1)*L];
        T[(2 * 0 + 1)*L] = _cdivf(z, y);
        /* Conditionally rescale the vector to prevent overflow. */
        w = cmagf(T[(2 * 0 + 1)*L]);
        if (w >= big) {
            exp = ilogbf(w);
            T[(2 * 0 + 1)*L] = cldexpf(T[(2 * 0 + 1)*L], -exp);
            T[(2 * 1 + 1)*L] = cldexpf(T[(2 * 1 + 1)*L], -exp);
        }
        T[(2 * 0 + 0)*L] = c1f;

        T++;
    }
} /* eigen_bksubst_2x2sf() */
#elif 0
{
    const xb_vecN_4xcf32 * restrict T_r;
          xb_vecN_4xcf32 * restrict T_w;

    const xb_vecN_4xcf32 _c1f = BBE_SELN_4XCF32I(BBE_CONSTN_4XCF32(0), BBE_CONSTN_4XCF32(1), BBE_SELI_INTERLEAVE_2_LO);    
    const xb_vecN_2xf32  _c0f = BBE_CONSTN_2XF32(0);
    const xb_vecN_2xf32 _EPS = EPS;
    const xb_vecN_2xf32 _small = FLT_MIN;
    const xb_vecN_2xf32 _big = 1e3f;
    const xb_vecNx16 exp_mask = BBE_MOVVINX16(BBE_FLOAT_MASK_EXP);

    xb_vecN_4xcf32 t00, t01, t11;
    xb_vecN_4xcf32 _y, _z;
    xb_vecN_2xf32 _w;
    xb_vecN_2xf32 _exp;
    vboolN_2 bez, ble;

    int k;

    for (k = 0; k < L/(BBE_SIMD_WIDTH/4); k++){
        T_r = T_w = (xb_vecN_4xcf32 *)T;

        t00 = BBE_LVN_4XCF32_I(T_r, 0);
        t11 = BBE_LVN_4XCF32_X(T_r, 3*L*sz_f32c);
        _y = BBE_SUBN_4XCF32(t11,t00);
        t11 = _c1f;
        BBE_SVN_4XCF32_X(t11, T_w, 3 * L*sz_f32c);

        _w = IT_ABS2N_4XCF32(_y);
        bez = BBE_OEQN_2XF32(_w,_c0f); // cabs2f(y) == 0
        _w = BBE_MAXN_2XF32(BBE_MULN_2XF32(_EPS, IT_CMAGN_4XCF32(t00)), _small);
        t00 = _c1f;
        BBE_SVN_4XCF32_I(t00, T_w, 0);
        _w = BBE_SELN_2XF32I(_c0f, _w, BBE_SELI_INTERLEAVE_2_EVEN);
        _y = BBE_MOVN_4XCF32T(BBE_MOVN_4XCF32_FROMN_2XF32(_w), _y, BBE_MOVN_4_FROMN_2(bez));

        t01 = BBE_LVN_4XCF32_X(T_r, 1 * L*sz_f32c);
        _z = t01;
        t01 = IT_DIVN_4XCF32(_z, _y, 0);
        BBE_SVN_4XCF32_X(t01, T_w, 1 * L*sz_f32c);
        /* Conditionally rescale the vector to prevent overflow. */
        _w = IT_CMAGN_4XCF32(t01);
        ble = BBE_OLEN_2XF32(_big, _w); //w >= big
        _exp = BBE_RECIP0N_2XF32(BBE_MOVN_2XF32_FROMNX16(BBE_ANDNX16(BBE_MOVNX16_FROMN_2XF32(_w), exp_mask)));
#if 0
        _w = BBE_MOVN_2XF32_FROMN_4XCF32(t01); BBE_MULN_2XF32T(_w, _w, _exp, ble);
        t01 = BBE_MOVN_4XCF32_FROMN_2XF32(_w);
        _w = BBE_MOVN_2XF32_FROMN_4XCF32(t11); BBE_MULN_2XF32T(_w, _w, _exp, ble);
        t11 = BBE_MOVN_4XCF32_FROMN_2XF32(_w);
        BBE_SVN_4XCF32_X(t01, T_w, 1 * L*sz_f32c);
        BBE_SVN_4XCF32_X(t11, T_w, 3 * L*sz_f32c);
#else
        //_w = BBE_MOVN_2XF32_FROMN_4XCF32(t01); BBE_MULN_2XF32T(_w, _w, _exp, ble);
        //t01 = BBE_MOVN_4XCF32_FROMN_2XF32(_w);
        _w = BBE_MOVN_2XF32_FROMN_4XCF32(t01); _w = BBE_MULN_2XF32(_w, _exp);
        t01 = BBE_MOVN_4XCF32_FROMN_2XF32(_w);
        //_w = BBE_MOVN_2XF32_FROMN_4XCF32(t11); BBE_MULN_2XF32T(_w, _w, _exp, ble);
        //t11 = BBE_MOVN_4XCF32_FROMN_2XF32(_w);
        t11 = BBE_MOVN_4XCF32_FROMN_2XF32(BBE_SELN_2XF32I(_c0f, _exp, BBE_SELI_INTERLEAVE_2_EVEN));

        BBE_SVN_4XCF32T_X(t01, T_w, 1 * L*sz_f32c, BBE_MOVN_4_FROMN_2(ble));
        BBE_SVN_4XCF32T_X(t11, T_w, 3 * L*sz_f32c, BBE_MOVN_4_FROMN_2(ble));
#endif
        T += (BBE_SIMD_WIDTH / 4);
    }
} /* eigen_bksubst_2x2sf() */
#else
{
    const xb_vecN_4xcf32 * restrict T_r;
          xb_vecN_4xcf32 * restrict T_w;

    const xb_vecN_4xcf32 _c1f = BBE_SELN_4XCF32I(BBE_CONSTN_4XCF32(0), BBE_CONSTN_4XCF32(1), BBE_SELI_INTERLEAVE_2_LO);
    const xb_vecN_2xf32  _c0f = BBE_CONSTN_2XF32(0);
    const xb_vecN_2xf32 _EPS = EPS;
    const xb_vecN_2xf32 _small = FLT_MIN;
    const xb_vecN_2xf32 _big = 1e3f;
    const xb_vecNx16 exp_mask = BBE_MOVVINX16(BBE_FLOAT_MASK_EXP);

    xb_vecN_4xcf32 t00, t01, t11;
    xb_vecN_4xcf32 _y, _z;
    xb_vecN_2xf32 _w;
    xb_vecN_2xf32 _exp;
    vboolN_2 bez, ble;

    int k;

    T_r = T_w = (xb_vecN_4xcf32 *)T;
    for (k = 0; k < L / (BBE_SIMD_WIDTH / 4); k++){
        t11 = BBE_LVN_4XCF32_X(T_r, 3 * L*sz_f32c);
        BBE_LVN_4XCF32_IP(t00, T_r, (BBE_SIMD_WIDTH / 4)*sz_f32c);
        _y = BBE_SUBN_4XCF32(t11, t00);
        _w = IT_ABS2N_4XCF32(_y);
        bez = BBE_OEQN_2XF32(_w, _c0f); // cabs2f(y) == 0
        _w = BBE_MAXN_2XF32(BBE_MULN_2XF32(_EPS, IT_CMAGN_4XCF32(t00)), _small);
        _w = BBE_SELN_2XF32I(_c0f, _w, BBE_SELI_INTERLEAVE_2_EVEN);
        _y = BBE_MOVN_4XCF32T(BBE_MOVN_4XCF32_FROMN_2XF32(_w), _y, BBE_MOVN_4_FROMN_2(bez));
        t11 = _c1f;
        BBE_SVN_4XCF32_X(t11, T_w, 3 * L*sz_f32c);
        t00 = _y;
        BBE_SVN_4XCF32_IP(t00, T_w, (BBE_SIMD_WIDTH / 4)*sz_f32c);
    }

    T_r = T_w = (xb_vecN_4xcf32 *)T;
    for (k = 0; k < L / (BBE_SIMD_WIDTH / 4); k++){
        BBE_LVN_4XCF32_XP(t00, T_r, 1 * L*sz_f32c);
        BBE_LVN_4XCF32_XP(t01, T_r, ((BBE_SIMD_WIDTH / 4) - 1 * L)*sz_f32c);
        _y = t00;  _z = t01;
        t00 = _c1f;  t01 = IT_DIVN_4XCF32(_z, _y, 0);
        BBE_SVN_4XCF32_X(t01, T_w, 1 * L*sz_f32c);
        BBE_SVN_4XCF32_IP(t00, T_w, (BBE_SIMD_WIDTH / 4)*sz_f32c);
    }

    /* Conditionally rescale the vector to prevent overflow. */
    T_r = T_w = (xb_vecN_4xcf32 *)&T[L];
    for (k = 0; k < L / (BBE_SIMD_WIDTH / 4); k++){
        BBE_LVN_4XCF32_IP(t01, T_r, (BBE_SIMD_WIDTH/4)*sz_f32c);
        _w = IT_CMAGN_4XCF32(t01);
        ble = BBE_OLEN_2XF32(_big, _w); //w >= big
        _exp = BBE_RECIP0N_2XF32(BBE_MOVN_2XF32_FROMNX16(BBE_ANDNX16(BBE_MOVNX16_FROMN_2XF32(_w), exp_mask)));
        _w = BBE_MOVN_2XF32_FROMN_4XCF32(t01);  _w = BBE_MULN_2XF32(_w, _exp);
        t01 = BBE_MOVN_4XCF32_FROMN_2XF32(_w);
        _w = BBE_MOVN_2XF32_FROMN_4XCF32(_c1f); _w = BBE_MULN_2XF32(_w, _exp);
        t11 = BBE_MOVN_4XCF32_FROMN_2XF32(_w);
        BBE_SVN_4XCF32T_X (t11, T_w, 2 * L*sz_f32c               , BBE_MOVN_4_FROMN_2(ble));
        BBE_SVN_4XCF32T_IP(t01, T_w, (BBE_SIMD_WIDTH / 4)*sz_f32c, BBE_MOVN_4_FROMN_2(ble));
    }
} /* eigen_bksubst_2x2sf() */
#endif

void eigen_tmulp_2x2sf(complex_float * restrict V, /* V[N*N][L]       */
                 const complex_float * restrict P, /* P[N*N][L]       */
                 const complex_float * restrict T, /* T[N*(N+1)/2][L] */
                       int L)
{
    /* A variant without exponent normalization. */
    const xb_vecN_4xcf32 * restrict P_r;
    const xb_vecN_4xcf32 * restrict T_r;
          xb_vecN_4xcf32 * V_w;

    xb_vecN_4xcf32 p00, p01, p10, p11;
    xb_vecN_4xcf32 v00, v01, v10, v11;
    xb_vecN_4xcf32 t00, t01, t11;
    xb_vecN_2xf32  z, n0, n1;

    int k;

    P_r = (xb_vecN_4xcf32 *)P;
    T_r = (xb_vecN_4xcf32 *)T;
    V_w = (xb_vecN_4xcf32 *)V;

    for (k = 0; k < L / (BBE_SIMD_WIDTH / 4); k++){
        BBE_LVN_4XCF32_XP(p00, P_r, 1 * L*sz_f32c);
        BBE_LVN_4XCF32_XP(p01, P_r, 1 * L*sz_f32c);
        BBE_LVN_4XCF32_XP(p10, P_r, 1 * L*sz_f32c);
        BBE_LVN_4XCF32_XP(p11, P_r, ((BBE_SIMD_WIDTH/4) - 3*L)*sz_f32c);

        BBE_LVN_4XCF32_XP(t00, T_r, 1 * L*sz_f32c);
        BBE_LVN_4XCF32_XP(t01, T_r, 2 * L*sz_f32c);
        BBE_LVN_4XCF32_XP(t11, T_r, ((BBE_SIMD_WIDTH / 4) - 3 * L)*sz_f32c);
        /* Compute P*T. */
        v00 = BBE_MULN_4XCF32(p00, t00); v10 = BBE_MULN_4XCF32(p10, t00);
        v01 = BBE_MULN_4XCF32(p00, t01); BBE_MULAN_4XCF32(v01, p01, t11);
        v11 = BBE_MULN_4XCF32(p10, t01); BBE_MULAN_4XCF32(v11, p11, t11);

        /* Compute the L2 norm. */
        z = BBE_MOVN_2XF32_FROMN_4XCF32(v00); n0 = BBE_MULN_2XF32(z, z);
        z = BBE_MOVN_2XF32_FROMN_4XCF32(v10); BBE_MULAN_2XF32(n0, z, z);
        n0 = BBE_ADDN_2XF32(n0, BBE_SHFLN_2XF32I(n0, BBE_SHFLI_SWAP_2));
        n0 = BBE_RSQRTN_2XF32(n0);

        z = BBE_MOVN_2XF32_FROMN_4XCF32(v01); n1 = BBE_MULN_2XF32(z, z);
        z = BBE_MOVN_2XF32_FROMN_4XCF32(v11); BBE_MULAN_2XF32(n1, z, z);
        n1 = BBE_ADDN_2XF32(n1, BBE_SHFLN_2XF32I(n1, BBE_SHFLI_SWAP_2));
        n1 = BBE_RSQRTN_2XF32(n1);

        /* Scale the column by its L2 norm. */
        z = BBE_MOVN_2XF32_FROMN_4XCF32(v00);
        v00 = BBE_MOVN_4XCF32_FROMN_2XF32(BBE_MULN_2XF32(z, n0));
        z = BBE_MOVN_2XF32_FROMN_4XCF32(v10);
        v10 = BBE_MOVN_4XCF32_FROMN_2XF32(BBE_MULN_2XF32(z, n0));

        z = BBE_MOVN_2XF32_FROMN_4XCF32(v01);
        v01 = BBE_MOVN_4XCF32_FROMN_2XF32(BBE_MULN_2XF32(z, n1));
        z = BBE_MOVN_2XF32_FROMN_4XCF32(v11);
        v11 = BBE_MOVN_4XCF32_FROMN_2XF32(BBE_MULN_2XF32(z, n1));

        BBE_SVN_4XCF32_XP(v00, V_w, 2 * L*sz_f32c);
        BBE_SVN_4XCF32_XP(v10, V_w, -1 * L*sz_f32c);
        BBE_SVN_4XCF32_XP(v01, V_w, 2 * L*sz_f32c);
        BBE_SVN_4XCF32_XP(v11, V_w, ((BBE_SIMD_WIDTH / 4) - 3 * L)*sz_f32c);
    }
}/* eigen_tmulp_2x2sf() */

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
void eigen2x2sf ( 
            void * pScr,
            complex_float * restrict e,
            complex_float * restrict V,
            complex_float * restrict A,
            int L )
{
  NASSERT_ALIGN( pScr, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( e   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( V   , 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( A   , 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/4)) );

  eigennxnsf(pScr,e,V,A,2,L);

} /* eigen2x2sf() */

size_t eigen2x2sf_getScratchSize ( int N, int L )
{
  NASSERT( 2==N );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/4)) );
  return ( eigennxnsf_getScratchSize(N,L) );
}
#else
void eigen2x2sf(
    void * pScr,
    complex_float * restrict e,
    complex_float * restrict V,
    complex_float * restrict A,
    int L)
{
    complex_float *Ps;
    //int p;

    NASSERT_ALIGN(pScr, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(e, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(V, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(A, 2 * BBE_SIMD_WIDTH);
    NASSERT(0 == (L % (BBE_SIMD_WIDTH / 4)));

    {
        void * p = pScr;
        Ps = (complex_float*)p; p = Ps + L * 2 * 2;
    }

    if (V) {
#if 0
        /* Initialize transformation matrices with 2x2 identity matrix. */
        const complex_float c0f = _makecomplexf(0.f, 0.f);
        const complex_float c1f = _makecomplexf(1.f, 0.f);
        for (p = 0; p<L; p++) {
            Ps[p + 0 * L] = c1f;
            Ps[p + 1 * L] = c0f;
            Ps[p + 2 * L] = c0f;
            Ps[p + 3 * L] = c1f;
        }
#endif
        /* Apply the QR algorithm and update the transformation matrix. */
        eigen_hqr_2x2sf(e, A, Ps, L, 1);
        /* Perfrom the backsubstitution to determine eigenvectors of
        * triangular form lying in Hb */
        eigen_bksubst_2x2sf(A, L);
#if 0
        /* Convert transformation matrices and eigenvectors of triangular forms
        * to stream order. */
        for (p = 0; p<L; p++){
            A[p + 0 * L] = A[p + 0 * L];
            A[p + 1 * L] = A[p + 1 * L];
            A[p + 2 * L] = A[p + 3 * L];
        }
#endif
        /* Left-multiply eigenvectors by transformation matrices and rescale them to
        * obtain eigenvectors for original input matrices. */
        eigen_tmulp_2x2sf(V, Ps, A, L);
    }
    else {
        eigen_hqr_2x2sf(e, A, Ps, L, 0);
    }
} /* eigen2x2sf() */

size_t eigen2x2sf_getScratchSize(int N, int L)
{
    NASSERT(2 == N);
    NASSERT(0 == (L % (BBE_SIMD_WIDTH / 4)));
    return (L*(2 * 2)*sz_f32c);  /* Ps: transformation matrix in strean order */
}
#endif

#else /* HAVE_VFPU */

DISCARD_FUN( void, eigen2x2sf, ( void * pScr,
                         complex_float * restrict e,
                         complex_float * restrict V,
                         complex_float * restrict A,
                         int L ) )

size_t eigen2x2sf_getScratchSize ( int N, int L ) 
{
  NASSERT(2==N);
  return (0);
}

#endif /* HAVE_VFPU */
