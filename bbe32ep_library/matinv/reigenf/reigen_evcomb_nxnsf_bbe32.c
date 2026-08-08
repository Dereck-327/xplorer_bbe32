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
    ombine eigenvectors of real matrices from re/im components.
    Real Data, Stream Order
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

#include <math.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"

#if !defined(COMPILER_MSVC)
#include <complex.h>
#endif

#if HAVE_VFPU

#define sz_f32    sizeof(float32_t)
#define sz_f32c   sizeof(complex_float)

#if 0
static complex_float _makecomplexf( float32_t re, float32_t im )
{
  union { float32_t r[2]; complex_float c; } u = {{re,im}};
  return (u.c);
}

static void r_evcomb_f( complex_float * restrict V,
                  const float32_t     * restrict U,
                  const complex_float * restrict e,
                  int N, int stride );
#endif

static void r_evcomb_fVec(complex_float * restrict V,
                    const float32_t     * restrict U,
                    const complex_float * restrict e,
                    int N, int stride);

/*
 * Combine eigenvectors of real matrices from re/im components. The function
 * exploits the special order of conjugate eigenvalues/eigenvectors by rhqr()
 * function. Eigenvectors are also rescaled so that L2 norm of a vector is 1.
 * Input:
 *   N          Matrix size
 *   L          Number of matrices
 *   U[N*N][L]  Real matrices of re/im components of eigenvectors
 *   e[N][L]    Complex eigenvalues
 * Output:
 *   V[N*N][L]  Complex matrices of eigenvectors
 * Restrictions:
 *   V,U,e      Must not overlap and must be aligned on 2*BBE_SIMD_WIDTH-byte
 *              boundary
 *   Variant functions may impose additional restrictions
 */

/* Real Data, Stream Order
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2 */
void reigen_evcomb_nxnsf ( complex_float * restrict V, /* V[N*N][L] */
                     const float32_t     * restrict U, /* U[N*N][L] */
                     const complex_float * restrict e, /* e[N][L]   */
                     int N, int L )
{
  int k;
  const int STEP = BBE_SIMD_WIDTH / 2;

  NASSERT_ALIGN( V, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( U, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( e, 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );
#if 0
  for ( k=0; k<L; k++ ) {
    r_evcomb_f(V,U,e,N,L);
    V++; U++; e++;
  }
#else
  for (k = 0; k<L/STEP; k++) {
      r_evcomb_fVec(V, U, e, N, L);
      V += STEP; U += STEP; e += STEP;
  }
#endif

} /* reigen_evcomb_nxnsf() */

#if 0
void r_evcomb_f( complex_float * restrict V,
           const float32_t     * restrict U,
           const complex_float * restrict e,
           int N, int stride )
{
  /*
   * MATLAB reference code:
   *
   *   V = zeros(N);
   *   for j=1:N
   *     if imag(e(j))>0
   *       % Combine re/im parts into complex vector.
   *       V(:,[j,j+1]) = U(:,[j,j+1])*[1,1;1j,-1j];
   *     elseif imag(ev(j))==0
   *       % Real vector
   *       V(:,j) = U(:,j);
   *     end
   *   end
   *   % Rescale vectors
   *   for j=1:N
   *     V(:,j) = V(:,j)/norm(V(:,j),2);
   *   end;
   */

  int i,j,exp=0;
  float32_t f,g,s=0.f;

  for ( j=0; j<N; j++ ) {
    if ( cimagf(e[j*stride])>0.f ) {
      s = 0.f;
      for ( i=0; i<N; i++ ) {
        f = fabsf(U[(i*N+j+0)*stride]);
        g = fabsf(U[(i*N+j+1)*stride]);
        if (s<f+g) s = f+g;
      }
      exp = ilogbf(s);
      s = 0.f;
      for ( i=0; i<N; i++ ) {
        f = ldexpf(U[(i*N+j+0)*stride], -exp);
        g = ldexpf(U[(i*N+j+1)*stride], -exp);
        s += f*f+g*g;
      }
      s = sqrtf(s);
      for ( i=0; i<N; i++ ) {
        f = ldexpf(U[(i*N+j+0)*stride], -exp);
        g = ldexpf(U[(i*N+j+1)*stride], -exp);
        V[(i*N+j)*stride] = _makecomplexf(f/s,g/s);
      }
    } else if ( cimagf(e[j*stride])<0.f ) {
      /* Use exp,s from the preceding iteration. */
      for ( i=0; i<N; i++ ) {
        NASSERT(j>0);
        f = ldexpf(U[(i*N+j-1)*stride], -exp);
        g = ldexpf(U[(i*N+j-0)*stride], -exp);
        V[(i*N+j)*stride] = _makecomplexf(f/s,-g/s);
      }
    } else {
      s = 0.f;
      for ( i=0; i<N; i++ ) {
        f = fabsf(U[(i*N+j)*stride]);
        if (s<f) s = f;
      }
      exp = ilogbf(s);
      s = 0.f;
      for ( i=0; i<N; i++ ) {
        f = ldexpf(U[(i*N+j)*stride], -exp);
        s += f*f;
      }
      s = sqrtf(s);
      for ( i=0; i<N; i++ ) {
        f = ldexpf(U[(i*N+j)*stride], -exp);
        V[(i*N+j)*stride] = _makecomplexf(f/s,0.f);
      }
    }
  }

} /* r_evcomb_f() */
#endif


void r_evcomb_fVec( complex_float * restrict V,
              const float32_t     * restrict U,
              const complex_float * restrict e,
              int N, int stride){
    /*
    * MATLAB reference code:
    *
    *   V = zeros(N);
    *   for j=1:N
    *     if imag(e(j))>0
    *       % Combine re/im parts into complex vector.
    *       V(:,[j,j+1]) = U(:,[j,j+1])*[1,1;1j,-1j];
    *     elseif imag(ev(j))==0
    *       % Real vector
    *       V(:,j) = U(:,j);
    *     end
    *   end
    *   % Rescale vectors
    *   for j=1:N
    *     V(:,j) = V(:,j)/norm(V(:,j),2);
    *   end;
    */

    int i, j;

#if 0
    int k;
    for (k = 0; k < BBE_SIMD_WIDTH / 2; k++) {
        float32_t f, g, s = 0.f;
        int exp = 0;

        for (j = 0; j<N; j++) {
            if (cimagf(e[j*stride])>0.f) {
                s = 0.f;
                for (i = 0; i < N; i++) {
                    f = fabsf(U[(i*N + j + 0)*stride]);
                    g = fabsf(U[(i*N + j + 1)*stride]);
                    if (s < f + g) s = f + g;
                }
                exp = ilogbf(s);
                s = 0.f;
                for (i = 0; i < N; i++) {
                    f = ldexpf(U[(i*N + j + 0)*stride], -exp);
                    g = ldexpf(U[(i*N + j + 1)*stride], -exp);
                    s += f*f + g*g;
                }
                s = sqrtf(s);
                for (i = 0; i < N; i++) {
                    f = ldexpf(U[(i*N + j + 0)*stride], -exp);
                    g = ldexpf(U[(i*N + j + 1)*stride], -exp);
                    V[(i*N + j)*stride] = _makecomplexf(f / s, g / s);
                }
            } else if (cimagf(e[j*stride]) < 0.f) {
                /* Use exp,s from the preceding iteration. */
                for (i = 0; i<N; i++) {
                    NASSERT(j>0);
                    f = ldexpf(U[(i*N + j - 1)*stride], -exp);
                    g = ldexpf(U[(i*N + j - 0)*stride], -exp);
                    V[(i*N + j)*stride] = _makecomplexf(f / s, -g / s);
                }
            } else {
                s = 0.f;
                for (i = 0; i < N; i++) {
                    f = fabsf(U[(i*N + j)*stride]);
                    if (s < f) s = f;
                }
                exp = ilogbf(s);
                s = 0.f;
                for (i = 0; i < N; i++) {
                    f = ldexpf(U[(i*N + j)*stride], -exp);
                    s += f*f;
                }
                s = sqrtf(s);
                for (i = 0; i < N; i++) {
                    f = ldexpf(U[(i*N + j)*stride], -exp);
                    V[(i*N + j)*stride] = _makecomplexf(f / s, 0.f);
                }
            }
        }
        V++;
        U++;
        e++;
    }
#elif 0
    const xb_vecN_2xf32  * U_r;
    const xb_vecN_4xcf32 * E_r;
          xb_vecN_4xcf32 * V_w;

    const xb_vecN_2xf32 c0f = BBE_CONSTN_2XF32(0);
    const xb_vecN_2x32v exp_mask = BBE_MOVN_2X32_FROMN_2XF32(BBE_MOVVIN_2XF32(BBE_FLOAT_MASK_EXP)); //0xFF800000
    const int step   = N*stride*sz_f32;
    const int step_c = N*stride*sz_f32c;

    xb_vecN_4xcf32 u0, u1;
    xb_vecN_2xf32  h0, h1, f, g, s, exp;
    xb_vecN_2x32v  v;
    vboolN_4 bmask_4_0 , bmask_4_1;
    vboolN_2 blez, bqez, bmoz, b0, b1;

    s   = BBE_ZERON_2XF32();
    exp = BBE_CONSTN_2XF32(1);
        
    E_r = (xb_vecN_4xcf32 *)e;

    for (j = 0; j < N; j++) {
        xb_vecN_2xf32  _s, _exp;

        u1 = BBE_LVN_4XCF32_I(E_r, (BBE_SIMD_WIDTH / 4)*sz_f32c);
        BBE_LVN_4XCF32_XP(u0, E_r, stride*sz_f32c);

        h0 = BBE_MOVN_2XF32_FROMN_4XCF32(u0);
        h1 = BBE_MOVN_2XF32_FROMN_4XCF32(u1);
        h0 = BBE_SELN_2XF32I(h1, h0, BBE_SELI_EXTRACT_2_OF_4_OFF_2);
        // cimagf(e[j*stride]) < 0.f
        blez = BBE_OLTN_2XF32(h0, c0f);
        // cimagf(e[j*stride]) > 0.f
        bmoz = BBE_OLTN_2XF32(c0f, h0);
  
        /*
         * cimagf(e[j*stride]) > 0.f
         */
        BBE_EXTRACTBN(b1, b0, BBE_MOVN_FROMN_2(bmoz));
        bmask_4_0 = BBE_MOVN_4_FROMN_2(b0);
        bmask_4_1 = BBE_MOVN_4_FROMN_2(b1);
        _s = BBE_ZERON_2XF32();
        U_r = (xb_vecN_2xf32 *)U;
        for (i = 0; i < N; i++) {
            g = BBE_LVN_2XF32_X(U_r, stride*sz_f32);
            BBE_LVN_2XF32_XP(f, U_r, step);
            f = BBE_ADDN_2XF32(BBE_ABSN_2XF32(f), BBE_ABSN_2XF32(g));
            _s = BBE_MOVN_2XF32T(f, _s, BBE_OLTN_2XF32(_s, f));
        }

        // Calculate scaling coefficient
        v = BBE_MOVN_2X32_FROMN_2XF32(_s);
        v = BBE_ANDN_2X32(v, exp_mask);
        _s = BBE_MOVN_2XF32_FROMN_2X32(v);
        _exp = BBE_RECIP0N_2XF32(_s);

        _s = BBE_ZERON_2XF32();
        U_r = (xb_vecN_2xf32 *)U;
        for (i = 0; i < N; i++) {
            g = BBE_LVN_2XF32_X(U_r, stride*sz_f32);
            BBE_LVN_2XF32_XP(f, U_r, step);
            f = BBE_MULN_2XF32(f, _exp);
            g = BBE_MULN_2XF32(g, _exp);
            BBE_MULAN_2XF32(_s, f, f);
            BBE_MULAN_2XF32(_s, g, g);
        }
        _s = BBE_RSQRTN_2XF32(_s);
        s   = BBE_MOVN_2XF32T(_s, s, bmoz);
        exp = BBE_MOVN_2XF32T(_exp, exp , bmoz);
        U_r = (xb_vecN_2xf32 *)U;
        V_w = (xb_vecN_4xcf32 *)V;
        for (i = 0; i < N; i++) {
            g = BBE_LVN_2XF32_X(U_r, stride*sz_f32);
            BBE_LVN_2XF32_XP(f, U_r, step);

            f = BBE_MULN_2XF32(f, _exp);
            g = BBE_MULN_2XF32(g, _exp);

            f = BBE_MULN_2XF32(f, _s);
            g = BBE_MULN_2XF32(g, _s);

            h1 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_HI);
            h0 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_LO);
            u0 = BBE_MOVN_4XCF32_FROMN_2XF32(h0);
            u1 = BBE_MOVN_4XCF32_FROMN_2XF32(h1);
                
            BBE_SVN_4XCF32T_I(u1, V_w, (BBE_SIMD_WIDTH / 4) * sz_f32c, bmask_4_1);
            BBE_SVN_4XCF32T_XP(u0, V_w, step_c , bmask_4_0);
        }

        /*
         * cimagf(e[j*stride]) < 0.f
         * Use exp,s from the preceding iteration.
         */
        BBE_EXTRACTBN(b1, b0, BBE_MOVN_FROMN_2(blez));
        bmask_4_0 = BBE_MOVN_4_FROMN_2(b0);
        bmask_4_1 = BBE_MOVN_4_FROMN_2(b1);

        U_r = (xb_vecN_2xf32 *)&U[-stride];
        V_w = (xb_vecN_4xcf32 *)V;
        for (i = 0; i < N; i++) {
            g = BBE_LVN_2XF32_X(U_r, stride*sz_f32);
            BBE_LVN_2XF32_XP(f, U_r, step);

            f = BBE_MULN_2XF32(f, exp);
            g = BBE_MULN_2XF32(g, exp);

            f = BBE_MULN_2XF32(f, s);
            g = BBE_MULN_2XF32(g, s);
            g = BBE_NEGN_2XF32(g);

            h1 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_HI);
            h0 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_LO);
            u0 = BBE_MOVN_4XCF32_FROMN_2XF32(h0);
            u1 = BBE_MOVN_4XCF32_FROMN_2XF32(h1);

            BBE_SVN_4XCF32T_I(u1, V_w, (BBE_SIMD_WIDTH / 4) * sz_f32c, bmask_4_1);
            BBE_SVN_4XCF32T_XP(u0, V_w, step_c, bmask_4_0);
        }

        /*
         * cimagf(e[j*stride]) == 0.f
         */
        bqez = BBE_NOTBN_2(BBE_ORBN_2(bmoz,blez));
        BBE_EXTRACTBN(b1, b0, BBE_MOVN_FROMN_2(bqez));
        bmask_4_0 = BBE_MOVN_4_FROMN_2(b0);
        bmask_4_1 = BBE_MOVN_4_FROMN_2(b1);

        _s = BBE_ZERON_2XF32();
        U_r = (xb_vecN_2xf32 *)U;
        for (i = 0; i < N; i++) {
            BBE_LVN_2XF32_XP(h0, U_r, step);
            h0 = BBE_ABSN_2XF32(h0);
            _s = BBE_MOVN_2XF32T(h0, _s, BBE_OLTN_2XF32(_s, h0));
        }
        // Calculate scaling coefficient
        v = BBE_MOVN_2X32_FROMN_2XF32(_s);
        v = BBE_ANDN_2X32(v, exp_mask);
        _s = BBE_MOVN_2XF32_FROMN_2X32(v);
        _exp = BBE_RECIP0N_2XF32(_s);

        U_r = (xb_vecN_2xf32 *)U;
        _s = BBE_ZERON_2XF32();
        for (i = 0; i < N; i++) {
            BBE_LVN_2XF32_XP(h0, U_r, step);
            h0 = BBE_MULN_2XF32(h0, _exp);
            BBE_MULAN_2XF32(_s, h0, h0);
        }

        U_r = (xb_vecN_2xf32 *)U;
        V_w = (xb_vecN_4xcf32 *)V;
        _s = BBE_RSQRTN_2XF32(_s);
        for (i = 0; i < N; i++) {
            BBE_LVN_2XF32_XP(h0, U_r, step);
            h0 = BBE_MULN_2XF32(h0, _exp);
            h0 = BBE_MULN_2XF32(h0, _s);

            h1 = BBE_SELN_2XF32I(c0f, h0, BBE_SELI_INTERLEAVE_2_HI);
            h0 = BBE_SELN_2XF32I(c0f, h0, BBE_SELI_INTERLEAVE_2_LO);

            u0 = BBE_MOVN_4XCF32_FROMN_2XF32(h0);
            u1 = BBE_MOVN_4XCF32_FROMN_2XF32(h1);
            BBE_SVN_4XCF32T_I(u1, V_w, (BBE_SIMD_WIDTH / 4) * sz_f32c, bmask_4_1);
            BBE_SVN_4XCF32T_XP(u0, V_w, step_c, bmask_4_0);
        }
        U += stride;
        V += stride;
    }
#elif 0
    const xb_vecN_2xf32  * U_r;
    const xb_vecN_4xcf32 * E_r;
    xb_vecN_4xcf32 * V_w;

    const xb_vecN_2xf32 c0f = BBE_CONSTN_2XF32(0);
    const xb_vecN_2x32v exp_mask = BBE_MOVN_2X32_FROMN_2XF32(BBE_MOVVIN_2XF32(BBE_FLOAT_MASK_EXP)); //0xFF800000
    const int step   = N*stride*sz_f32;
    const int step_c = N*stride*sz_f32c;

    xb_vecN_4xcf32 u0, u1;
    xb_vecN_2xf32  h0, h1, f, g, s, exp;
    xb_vecN_2x32v  v;
    vboolN_4 bmask_4_0 , bmask_4_1;
    vboolN_2 blez, bmoz, b0, b1;

    s   = BBE_ZERON_2XF32();
    exp = BBE_CONSTN_2XF32(1);

    E_r = (xb_vecN_4xcf32 *)e;

    for (j = 0; j < N; j++) {
        xb_vecN_2xf32  _s, _exp;

        u1 = BBE_LVN_4XCF32_I(E_r, (BBE_SIMD_WIDTH / 4)*sz_f32c);
        BBE_LVN_4XCF32_XP(u0, E_r, stride*sz_f32c);

        h0 = BBE_MOVN_2XF32_FROMN_4XCF32(u0);
        h1 = BBE_MOVN_2XF32_FROMN_4XCF32(u1);
        h0 = BBE_SELN_2XF32I(h1, h0, BBE_SELI_EXTRACT_2_OF_4_OFF_2);
        // cimagf(e[j*stride]) < 0.f
        blez = BBE_OLTN_2XF32(h0, c0f);
        // cimagf(e[j*stride]) > 0.f
        bmoz = BBE_OLTN_2XF32(c0f, h0);

        /*
        * cimagf(e[j*stride]) >= 0.f or NaN
        */        
        _s = BBE_ZERON_2XF32();
        U_r = (xb_vecN_2xf32 *)U;
        for (i = 0; i < N; i++) {
            g = BBE_LVN_2XF32T_X(U_r, stride*sz_f32, bmoz);
            BBE_LVN_2XF32_XP(f, U_r, step);
            f = BBE_ADDN_2XF32(BBE_ABSN_2XF32(f), BBE_ABSN_2XF32(g));
            _s = BBE_MOVN_2XF32T(f, _s, BBE_OLTN_2XF32(_s, f));
        }

        // Calculate scaling coefficient
        v = BBE_MOVN_2X32_FROMN_2XF32(_s);
        v = BBE_ANDN_2X32(v, exp_mask);
        _s = BBE_MOVN_2XF32_FROMN_2X32(v);
        _exp = BBE_RECIP0N_2XF32(_s);

        _s = BBE_ZERON_2XF32();
        U_r = (xb_vecN_2xf32 *)U;
        for (i = 0; i < N; i++) {
            g = BBE_LVN_2XF32T_X(U_r, stride*sz_f32, bmoz);
            BBE_LVN_2XF32_XP(f, U_r, step);
            f = BBE_MULN_2XF32(f, _exp);
            g = BBE_MULN_2XF32(g, _exp);
            BBE_MULAN_2XF32(_s, f, f);
            BBE_MULAN_2XF32(_s, g, g); // BBE_MULAN_2XF32T(_s, g, g, bmoz);
        }
        _s = BBE_RSQRTN_2XF32(_s);
        s = BBE_MOVN_2XF32T(s, _s, blez);
        exp = BBE_MOVN_2XF32T(exp, _exp, blez);

        U_r = (xb_vecN_2xf32 *)U;
        V_w = (xb_vecN_4xcf32 *)V;

        BBE_EXTRACTBN(b1, b0, BBE_MOVN_FROMN_2(blez));
        bmask_4_0 = BBE_MOVN_4_FROMN_2(b0);
        bmask_4_1 = BBE_MOVN_4_FROMN_2(b1);
        for (i = 0; i < N; i++) {
            g = BBE_LVN_2XF32T_X(U_r, stride*sz_f32, bmoz);
            BBE_LVN_2XF32_XP(f, U_r, step);

            f = BBE_MULN_2XF32(f, exp);
            g = BBE_MULN_2XF32(g, exp); // BBE_MULN_2XF32T(g, g, _exp, bmoz);

            f = BBE_MULN_2XF32(f, s);
            g = BBE_MULN_2XF32(g, s); //  BBE_MULN_2XF32T(g, g, _s, bmoz);

            h1 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_HI);
            h0 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_LO);
            u0 = BBE_MOVN_4XCF32_FROMN_2XF32(h0);
            u1 = BBE_MOVN_4XCF32_FROMN_2XF32(h1);

            BBE_SVN_4XCF32F_I(u1, V_w, (BBE_SIMD_WIDTH / 4) * sz_f32c, bmask_4_1);
            BBE_SVN_4XCF32F_XP(u0, V_w, step_c , bmask_4_0);
        }    

        /*
        * cimagf(e[j*stride]) < 0.f
        * Use exp,s from the preceding iteration.
        */
        BBE_EXTRACTBN(b1, b0, BBE_MOVN_FROMN_2(blez));
        bmask_4_0 = BBE_MOVN_4_FROMN_2(b0);
        bmask_4_1 = BBE_MOVN_4_FROMN_2(b1);

        U_r = (xb_vecN_2xf32 *)&U[-stride];
        V_w = (xb_vecN_4xcf32 *)V;
        for (i = 0; i < N; i++) {
            g = BBE_LVN_2XF32_X(U_r, stride*sz_f32);
            BBE_LVN_2XF32_XP(f, U_r, step);

            f = BBE_MULN_2XF32(f, exp);
            g = BBE_MULN_2XF32(g, exp);

            f = BBE_MULN_2XF32(f, s);
            g = BBE_MULN_2XF32(g, s);
            g = BBE_NEGN_2XF32(g);

            h1 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_HI);
            h0 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_LO);
            u0 = BBE_MOVN_4XCF32_FROMN_2XF32(h0);
            u1 = BBE_MOVN_4XCF32_FROMN_2XF32(h1);

            BBE_SVN_4XCF32T_I(u1, V_w, (BBE_SIMD_WIDTH / 4) * sz_f32c, bmask_4_1);
            BBE_SVN_4XCF32T_XP(u0, V_w, step_c, bmask_4_0);
        }

        U += stride;
        V += stride;
    }    
#else
    const xb_vecN_2xf32  * U_r;
    const xb_vecN_4xcf32 * E_r;
          xb_vecN_4xcf32 * V_w;

    const xb_vecN_2xf32 c0f = BBE_CONSTN_2XF32(0);
    const xb_vecN_2x32v exp_mask = BBE_MOVN_2X32_FROMN_2XF32(BBE_MOVVIN_2XF32(BBE_FLOAT_MASK_EXP)); //0xFF800000
    const int step   = N*stride*sz_f32;
    const int step_c = N*stride*sz_f32c;

    xb_vecN_4xcf32 u0, u1;
    xb_vecN_2xf32  h0, h1, f, _f, g, s, exp;
    xb_vecN_2x32v  v;
    vboolN_2 blez, bmoz;

    s   = BBE_ZERON_2XF32();
    exp = BBE_CONSTN_2XF32(1);
        
    E_r = (xb_vecN_4xcf32 *)e;

    for (j = 0; j < N; j++) {
        xb_vecN_2xf32  _s, _exp;

        u1 = BBE_LVN_4XCF32_I(E_r, (BBE_SIMD_WIDTH / 4)*sz_f32c);
        BBE_LVN_4XCF32_XP(u0, E_r, stride*sz_f32c);

        h0 = BBE_MOVN_2XF32_FROMN_4XCF32(u0);
        h1 = BBE_MOVN_2XF32_FROMN_4XCF32(u1);
        h0 = BBE_SELN_2XF32I(h1, h0, BBE_SELI_EXTRACT_2_OF_4_OFF_2);
        // cimagf(e[j*stride]) < 0.f
        blez = BBE_OLTN_2XF32(h0, c0f);
        // cimagf(e[j*stride]) > 0.f
        bmoz = BBE_OLTN_2XF32(c0f, h0);
  
      
        _s = BBE_ZERON_2XF32();
        U_r = (xb_vecN_2xf32 *)U;
        for (i = 0; i < N; i++) {
            g = BBE_LVN_2XF32T_X(U_r, stride*sz_f32, bmoz);
            BBE_LVN_2XF32_XP(f, U_r, step);
            f = BBE_ADDN_2XF32(BBE_ABSN_2XF32(f), BBE_ABSN_2XF32(g));
            _s = BBE_MOVN_2XF32T(f, _s, BBE_OLTN_2XF32(_s, f));
        }

        // Calculate scaling coefficient
        v = BBE_MOVN_2X32_FROMN_2XF32(_s);
        v = BBE_ANDN_2X32(v, exp_mask);
        _exp = BBE_RECIP0N_2XF32(BBE_MOVN_2XF32_FROMN_2X32(v));
        exp = BBE_MOVN_2XF32T(exp, _exp, blez);
#if 1
        _s = BBE_ZERON_2XF32();
        U_r = (xb_vecN_2xf32 *)U;
        for (i = 0; i < N; i++) {
            g = BBE_LVN_2XF32T_X(U_r, stride*sz_f32, bmoz);
            BBE_LVN_2XF32_XP(f, U_r, step);
            f = BBE_MULN_2XF32(f, exp);
            g = BBE_MULN_2XF32(g, exp);
            BBE_MULAN_2XF32(_s, f, f);
            BBE_MULAN_2XF32(_s, g, g); // BBE_MULAN_2XF32T(_s, g, g, bmoz);
        }
#else
        _s = BBE_ZERON_2XF32();
        U_r = (xb_vecN_2xf32 *)U;
        for (i = 0; i < N; i++) {
            g = BBE_LVN_2XF32_X(U_r, stride*sz_f32);
            BBE_LVN_2XF32_XP(f, U_r, step);
            g = BBE_MOVN_2XF32T(g, c0f, bmoz);
            f = BBE_MULN_2XF32(f, exp);
            g = BBE_MULN_2XF32(g, exp);
            BBE_MULAN_2XF32(_s, f, f);
            BBE_MULAN_2XF32(_s, g, g); // BBE_MULAN_2XF32T(_s, g, g, bmoz);
        }
#endif

        _s = BBE_RSQRTN_2XF32(_s);
        s = BBE_MOVN_2XF32T(s, _s, blez);

        U_r = (xb_vecN_2xf32 *)&U[-stride];
        V_w = (xb_vecN_4xcf32 *)V;

        for (i = 0; i < N; i++) {
            _f = BBE_LVN_2XF32_X(U_r, stride*sz_f32);       // U[(i*N + j + 0)*stride]
            g = BBE_LVN_2XF32T_X(U_r, 2*stride*sz_f32, bmoz); // U[(i*N + j + 1)*stride]
            BBE_LVN_2XF32T_XP(f, U_r, step, blez); // U[(i*N + j - 1)*stride]
            g = BBE_MOVN_2XF32T(_f, g, blez);
            f = BBE_MOVN_2XF32T(f, _f, blez);

            f = BBE_MULN_2XF32(f, exp);
            g = BBE_MULN_2XF32(g, exp); // BBE_MULN_2XF32T(g, g, _exp, bmoz);

            f = BBE_MULN_2XF32(f, s);
            g = BBE_MULN_2XF32(g, s); //  BBE_MULN_2XF32T(g, g, _s, bmoz);
            BBE_NEGN_2XF32T(g, g, blez); // if cimagf(e[j*stride]) < 0.f

            h1 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_HI);
            h0 = BBE_SELN_2XF32I(g, f, BBE_SELI_INTERLEAVE_2_LO);
            u0 = BBE_MOVN_4XCF32_FROMN_2XF32(h0);
            u1 = BBE_MOVN_4XCF32_FROMN_2XF32(h1);
                
            BBE_SVN_4XCF32_I(u1, V_w, (BBE_SIMD_WIDTH / 4) * sz_f32c);
            BBE_SVN_4XCF32_XP(u0, V_w, step_c);
        }
        U += stride;
        V += stride;
    }    
#endif
} /* r_evcomb_f() */

#endif
