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
    Left-multiply real square upper-Hessenberg matrices H by orthogonal 
    matrices P and store resulting matrices to U: U <- P*H. 
    Real Data, Stream Order
    C code optimized for BBE32 with VFPU
    IntegrIT, 2006-2017
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility declarations. */
#include "common.h"
/* Eigenvalues and eigenvectors common declarations. */
#include "eigen_common.h"

#if HAVE_VFPU

/* Index of (i,j)-th element of an NxN upper-Hessenberg matrix stored
 * in compact packed format. Compactness implies that zeros below the
 * first subdiagonal aren't actually stored in memory. */
#define HIDX(i,j)   ( (i)*(N) + (i)*(1-(i))/2 + (j) )

#define sz_f32    sizeof(float32_t)

/*
 * Left-multiply real square upper-Hessenberg matrices H by orthogonal 
 * matrices P and store resulting matrices to U: U <- P*H. 
 * Note:
 *   TBD If balancing is implemented, the complexity of this functions may
 *       be reduced by taking into account the block diagonal structure of
 *       input matrices. See the MATLAB reference.
 * Input:
 *   N      Matrix size
 *   L      Number of matrices
 *   P[]    Orthogonal transformation matrices
 *   H[]    Upper-Hessenberg matrices. Zero elements below the first 
 *          subdiagonal are not stored
 * Output
 *   U[]    Resulting matrices
 * Restrictions:
 *   U,P,H  Must not overlap and must be aligned on 2*BBE_SIMD_WIDTH-byte
 *          boundary
 *   Variant functions may impose additional restrictions
 */

/* Real Data, Stream Order
 * Restrictions:
 *   L  Must be a multiple of BBE_SIMD_WIDTH/2 */
void reigen_hmulp_nxnsf ( float32_t * restrict U, /* U[N*N][L]         */
                    const float32_t * restrict P, /* P[N*N][L]         */
                    const float32_t * restrict H, /* H[N*(N+3)/2-1][L] */
                    int N, int L )
#if 0
{
  int i,j,k,p,J;
  float32_t s;
  NASSERT_ALIGN( U, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( P, 2*BBE_SIMD_WIDTH );
  NASSERT_ALIGN( H, 2*BBE_SIMD_WIDTH );
  NASSERT( 0==(L%(BBE_SIMD_WIDTH/2)) );

  for ( k=0; k<L; k++ ) {
    for ( i=0; i<N; i++ ) {
      for ( j=0; j<N; j++ ) {
        s = 0.f; J = ( j<N-1 ? j+1 : N-1 );
        for ( p=0; p<=J; p++ )
          s += P[(i*N+p)*L]*H[HIDX(p,j)*L];
        U[(i*N+j)*L] = s;
      }
    }
    U++; P++; H++;
  }


} /* reigen_hmulp_nxnsf() */
#else
{
    int i, j, k, p, step;
    
    const xb_vecN_2xf32 * P_r;
    const xb_vecN_2xf32 * H_r;
          xb_vecN_2xf32 * U_w;

    NASSERT_ALIGN(U, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(P, 2 * BBE_SIMD_WIDTH);
    NASSERT_ALIGN(H, 2 * BBE_SIMD_WIDTH);
    NASSERT(0 == (L % (BBE_SIMD_WIDTH / 2)));
#if 0
    for (k = 0; k<L/(BBE_SIMD_WIDTH / 2); k++) {
        for (i = 0; i<N; i++) {
            U_w = (xb_vecN_2xf32*)&U[i*N*L];
            
            for (j = 0; j<N; j++) {
                int J;
                int step = (N - 1)*L;
                xb_vecN_2xf32 u, v;
                xb_vecN_2xf32 s;

                P_r = (xb_vecN_2xf32*)&P[i*N*L];
                H_r = (xb_vecN_2xf32*)&H[j*L];

                BBE_LVN_2XF32_XP(u, P_r, L*sz_f32);
                BBE_LVN_2XF32_XP(v, H_r, L*N*sz_f32);
                s = BBE_MULN_2XF32(u, v);

                J = (j<N - 1 ? j : N - 2);

                for (p = 0; p <= J; p++) {
                    BBE_LVN_2XF32_XP(u, P_r, L*sz_f32);
                    BBE_LVN_2XF32_XP(v, H_r, step*sz_f32);
                    BBE_MULAN_2XF32(s, u, v);
                    step-=L;
                }
                BBE_SVN_2XF32_XP(s, U_w, L*sz_f32);
            }
        }
        U += (BBE_SIMD_WIDTH / 2); P += (BBE_SIMD_WIDTH / 2); H += (BBE_SIMD_WIDTH / 2);
    }
#else 
    for (k = 0; k<L / (BBE_SIMD_WIDTH / 2); k++) {
        for (i = 0; i<N; i++) {
            xb_vecN_2xf32 u, v;
            xb_vecN_2xf32 s;
            U_w = (xb_vecN_2xf32*)&U[i*N*L];
            for (j = 0; j<N - 1; j++) {
                step = (N - 1)*L;

                P_r = (xb_vecN_2xf32*)&P[i*N*L];
                H_r = (xb_vecN_2xf32*)&H[j*L];

                BBE_LVN_2XF32_XP(u, P_r, L*sz_f32);
                BBE_LVN_2XF32_XP(v, H_r, L*N*sz_f32);
                s = BBE_MULN_2XF32(u, v);

                for (p = 0; p <= j; p++) {
                    BBE_LVN_2XF32_XP(u, P_r, L*sz_f32);
                    BBE_LVN_2XF32_XP(v, H_r, step*sz_f32);
                    BBE_MULAN_2XF32(s, u, v);
                    step -= L;
                }
                BBE_SVN_2XF32_XP(s, U_w, L*sz_f32);
            }

            step = (N - 1)*L;
            P_r = (xb_vecN_2xf32*)&P[i*N*L];
            H_r = (xb_vecN_2xf32*)&H[j*L];

            BBE_LVN_2XF32_XP(u, P_r, L*sz_f32);
            BBE_LVN_2XF32_XP(v, H_r, L*N*sz_f32);
            s = BBE_MULN_2XF32(u, v);

            for (p = 0; p <= (N-2); p++) {
                BBE_LVN_2XF32_XP(u, P_r, L*sz_f32);
                BBE_LVN_2XF32_XP(v, H_r, step*sz_f32);
                BBE_MULAN_2XF32(s, u, v);
                step -= L;
            }
            BBE_SVN_2XF32_XP(s, U_w, L*sz_f32);
        }
        U += (BBE_SIMD_WIDTH / 2); P += (BBE_SIMD_WIDTH / 2); H += (BBE_SIMD_WIDTH / 2);
    }
#endif
} /* reigen_hmulp_nxnsf() */
#endif

#endif /* HAVE_VFPU */
