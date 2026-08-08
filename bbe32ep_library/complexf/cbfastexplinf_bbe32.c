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
  NatureDSP_Baseband library. Complex Math functions
    Complex Exponential with Linearly Evolving Phase
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Complex Math Functions. */
#include "NatureDSP_Baseband_complex.h"

#if !(HAVE_VFPU)
DISCARD_FUN(float32_t, cbfastexplinf, ( complex_float   * restrict z,  float32_t phase0, float32_t delta, int N ))
#else
/* Value of 2/pi, 4/pi, etc. */
#include "inv2pif_tbl.h"
/* sine/cosine approximation polynomial coeffs. */
#include "sinf_tbl.h"

/*-------------------------------------------------------------------------
Complex Exponential with Linearly Evolving Phase

Description: These functions compute a series of complex numbers lying on
the unit circle, with the phase angle being incremented for each successive
complex number.

Representation:
cbexplin,           16-bit signed fixed-point format
cbexplin_fast       Initial phase and phase increment (phase0 and delta) are 
                    Q15 angular values normalized by pi.
                    Output data z[] format is Q15.
cbexplinf,          IEEE-754 Std. single precision floating-point format
cexplinf            Initial phase and phase increment (phase0 and delta) are
                    in radians. See Note 3.

Accuracy:
cbexplin            Absolute error for re/im components does not exceed
                    3.7e-4 (12 in Q15)
cbexplin_fast       Absolute error for re/im components for n-th result is
                    ~(floor(n/8)+1)*3.7e-4, n=0..N-1
cbexplinf,cexplinf  2 ULP 
cbfastexplinf       3 ULP

Notes:
1. cbexplin() function computes each complex value from a full-precision
   linearly evolving angle, thus the function is able to generate long data
   sequences with no risk of phase errors accumulation.
2. cbexplin_fast() function is approximately twice faster but less accurate 
   because it performs successive rotations by small angle on each iteration.
   Avoid using this function for long sequences or tiny phase increments.
3. For the non-fast floating-point functions (cbexplinf, cexplinf) , the 
   phase (initial and updated) should belong to the range [-102940.0, 
   102940.0], otherwise the respective resultis 0+0j.
4. Scalar function cexplinf uses the initial phase angle (phase0) to compute 
   a single complex number, and returns the once incremented phase value. 

Input domain for 'fast' version (cbfastexplinf)
|phase|<804.2477
for all initial and updated phases 
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
phase0  Initial phase angle (phase of the first complex number in a series)
delta   Phase increment step
N       Length of output vector
Output:
z[N]    Complex numbers
Returned Value:
        Updated phase, phase0+N*delta

Restrictions:
z       Aligned on 32-byte boundary
z       Must not overlap
N       Multiple of 16 (cbexplln, cbexplln_fast) or 8 (cbexplinf,cbfastexplinf)
-------------------------------------------------------------------------*/

float32_t cbfastexplinf ( complex_float   * restrict z,  float32_t phase0, float32_t delta, int N )
{
  const xb_vecN_2xf32 * restrict PH0;
  const xb_vecN_2xf32 * restrict PH1;
  const xb_vecN_2xf32 * restrict PH2;
        xb_vecN_2xf32 * restrict Z;
  const xtfloat       * restrict TBL_PI2;
  const xtfloat       * restrict TBL_SIN;
  const xtfloat       * restrict TBL_COS;

  xb_vecN_2xf32 ph, ph_, phInc0, phInc1,
                zoutre, zoutim,
                zout1, zout2;
  xb_vecN_2xf32 p2, zsn, zcs;
  xb_vecN_2xf32 sn0, sn1, sn2;
  xb_vecN_2xf32 cs0, cs1, cs2;
  xb_vecN_2xf32 jf, dx;
  vboolN_2 bdom, b_cs, b_sx, b_scs, b_ssn;
  vboolN b_t;
  xb_vecN_2x32v ji, sn;
  xb_vecNx16 tmp0, ji16, ss, sc, neg;
  xb_vecN_2xf32 pi2fc0, pi2fc1, pi2fc2;
  float32_t p;
  int n;

  NASSERT_ALIGN(z, (2*BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH/2) == 0);
  NASSERT_ALIGN(polysinf_tbl, (2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(polycosf_tbl, (2*BBE_SIMD_WIDTH));
  if (N <= 0) return phase0;

  TBL_PI2 = (const xtfloat *)pi2fc;
  TBL_SIN = (const xtfloat *)polysinf_tbl;
  TBL_COS = (const xtfloat *)polycosf_tbl;

  /*
   * First stage:
   * compute phase for all N complex values and save it to the output.
   */

  Z = (xb_vecN_2xf32 *)z;

  /* Initialize phase */
  tmp0 = BBE_SELNX16I(BBE_ZERONX16(), BBE_SEQNX16(), BBE_SELI_INTERLEAVE_1_LO);
  sn = BBE_MOVN_2X32_FROMNX16(tmp0);
  phInc0 = BBE_FLOATN_2X32(sn, 0);
  phInc1 = BBE_ADDN_2XF32(phInc0, 8.0f);

  b_sx = BBE_LTRN_2I(1);
  b_sx = BBE_NOTBN_2(b_sx);
  ph = BBE_CONSTN_2XF32(0);
  BBE_MULN_2XF32T(ph, phInc0, delta, b_sx);
  ph_ = BBE_MULN_2XF32(phInc1, delta);
  ph  = BBE_ADDN_2XF32(ph, phase0);
  ph_ = BBE_ADDN_2XF32(ph_, phase0); 

  for (n = 0; n<(N>>(LOG2_BBE_SIMD_WIDTH)); n++)
  {
    BBE_SVN_2XF32_IP(ph , Z, 4*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(ph_, Z, 4*BBE_SIMD_WIDTH);

    /* Increment phase */
    phInc0 = BBE_ADDN_2XF32(phInc0, 16.0f);
    ph = BBE_MULN_2XF32(phInc0, delta);
    ph = BBE_ADDN_2XF32(ph, phase0); 
    phInc1 = BBE_ADDN_2XF32(phInc1, 16.0f);
    ph_ = BBE_MULN_2XF32(phInc1, delta);
    ph_ = BBE_ADDN_2XF32(ph_, phase0); 
  }
  if (N&8)
  {
    BBE_SVN_2XF32_IP(ph, Z, 4*BBE_SIMD_WIDTH);
    ph = ph_;
  }
  if (N&4)
  {
    BBE_SVN_2XF32_I(ph, Z, 0);
  }

  __Pragma("no_reorder");

  /*
   * Second stage:
   * compute complex values from previously computed phases.
   */
  PH0 = (const xb_vecN_2xf32 *)z;
  PH1 = (const xb_vecN_2xf32 *)z;
  PH2 = (const xb_vecN_2xf32 *)z;
  Z   = (xb_vecN_2xf32 *)z;

  for (n = 0; n<(N>>(LOG2_BBE_SIMD_WIDTH-1)); n++)
  {
    BBE_LVN_2XF32_IP(ph, PH0, 4*BBE_SIMD_WIDTH);

    /*
     * Argument reduction.
     */

    jf = BBE_MULN_2XF32(ph, inv2pif.f);
    jf = BBE_FIROUNDN_2XF32(jf);
    ji = BBE_TRUNCN_2XF32(jf, 0);
    ji16 = BBE_MOVNX16_FROMN_2X32(ji);

    ss = ji16;
    sc = BBE_ADDNX16(ji16, 1);
    neg = BBE_SELNX16I(ss, sc, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    neg = BBE_SLLINX16(neg, 14);
    tmp0 = BBE_SLLINX16(neg, 1);
    b_t = BBE_LTNX16(neg, BBE_ZERONX16());
    BBE_EXTRACTBN(b_ssn, b_scs, b_t);
    b_t = BBE_LTNX16(tmp0, BBE_ZERONX16());
    BBE_EXTRACTBN(b_cs, b_sx, b_t);

    pi2fc0 = BBE_LSN_2XF32_I(TBL_PI2, 0*sizeof(float32_t)); pi2fc0 = BBE_REPN_2XF32(pi2fc0, 0);
    pi2fc1 = BBE_LSN_2XF32_I(TBL_PI2, 1*sizeof(float32_t)); pi2fc1 = BBE_REPN_2XF32(pi2fc1, 0);
    pi2fc2 = BBE_LSN_2XF32_I(TBL_PI2, 2*sizeof(float32_t)); pi2fc2 = BBE_REPN_2XF32(pi2fc2, 0);

    BBE_LVN_2XF32_IP(ph, PH1, 4*BBE_SIMD_WIDTH);
    BBE_MULSN_2XF32(ph, jf, pi2fc0);
    BBE_MULSN_2XF32(ph, jf, pi2fc1);
    //BBE_MULSN_2XF32(ph, jf, pi2fc2);
    dx = BBE_MULN_2XF32(jf, pi2fc2);
    dx = BBE_NEGN_2XF32(dx);

    /*
    * Compute sine/cosine via minmax polynomial.
    */

    sn0 = BBE_LSN_2XF32_I(TBL_SIN, 0*sizeof(float32_t)); sn0 = BBE_REPN_2XF32(sn0, 0);
    sn1 = BBE_LSN_2XF32_I(TBL_SIN, 1*sizeof(float32_t)); sn1 = BBE_REPN_2XF32(sn1, 0);
    sn2 = BBE_LSN_2XF32_I(TBL_SIN, 2*sizeof(float32_t)); sn2 = BBE_REPN_2XF32(sn2, 0);

    cs0 = BBE_LSN_2XF32_I(TBL_COS, 0*sizeof(float32_t)); cs0 = BBE_REPN_2XF32(cs0, 0);
    cs1 = BBE_LSN_2XF32_I(TBL_COS, 1*sizeof(float32_t)); cs1 = BBE_REPN_2XF32(cs1, 0);
    cs2 = BBE_LSN_2XF32_I(TBL_COS, 2*sizeof(float32_t)); cs2 = BBE_REPN_2XF32(cs2, 0);

    p2 = BBE_MULN_2XF32(ph, ph);

    /* Compute polynomials by Horner's method. */
    zsn = BBE_CONSTN_2XF32(1);
    BBE_MULAN_2XF32(sn1, sn0, p2);
    BBE_MULAN_2XF32(sn2, sn1, p2);
    BBE_MULAN_2XF32(zsn, sn2, p2);

    BBE_MULAN_2XF32(dx, zsn, ph);
    zsn = dx;
    
    zcs = BBE_CONSTN_2XF32(1);
    BBE_MULAN_2XF32(cs1, cs0, p2);
    BBE_MULAN_2XF32(cs2, cs1, p2);
    BBE_MULAN_2XF32(zcs, cs2, p2);

    /* Select sine/cosine */
    zoutim = BBE_MOVN_2XF32T(zcs, zsn, b_cs);
    zoutre = BBE_MOVN_2XF32T(zsn, zcs, b_cs);

    /* restore the sign */
    BBE_NEGN_2XF32T(zoutim, zoutim, b_ssn);
    BBE_NEGN_2XF32T(zoutre, zoutre, b_scs);

    /* Check for a domain violation. */
    BBE_LVN_2XF32_IP(ph, PH2, 4*BBE_SIMD_WIDTH);
    bdom = BBE_UNN_2XF32(p2, p2);
    bdom = BBE_ANDNOTBN_2(bdom, BBE_UNN_2XF32(ph, ph));

    BBE_CONSTN_2XF32T(zoutim, 0, bdom);
    BBE_CONSTN_2XF32T(zoutre, 0, bdom);

    BBE_DSELN_2XF32I(zout2, zout1, zoutim, zoutre, BBE_DSELI_INTERLEAVE_2);

    BBE_SVN_2XF32_IP(zout1, Z, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(zout2, Z, 2*BBE_SIMD_WIDTH);
  }

#if 1
  /* Process last 4 values */
  if (N&4)
  {
    ph = BBE_LVN_2XF32_I(PH2, 0);

    /*
     * Argument reduction.
     */

    jf = BBE_MULN_2XF32(ph, inv2pif.f);
    jf = BBE_FIROUNDN_2XF32(jf);
    ji = BBE_TRUNCN_2XF32(jf, 0);
    ji16 = BBE_MOVNX16_FROMN_2X32(ji);

    ss = ji16;
    sc = BBE_ADDNX16(ji16, 1);
    neg = BBE_SELNX16I(ss, sc, BBE_SELI_EXTRACT_1_OF_2_OFF_0);
    neg = BBE_SLLINX16(neg, 14);
    tmp0 = BBE_SLLINX16(neg, 1);
    b_t = BBE_LTNX16(neg, BBE_ZERONX16());
    BBE_EXTRACTBN(b_ssn, b_scs, b_t);
    b_t = BBE_LTNX16(tmp0, BBE_ZERONX16());
    BBE_EXTRACTBN(b_cs, b_sx, b_t);

    pi2fc0 = BBE_LSN_2XF32_I(TBL_PI2, 0*sizeof(float32_t)); pi2fc0 = BBE_REPN_2XF32(pi2fc0, 0);
    pi2fc1 = BBE_LSN_2XF32_I(TBL_PI2, 1*sizeof(float32_t)); pi2fc1 = BBE_REPN_2XF32(pi2fc1, 0);
    pi2fc2 = BBE_LSN_2XF32_I(TBL_PI2, 2*sizeof(float32_t)); pi2fc2 = BBE_REPN_2XF32(pi2fc2, 0);

    BBE_MULSN_2XF32(ph, jf, pi2fc0);
    BBE_MULSN_2XF32(ph, jf, pi2fc1);
    //BBE_MULSN_2XF32(ph, jf, pi2fc2);
    dx = BBE_MULN_2XF32(jf, pi2fc2);
    dx = BBE_NEGN_2XF32(dx);

    /*
    * Compute sine/cosine via minmax polynomial.
    */

    sn0 = BBE_LSN_2XF32_I(TBL_SIN, 0*sizeof(float32_t)); sn0 = BBE_REPN_2XF32(sn0, 0);
    sn1 = BBE_LSN_2XF32_I(TBL_SIN, 1*sizeof(float32_t)); sn1 = BBE_REPN_2XF32(sn1, 0);
    sn2 = BBE_LSN_2XF32_I(TBL_SIN, 2*sizeof(float32_t)); sn2 = BBE_REPN_2XF32(sn2, 0);

    cs0 = BBE_LSN_2XF32_I(TBL_COS, 0*sizeof(float32_t)); cs0 = BBE_REPN_2XF32(cs0, 0);
    cs1 = BBE_LSN_2XF32_I(TBL_COS, 1*sizeof(float32_t)); cs1 = BBE_REPN_2XF32(cs1, 0);
    cs2 = BBE_LSN_2XF32_I(TBL_COS, 2*sizeof(float32_t)); cs2 = BBE_REPN_2XF32(cs2, 0);

    p2 = BBE_MULN_2XF32(ph, ph);

    /* Compute polynomials by Horner's method. */
    zsn = BBE_CONSTN_2XF32(1);
    BBE_MULAN_2XF32(sn1, sn0, p2);
    BBE_MULAN_2XF32(sn2, sn1, p2);
    BBE_MULAN_2XF32(zsn, sn2, p2);

    BBE_MULAN_2XF32(dx, zsn, ph);
    zsn = dx;
    
    zcs = BBE_CONSTN_2XF32(1);
    BBE_MULAN_2XF32(cs1, cs0, p2);
    BBE_MULAN_2XF32(cs2, cs1, p2);
    BBE_MULAN_2XF32(zcs, cs2, p2);

    /* Select sine/cosine */
    zoutim = BBE_MOVN_2XF32T(zcs, zsn, b_cs);
    zoutre = BBE_MOVN_2XF32T(zsn, zcs, b_cs);

    /* restore the sign */
    BBE_NEGN_2XF32T(zoutim, zoutim, b_ssn);
    BBE_NEGN_2XF32T(zoutre, zoutre, b_scs);

    /* Check for a domain violation. */
    ph = BBE_LVN_2XF32_I(PH1, 0);
    bdom = BBE_UNN_2XF32(p2, p2);
    bdom = BBE_ANDNOTBN_2(bdom, BBE_UNN_2XF32(ph, ph));

    BBE_CONSTN_2XF32T(zoutim, 0, bdom);
    BBE_CONSTN_2XF32T(zoutre, 0, bdom);

    zout1 = BBE_SELN_2XF32I(zoutim, zoutre, BBE_SELI_INTERLEAVE_2_LO);

    BBE_SVN_2XF32_I(zout1, Z, 0);
  }
#endif

  /* Compute updated phase */
  p = XT_MUL_S((float32_t)N, delta);
  return (float32_t)XT_ADD_S(phase0, p);
} /* cbfastexplinf() */

#endif/* !HAVE_VFPU */
