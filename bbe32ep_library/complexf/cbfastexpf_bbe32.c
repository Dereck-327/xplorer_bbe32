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
    Complex Exponential
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
DISCARD_FUN(void, cbfastexpf, ( complex_float * restrict z, const float32_t * restrict phase, int N ))
#else
/* Value of 2/pi, 4/pi, etc. */
#include "inv2pif_tbl.h"
/* sine/cosine approximation polynomial coeffs. */
#include "sinf_tbl.h"

/*-------------------------------------------------------------------------
Complex Exponential

Description: These functions compute complex number lying on the unit 
circle from input phase data.

Representation:
cbexp         16-bit signed fixed-point format
              Input data phase[] are Q15 angular values normalized by pi.
              Output data z[] format is Q(15-sh), where sh is the shift
              control argument
cbexpf,sexpf  IEEE-754 Std. single precision floating-point format
              Input phase data phase[] are in radians. See the Note.
              
Accuracy:
3.7e-4 (12 in Q15) for cbexp
2 ULP for cbexpf, sexpf
3 ULP for cbfastexpf

Note for (non-fast versions):
For the floating-point functions, input phase data should belong to the
range [-102940.0, 102940.0], otherwise the respective result is 0+0j.

Input domain for 'fast' version cfastbexpf()
|phase|<804.2477
The output value is not defined outside of this range or accuracy is degraded.

Parameters:
Input:
phase[N]  Phase values
sh        (cbexp only) left bit shift amount, [0..15]
N         Length of vectors
Output:
z[N]      Results

Restrictions:
z,phase   Aligned on 32-byte boundary
z,phase   Must not overlap
N         Multiple of 16 (cbexp) or 8 (cbexpf,cbfastexpf)
-------------------------------------------------------------------------*/

void cbfastexpf ( complex_float   * restrict z, const float32_t * restrict phase, int N )
{
#if 1
  const xb_vecN_2xf32 * restrict PH0;
  const xb_vecN_2xf32 * restrict PH1;
  const xb_vecN_2xf32 * restrict PH2;
        xb_vecN_2xf32 * restrict Z;
  const xtfloat       * restrict TBL_PI2;
  const xtfloat       * restrict TBL_SIN;
  const xtfloat       * restrict TBL_COS;

  xb_vecN_2xf32 phin,
                zoutre, zoutim,
                zout1, zout2;
  xb_vecN_2xf32 dx, p2, zsn, zcs;
  xb_vecN_2xf32 sn0, sn1, sn2;
  xb_vecN_2xf32 cs0, cs1, cs2;
  xb_vecN_2xf32 jf;
  vboolN_2 bdom, b_cs, b_sx, b_scs, b_ssn;
  vboolN b_t;
  xb_vecN_2x32v ji;
  xb_vecNx16 tmp0, ji16, ss, sc, neg;
  xb_vecN_2xf32 pi2fc0, pi2fc1, pi2fc2;
  int n;

  NASSERT_ALIGN(z,     (2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(phase, (2*BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH/2) == 0);
  NASSERT_ALIGN(polysinf_tbl, (2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(polycosf_tbl, (2*BBE_SIMD_WIDTH));

  PH0 = (const xb_vecN_2xf32 *)phase;
  PH1 = (const xb_vecN_2xf32 *)phase;
  PH2 = (const xb_vecN_2xf32 *)phase;
  Z = (xb_vecN_2xf32 *)z;
  TBL_PI2 = (const xtfloat *)pi2fc;
  TBL_SIN = (const xtfloat *)polysinf_tbl;
  TBL_COS = (const xtfloat *)polycosf_tbl;

  for (n = 0; n<(N>>(LOG2_BBE_SIMD_WIDTH-1)); n++)
  {
    BBE_LVN_2XF32_IP(phin, PH0, 2*BBE_SIMD_WIDTH);

    /*
     * Argument reduction.
     */

    jf = BBE_MULN_2XF32(phin, inv2pif.f);
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

    BBE_LVN_2XF32_IP(phin, PH1, 2*BBE_SIMD_WIDTH);
    BBE_MULSN_2XF32(phin, jf, pi2fc0);
    BBE_MULSN_2XF32(phin, jf, pi2fc1);
    //BBE_MULSN_2XF32(phin, jf, pi2fc2);
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

    p2 = BBE_MULN_2XF32(phin, phin);

    /* Compute polynomials by Horner's method. */
    zsn = BBE_CONSTN_2XF32(1);
    BBE_MULAN_2XF32(sn1, sn0, p2);
    BBE_MULAN_2XF32(sn2, sn1, p2);
    BBE_MULAN_2XF32(zsn, sn2, p2);

    BBE_MULAN_2XF32(dx, zsn, phin);
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
    BBE_LVN_2XF32_IP(phin, PH2, 2*BBE_SIMD_WIDTH);
    bdom = BBE_UNN_2XF32(p2, p2);
    bdom = BBE_ANDNOTBN_2(bdom, BBE_UNN_2XF32(phin, phin));

    BBE_CONSTN_2XF32T(zoutim, 0, bdom);
    BBE_CONSTN_2XF32T(zoutre, 0, bdom);

    BBE_DSELN_2XF32I(zout2, zout1, zoutim, zoutre, BBE_DSELI_INTERLEAVE_2);

    BBE_SVN_2XF32_IP(zout1, Z, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(zout2, Z, 2*BBE_SIMD_WIDTH);
  }

#else
  const xb_vecN_2xf32 * restrict PH0;
  const xb_vecN_2xf32 * restrict PH1;
  const xb_vecN_2xf32 * restrict Zld;
        xb_vecN_2xf32 * restrict Zst;
  const xtfloat       * restrict TBL_PI2;
  const xtfloat       * restrict TBL_SIN;
  const xtfloat       * restrict TBL_COS;

 xb_vecN_2xf32 phin,
               zoutre, zoutim,
               zout1, zout2;
  xb_vecN_2xf32 p2, zsn, zcs;
  xb_vecN_2xf32 sn0, sn1, sn2;
  xb_vecN_2xf32 cs0, cs1, cs2;
  xb_vecN_2xf32 jf, dx;
  vboolN_2 bdom, b_cs, b_sx, b_scs, b_ssn;
  vboolN b_t;
  xb_vecN_2x32v ji;
  xb_vecNx16 tmp0, ji16, ss, sc, neg;
  xb_vecN_2xf32 pi2fc0, pi2fc1, pi2fc2;
  int n;

  NASSERT_ALIGN(z,     (2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(phase, (2*BBE_SIMD_WIDTH));
  NASSERT(N % (BBE_SIMD_WIDTH/2) == 0);
  NASSERT_ALIGN(polysinf_tbl, (2*BBE_SIMD_WIDTH));
  NASSERT_ALIGN(polycosf_tbl, (2*BBE_SIMD_WIDTH));

  TBL_PI2 = (const xtfloat *)pi2fc;
  TBL_SIN = (const xtfloat *)polysinf_tbl;
  TBL_COS = (const xtfloat *)polycosf_tbl;

  /*
   * First stage:
   * make argument range reduction, calculate the sine approximation polynomial;
   * save reduced argument and sine to the output.
   */
  PH0 = (const xb_vecN_2xf32 *)phase;
  PH1 = (const xb_vecN_2xf32 *)phase;
  Zst = (      xb_vecN_2xf32 *)z;

  for (n = 0; n<(N>>(LOG2_BBE_SIMD_WIDTH-1)); n++)
  {
    BBE_LVN_2XF32_IP(phin, PH0, 2*BBE_SIMD_WIDTH);

    /*
     * Argument reduction.
     */

    jf = BBE_MULN_2XF32(phin, inv2pif.f);
    jf = BBE_FIROUNDN_2XF32(jf);

    pi2fc0 = BBE_LSN_2XF32_I(TBL_PI2, 0*sizeof(float32_t)); pi2fc0 = BBE_REPN_2XF32(pi2fc0, 0);
    pi2fc1 = BBE_LSN_2XF32_I(TBL_PI2, 1*sizeof(float32_t)); pi2fc1 = BBE_REPN_2XF32(pi2fc1, 0);
    pi2fc2 = BBE_LSN_2XF32_I(TBL_PI2, 2*sizeof(float32_t)); pi2fc2 = BBE_REPN_2XF32(pi2fc2, 0);

    BBE_LVN_2XF32_IP(phin, PH1, 2*BBE_SIMD_WIDTH);
    BBE_MULSN_2XF32(phin, jf, pi2fc0);
    BBE_MULSN_2XF32(phin, jf, pi2fc1);
    //BBE_MULSN_2XF32(phin, jf, pi2fc2);
    dx = BBE_MULN_2XF32(jf, pi2fc2);
    dx = BBE_NEGN_2XF32(dx);

    /*
     * Compute sine via minmax polynomial.
     */

    sn0 = BBE_LSN_2XF32_I(TBL_SIN, 0*sizeof(float32_t)); sn0 = BBE_REPN_2XF32(sn0, 0);
    sn1 = BBE_LSN_2XF32_I(TBL_SIN, 1*sizeof(float32_t)); sn1 = BBE_REPN_2XF32(sn1, 0);
    sn2 = BBE_LSN_2XF32_I(TBL_SIN, 2*sizeof(float32_t)); sn2 = BBE_REPN_2XF32(sn2, 0);

    p2 = BBE_MULN_2XF32(phin, phin);

    /* Compute polynomials by Horner's method. */
    zsn = BBE_CONSTN_2XF32(1);
    BBE_MULAN_2XF32(sn1, sn0, p2);
    BBE_MULAN_2XF32(sn2, sn1, p2);
    BBE_MULAN_2XF32(zsn, sn2, p2);

    BBE_MULAN_2XF32(dx, zsn, phin);
    zsn = dx;

    BBE_SVN_2XF32_IP(p2 , Zst, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(zsn, Zst, 2*BBE_SIMD_WIDTH);
  }

  /*
   * Second stage:
   * calculate the cosine approximation polynomial using reduced argument;
   * restore range and sign of results, check for non-domain input data.
   */
  PH0 = (const xb_vecN_2xf32 *)phase;
  Zld = (const xb_vecN_2xf32 *)z;
  Zst = (      xb_vecN_2xf32 *)z;

  for (n = 0; n<(N>>(LOG2_BBE_SIMD_WIDTH-1)); n++)
  {
    /* Load previously computed reduced argument and sine */
    BBE_LVN_2XF32_IP(p2 , Zld, 2*BBE_SIMD_WIDTH);
    BBE_LVN_2XF32_IP(zsn, Zld, 2*BBE_SIMD_WIDTH);

    /*
     * Compute cosine via minmax polynomial.
     */

    cs0 = BBE_LSN_2XF32_I(TBL_COS, 0*sizeof(float32_t)); cs0 = BBE_REPN_2XF32(cs0, 0);
    cs1 = BBE_LSN_2XF32_I(TBL_COS, 1*sizeof(float32_t)); cs1 = BBE_REPN_2XF32(cs1, 0);
    cs2 = BBE_LSN_2XF32_I(TBL_COS, 2*sizeof(float32_t)); cs2 = BBE_REPN_2XF32(cs2, 0);

    /* Compute polynomials by Horner's method. */
    BBE_MULAN_2XF32(cs1, cs0, p2);
    BBE_MULAN_2XF32(cs2, cs1, p2);
    zcs = BBE_CONSTN_2XF32(1);
    BBE_MULAN_2XF32(zcs, cs2, p2);

    /* Determine the quadrant of input phase */
    BBE_LVN_2XF32_IP(phin, PH0, 2*BBE_SIMD_WIDTH);
    jf = BBE_MULN_2XF32(phin, inv2pif.f);
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

    /* Check for a domain violation. */
    bdom = BBE_UNN_2XF32(p2, p2);
    bdom = BBE_ANDNOTBN_2(bdom, BBE_UNN_2XF32(phin, phin));

    /* restore the sign */
    BBE_NEGN_2XF32T(zsn, zsn, b_scs);
    BBE_NEGN_2XF32T(zcs, zcs, b_ssn);

    /* Select sine/cosine */
    zoutim = BBE_MOVN_2XF32T(zcs, zsn, b_cs);
    zoutre = BBE_MOVN_2XF32T(zsn, zcs, b_cs);

    BBE_CONSTN_2XF32T(zoutim, 0, bdom);
    BBE_CONSTN_2XF32T(zoutre, 0, bdom);

    BBE_DSELN_2XF32I(zout2, zout1, zoutim, zoutre, BBE_DSELI_INTERLEAVE_2);

    BBE_SVN_2XF32_IP(zout1, Zst, 2*BBE_SIMD_WIDTH);
    BBE_SVN_2XF32_IP(zout2, Zst, 2*BBE_SIMD_WIDTH);
  }

#endif
} /* cbfastexpf() */

#endif/* !HAVE_VFPU */
