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
#include "NatureDSP_types.h"
#include "common.h"
#include "lognTbl.h"

/* table for natural logariphm
   tablated coefficients for 2-nd order  polynomials 
*/

// Taylor series expansion for the natural logarithm.
// MATLAB code:
// l2_x = 0.5+(1:2:31)/64;
// le_0_q11 = floor(0.5+2^11*log(l2_x));
// le_1_q14 = floor(0.5+(2^14*1./l2_x*1/64)*2);
// le_2_q17 = floor(0.5-(2^17*1./(l2_x.^2)*1/2*1/64^2)*4);
//

const int16_t ALIGN(32) logNatTbl[3 * 16] =
{
  -1357,-1236,-1122,-1014,-912,-814,-721, -632,-547,-465,-386,-310,-237,-167,-98,-32, /* Order 0 derivative, Q11 */
    993,  936,  886,  840, 799, 762, 728,  697, 669, 643, 618, 596, 575, 555,537,520, /* Order 1 derivative, Q14 */
   -241, -214, -191, -172,-156,-142,-129, -119,-109,-101, -93, -87, -81, -75,-70,-66  /* Order 2 derivative, Q17 */
};
