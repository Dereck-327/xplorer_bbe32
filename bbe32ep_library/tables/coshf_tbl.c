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
    tables for coshf(x) approximation
*/

#include "NatureDSP_types.h"
#include "coshf_tbl.h"
#include "common.h"

/* Maximum magnitude of input value such that coshf does not overflow. */
const union ufloat32uint32 coshf_maxarg = { 0x42b2d4fc }; /* 89.41598510742188 */

#if COSHF_ALG==0
/* polynomial coefficients for cosh(x) in range +/-0.25 

   x=(-0.25:pow2(2,-16):0.25); y=cosh(x);
   polyfit(x,y,4)
   p=p(1:2:end); p(end)=[];
*/
const union ufloat32uint32 ALIGN(32) polycoshf_tbl[]=
{
        {0x3d2b26f4}, /*4.1785194725e-002*/
        {0x3effffad}, /*4.9999752938e-001*/
};
#endif
