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
#include "log2Tbl.h"

/* table for log2 
   tablated coefficients for 2-nd order  polynomials 
*/

const int16_t ALIGN(32) log2Tbl[3*16]=
{
  -1957,-1783,-1619,-1463,-1316,-1175,-1041, -912,-789,-671,-557,-448,-342,-240,-142,-47, /* Q11 */
   1433, 1351, 1278, 1212, 1153, 1099, 1051, 1006, 965, 927, 892, 860, 829, 801, 775,750, /* Q15 */
   -347, -309, -276, -249, -225, -205, -187, -171,-158,-145,-135,-125,-116,-109,-102,-95  /* Q19 */
};
