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
  NatureDSP_Baseband library. Communications
    Convolution encoder
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Basic operations for the reference code. */
#include "NatureDSP_Math.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_comm.h"

#if !(HAVE_INTLV && HAVE_LFSR && 1)
DISCARD_FUN(size_t, convenc_alloc,())
DISCARD_FUN(convenc_handle_t, convenc_init,(void * objmem, 
                              int K, int R, const int16_t * restrict poly))
#else
/*---------------------------------------------------------------------------
Convolution encoder
Description: Allows to encode data with arbitrary selected polynomials with constraint length 
from 3 to 9 and the code rate up to 4.
---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------
Object allocation
Input:
none
Output:
none
Returns: size of memory in bytes to be allocated
---------------------------------------------------------------------------*/

#define MAGIC 0x756bb623
typedef struct
{
    uint32_t magic;       // Instance pointer validation number
    uint16_t  polyMask[4];// polynomials 
    int16_t  R;           // code rate (2…4)
    int16_t  K;           // constraint length (3…9)
    int16_t  state;       // state of encoder
}
tConvenc_;

size_t convenc_alloc ()
{
    return sizeof(tConvenc_) + 7;
} /* convenc_alloc() */

/*---------------------------------------------------------------------------
Object initialization

Input:
objmem  Allocated memory block
K       Constraint length (3-9)
R       Code rate (2-4)
poly[R] Polynomials in octal form
Output:
None
Returns: handle to the object
---------------------------------------------------------------------------*/

convenc_handle_t convenc_init ( void * objmem, 
                                int K, int R,
                                const int16_t * restrict poly )
{
    tConvenc_ *pEnc;
    uintptr_t a;
    a = (uintptr_t)objmem;
    a = (a + 7)&(~7);
    pEnc = (tConvenc_ *)a;
    pEnc->magic = MAGIC;
    if (R <= 2) R = 2;
    if (R >= 4) R = 4;
    pEnc->R = R;
    pEnc->K = K;
    pEnc->state = 0;
    // convert octals to masks
    {
        int32_t r, nsa;

        for (r = 0; r<R; r++)
        {
            unsigned int mask = 0;
            int16_t d, d0, d1, d2;
            d = poly[r];
            d2 = (int16_t)(L_mul_ss(d, 20972) >> 21);//(d/100);
            d = (int16_t)(d - L_mul_ss(d2, 100));
            d1 = (int16_t)(L_mul_ss(d, 26215) >> 18);//d/10;
            d0 = (int16_t)(d - L_mul_ss(d1, 10));
            mask = (d0 & 7) | ((d1 & 7) << 3) | ((d2 & 3) << 6);
            mask = ((uint16_t)mask | ((1 << (K - 1)) + 1));
            mask = mask&((1 << K) - 1);
            nsa = XT_NSA(mask);
            nsa = nsa - 16;
            mask = XT_SSL_SLL(mask, nsa);
            pEnc->polyMask[r] = (uint16_t)mask;
        }
    }
    return (convenc_handle_t)pEnc;
} /* convenc_init() */
#endif
