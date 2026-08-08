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

#ifndef __BASEOP_H__
#define __BASEOP_H__

/* ------------------------------------------------------------------------ */
/*          Basic operations -- needed for reference code only              */
/* ------------------------------------------------------------------------ */

#include "NatureDSP_types.h"

typedef struct tagInt128_s
{
    uint64_t lo;
    int64_t hi;
} int128_t;


/* saturate Q31 to 32 bit */
inline_ int32_t satQ31(int64_t x)   {
    if (x > (int64_t)0x000000007FFFFFFFLL) x = 0x7FFFFFFF;
    if (x < (int64_t)0xFFFFFFFF80000000LL) x = 0x80000000;
    return (int32_t)x;
}

/* saturate 64 bit to 40 bit */
inline_ int32_t sat40(int64_t x)   {
    if (x > (int64_t)0x0000007FFFFFFFFFLL) x = 0x7FFFffFFFF;
    if (x < (int64_t)0xFFFFFF8000000000LL) x = 0x8000000000;
    return (int32_t)x;
}

/* saturate double to 64 bit */
inline_ int64_t satQ63(double x)   {
    int64_t z = (int64_t)x;
    if (x >= ((double)0x7FFFFFFFFFFFFFFFLL)) z = (int64_t)0x7FFFFFFFFFFFFFFFLL;
    if (x < (-(double)0x7FFFFFFFFFFFFFFFLL)) z = (int64_t)0x8000000000000000UL;
    return z;
}

// saturate 40-bit register to 20-bit
inline_ int32_t L_sature_20(int64_t A)
{
    A <<= 24;
    A >>= 24;
    if (A < (-1LL << 19)) A = (-1LL << 19);
    if (A > (1LL << 19) - 1) A = (1LL << 19) - 1;
    return (int32_t)A;
}

/* saturate 40-bit register to 32-bit */
inline_ int32_t L_sature_40(int64_t A)
{
    A <<= 24;
    A >>= 24;
    if (A > MAX_INT32) A = MAX_INT32;
    if (A < MIN_INT32) A = MIN_INT32;
    return (int32_t)A;
}

/* saturate 40-bit register to 16-bit */
inline_ int16_t S_sature_40(int64_t A)
{
    A <<= 24;
    A >>= 24;
    if (A > MAX_INT16) A = MAX_INT16;
    if (A < MIN_INT16) A = MIN_INT16;
    return (int16_t)A;
}

/* saturate 64-bit register to 32-bit */
inline_ int32_t L_sature_64(int64_t A)
{
    if (A > MAX_INT32) A = MAX_INT32;
    if (A < MIN_INT32) A = MIN_INT32;
    return (int32_t)A;
}

/* saturate 64-bit register to 16-bit */
inline_ int16_t S_sature_64(int64_t A)
{
    if (A > MAX_INT16) A = MAX_INT16;
    if (A < MIN_INT16) A = MIN_INT16;
    return (int16_t)A;
}

/* truncate 32-bit register to 20-bit */
inline_ int32_t L_trunc_20(int32_t A)
{
    A <<= 12;
    A >>= 12;
    return A;
}

/* truncate 32-bit register to 18-bit */
inline_ int32_t L_trunc_18(int32_t A)
{
    A <<= 14;
    A >>= 14;
    return A;
}

/* 18x18-bit multiplication with 36-bit result  */
inline_ int64_t LL_mul_18(int32_t x, int32_t y)
{
    return ((int64_t)((x << 14) >> 14) *
        ((y << 14) >> 14));
}

/* 18x18-bit fractional multiplication with asymmetric rounding,
   20-bit result: Q4.15 <- Q2.15*Q2.15 - 15 */
inline_ int32_t L_mul_18_packq(int32_t x, int32_t y)
{
    int64_t z;

    z = (int64_t)((x << 14) >> 14)*((y << 14) >> 14);

    return (int32_t)((z + (1LL << 14)) >> 15);
}

/* 16x16-bit multiplication with 32-bit result in a 40-bit register */
inline_ int64_t LL_mul_16(int16_t x, int16_t y)
{
    return ((int64_t)x * y);
}

/* fractional multiple Q31xQ15 forming Q16.47 result */
inline_ int64_t mpy_Q31Q15(int32_t x, int16_t y)
{
    return (((int64_t)x*y) << 1);
}
inline_ int64_t mul_Q31Q15(int32_t x, int16_t y)
{
    return (((int64_t)x*y));
}

/* Fractional mupltiplication Q31xQ31 forming Q63 result */
inline_ int64_t mpy_Q31Q31(int32_t x, int32_t y)
{
    int64_t z = (int64_t)x*y;

    if (x == (int32_t)0x80000000 && y == (int32_t)0x80000000)
    {
        return (0x7fffffffffffffffLL);
    }

    return (z << 1);
}

/* 64-bit addition with saturation. */
inline_ int64_t add_64(int64_t x, int64_t y)
{
    int64_t z;

    z = x + y;

    if ((((x^y) & MIN_INT64) == 0LL)
        &&
        (((z^x) & MIN_INT64) != 0LL))
    {
        z = (x < 0 ? MIN_INT64 : MAX_INT64);
    }

    return z;
}

/* 64-bit subtraction with saturation. */
inline_ int64_t sub_64(int64_t x, int64_t y)
{
    int64_t z;

    z = x - y;

    if ((((x^y) & MIN_INT64) != 0LL)
        &&
        (((z^x) & MIN_INT64) != 0LL))
    {
        z = (x < 0 ? MIN_INT64 : MAX_INT64);
    }

    return z;
}

/* Reverse the bit order for bit positions 0..(nbits-1) */
inline_ int bitrev(int x, int nbits)
{
    int y, n;

    for (y = 0, n = 1; n <= nbits; n++, x >>= 1)
    {
        y |= ((x & 1) << (nbits - n));
    }

    return (y);

} /* bitrev() */

/* 16-bit fractional complex multiplication:
 * CQ15 <- CQ15*CQ15 - 15  with rounding, saturation */
inline_ complex_fract16 mpy_CQ15CQ15(complex_fract16 x, complex_fract16 y)
{
    complex_fract16 z;

    z.s.re = S_sature_64(((int64_t)x.s.re*y.s.re - (int64_t)x.s.im*y.s.im + (1 << 14)) >> 15);
    z.s.im = S_sature_64(((int64_t)x.s.im*y.s.re + (int64_t)x.s.re*y.s.im + (1 << 14)) >> 15);

    return (z);
}

/* 32-bit fractional complex multiplication:
 * CQ31 <- CQ31*CQ31 - 31  with rounding, saturation */
inline_ complex_fract32 mpy_CQ31CQ31(complex_fract32 x, complex_fract32 y)
{
    int64_t z_re, z_im;
    complex_fract32 z;

    z_re = sub_64((int64_t)x.s.re*y.s.re, (int64_t)x.s.im*y.s.im);
    z_im = add_64((int64_t)x.s.im*y.s.re, (int64_t)x.s.re*y.s.im);

    z.s.re = L_sature_64(add_64(z_re, (1LL << 30)) >> 31);
    z.s.im = L_sature_64(add_64(z_im, (1LL << 30)) >> 31);

    return (z);

} /* mpy_CQ31CQ31() */

/* truncate Q63 to Q31 */
inline_ int32_t roundQ63_Q31(int64_t x)   { return (int32_t)(add_64(x, ((int64_t)1) << 31) >> 32); }

/* extract from "NatureDSP_Math.h" */
/* shift left with saturation: Qx->Q(x+t)         */
inline_ int32_t L_shl_l(int32_t x, int16_t t)
{
    int32_t a;
    if (t > 0)
    {
        a = (int32_t)(x >> (31 - t));
        if ((a == 0) || (a == -1))
            x <<= t;
        else
            x = (a > 0) ? MAX_INT32 : MIN_INT32;
    }
    else
    {
        if (t < -31)    t = -31;
        x >>= (-t);
    }
    return x;
}

/* getting exponent:exp(Q31), return 31 if x==0                  */
inline_ int16_t S_exp0_l(int32_t x)
{
    int16_t z = 0;
    if (x == 0)  return 0x1F;
    while ((int32_t)(x ^ (x << 1)) > 0) //MSB != MSB-1
    {
        x <<= 1;
        z++;
    }
    return z;
}

/* getting exponent:exp(Q31), return  0 if x==0                  */
inline_ int16_t S_exp_l(int32_t x)
{
    int16_t z = 0;
    if (x == 0)    return 0;
    while ((int32_t)(x ^ (x << 1)) > 0) /* MSB != MSB-1 */
    {
        x <<= 1;
        z++;
    }
    return z;
}

/* misc conversions  saturation: Q16.15->Q15        */
inline_ int16_t S_sature_l(int32_t x)
{
    return (int16_t)(L_shl_l(x, 16) >> 16);
}

/* misc conversions  saturation: Q31->Q15 (extract higher 16 bits) */
inline_ int16_t S_extract_l(int32_t x)
{
    return (int16_t)(x >> 16);
}

/* misc conversions  saturation: Q15->Q31 (shift left by 16) */
inline_ int32_t L_deposit_s(int16_t x)
{
    return ((int32_t)x) << 16;
}

/* addition with saturation: Q31+Q31->Q31        */
inline_ int32_t L_add_ll(int32_t x, int32_t y)
{
    int32_t z;

    z = x + y;

    if ((((x ^ y) & MIN_INT32) == 0L)
        && ((z ^ x) & MIN_INT32))
        z = (x < 0) ? MIN_INT32 : MAX_INT32;
    return z;
}

/* addition with saturation: Q15+Q15->Q31        */
inline_ int32_t L_add_ss(int16_t x, int16_t y)
{
    return L_add_ll(L_deposit_s(x), L_deposit_s(y));
}

/* addition with saturation: Q15+Q15->Q15        */
inline_ int16_t S_add_ss(int16_t x, int16_t y)
{
    return S_sature_l(((int32_t)x) + y);
}

/* subtraction with saturation  Q31-Q31->Q31      */
inline_ int32_t L_sub_ll(int32_t x, int32_t y)
{
    int32_t z;

    z = x - y;
    if (((x ^ y) & MIN_INT32) &&
        ((z ^ x) & MIN_INT32))
    {
        z = (x < 0) ? MIN_INT32 : MAX_INT32;
    }
    return z;
}

/* subtraction with saturation  Q15-Q15->Q31      */
inline_ int32_t L_sub_ss(int16_t x, int16_t y)
{
    return L_sub_ll(L_deposit_s(x), L_deposit_s(y));
}

/* subtraction with saturation  Q15-Q15->Q15      */
inline_ int16_t S_sub_ss(int16_t x, int16_t y)
{
    return S_sature_l(((int32_t)x) - y);
}

/* negation -Q15->Q15                            */
inline_ int32_t L_neg_l(int32_t x)
{
    return (x == MIN_INT32) ? MAX_INT32 : -x;
}

/* negation -Q15->Q15                            */
inline_ int16_t S_neg_s(int16_t x)
{
    return (x == (int16_t)MIN_INT16) ? MAX_INT16 : -x;
}

/* absolute value with saturation  Q15->Q15   */
inline_ int16_t S_abs_s(int16_t x)
{
    if (x == (int16_t)MIN_INT16) return MAX_INT16;
    return (x < 0) ? -x : x;
}

/* absolute value with saturation  Q31->Q31  */
inline_ int32_t L_abs_l(int32_t x)
{
    if (x == MIN_INT32) return MAX_INT32;
    return (x < 0) ? -x : x;
}



/* rounding (adding a 0x8000) with saturation:
Q31->Q31 - add 0x8000 with saturation                             */
inline_ int32_t L_round_l(int32_t x)
{
    return (L_add_ll(x, (int32_t)0x8000U)) & 0xFFFF0000;
}

/* rounding (adding a 0x8000) with saturation:
Q31->Q15 - add 0x8000 with saturation and shift left by 16       */
inline_ int16_t S_round_l(int32_t x)
{
    return (int16_t)(L_add_ll(x, (int32_t)0x8000U) >> 16);
}



/* getting exponent: exp(Q15), return  0 if x==0                 */
inline_ int16_t S_exp_s(int16_t x)
{
    return S_exp_l(L_deposit_s(x));
}

/* getting exponent: exp(Q15), return 15 if x==0                 */
inline_ int16_t S_exp0_s(int16_t x)
{
    return S_exp0_l((int32_t)x) - 16;
}


/* shift left with saturation: Qx->Q(x+t-16)      */
inline_ int32_t L_shl_s(int16_t x, int16_t t)
{
    return L_shl_l(L_deposit_s(x), t);
}

/* shift left with saturation: Qx->Q(x+t)         */
inline_ int16_t S_shl_s(int16_t x, int16_t t)
{
    return (int16_t)(L_shl_l(L_deposit_s(x), t) >> 16);
}

/* shift right with saturation: Qx->Q(x-t)        */
inline_ int32_t L_shr_l(int32_t x, int16_t t)
{
    return L_shl_l(x, (int16_t)-t);
}

/* shift right with saturation: Qx->Q(x-t-16)     */
inline_ int32_t L_shr_s(int16_t x, int16_t t)
{
    return L_shl_s(x, (int16_t)-t);
}

/* shift right with saturation: Qx->Q(x-t)        */
inline_ int16_t S_shr_s(int16_t x, int16_t t)
{
    return S_shl_s(x, (int16_t)-t);
}

/* fractional multiplication: Q15xQ15->Q31         */
inline_ int32_t	L_mpy_ss(int16_t x, int16_t y)
{
    int32_t z;
    z = (int32_t)x * (int32_t)y;
    if ((x != (int16_t)MIN_INT16) || (y != (int16_t)MIN_INT16))
        return (z << 1);
    return MAX_INT32;
}

/* fractional multiplication: Q15xQ15->Q15         */
inline_ int16_t	S_mpy_ss(int16_t x, int16_t y)
{
    return (int16_t)(L_mpy_ss(x, y) >> 16);
}

/* fractional multiplication: Q31xQ15->Q31         */
inline_ int32_t	L_mpy_ls(int32_t x, int16_t y)
{
    int32_t z;
    int16_t t;

    t = ((int16_t)(x >> 1)) & 0x7FFF;
    z = L_mpy_ss(y, t) >> 15;
    z = L_add_ll(z, L_mpy_ss(y, S_extract_l(x)));
    return z;
}


/* multiplication : no saturation and
 overflow control                     */
inline_ int32_t	L_mul_ss(int16_t x, int16_t y) { return ((int32_t)x)*y; }
inline_ int16_t	S_mul_ss(int16_t x, int16_t y)	{ return x*y; }
inline_ int32_t	L_mul_ls(int32_t x, int16_t y) { return x*y; }
inline_ int32_t	L_mul_ll(int32_t x, int32_t y) { return x*y; }

/* MAC operation z=z+x*y: Q31+Q15xQ15->Q31            */
inline_ int32_t	L_mac_ss(int32_t z, int16_t x, int16_t y)
{
    return L_add_ll(z, L_mpy_ss(x, y));
}
/* MAC operation z=z+x*y: Q31+Q15xQ15->Q31            */
inline_ int32_t	L_mac_ls(int32_t z, int32_t x, int16_t y)
{
    return L_add_ll(z, L_mpy_ls(x, y));
}


/* MAS operation z=z-x*y: Q31-Q15xQ15->Q31            */
inline_ int32_t	L_mas_ss(int32_t z, int16_t x, int16_t y)
{
    return L_sub_ll(z, L_mpy_ss(x, y));
}
/* MAS operation z=z-x*y: Q31-Q15xQ15->Q31            */
inline_ int32_t	L_mas_ls(int32_t z, int32_t x, int16_t y)
{
    return L_sub_ll(z, L_mpy_ls(x, y));
}

/* fractional multiplication: Q31xQ31->Q31         */
inline_ int32_t	L_mpy_ll(int32_t x, int32_t y)
{
    int16_t xl, yl;
    int16_t xh, yh;
    int32_t  Acc;

    xl = (((int16_t)(x)) >> 1) & 0x7FFF;
    yl = (((int16_t)(y)) >> 1) & 0x7FFF;
    xh = S_extract_l(x);
    yh = S_extract_l(y);

    Acc = (L_mpy_ss(xh, yl) >> 1) +
        (L_mpy_ss(yh, xl) >> 1) +
        (L_mpy_ss(xl, yl) >> 16);
    Acc = Acc >> 14;
    Acc += L_mpy_ss(xh, yh);
    return Acc;
}

/* addition without saturation */
inline_ int128_t int128_add(int128_t a, int128_t b)
{
    int128_t r = { a.lo + b.lo, a.hi + b.hi };
    r.hi += (r.lo < a.lo) || (r.lo < b.lo);
    return r;
}

/* subtraction without saturation */
inline_ int128_t int128_sub(int128_t a, int128_t b)
{
    int128_t r = { a.lo - b.lo, a.hi - b.hi };
    if (r.lo > a.lo) r.hi--;
    return r;
}

/* Multiply two 32-bit values  */
inline_ int128_t int128_mul_ll(int32_t x, int32_t y)
{
    int64_t t;
    int128_t r;
    t = (int64_t)x*y;
    r.lo = (uint64_t)t;
    r.hi = (t >> 63);
    return r;
}

/* Multiply two 32-bit values and add the result to 128-bit accumulator, without saturation */
inline_ int128_t int128_mac_ll(int128_t z, int32_t x, int32_t y)
{
    int64_t t;
    int128_t p;
    t = (int64_t)x*y;
    p.lo = (uint64_t)t;
    p.hi = (t >> 63);
    return int128_add(z, p);
}

/* Multiply two 32-bit values and subtract the result from 128-bit accumulator, without saturation */
inline_ int128_t int128_mas_ll(int128_t z, int32_t x, int32_t y)
{
    int64_t t;
    int128_t p;
    t = (int64_t)x*y;
    p.lo = (uint64_t)t;
    p.hi = (t >> 63);
    return int128_sub(z, p);
}

/*  shift right  */
inline_ int128_t int128_rsh(int128_t x, int sh)
{
    int128_t y;
    int64_t tmp;
    ASSERT(sh >= 0);
    ASSERT(sh < 64);
    y.hi = 0; y.lo = ((int64_t)1 << sh) >> 1;
    x = int128_add(x, y);
    tmp = x.hi&(((int64_t)1 << sh) - 1);
    y.hi = x.hi >> sh;
    x.lo = x.lo >> sh;
    y.lo = x.lo | (tmp << (64 - sh));
    return y;
}


inline_ int128_t int128_set64(int64_t x)
{
    int128_t a;
    a.hi = (int64_t)(x >> 63);
    a.lo = (uint64_t)x;
    return a;
}

inline_ int128_t int128_lsh(int128_t x, int sh)
{
    int128_t y;
    y = x;
    ASSERT(sh >= 0);
    ASSERT(sh < 64);
    y.hi = x.hi << sh;
    y.lo = y.lo << sh;
    y.hi |= (sh == 0) ? 0 : (x.lo >> (64 - sh));
    return y;
}

/* returns exponent of 80-bit number. For 0 returns 79 */
inline_ int16_t S_exp0_80(int128_t x)
{
    int32_t t;
    int16_t z = 0;
    if (x.lo == 0 && x.hi == 0)  return 79;
    while (1)
    {
        t = (int32_t)x.hi;
        if (((t ^ (t << 1)) & 0x8000) != 0) break;
        x.hi = (x.hi << 1) | (x.lo >> 63);
        x.lo <<= 1;
        z++;
    }
    return z;
}

/* returns exponent of 40-bit number. For 0 returns 39 */
inline_ int16_t S_exp0_40(int64_t x)
{
    int32_t t;
    int16_t z = 0;
    if (x == 0)  return 39;
    while (1)
    {
        t = (int32_t)(x >> 32);
        if (((t ^ (t << 1)) & 0x80) != 0) break;
        x <<= 1;
        z++;
    }
    return z;
}

/* returns exponent of 64-bit number. For 0 returns 63 */
inline_ int16_t S_exp0_64(int64_t x)
{
    int32_t t;
    int16_t z = 0;
    if (x == 0)  return 63;
    while (1)
    {
        t = (int32_t)(x >> 32);
        if (((t ^ (t << 1)) & 0x80000000) != 0) break;
        x <<= 1;
        z++;
    }
    return z;
}

/* multiply int32 number by lower 64 bit number from y */
inline_ int128_t int128_mul_32x64(int32_t x, int64_t y)
{
    int128_t a, b;
    int64_t t;
    t = x*(y >> 32);                  a = int128_set64(t);
    t = (int64_t)x*((uint32_t)y);   b = int128_set64(t);
    a = int128_lsh(a, 32);
    a = int128_add(a, b);
    return a;
}

/*  2-way shift left
    shifts greater than 80 will output zero
    shifts less than -80 will output the sign
    */
inline_ int128_t int128_sla(int128_t x, int sh)
{
    int128_t y;
    if (sh == 0) return x;
    if (sh >= 80)
    {
        y.lo = 0; y.hi = 0; return y;
    }
    if (sh <= -80)
    {
        y.hi = x.hi >> 63;
        y.lo = y.hi; return y;
    }

    if (sh < 0)
    {    /* shift right for negative inputs */
        sh = -sh;
        if (sh >= 64)
        {
            sh -= 64;
            y.lo = x.hi >> sh;
            y.hi = x.hi >> 63;
        }
        else
        {
            y.lo = x.lo >> sh;
            y.hi = x.hi >> sh;
            y.lo |= x.hi << (64 - sh);
        }
    }
    else
    {
        /* shift left for positive inputs */
        if (sh >= 64)
        {
            sh -= 64;
            y.hi = x.lo << sh;
            y.lo = 0;
        }
        else
        {
            y.hi = x.hi << sh;
            y.lo = x.lo << sh;
            y.hi |= x.lo >> (64 - sh);
        }
    }
    return y;
}

/* take sign of floating-point number */
inline_ int takesignf(float32_t x)
{
  union
  {
    float32_t f;
    uint32_t  u;
  }y;
  y.f = x;
  return (int)(y.u >> 31);
}
/* set sign of floating-point number */
inline_ float32_t setsignf(float32_t x, int sign)
{
  union
  {
    float32_t f;
    uint32_t  u;
  }y;
  y.f = x;
  y.u &= 0x7fffffff;
  sign = sign ? 1 : 0;
  y.u |= ((uint32_t)sign) << 31;
  return y.f;
}
/* change sign of floating-point number */
inline_ float32_t changesignf(float32_t x, int sign)
{
  union
  {
    float32_t f;
    uint32_t  u;
  }y;
  y.f = x;
  sign = sign ? 1 : 0;
  y.u ^= ((uint32_t)sign) << 31;
  return y.f;
}

/* return non-zero if argument is a signalling NaN */
inline_ int is_snanf(float32_t x)
{
  union
  {
    float32_t f;
    uint32_t  u;
  }y;
  y.f = x;
  return ((y.u & 0x7f800000) == 0x7f800000 &&
    (y.u & 0x00400000) == 0 &&
    (y.u & 0x003fffff) != 0);
}
#endif /* __BASEOP_H__ */
