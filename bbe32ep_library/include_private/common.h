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

#ifndef __COMMON_H__
#define __COMMON_H__

#if defined COMPILER_XTENSA
  #include <xtensa/config/core-isa.h>
  #include <xtensa/tie/xt_core.h>
  #include <xtensa/tie/xt_misc.h>
  #include <xtensa/tie/xt_bben.h>
  #include <xtensa/tie/xt_density.h>
  #include <xtensa/tie/xt_bben_verification.h>
  #ifndef BBE_SIMD_WIDTH
    #define BBE_SIMD_WIDTH XCHAL_BBEN_SIMD_WIDTH
  #endif
#else

/* the code below causes inclusion of file "cstub-"XTENSA_CORE"-c.h" */
#define PPCAT_NX(A, B) A-B
#define PPCAT(A, B) PPCAT_NX(A, B)      /* Concatenate preprocessor tokens A and B after macro-expanding them. */
#define STRINGIZE_NX(A) #A              /* Turn A into a string literal without expanding macro definitions */
#define STRINGIZE(A) STRINGIZE_NX(A)    /*  Turn A into a string literal after macro-expanding it. */
//#include STRINGIZE(PPCAT(PPCAT(cstub,XTENSA_CORE),c.h))
#include "xtensa/tie/xt_bben.h"
#include "xtensa/config/core-isa.h"

#define BBE_SIMD_WIDTH XCHAL_BBEN_SIMD_WIDTH

#endif

#ifdef _MSC_VER
    #define ALIGN(x)    _declspec(align(x)) 
#else
    #define ALIGN(x)    __attribute__((aligned(x))) 
#endif

#ifdef COMPILER_XTENSA
  #define ATTRIBUTE_ALWAYS_INLINE __attribute__((always_inline))
  #define ATTRIBUTE_NEVER_INLINE  __attribute__((noinline))
#else
  #define ATTRIBUTE_ALWAYS_INLINE
  #define ATTRIBUTE_NEVER_INLINE
#endif

//-----------------------------------------------------
// log2(BBE_SIMD_WIDTH)
//-----------------------------------------------------
#if BBE_SIMD_WIDTH==8
#define LOG2_BBE_SIMD_WIDTH 3
#define ALIGN_SIMD ALIGN(16)
#elif BBE_SIMD_WIDTH==16
#define LOG2_BBE_SIMD_WIDTH 4
#define ALIGN_SIMD ALIGN(32)
#elif BBE_SIMD_WIDTH==32
#define LOG2_BBE_SIMD_WIDTH 5
#define ALIGN_SIMD ALIGN(64)
#else
#error unsupported BBE_SIMD_WIDTH 
#endif
//-----------------------------------------------------
// some C++ support
//-----------------------------------------------------
#ifdef __cplusplus

typedef struct { int16_t x[BBE_SIMD_WIDTH]; } ALIGN(32) _xb_vecNx16;
typedef struct { int16_t s[BBE_SIMD_WIDTH]; } ALIGN(32) _vselN;
typedef struct { int16_t v[BBE_SIMD_WIDTH]; } ALIGN(32) _vsaN;
#if defined (COMPILER_MSVC)
#define _X(ptr_xb_vecNx16) (*(const xb_vecNx16*)&(ptr_xb_vecNx16))
#define _V(ptr_vsaN)       (*(const vsaN*) (&(ptr_vsaN)))
#define _S(ptr_vselN)      (*(const vselN*)(&(ptr_vselN)))
#elif defined (COMPILER_GNU)
inline_ const xb_vecNx16* convert_dummy_xb_vecNx16_ptr(const void* ptr) { return (const xb_vecNx16*) ptr; }
inline_ const vsaN*       convert_dummy_vsaN_ptr      (const void* ptr) { return (const vsaN*)       ptr; }
inline_ const vselN*      convert_dummy_vselN_ptr     (const void* ptr) { return (const vselN*)      ptr; }
#define _X(ptr_xb_vecNx16) (*convert_dummy_xb_vecNx16_ptr(&(ptr_xb_vecNx16)))
#define _V(ptr_vsaN)       (*convert_dummy_vsaN_ptr      (&(ptr_vsaN)))
#define _S(ptr_vselN)      (*convert_dummy_vselN_ptr     (&(ptr_vselN)))
#else
#error unsupported compiler
#endif

#else
// in C we may simply avoid stupid conversions!!!
typedef xb_vecNx16 _xb_vecNx16;
typedef vselN      _vselN;
typedef vsaN       _vsaN;
#define _X(ptr_xb_vecNx16) (ptr_xb_vecNx16)
#define _V(ptr_vsaN)       (ptr_vsaN)
#define _S(ptr_vselN)      (ptr_vselN)
#endif

// special XCC type casting of pointers
#ifdef __cplusplus
#define castxcc(type_,ptr)  (ptr)
#else
#define castxcc(type_,ptr)  (type_ *)(ptr)
#endif

//-----------------------------------------------------
// C99 pragma wrapper
//-----------------------------------------------------

#ifdef COMPILER_XTENSA
#define __Pragma(a) _Pragma(a)
#else
#define __Pragma(a)
#endif

//-----------------------------------------------------
// Patches for malfuntioning instruction intrinsics.
//-----------------------------------------------------

#ifndef BBE_LV4X16_I
#define BBE_LV4X16_I(ars, i) BBE_MOVNX16_FROMN_4X64(BBE_LSN_4X64_I((const long long *)ars, i));
#endif
#ifndef BBE_SV4X16_I
#define BBE_SV4X16_I(vr, ars, i) BBE_SSN_4X64_I(BBE_MOVN_4X64_FROMNX16(vr),(long long *)ars,i);
#endif

/*without advanced precision*/
#if XCHAL_HAVE_BBEN_ADVPRECISION==0
  

  #ifndef BBE_SUBSAVSN
  #define BBE_SUBSAVSN(ar, shft) __BBE_SUBSAVSN(ar, &shft)
  inline_ vsaN __BBE_SUBSAVSN(unsigned int a, const void* pb)
  {
    vsaN b, c;
    xb_vecNx16 tmp;
    b = *(const vsaN*)pb;
    tmp = BBE_SUBNX16(BBE_MOVVA16(a), BBE_MOVVVS(b));
    tmp = BBE_MAXNX16(tmp, -64);
    tmp = BBE_MINNX16(tmp, 63);
    c = BBE_MOVVSV(tmp, 0);
    return c;
  }
#endif

#ifndef BBE_PACKNVNX40
#define BBE_PACKNVNX40(w, shft) __BBE_PACKNVNX40(&w, &shft)  
#endif // !BBE_PACKNVNX40
  inline_ xb_vecNx16 __BBE_PACKNVNX40(const void* pw, const void* ps)
  {

    xb_vecNx40 w;
    xb_vecNx16 z;
    vsaN shft;
    w = *(const xb_vecNx40*)pw;
    shft = *(const vsaN*)ps;
    shft = BBE_SUBSAVSN(15, shft);
    w = BBE_RNDADJNX40(w, shft);
    z = BBE_PACKVNX40(w, shft);
    return z;
  }

  #define BBE_UNPKNVNX16(x, shft) __BBE_UNPKNVNX16(&x, &shft)
  inline_ xb_vecNx40 __BBE_UNPKNVNX16(const void* px, const void* ps)
  {

    xb_vecNx40 w;
    xb_vecNx16 x;
    vsaN shft;
    xb_vecNx16 _0x8000 = 0x8000;
    x = *(const xb_vecNx16*)px;
    shft = *(const vsaN*)ps;
    w = BBE_MULUSNX16(_0x8000,x);
    w = BBE_SRANX40(w, shft);
    return w;
  }

  #define BBE_PACKHNX40(w) __BBE_PACKHNX40(&w)
  inline_ xb_vecNx16 __BBE_PACKHNX40(const void* pw)
  {
    xb_vecNx40 w;
    xb_vecNx16 z;
    vsaN shft = BBE_MOVVSA32(1);
    w = *(const xb_vecNx40*)pw;
    w = BBE_SRANX40(w, BBE_MOVVSA32(23));
    w = BBE_RNDADJNX40(w, shft);
    z = BBE_PACKVNX40(w, shft);
    return z;
  }

  #define BBE_ADDSVSN(a, b) __BBE_ADDSVSN(&a, &b)

  inline_ vsaN __BBE_ADDSVSN(const void* pa, const void* pb)
  {

    xb_vecNx16 tmp;
    vsaN a, b, c;
    a = *(const vsaN*)pa;
    b = *(const vsaN*)pb;
    tmp = BBE_ADDNX16(BBE_MOVVVS(a), BBE_MOVVVS(b));
    tmp = BBE_MAXNX16(tmp, -64);
    tmp = BBE_MINNX16(tmp, 63);
    c = BBE_MOVVSV(tmp, 0);
    return c;
  }

  #define BBE_SUBSVSN(a, b) __BBE_SUBSVSN(&a, &b)

  inline_ vsaN __BBE_SUBSVSN(const void* pa, const void* pb)
  {

    xb_vecNx16 tmp;
    vsaN a, b, c;
    a = *(const vsaN*)pa;
    b = *(const vsaN*)pb;
    tmp = BBE_SUBNX16(BBE_MOVVVS(a), BBE_MOVVVS(b));
    tmp = BBE_MAXNX16(tmp, -64);
    tmp = BBE_MINNX16(tmp, 63);
    c = BBE_MOVVSV(tmp, 0);
    return c;
  }

  #define BBE_MINVSN(a, b) __BBE_MINVSN(&a, &b)

  inline_ vsaN __BBE_MINVSN(const void* pa, const void* pb)
  {

    xb_vecNx16 tmp;
    vsaN a, b, c;
    a = *(const vsaN*)pa;
    b = *(const vsaN*)pb;
    tmp = BBE_MINNX16(BBE_MOVVVS(a), BBE_MOVVVS(b));
    c = BBE_MOVVSV(tmp, 0);
    return c;
  }

#endif

//-----------------------------------------------------
// API annotation helper macro
//-----------------------------------------------------

#ifdef COMPILER_XTENSA
#define ANNOTATE_ATTR __attribute__ ((section (".rodata.NatureDSP_Baseband_annotation")))
#else
#define ANNOTATE_ATTR
#endif

#define ANNOTATE_FUN_REF(fun)   NatureDSP_Baseband_annotation_##fun

#ifdef __cplusplus
#define ANNOTATE_FUN(fun,text) \
  extern "C" const char ANNOTATE_ATTR ANNOTATE_FUN_REF(fun)[] = (text)
#else
#define ANNOTATE_FUN(fun,text) \
  const char ANNOTATE_ATTR ANNOTATE_FUN_REF(fun)[] = (text)
#endif

//-----------------------------------------------------
// Conditionalization support
//-----------------------------------------------------
/* place DISCARD_FUN(retval_type,name) instead of function definition for functions
   to be discarded from the executable 
   THIS WORKS only for external library functions declared as extern "C" and
   not supported for internal references without "C" qualifier!
*/
#ifdef COMPILER_MSVC
#pragma section( "$DISCARDED_FUNCTIONS" , execute, discard )
#pragma section( "$$$$$$$$$$" , execute, discard )
#define DISCARD_FUN(retval_type,name,arglist) __pragma (alloc_text( "$DISCARDED_FUNCTIONS",name))\
__pragma(section( "$DISCARDED_FUNCTIONS" , execute, discard ))\
__pragma (warning(push))\
__pragma (warning( disable : 4026 4716))\
retval_type name arglist {}\
__pragma (warning(pop))
#endif

#if defined (COMPILER_GNU)
#define F_UNDERSCORE " "
#define DISCARD_FUN(retval_type,name,arglist)    \
__asm__                        \
(                              \
".section unused_section\n"    \
".globl " F_UNDERSCORE STRINGIZE(name) "\n" \
".type "F_UNDERSCORE STRINGIZE(name)", @function \n"\
F_UNDERSCORE STRINGIZE(name) ":\n"          \
".text"                        \
);
#endif

#if defined(COMPILER_XTENSA)
#define DISCARD_FUN(retval_type,name,arglist) \
__attribute__ ((section ("/DISCARD/"))) \
retval_type name arglist \
{  }
#endif

/*------ LIST OF DEFINES DEPENDING ON ISA OPTIONS ------*/

/*
presence of:
BBE_ADDPNX16RCU
BBE_ADDPNX16RCU
BBE_ADDPNX16RCUMBC
BBE_ADDPNX16RCUMBCIAD
BBE_ADDPNX16RRU
BBE_ADDPNX16RRUMBC
BBE_ADDPNX16RRUMBCIAD
*/
#define HAVE_ADDPN XCHAL_HAVE_BBEN_SYM_FIR_STATE

/*
BBE_ADDSAVSN
BBE_VSA_SHFLI_DUPLICATE_1_EVEN
BBE_VSA_SHFLI_DUPLICATE_1_ODD
BBE_SHFLVSNI
BBE_SUBSAVSN
BBE_SUBSR1SAVSN
*/
#define HAVE_VSAMATH (XCHAL_HAVE_BBEN_ADVPRECISION || XCHAL_HAVE_BBEN_ADVPRECRSQRT || XCHAL_HAVE_BBEN_ADVPRECRRECIP)

/*
BBE_ADDSVSN
BBE_SUBSVSN
*/
#define HAVE_VSAMATHEX (XCHAL_HAVE_BBEN_ADVPRECISION)

/*
BBE_NSAENX40
*/
#if ((XCHAL_HAVE_BBEN_ADVPRECISION) || (XCHAL_HAVE_BBEN_ADVPRECRSQRT))
#define HAVE_NSAENX40 1
#else
#define HAVE_NSAENX40 0
#endif

/*
BBE_INTLVNX16X1H
BBE_INTLVNX16X1L
*/
#if (defined(BBE_INTLVNX16X1H) || defined(BBE_INTLVNX16X1L))
#define HAVE_INTLV 1
#else
#define HAVE_INTLV 0
#endif

/*
BBE_DIPACKxxx
*/
#if (defined(BBE_DIPACKQNX40C) || defined(BBE_DIPACKSNX40C) || defined(BBE_DIPACKLNX40C))
#define HAVE_DIPACK 1
#else
#define HAVE_DIPACK 0
#endif

/*
BBE_MULUUNX16PACKH
*/
#if (defined(BBE_MULUUNX16PACKH) )
#define HAVE_MULUUNX16PACKH 1
#else
#define HAVE_MULUUNX16PACKH 0
#endif

/*
BBE_BMUL32A
BBE_CC64
BBE_MOVABMULACC
BBE_MOVBMULACCA
BBE_MOVBMULSTATEV
*/
#if ((XCHAL_HAVE_BBEN_LFSR) )
#define HAVE_LFSR 1
#else
#define HAVE_LFSR 0
#endif

/*
BBE_DIVNX16S_3STEPN_0
BBE_DIVNX16S_3STEPN_1
BBE_DIVNX16S_4STEP_0
BBE_DIVNX16S_4STEP_1
BBE_DIVNX16U_4STEP_0
BBE_DIVNX16U_4STEPN_0
BBE_DIVNX32S_5STEP0_0
BBE_DIVNX32S_5STEP0_1
BBE_DIVNX32U
BBE_DIVNX32U_4STEP0_0
BBE_QUONX32
BBE_MOVVREM
*/
#if ((XCHAL_HAVE_BBEN_VECDIVIDE) )
#define HAVE_DIV 1
#else
#define HAVE_DIV 0
#endif

/*
BBE_DSPR1DNX16CSF4
BBE_DSPR1DNX16CSF8
*/
#if ((XCHAL_HAVE_BBEN_1D_DESPREAD_OPS) )
#define HAVE_DSPR 1
#else
#define HAVE_DSPR 0
#endif

/*
BBE_FFTADD4SABN_2XCQ15
BBE_FFTADD4SABNX16
BBE_FFTADD4SCDN_2XCQ15
BBE_FFTADD4SCDNX16
BBE_FFTADDSSRN_2XCQ15
BBE_FFTADDSSRNX16
BBE_FFTAVGN_2XCQ15SB
BBE_FFTAVGNX16SB
BBE_FFTSRAN_2XCQ15
BBE_FFTSRANX16
BBE_FFTSUBSSRN_2XCQ15
BBE_FFTSUBSSRNX16
BBE_FFTWMODE
BBE_SALIGNVRNX16_XP
BBE_SAVRN_2XCQ15_XP
BBE_SVINTLARNX16_XP
BBE_SVRN_2XCQ15_IP
BBE_SVRN_2XCQ15_X
BBE_SVRN_2XCQ15_XP
BBE_SVRNX16_IP
BBE_SVRNX16_XP
BBE_MOVUVR
BBE_RANGENX16
BBE_WRANGE
*/
#if ((XCHAL_HAVE_BBEN_VECTORFFT) )
#define HAVE_FFT 1
#else
#define HAVE_FFT 0
#endif

/*
BBE_FPRECIPNX16_0
BBE_FPRECIPNX16_1
*/
#if ((XCHAL_HAVE_BBEN_FASTRECIP) )
#define HAVE_FPRECIP 1
#else
#define HAVE_FPRECIP 0
#endif

/*
BBE_FPRSQRTNX16_0
BBE_FPRSQRTNX16_1
*/
#if ((XCHAL_HAVE_BBEN_FASTRSQRT ) )
#define HAVE_FPRSQRT 1
#else
#define HAVE_FPRSQRT 0
#endif

/*
BBE_LVA_IP
BBE_LVA_XP
BBE_LVB_IP
BBE_LVB_XP
BBE_LVC_IP
BBE_LVC_XP
BBE_LVD_IP
BBE_LVD_XP
BBE_MOVSAV
BBE_MOVSBV
BBE_MOVSCV
BBE_MOVSDV
*/
#if ((XCHAL_HAVE_BBEN_VECTORFFT) || (XCHAL_HAVE_BBEN_SYM_FIR_STATE  ) )
#define HAVE_ABCD 1
#else
#define HAVE_ABCD 0
#endif

/*
BBE_MULANX16PC_0
BBE_MULANX16PC_1
BBE_MULNX16PC_0
BBE_MULNX16PC_1
BBE_MULRNX16PC_0
*/
#if (defined(BBE_MULANX16PC_0) || defined(BBE_MULANX16PC_1  ) || defined(BBE_MULNX16PC_0) || defined(BBE_MULNX16PC_1) || defined(BBE_MULRNX16PC_0) )
#define HAVE_MULPC 1
#else
#define HAVE_MULPC 0
#endif

/*
BBE_RECIPLUNX40_0
BBE_RECIPLUNX40_1
*/
#if ((XCHAL_HAVE_BBEN_ADVPRECISION)|| (XCHAL_HAVE_BBEN_ADVPRECRRECIP))
#define HAVE_RECIP 1
#else
#define HAVE_RECIP 0
#endif

/*
BBE_RSQRTLUNX40_0
BBE_RSQRTLUNX40_1
*/
#if ((XCHAL_HAVE_BBEN_ADVPRECRSQRT) )
#define HAVE_RSQRT 1
#else
#define HAVE_RSQRT 0
#endif

/*
BBE_SELI_MMC3X3X3X1_OFFSET_M1_0_STEP_1
BBE_SELI_MMC3X3X3X1_OFFSET_M1_0_STEP_2
BBE_SELI_MMC3X3X3X1_OFFSET_M1_0_STEP_3
BBE_SELI_MMC3X3X3X3_OFFSET_M2_0_STEP_3
BBE_SHFLI_MMC1X4X4X4_M1_STEP_1_HIGH_HALF
BBE_SHFLI_MMC1X4X4X4_M1_STEP_1_LOW_HALF
BBE_SHFLI_MMC1X4X4X4_M1_STEP_2_HIGH_HALF
BBE_SHFLI_MMC1X4X4X4_M1_STEP_2_LOW_HALF
BBE_SHFLI_MMC2X2X2X2_M1_STEP_1
BBE_SHFLI_MMC2X2X2X2_M1_STEP_2
BBE_SHFLI_MMC2X2X2X2_M2_STEP_1
BBE_SHFLI_MMC2X2X2X2_M2_STEP_2
BBE_SHFLI_MMC3X3X3X3_OFFSET_M1_0_STEP_1
BBE_SHFLI_MMC3X3X3X3_OFFSET_M1_0_STEP_2
BBE_SHFLI_MMC3X3X3X3_OFFSET_M2_0_STEP_1
BBE_SHFLI_MMC3X3X3X3_OFFSET_M2_0_STEP_2
BBE_SHFLI_MMC4X4X4X1_M1_STEP_1_HIGH_HALF
BBE_SHFLI_MMC4X4X4X1_M1_STEP_1_LOW_HALF
BBE_SHFLI_MMC4X4X4X1_M1_STEP_2_HIGH_HALF
BBE_SHFLI_MMC4X4X4X1_M1_STEP_2_LOW_HALF
BBE_SHFLI_MMC4X4X4X4_M1_STEP_1
BBE_SHFLI_MMC4X4X4X4_M1_STEP_2
BBE_SHFLI_MMC4X4X4X4_M1_STEP_3
BBE_SHFLI_MMC4X4X4X4_M1_STEP_4
BBE_SHFLI_MMC4X4X4X4_M2_STEP_1
BBE_SHFLI_MMC4X4X4X4_M2_STEP_2
BBE_SHFLI_REP_2
BBE_SHFLI_REP_2X4_OFFSET_0
BBE_SHFLI_REP_2X4_OFFSET_1
BBE_SHFLI_REP_2X4_OFFSET_2
BBE_SHFLI_REP_2X4_OFFSET_3
BBE_SHFLNX16I_F5S2
*/
#if (defined(BBE_SELI_MMC3X3X3X1_OFFSET_M1_0_STEP_1) || \
        defined(BBE_SELI_MMC3X3X3X1_OFFSET_M1_0_STEP_2) || \
        defined(BBE_SELI_MMC3X3X3X1_OFFSET_M1_0_STEP_3) || \
        defined(BBE_SELI_MMC3X3X3X3_OFFSET_M2_0_STEP_3) || \
        defined(BBE_SHFLI_MMC1X4X4X4_M1_STEP_1_HIGH_HALF) || \
        defined(BBE_SHFLI_MMC1X4X4X4_M1_STEP_1_LOW_HALF) || \
        defined(BBE_SHFLI_MMC1X4X4X4_M1_STEP_2_HIGH_HALF) || \
        defined(BBE_SHFLI_MMC1X4X4X4_M1_STEP_2_LOW_HALF) || \
        defined(BBE_SHFLI_MMC2X2X2X2_M1_STEP_1) || \
        defined(BBE_SHFLI_MMC2X2X2X2_M1_STEP_2) || \
        defined(BBE_SHFLI_MMC2X2X2X2_M2_STEP_1) || \
        defined(BBE_SHFLI_MMC2X2X2X2_M2_STEP_2) || \
        defined(BBE_SHFLI_MMC3X3X3X3_OFFSET_M1_0_STEP_1) || \
        defined(BBE_SHFLI_MMC3X3X3X3_OFFSET_M1_0_STEP_2) || \
        defined(BBE_SHFLI_MMC3X3X3X3_OFFSET_M2_0_STEP_1) || \
        defined(BBE_SHFLI_MMC3X3X3X3_OFFSET_M2_0_STEP_2) || \
        defined(BBE_SHFLI_MMC4X4X4X1_M1_STEP_1_HIGH_HALF) ||\
        defined(BBE_SHFLI_MMC4X4X4X1_M1_STEP_1_LOW_HALF) || \
        defined(BBE_SHFLI_MMC4X4X4X1_M1_STEP_2_HIGH_HALF) || \
        defined(BBE_SHFLI_MMC4X4X4X1_M1_STEP_2_LOW_HALF) || \
        defined(BBE_SHFLI_MMC4X4X4X4_M1_STEP_1) || \
        defined(BBE_SHFLI_MMC4X4X4X4_M1_STEP_2) || \
        defined(BBE_SHFLI_MMC4X4X4X4_M1_STEP_3) || \
        defined(BBE_SHFLI_MMC4X4X4X4_M1_STEP_4) || \
        defined(BBE_SHFLI_MMC4X4X4X4_M2_STEP_1) || \
        defined(BBE_SHFLI_MMC4X4X4X4_M2_STEP_2) || \
        defined(BBE_SHFLI_REP_2) || \
        defined(BBE_SHFLI_REP_2X4_OFFSET_0) || \
        defined(BBE_SHFLI_REP_2X4_OFFSET_1) || \
        defined(BBE_SHFLI_REP_2X4_OFFSET_2) || \
        defined(BBE_SHFLI_REP_2X4_OFFSET_3) || \
        defined(BBE_SHFLNX16I_F5S2))
#define HAVE_PACKEDMUL 1
#else
#define HAVE_PACKEDMUL 0
#endif

/*
  BBE_SDMAP16QAMNX16C
  BBE_SDMAP256QAMNX16C
  BBE_SDMAP64QAMNX16C
  BBE_SDMAPQPSKNX16C
*/
#if ( (XCHAL_HAVE_BBEN_EP_3GPP_SOFT_BIT_DEMAP) )
#define HAVE_SDMAP 1
#else
#define HAVE_SDMAP 0
#endif

/*
BBE_MOVSVWXH
BBE_MOVSVWXL
BBE_MULMNX16
BBE_SRAIWADDMNX40
*/
#if ((XCHAL_HAVE_BBEN_ADVPRECISION) )
#define HAVE_ADVPMUL 1
#else
#define HAVE_ADVPMUL 0
#endif

/* all single precision floating point instructions */
#if ( (XCHAL_HAVE_BBENEP_SP_VFPU) )
#define HAVE_VFPU 1
#else
#define HAVE_VFPU 0
#endif
#ifdef __cplusplus
#define externC extern "C" 
#else
#define externC extern 
#endif
#if defined(COMPILER_MSVC)

/* make complex_float from xb_vecN_4xcf32 */
#define MOV_COMPLEXFLOAT_FROM_N_4XCF32(a,b) { complex_float ALIGN(8) tmp; BBE_SSN_4XCF32_I (b,(xtcomplexfloat*)&tmp,0); a=tmp; }
#else

#if defined (COMPILER_XTENSA)
/* make complex_float from vecM2xcf32 */
#define MOV_COMPLEXFLOAT_FROM_N_4XCF32(a,b) {a=b;}

#else
/* NOTE: gcc is not able to convert xtcomplexfloat<->complex float
so we have to use load/store instructions here
*/
#define MOV_COMPLEXFLOAT_FROM_N_4XCF32(a,b) { complex_float ALIGN(8) tmp; BBE_SSN_4XCF32_I (b,(xtcomplexfloat*)&tmp,0); a=tmp; }

#endif
#endif

/* patch for RG2016.4 - missing protos */
#if HAVE_VFPU

/* Move a complex_Float value. */
inline_ void MOV_COMPLEX_FLOAT( complex_float * dst, const complex_float * src ) { *dst = *src; }

/* Convert a complex_float value to xtcomplexfloat type. */
inline_ xtcomplexfloat MOV_CF32_FROM_COMPLEX_FLOAT( complex_float a ) { 
  xtcomplexfloat b;  MOV_COMPLEX_FLOAT((complex_float*)&b, &a); return (b); }

/* Convert an xtcomplexfloat value to complex_float type. */
#define MOV_COMPLEX_FLOAT_FROM_CF32(a)    *__MOV_COMPLEX_FLOAT_FROM_CF32(&(a))
inline_ complex_float * __MOV_COMPLEX_FLOAT_FROM_CF32( const xtcomplexfloat * p ) { return (complex_float*)p; }

#ifndef BBE_MOVN_4XCF32T
#define BBE_MOVN_4XCF32T(x,y,b) (BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_4XCF32(x),BBE_MOVNX16_FROMN_4XCF32(y),BBE_MOVN_FROMN_4(b))))
#endif //BBE_MOVN_4XCF32T

#ifndef BBE_MOVN_2XF32T
#define BBE_MOVN_2XF32T(x,y,b) (BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16T(BBE_MOVNX16_FROMN_2XF32(x),BBE_MOVNX16_FROMN_2XF32(y),BBE_MOVN_FROMN_2(b))))
#endif //BBE_MOVN_2XF32T

#ifndef BBE_MOVN_2XF32_FROMN_4XCF32
#define BBE_MOVN_2XF32_FROMN_4XCF32(x) BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_4XCF32(x))
#endif //BBE_MOVN_2XF32_FROMN_4XCF32

#ifndef BBE_MOVN_4XCF32_FROMN_2XF32
#define BBE_MOVN_4XCF32_FROMN_2XF32(x) BBE_MOVN_4XCF32_FROMNX16(BBE_MOVNX16_FROMN_2XF32(x))
#endif //BBE_MOVN_4XCF32_FROMN_2XF32

#ifndef BBE_SELN_2XF32
#define BBE_SELN_2XF32(b,c,d) BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMN_2XC16(BBE_SELN_2XC16(BBE_MOVN_2XC16_FROMNX16(BBE_MOVNX16_FROMN_2XF32(b)),BBE_MOVN_2XC16_FROMNX16(BBE_MOVNX16_FROMN_2XF32(c)),d)))
#endif //BBE_SELN_2XF32

#ifndef BBE_LTRN_4I
#define BBE_LTRN_4I(a) __BBE_LTRN_4I(BBE_MOVN_FROMN_2(BBE_LTRN_2I(a)))
inline_ vboolN_4 __BBE_LTRN_4I(vboolN b)
{
    vboolN c0,c1;
    BBE_EXTRACTB(c1,c0,b);
    return BBE_MOVN_4_FROMN(c0);
}
#endif // BBE_LTRN_4I

#ifndef BBE_LTRN_4
#define BBE_LTRN_4(a) __BBE_LTRN_4(BBE_MOVN_FROMN_2(BBE_LTRN_2(a)))
inline_ vboolN_4 __BBE_LTRN_4(vboolN b)
{
    vboolN c0,c1;
    BBE_EXTRACTB(c1,c0,b);
    return BBE_MOVN_4_FROMN(c0);
}
#endif // BBE_LTRN_4

#ifndef BBE_SELN_4XCF32
#define BBE_SELN_4XCF32(b,c,d) __BBE_SELN_4XCF32(&(b),&(c),&(d))
inline_ xb_vecN_4xcf32 __BBE_SELN_4XCF32(const void* pb,const void* pc,const void* pd)
{
    xb_vecN_4xcf32 res;
    xb_vecNx16 a,b,c;
    vselN v;
    b=*(xb_vecNx16*)pb;
    c=*(xb_vecNx16*)pc;
    v=*(vselN*)pd;
    a=BBE_SELNX16(b,c,v);
    res=BBE_MOVN_4XCF32_FROMNX16(a);
    return res;
}
#endif //BBE_SELN_4XCF32

#ifndef BBE_MOVF32T
#define BBE_MOVF32T(a,b,c)    __BBE_MOVF32T(&(a),&(b),&(c))
inline_ xtfloat __BBE_MOVF32T( const void * pa, const void *pb, const void * pc )
{
  xb_vecN_2xf32 a,b,d;
  vboolN_2 c;
  a = BBE_MOVN_2XF32_FROMF32(*(xtfloat*)pa);
  b = BBE_MOVN_2XF32_FROMF32(*(xtfloat*)pb);
  c = BBE_MOVN_2_FROMN(BBE_EXT0IB(BBE_MOVN_FROM1(*(vbool1*)pc),2));
  d = BBE_MOVN_2XF32T(a,b,c);
  return (BBE_MOVF32_FROMN_2XF32(d));
}
#endif  // BBE_MOVF32T

#ifndef BBE_MOVCF32T
#define BBE_MOVCF32T(a,b,c)   __BBE_MOVCF32T(&(a),&(b),&(c))
inline_ xtcomplexfloat __BBE_MOVCF32T( const void * pa, const void * pb, const void * pc )
{
  xb_vecN_4xcf32 a,b,d;
  vboolN_4 c;
  a = BBE_MOVN_4XCF32_FROMCF32(*(xtcomplexfloat*)pa);
  b = BBE_MOVN_4XCF32_FROMCF32(*(xtcomplexfloat*)pb);
  c = BBE_MOVN_4_FROMN(BBE_EXT0IB(BBE_MOVN_FROM1(*(vbool1*)pc),4));
  d = BBE_MOVN_4XCF32T(a,b,c);
  return (BBE_MOVCF32_FROMN_4XCF32(d));
}
#endif // BBE_MOVCF32T

#ifndef BBE_MOVN_2XF32_FROMCF32
#define BBE_MOVN_2XF32_FROMCF32(x)    BBE_MOVN_2XF32_FROMNX16(BBE_MOVNX16_FROMCF32(x))
#endif // BBE_MOVN_2XF32_FROMCF32

#ifndef BBE_MOVCF32_FROMN_2XF32
#define BBE_MOVCF32_FROMN_2XF32(x)    BBE_MOVCF32_FROMN_2X32(BBE_MOVN_2X32_FROMN_2XF32(x))
#endif // BBE_MOVCF32_FROMN_2XF32

#ifndef BBE_MOVF32_FROMNX16
#define BBE_MOVF32_FROMNX16(x)    BBE_MOVF32_FROMN_2XF32(BBE_MOVN_2XF32_FROMNX16(x))
#endif // BBE_MOVF32_FROMNX16

#ifndef BBE_MOVVSELN_2_FROMVSELN
#define BBE_MOVVSELN_2_FROMVSELN(x)   vsaN_2C_rtor_vselN_2(vsaN_rtor_vsaN_2C(vselN_rtor_vsaN(x)))
#endif

#ifndef BBE_MOVVSELN_FROMVSELN_2
#define BBE_MOVVSELN_FROMVSELN_2(x)   vsaN_rtor_vselN(vsaN_2C_rtor_vsaN(vselN_2_rtor_vsaN_2C(x)))
#endif

#ifndef BBE_MOVN_4_FROMN_2
#define BBE_MOVN_4_FROMN_2(b)   BBE_MOVN_4_FROMN(BBE_MOVN_FROMN_2(b))
#endif

#ifndef BBE_SHFLN_2XF32
#define BBE_SHFLN_2XF32(a,b)    BBE_MOVN_2XF32_FROMNX16(BBE_SHFLNX16(BBE_MOVNX16_FROMN_2XF32(a),BBE_MOVVSELN_FROMVSELN_2(b)))
#endif

#ifndef BBE_SAVN_2XF32POS_FP
#define BBE_SAVN_2XF32POS_FP(va,p)    BBE_SAVN_2XC16POS_FP((va),(xb_vecN_2xc16*)(p))
#endif // BBE_SAVN_2XF32POS_FP

#ifndef BBE_NEGCF32
#define BBE_NEGCF32(x)    BBE_MOVCF32_FROMN_2XF32(BBE_NEGN_2XF32(BBE_MOVN_2XF32_FROMCF32(x)))
#endif // BBE_NEGCF32

#endif
#define MAX_ALLOCA_SZ 2048

#endif // __COMMON_H__
