/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs (“Cadence    */
/* Libraries”) are the copyrighted works of Cadence Design Systems Inc.     */
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
/*          Copyright (C) 2009-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */

#ifndef __NATUREDSP_BASEBAND_ID_H
#define __NATUREDSP_BASEBAND_ID_H
   
/* library identification */
#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
  Identification Routines
  NatureDSP_Baseband_get_library_version     Library Version Request
  NatureDSP_Baseband_get_library_api_version Library API Version Request
===========================================================================*/

/*-------------------------------------------------------------------------
Returns library version string
Input:
    none
Output:
     version_string   Pre-allocated buffer for version string.
Restrictions:
     version_string must points to a buffer large enough to hold up to
     30 characters.
-------------------------------------------------------------------------*/

void NatureDSP_Baseband_get_library_version(char *version_string);

/*-------------------------------------------------------------------------
Returns library API version string
Input:
    none
Output:
     version_string   Pre-allocated buffer for API version string.
Restrictions:
     version_string must points to a buffer large enough to hold up to
     30 characters.
-------------------------------------------------------------------------*/

void NatureDSP_Baseband_get_library_api_version(char *version_string);

/*-------------------------------------------------------------------------
Returns non-zero if given function (by its address) is supported by
specific processor capabilities
Input:
     fun    one of function from the list above
Output:
none
NOTE:
in gcc/xcc environment, calls of this function is not neccessary - if 
function pointer is non-zero it means it is supported. VisualStudio linker 
does not support section removal so this function might be used for 
running library under MSVC environment
-------------------------------------------------------------------------*/

typedef void(*NatureDSP_Baseband_funptr)();
int __NatureDSP_Baseband_isPresent(NatureDSP_Baseband_funptr fun);
#define NatureDSP_Baseband_isPresent(fun) __NatureDSP_Baseband_isPresent((NatureDSP_Baseband_funptr)fun)

#ifdef __cplusplus
}
#endif

#endif /* __NATUREDSP_BASEBAND_ID_H */
