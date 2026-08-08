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
    support for conditionalization: detect if some function is present or 
    not
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_id.h"
#include "common.h"

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

#ifdef COMPILER_MSVC
/*  ------------------------- MSVC code ------------------------- */

#if (_MSC_VER<1500) // (Visual Studio 2005 and earlier)
unsigned char _BitScanReverse(
   unsigned long * Index,
   unsigned long Mask
)
{
    int idx=0;
    if(Mask==0) return 0;
    while(Mask & 1)
    {
        idx++;
        Mask>>=1;
    }
    return idx;
}
#endif  


/* 
    dummy function that does nothing. It just marks the 
    memory region where discarding library functions 
    begins from
*/
extern  "C" void ____begin_discarding____ (void) ;
#pragma alloc_text( "$$$$$$$$$$", ____begin_discarding____ )
static void ____begin_discarding____ (void) {}

int __NatureDSP_Baseband_isPresent(NatureDSP_Baseband_funptr fun)
{
    return ((uintptr_t)fun)<((uintptr_t)____begin_discarding____) && (fun!=NULL);
}

#endif
/*  ------------------------- GCC code ------------------------- */
#ifdef COMPILER_GNU
int __NatureDSP_Baseband_isPresent(NatureDSP_Baseband_funptr fun)
{
    return fun!=NULL;
}
#endif

#ifdef COMPILER_XTENSA
/*  ------------------------- XCC code ------------------------- */
int __NatureDSP_Baseband_isPresent(NatureDSP_Baseband_funptr fun)
{
    return fun!=NULL;
}
#endif
