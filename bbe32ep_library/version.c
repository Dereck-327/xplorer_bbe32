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
#include "NatureDSP_Baseband_id.h"
#include "common.h"
#include <string.h>

static const char NatureDSP_Library_version[]     = "4.0.0";
static const char NatureDSP_Library_API_version[] = "1.54";

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
void NatureDSP_Baseband_get_library_version(char *version_string)
{
  NASSERT(sizeof(NatureDSP_Library_version)<=30);
  strncpy (version_string, NatureDSP_Library_version, 30);
}

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
void NatureDSP_Baseband_get_library_api_version(char *version_string)
{
  NASSERT(sizeof(NatureDSP_Library_API_version)<=30);
  strncpy (version_string, NatureDSP_Library_API_version, 30);
}
