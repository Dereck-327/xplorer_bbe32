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
    NatureDSP_Baseband library. Direct Matrix Inversion
    Direct inversion of 3x3 floating point matrices 
    C code optimized for BBE32EP with VFPU
    IntegrIT, 2006-2017
*/
#ifndef MATINV3X3TBL_H__
#define MATINV3X3TBL_H__

#define MATINV3X3_PERM_TBL_SIZE 8

extern const int16_t ALIGN(32) matinv3x3_fwd_perm_tbl[9*MATINV3X3_PERM_TBL_SIZE];
extern const int16_t ALIGN(32) matinv3x3_bkw_perm_tbl[9*MATINV3X3_PERM_TBL_SIZE];

// for ported variants
extern const int16_t ALIGN(32) matinv3x3_bkw_perm_tbl_bbe32[MATINV3X3_PERM_TBL_SIZE*32];
extern const int16_t ALIGN(32) matinv3x3_fwd_perm_tbl_bbe32[MATINV3X3_PERM_TBL_SIZE*32];
extern const int16_t ALIGN(32) matinv3x3sf_searchTbl[16*5];
#endif
