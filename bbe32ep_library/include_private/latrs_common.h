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
#ifndef LATRS_COMMON__
#define LATRS_COMMON__

/* Portable data types. */
#include "NatureDSP_types.h"

/* Lattice filter processing function. */
typedef void (proc_fxn_t)( int16_t * restrict r,
                           int16_t * restrict d,
                     const int16_t *          x,
                     const int16_t *          coef,
                           int16_t            gain,
                           int N, int L );
typedef proc_fxn_t * proc_fxn_ptr_t;


/* Lattice real block IIR processing function, Fast fixed-point implementation. */
proc_fxn_t latrs_sp_proc1;
proc_fxn_t latrs_sp_proc2;
proc_fxn_t latrs_sp_proc3;
proc_fxn_t latrs_sp_proc4;
proc_fxn_t latrs_sp_proc5;
proc_fxn_t latrs_sp_proc6;
proc_fxn_t latrs_sp_proc7;
proc_fxn_t latrs_sp_proc8;

/* Lattice real block IIR processing function, Low Noise fixed-point implementation. */
proc_fxn_t latrs_dp_proc1;
proc_fxn_t latrs_dp_proc2;
proc_fxn_t latrs_dp_proc3;
proc_fxn_t latrs_dp_proc4;
proc_fxn_t latrs_dp_proc5;
proc_fxn_t latrs_dp_proc6;
proc_fxn_t latrs_dp_proc7;
proc_fxn_t latrs_dp_proc8;

#endif // LATRS_COMMON__
