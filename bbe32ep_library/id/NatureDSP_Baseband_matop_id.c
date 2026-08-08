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
 * NatureDSP_Baseband Library API
 * Matrix Operations
 * Annotations
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matop.h"
#include "common.h"

ANNOTATE_FUN(matmul2x2n  ,"Matrix multiply (block order real 16-bit data)");
ANNOTATE_FUN(matmul3x3n  ,"Matrix multiply (block order real 16-bit data)"); 
ANNOTATE_FUN(matmul4x4n  ,"Matrix multiply (block order real 16-bit data)");
ANNOTATE_FUN(matmul8x8n  ,"Matrix multiply (block order real 16-bit data)");
ANNOTATE_FUN(matmul16x16n,"Matrix multiply (block order real 16-bit data)");
ANNOTATE_FUN(matmulnxnn  ,"Matrix multiply (block order real 16-bit data)");
ANNOTATE_FUN(matmulnxmn  ,"Matrix multiply (block order real 16-bit data)");
ANNOTATE_FUN(matmul2x2nf ,"Matrix multiply (block order real floating-point data)");
ANNOTATE_FUN(matmul3x3nf ,"Matrix multiply (block order real floating-point data)");
ANNOTATE_FUN(matmul4x4nf ,"Matrix multiply (block order real floating-point data)");
ANNOTATE_FUN(matmul8x8nf ,"Matrix multiply (block order real floating-point data)");
ANNOTATE_FUN(matmulnxmnf ,"Matrix multiply (block order real floating-point data)");

ANNOTATE_FUN(matvmul2x2n  ,"Matrix-vector multiply (block order real 16-bit data)");
ANNOTATE_FUN(matvmul3x3n  ,"Matrix-vector multiply (block order real 16-bit data)");
ANNOTATE_FUN(matvmul4x4n  ,"Matrix-vector multiply (block order real 16-bit data)");
ANNOTATE_FUN(matvmul8x8n  ,"Matrix-vector multiply (block order real 16-bit data)");
ANNOTATE_FUN(matvmul16x16n,"Matrix-vector multiply (block order real 16-bit data)");
ANNOTATE_FUN(matvmulnxnn  ,"Matrix-vector multiply (block order real 16-bit data)");
ANNOTATE_FUN(matvmulnxmn  ,"Matrix-vector multiply (block order real 16-bit data)");
ANNOTATE_FUN(matvmul2x2nf ,"Matrix-vector multiply (block order real floating-point data)");
ANNOTATE_FUN(matvmul3x3nf ,"Matrix-vector multiply (block order real floating-point data)");
ANNOTATE_FUN(matvmul4x4nf ,"Matrix-vector multiply (block order real floating-point data)");
ANNOTATE_FUN(matvmul8x8nf ,"Matrix-vector multiply (block order real floating-point data)");
ANNOTATE_FUN(matvmulnxmnf ,"Matrix-vector multiply (block order real floating-point data)");

ANNOTATE_FUN(matmul2x2s   ,"Matrix multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matmul3x3s   ,"Matrix multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matmul4x4s   ,"Matrix multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matmul8x8s   ,"Matrix multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matmul16x16s ,"Matrix multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matmulnxns   ,"Matrix multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matmulnxms   ,"Matrix multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matmul2x2sf  ,"Matrix multiply (streaming order real floating-point data)");
ANNOTATE_FUN(matmul3x3sf  ,"Matrix multiply (streaming order real floating-point data)");
ANNOTATE_FUN(matmul4x4sf  ,"Matrix multiply (streaming order real floating-point data)");
ANNOTATE_FUN(matmul8x8sf  ,"Matrix multiply (streaming order real floating-point data)");
ANNOTATE_FUN(matmulnxmsf  ,"Matrix multiply (streaming order real floating-point data)");

ANNOTATE_FUN(matvmul2x2s   ,"Matrix-vector multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matvmul3x3s   ,"Matrix-vector multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matvmul4x4s   ,"Matrix-vector multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matvmul8x8s   ,"Matrix-vector multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matvmul16x16s ,"Matrix-vector multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matvmulnxns   ,"Matrix-vector multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matvmulnxms   ,"Matrix-vector multiply (streaming order real 16-bit data)");
ANNOTATE_FUN(matvmul2x2sf  ,"Matrix-vector multiply (streaming order real floating-point data)");
ANNOTATE_FUN(matvmul3x3sf  ,"Matrix-vector multiply (streaming order real floating-point data)");
ANNOTATE_FUN(matvmul4x4sf  ,"Matrix-vector multiply (streaming order real floating-point data)");
ANNOTATE_FUN(matvmul8x8sf  ,"Matrix-vector multiply (streaming order real floating-point data)");
ANNOTATE_FUN(matvmulnxmsf  ,"Matrix-vector multiply (streaming order real floating-point data)");

ANNOTATE_FUN(cmatmul2x2n   ,"Matrix multiply (block order complex 16-bit data)");
ANNOTATE_FUN(cmatmul3x3n   ,"Matrix multiply (block order complex 16-bit data)"); 
ANNOTATE_FUN(cmatmul4x4n   ,"Matrix multiply (block order complex 16-bit data)");
ANNOTATE_FUN(cmatmul8x8n   ,"Matrix multiply (block order complex 16-bit data)");
ANNOTATE_FUN(cmatmul16x16n ,"Matrix multiply (block order complex 16-bit data)");
ANNOTATE_FUN(cmatmulnxnn   ,"Matrix multiply (block order complex 16-bit data)");
ANNOTATE_FUN(cmatmulnxmn   ,"Matrix multiply (block order complex 16-bit data)");
ANNOTATE_FUN(cmatmul2x2nf  ,"Matrix multiply (block order complex floating-point data)");
ANNOTATE_FUN(cmatmul3x3nf  ,"Matrix multiply (block order complex floating-point data)");
ANNOTATE_FUN(cmatmul4x4nf  ,"Matrix multiply (block order complex floating-point data)");
ANNOTATE_FUN(cmatmul8x8nf  ,"Matrix multiply (block order complex floating-point data)");
ANNOTATE_FUN(cmatmulnxmnf  ,"Matrix multiply (block order complex floating-point data)");

ANNOTATE_FUN(cmatvmul2x2n  ,"Matrix-vector multiply (block order complex 16-bit data)");
ANNOTATE_FUN(cmatvmul3x3n  ,"Matrix-vector multiply (block order complex 16-bit data)");
ANNOTATE_FUN(cmatvmul4x4n  ,"Matrix-vector multiply (block order complex 16-bit data)");
ANNOTATE_FUN(cmatvmul8x8n  ,"Matrix-vector multiply (block order complex 16-bit data)");
ANNOTATE_FUN(cmatvmul16x16n,"Matrix-vector multiply (block order complex 16-bit data)");
ANNOTATE_FUN(cmatvmulnxnn  ,"Matrix-vector multiply (block order complex 16-bit data)");
ANNOTATE_FUN(cmatvmulnxmn  ,"Matrix-vector multiply (block order complex 16-bit data)");
ANNOTATE_FUN(cmatvmul2x2nf ,"Matrix-vector multiply (block order complex floating-point data)");
ANNOTATE_FUN(cmatvmul3x3nf ,"Matrix-vector multiply (block order complex floating-point data)");
ANNOTATE_FUN(cmatvmul4x4nf ,"Matrix-vector multiply (block order complex floating-point data)");
ANNOTATE_FUN(cmatvmul8x8nf ,"Matrix-vector multiply (block order complex floating-point data)");
ANNOTATE_FUN(cmatvmulnxmnf ,"Matrix-vector multiply (block order complex floating-point data)");

ANNOTATE_FUN(cmatmul2x2s   ,"Matrix multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatmul3x3s   ,"Matrix multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatmul4x4s   ,"Matrix multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatmul8x8s   ,"Matrix multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatmul16x16s ,"Matrix multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatmulnxns   ,"Matrix multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatmulnxms   ,"Matrix multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatmul2x2sf  ,"Matrix multiply (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatmul3x3sf  ,"Matrix multiply (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatmul4x4sf  ,"Matrix multiply (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatmul8x8sf  ,"Matrix multiply (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatmulnxmsf  ,"Matrix multiply (streaming order complex floating-point data)");

ANNOTATE_FUN(cmatvmul2x2s  ,"Matrix-vector multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatvmul3x3s  ,"Matrix-vector multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatvmul4x4s  ,"Matrix-vector multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatvmul8x8s  ,"Matrix-vector multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatvmul16x16s,"Matrix-vector multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatvmulnxns  ,"Matrix-vector multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatvmulnxms  ,"Matrix-vector multiply (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatvmul2x2sf ,"Matrix-vector multiply (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatvmul3x3sf ,"Matrix-vector multiply (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatvmul4x4sf ,"Matrix-vector multiply (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatvmul8x8sf ,"Matrix-vector multiply (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatvmulnxmsf ,"Matrix-vector multiply (streaming order complex floating-point data)");

ANNOTATE_FUN(rcmatmul2x2n  ,"Matrix multiply (block order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatmul3x3n  ,"Matrix multiply (block order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatmul4x4n  ,"Matrix multiply (block order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatmul8x8n  ,"Matrix multiply (block order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatmul16x16n,"Matrix multiply (block order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatmulnxnn  ,"Matrix multiply (block order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatmulnxmn  ,"Matrix multiply (block order real/complex 16-bit data)");

ANNOTATE_FUN(rcmatvmul2x2n  ,"Matrix-vector multiply (block order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatvmul3x3n  ,"Matrix-vector multiply (block order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatvmul4x4n  ,"Matrix-vector multiply (block order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatvmul8x8n  ,"Matrix-vector multiply (block order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatvmul16x16n,"Matrix-vector multiply (block order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatvmulnxnn  ,"Matrix-vector multiply (block order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatvmulnxmn  ,"Matrix-vector multiply (block order real/complex 16-bit data)");

ANNOTATE_FUN(rcmatmul2x2s   ,"Matrix multiply (streaming order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatmul3x3s   ,"Matrix multiply (streaming order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatmul4x4s   ,"Matrix multiply (streaming order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatmul8x8s   ,"Matrix multiply (streaming order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatmul16x16s ,"Matrix multiply (streaming order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatmulnxns   ,"Matrix multiply (streaming order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatmulnxms   ,"Matrix multiply (streaming order real/complex 16-bit data)");

ANNOTATE_FUN(rcmatvmul2x2s  ,"Matrix-vector multiply (streaming order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatvmul3x3s  ,"Matrix-vector multiply (streaming order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatvmul4x4s  ,"Matrix-vector multiply (streaming order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatvmul8x8s  ,"Matrix-vector multiply (streaming order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatvmul16x16s,"Matrix-vector multiply (streaming order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatvmulnxns  ,"Matrix-vector multiply (streaming order real/complex 16-bit data)");
ANNOTATE_FUN(rcmatvmulnxms  ,"Matrix-vector multiply (streaming order real/complex 16-bit data)");

ANNOTATE_FUN(mattran2x2n    ,"Matrix transpose (block order real 16-bit data)");
ANNOTATE_FUN(mattran4x4n    ,"Matrix transpose (block order real 16-bit data)");
ANNOTATE_FUN(mattran8x8n    ,"Matrix transpose (block order real 16-bit data)");
ANNOTATE_FUN(mattran16x16n  ,"Matrix transpose (block order real 16-bit data)");
ANNOTATE_FUN(mattrannxnn    ,"Matrix transpose (block order real 16-bit data)");
ANNOTATE_FUN(mattrannxmn    ,"Matrix transpose (block order real 16-bit data)");
ANNOTATE_FUN(mattran2x2nf   ,"Matrix transpose (block order real floating-point data)");
ANNOTATE_FUN(mattran4x4nf   ,"Matrix transpose (block order real floating-point data)");
ANNOTATE_FUN(mattran8x8nf   ,"Matrix transpose (block order real floating-point data)");
ANNOTATE_FUN(mattrannxmnf   ,"Matrix transpose (block order real floating-point data)");

ANNOTATE_FUN(mattran2x2s    ,"Matrix transpose (streaming order real 16-bit data)");
ANNOTATE_FUN(mattran4x4s    ,"Matrix transpose (streaming order real 16-bit data)");
ANNOTATE_FUN(mattran8x8s    ,"Matrix transpose (streaming order real 16-bit data)");
ANNOTATE_FUN(mattran16x16s  ,"Matrix transpose (streaming order real 16-bit data)");
ANNOTATE_FUN(mattrannxns    ,"Matrix transpose (streaming order real 16-bit data)");
ANNOTATE_FUN(mattrannxms    ,"Matrix transpose (streaming order real 16-bit data)");
ANNOTATE_FUN(mattran2x2sf   ,"Matrix transpose (streaming order real floating-point data)");
ANNOTATE_FUN(mattran4x4sf   ,"Matrix transpose (streaming order real floating-point data)");
ANNOTATE_FUN(mattran8x8sf   ,"Matrix transpose (streaming order real floating-point data)");
ANNOTATE_FUN(mattrannxmsf   ,"Matrix transpose (streaming order real floating-point data)");

ANNOTATE_FUN(cmattran2x2n   ,"Matrix transpose (block order complex 16-bit data)");
ANNOTATE_FUN(cmattran4x4n   ,"Matrix transpose (block order complex 16-bit data)");
ANNOTATE_FUN(cmattran8x8n   ,"Matrix transpose (block order complex 16-bit data)");
ANNOTATE_FUN(cmattran16x16n ,"Matrix transpose (block order complex 16-bit data)");
ANNOTATE_FUN(cmattrannxnn   ,"Matrix transpose (block order complex 16-bit data)");
ANNOTATE_FUN(cmattrannxmn   ,"Matrix transpose (block order complex 16-bit data)");
ANNOTATE_FUN(cmattran2x2nf  ,"Matrix transpose (block order complex floating-point data)");
ANNOTATE_FUN(cmattran4x4nf  ,"Matrix transpose (block order complex floating-point data)");
ANNOTATE_FUN(cmattran8x8nf  ,"Matrix transpose (block order complex floating-point data)");
ANNOTATE_FUN(cmattrannxmnf  ,"Matrix transpose (block order complex floating-point data)");

ANNOTATE_FUN(cmattran2x2s   ,"Matrix transpose (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmattran4x4s   ,"Matrix transpose (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmattran8x8s   ,"Matrix transpose (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmattran16x16s ,"Matrix transpose (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmattrannxns   ,"Matrix transpose (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmattrannxms   ,"Matrix transpose (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmattran2x2sf  ,"Matrix transpose (streaming order complex floating-point data)");
ANNOTATE_FUN(cmattran4x4sf  ,"Matrix transpose (streaming order complex floating-point data)");
ANNOTATE_FUN(cmattran8x8sf  ,"Matrix transpose (streaming order complex floating-point data)");
ANNOTATE_FUN(cmattrannxmsf  ,"Matrix transpose (streaming order complex floating-point data)");

ANNOTATE_FUN(cmatherm2x2n   ,"Matrix Hermitian (block order complex 16-bit data)");
ANNOTATE_FUN(cmatherm4x4n   ,"Matrix Hermitian (block order complex 16-bit data)");
ANNOTATE_FUN(cmatherm8x8n   ,"Matrix Hermitian (block order complex 16-bit data)");
ANNOTATE_FUN(cmatherm16x16n ,"Matrix Hermitian (block order complex 16-bit data)");
ANNOTATE_FUN(cmathermnxnn   ,"Matrix Hermitian (block order complex 16-bit data)");
ANNOTATE_FUN(cmathermnxmn   ,"Matrix Hermitian (block order complex 16-bit data)");
ANNOTATE_FUN(cmatherm2x2nf  ,"Matrix Hermitian (block order complex floating-point data)");
ANNOTATE_FUN(cmatherm4x4nf  ,"Matrix Hermitian (block order complex floating-point data)");
ANNOTATE_FUN(cmatherm8x8nf  ,"Matrix Hermitian (block order complex floating-point data)");
ANNOTATE_FUN(cmathermnxmnf  ,"Matrix Hermitian (block order complex floating-point data)");

ANNOTATE_FUN(cmatherm2x2s   ,"Matrix Hermitian (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatherm4x4s   ,"Matrix Hermitian (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatherm8x8s   ,"Matrix Hermitian (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatherm16x16s ,"Matrix Hermitian (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmathermnxns   ,"Matrix Hermitian (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmathermnxms   ,"Matrix Hermitian (streaming order complex 16-bit data)");
ANNOTATE_FUN(cmatherm2x2sf  ,"Matrix Hermitian (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatherm4x4sf  ,"Matrix Hermitian (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatherm8x8sf  ,"Matrix Hermitian (streaming order complex floating-point data)");
ANNOTATE_FUN(cmathermnxmsf  ,"Matrix Hermitian (streaming order complex floating-point data)");

ANNOTATE_FUN(rsb2x2  , "Stream to block order conversion (real 16-bit data)");
ANNOTATE_FUN(rsb3x3  , "Stream to block order conversion (real 16-bit data)");
ANNOTATE_FUN(rsb4x4  , "Stream to block order conversion (real 16-bit data)");
ANNOTATE_FUN(rsb8x8  , "Stream to block order conversion (real 16-bit data)");
ANNOTATE_FUN(rsb2x1  , "Stream to block order conversion (real 16-bit data)");
ANNOTATE_FUN(rsb4x1  , "Stream to block order conversion (real 16-bit data)");
ANNOTATE_FUN(rsb8x1  , "Stream to block order conversion (real 16-bit data)");
ANNOTATE_FUN(rsbmxn  , "Stream to block order conversion (real 16-bit data)");
ANNOTATE_FUN(rsb2x2f , "Stream to block order conversion (real floating-point data)");
ANNOTATE_FUN(rsb3x3f , "Stream to block order conversion (real floating-point data)");
ANNOTATE_FUN(rsb4x4f , "Stream to block order conversion (real floating-point data)");
ANNOTATE_FUN(rsb8x8f , "Stream to block order conversion (real floating-point data)");
ANNOTATE_FUN(rsb2x1f , "Stream to block order conversion (real floating-point data)");
ANNOTATE_FUN(rsb4x1f , "Stream to block order conversion (real floating-point data)");
ANNOTATE_FUN(rsb8x1f , "Stream to block order conversion (real floating-point data)");
ANNOTATE_FUN(rsbmxnf , "Stream to block order conversion (real floating-point data)");

ANNOTATE_FUN(csb2x2  , "Stream to block order conversion (complex 16-bit data)");
ANNOTATE_FUN(csb3x3  , "Stream to block order conversion (complex 16-bit data)");
ANNOTATE_FUN(csb4x4  , "Stream to block order conversion (complex 16-bit data)");
ANNOTATE_FUN(csb8x8  , "Stream to block order conversion (complex 16-bit data)");
ANNOTATE_FUN(csb2x1  , "Stream to block order conversion (complex 16-bit data)");
ANNOTATE_FUN(csb4x1  , "Stream to block order conversion (complex 16-bit data)");
ANNOTATE_FUN(csb8x1  , "Stream to block order conversion (complex 16-bit data)");
ANNOTATE_FUN(csbmxn  , "Stream to block order conversion (complex 16-bit data)");
ANNOTATE_FUN(csb2x2f , "Stream to block order conversion (complex floating-point data)");
ANNOTATE_FUN(csb3x3f , "Stream to block order conversion (complex floating-point data)");
ANNOTATE_FUN(csb4x4f , "Stream to block order conversion (complex floating-point data)");
ANNOTATE_FUN(csb8x8f , "Stream to block order conversion (complex floating-point data)");
ANNOTATE_FUN(csb2x1f , "Stream to block order conversion (complex floating-point data)");
ANNOTATE_FUN(csb4x1f , "Stream to block order conversion (complex floating-point data)");
ANNOTATE_FUN(csb8x1f , "Stream to block order conversion (complex floating-point data)");
ANNOTATE_FUN(csbmxnf , "Stream to block order conversion (complex floating-point data)");

ANNOTATE_FUN(rbs2x2  , "Block to stream format conversion (real 16-bit data)");
ANNOTATE_FUN(rbs3x3  , "Block to stream format conversion (real 16-bit data)");
ANNOTATE_FUN(rbs4x4  , "Block to stream format conversion (real 16-bit data)");
ANNOTATE_FUN(rbs8x8  , "Block to stream format conversion (real 16-bit data)");
ANNOTATE_FUN(rbs2x1  , "Block to stream format conversion (real 16-bit data)");
ANNOTATE_FUN(rbs4x1  , "Block to stream format conversion (real 16-bit data)");
ANNOTATE_FUN(rbs8x1  , "Block to stream format conversion (real 16-bit data)");
ANNOTATE_FUN(rbsmxn  , "Block to stream format conversion (real 16-bit data)");
ANNOTATE_FUN(rbs2x2f , "Block to stream format conversion (real floating-point data)");
ANNOTATE_FUN(rbs3x3f , "Block to stream format conversion (real floating-point data)");
ANNOTATE_FUN(rbs4x4f , "Block to stream format conversion (real floating-point data)");
ANNOTATE_FUN(rbs8x8f , "Block to stream format conversion (real floating-point data)");
ANNOTATE_FUN(rbs2x1f , "Block to stream format conversion (real floating-point data)");
ANNOTATE_FUN(rbs4x1f , "Block to stream format conversion (real floating-point data)");
ANNOTATE_FUN(rbs8x1f , "Block to stream format conversion (real floating-point data)");
ANNOTATE_FUN(rbsmxnf , "Block to stream format conversion (real floating-point data)");

ANNOTATE_FUN(cbs2x2  , "Block to stream format conversion (complex 16-bit data)");
ANNOTATE_FUN(cbs3x3  , "Block to stream format conversion (complex 16-bit data)");
ANNOTATE_FUN(cbs4x4  , "Block to stream format conversion (complex 16-bit data)");
ANNOTATE_FUN(cbs8x8  , "Block to stream format conversion (complex 16-bit data)");
ANNOTATE_FUN(cbs2x1  , "Block to stream format conversion (complex 16-bit data)");
ANNOTATE_FUN(cbs4x1  , "Block to stream format conversion (complex 16-bit data)");
ANNOTATE_FUN(cbs8x1  , "Block to stream format conversion (complex 16-bit data)");
ANNOTATE_FUN(cbsmxn  , "Block to stream format conversion (complex 16-bit data)");
ANNOTATE_FUN(cbs2x2f , "Block to stream format conversion (complex floating-point data)");
ANNOTATE_FUN(cbs3x3f , "Block to stream format conversion (complex floating-point data)");
ANNOTATE_FUN(cbs4x4f , "Block to stream format conversion (complex floating-point data)");
ANNOTATE_FUN(cbs8x8f , "Block to stream format conversion (complex floating-point data)");
ANNOTATE_FUN(cbs2x1f , "Block to stream format conversion (complex floating-point data)");
ANNOTATE_FUN(cbs4x1f , "Block to stream format conversion (complex floating-point data)");
ANNOTATE_FUN(cbs8x1f , "Block to stream format conversion (complex floating-point data)");
ANNOTATE_FUN(cbsmxnf , "Block to stream format conversion (complex floating-point data)");

ANNOTATE_FUN(cinterleave2  ,"Interleaving (complex 16-bit data)");
ANNOTATE_FUN(cinterleave3  ,"Interleaving (complex 16-bit data)");
ANNOTATE_FUN(cinterleave4  ,"Interleaving (complex 16-bit data)");
ANNOTATE_FUN(cinterleave8  ,"Interleaving (complex 16-bit data)");
ANNOTATE_FUN(cinterleave2f ,"Interleaving (complex floating-point data)");
ANNOTATE_FUN(cinterleave3f ,"Interleaving (complex floating-point data)");
ANNOTATE_FUN(cinterleave4f ,"Interleaving (complex floating-point data)");
ANNOTATE_FUN(cinterleavemf ,"Interleaving (complex floating-point data)");

ANNOTATE_FUN(rinterleave2  ,"Interleaving (real 16-bit data)");
ANNOTATE_FUN(rinterleave3  ,"Interleaving (real 16-bit data)");
ANNOTATE_FUN(rinterleave4  ,"Interleaving (real 16-bit data)");
ANNOTATE_FUN(rinterleave8  ,"Interleaving (real 16-bit data)");
ANNOTATE_FUN(rinterleave2f ,"Interleaving (real floating-point data)");
ANNOTATE_FUN(rinterleave3f ,"Interleaving (real floating-point data)");
ANNOTATE_FUN(rinterleave4f ,"Interleaving (real floating-point data)");
ANNOTATE_FUN(rinterleavemf ,"Interleaving (real floating-point data)");

ANNOTATE_FUN(cdeinterleave2  ,"Deinterleaving (complex 16-bit data)");
ANNOTATE_FUN(cdeinterleave3  ,"Deinterleaving (complex 16-bit data)");
ANNOTATE_FUN(cdeinterleave4  ,"Deinterleaving (complex 16-bit data)");
ANNOTATE_FUN(cdeinterleave8  ,"Deinterleaving (complex 16-bit data)");
ANNOTATE_FUN(cdeinterleave2f ,"Deinterleaving (complex floating-point data)");
ANNOTATE_FUN(cdeinterleave3f ,"Deinterleaving (complex floating-point data)");
ANNOTATE_FUN(cdeinterleave4f ,"Deinterleaving (complex floating-point data)");
ANNOTATE_FUN(cdeinterleavemf ,"Deinterleaving (complex floating-point data)");

ANNOTATE_FUN(rdeinterleave2  ,"Deinterleaving (real 16-bit data)");
ANNOTATE_FUN(rdeinterleave3  ,"Deinterleaving (real 16-bit data)");
ANNOTATE_FUN(rdeinterleave4  ,"Deinterleaving (real 16-bit data)");
ANNOTATE_FUN(rdeinterleave8  ,"Deinterleaving (real 16-bit data)");
ANNOTATE_FUN(rdeinterleave2f ,"Deinterleaving (real floating-point data)");
ANNOTATE_FUN(rdeinterleave3f ,"Deinterleaving (real floating-point data)");
ANNOTATE_FUN(rdeinterleave4f ,"Deinterleaving (real floating-point data)");
ANNOTATE_FUN(rdeinterleavemf ,"Deinterleaving (real floating-point data)");
