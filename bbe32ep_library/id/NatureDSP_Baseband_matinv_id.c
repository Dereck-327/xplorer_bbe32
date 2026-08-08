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
 * Matrix Decomposition and Inversion Functions
 * Annotations
 */
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_matinv.h"
#include "common.h"


ANNOTATE_FUN(cmatinv2x2n ,"Direct matrix inversion (block order complex 16-bit data)");
ANNOTATE_FUN(cmatinv4x4n ,"Direct matrix inversion (block order complex 16-bit data)");
ANNOTATE_FUN(cmatinv2x2nf,"Direct matrix inversion (block order complex floating-point data)");
ANNOTATE_FUN(cmatinv3x3nf,"Direct matrix inversion (block order complex floating-point data)");
ANNOTATE_FUN(cmatinv4x4nf,"Direct matrix inversion (block order complex floating-point data)");
ANNOTATE_FUN(cmatinv2x2sf,"Direct matrix inversion (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatinv3x3sf,"Direct matrix inversion (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatinv4x4sf,"Direct matrix inversion (streaming order complex floating-point data)");
ANNOTATE_FUN(matinv2x2nf,"Direct matrix inversion (block order real floating-point data)");
ANNOTATE_FUN(matinv3x3nf,"Direct matrix inversion (block order real floating-point data)");
ANNOTATE_FUN(matinv4x4nf,"Direct matrix inversion (block order real floating-point data)");
ANNOTATE_FUN(matinv2x2sf,"Direct matrix inversion (streaming order real floating-point data)");
ANNOTATE_FUN(matinv3x3sf,"Direct matrix inversion (streaming order real floating-point data)");
ANNOTATE_FUN(matinv4x4sf,"Direct matrix inversion (streaming order real floating-point data)");
ANNOTATE_FUN(chermmatinv2x2n,"Direct Hermitian matrix inversion (block order complex 16-bit data)");
ANNOTATE_FUN(chermmatinv4x4n,"Direct Hermitian matrix inversion (block order complex 16-bit data)");

ANNOTATE_FUN(matinvgj8x8nf   ,"Gauss-Jordan inversion (block order real floating-point data)");
ANNOTATE_FUN(matinvgj16x16nf ,"Gauss-Jordan inversion (block order real floating-point data)");
ANNOTATE_FUN(matinvgjnxnnf   ,"Gauss-Jordan inversion (block order real floating-point data)");
ANNOTATE_FUN(cmatinvgj8x8nf  ,"Gauss-Jordan inversion (block order complex floating-point data)");
ANNOTATE_FUN(cmatinvgj16x16nf,"Gauss-Jordan inversion (block order complex floating-point data)");
ANNOTATE_FUN(cmatinvgjnxnnf  ,"Gauss-Jordan inversion (block order complex floating-point data)");

ANNOTATE_FUN(matinvgj2x2sf  ,"Gauss-Jordan inversion (streaming order real floating-point data)");
ANNOTATE_FUN(matinvgj3x3sf  ,"Gauss-Jordan inversion (streaming order real floating-point data)");
ANNOTATE_FUN(matinvgj4x4sf  ,"Gauss-Jordan inversion (streaming order real floating-point data)");
ANNOTATE_FUN(matinvgjnxnsf  ,"Gauss-Jordan inversion (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatinvgj2x2sf ,"Gauss-Jordan inversion (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatinvgj3x3sf ,"Gauss-Jordan inversion (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatinvgj4x4sf ,"Gauss-Jordan inversion (streaming order complex floating-point data)");
ANNOTATE_FUN(cmatinvgjnxnsf ,"Gauss-Jordan inversion (streaming order complex floating-point data)");

ANNOTATE_FUN(lu8x8nf   ,"LU decomposition (block order real floating-point data)");
ANNOTATE_FUN(lu16x16nf ,"LU decomposition (block order real floating-point data)");
ANNOTATE_FUN(lunxnnf   ,"LU decomposition (block order real floating-point data)");
ANNOTATE_FUN(clu8x8nf  ,"LU decomposition (block order complex floating-point data)");
ANNOTATE_FUN(clu16x16nf,"LU decomposition (block order complex floating-point data)");
ANNOTATE_FUN(clunxnnf  ,"LU decomposition (block order complex floating-point data)");
ANNOTATE_FUN(lu2x2sf   ,"LU decomposition (streaming order real floating-point data)");
ANNOTATE_FUN(lu3x3sf   ,"LU decomposition (streaming order real floating-point data)");
ANNOTATE_FUN(lu4x4sf   ,"LU decomposition (streaming order real floating-point data)");
ANNOTATE_FUN(lunxnsf   ,"LU decomposition (streaming order real floating-point data)");
ANNOTATE_FUN(clu2x2sf  ,"LU decomposition (streaming order complex floating-point data)");
ANNOTATE_FUN(clu3x3sf  ,"LU decomposition (streaming order complex floating-point data)");
ANNOTATE_FUN(clu4x4sf  ,"LU decomposition (streaming order complex floating-point data)");
ANNOTATE_FUN(clunxnsf  ,"LU decomposition (streaming order complex floating-point data)");

ANNOTATE_FUN(det8x8nf   ,"Determinant (block order real floating-point data)");
ANNOTATE_FUN(det16x16nf ,"Determinant (block order real floating-point data)");
ANNOTATE_FUN(detnxnnf   ,"Determinant (block order real floating-point data)");
ANNOTATE_FUN(cdet8x8nf  ,"Determinant (block order complex floating-point data)");
ANNOTATE_FUN(cdet16x16nf,"Determinant (block order complex floating-point data)");
ANNOTATE_FUN(cdetnxnnf  ,"Determinant (block order complex floating-point data)");

ANNOTATE_FUN(det2x2sf   ,"Determinant (streaming order real floating-point data)");
ANNOTATE_FUN(det3x3sf   ,"Determinant (streaming order real floating-point data)");
ANNOTATE_FUN(det4x4sf   ,"Determinant (streaming order real floating-point data)");
ANNOTATE_FUN(detnxnsf   ,"Determinant (streaming order real floating-point data)");
ANNOTATE_FUN(cdet2x2sf  ,"Determinant (streaming order real floating-point data)");
ANNOTATE_FUN(cdet3x3sf  ,"Determinant (streaming order real floating-point data)");
ANNOTATE_FUN(cdet4x4sf  ,"Determinant (streaming order real floating-point data)");
ANNOTATE_FUN(cdetnxnsf  ,"Determinant (streaming order real floating-point data)");

ANNOTATE_FUN(chol2x2s,"Cholesky decomposition (streaming order complex 16-bit data)");  
ANNOTATE_FUN(chol4x4s,"Cholesky decomposition (streaming order complex 16-bit data)");  
ANNOTATE_FUN(chol8x8s,"Cholesky decomposition (streaming order complex 16-bit data)");  
ANNOTATE_FUN(cholmxns,"Cholesky decomposition (streaming order complex 16-bit data)");  
ANNOTATE_FUN(cholfwd2x2x1s,"Cholesky forward substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cholfwd4x4x1s,"Cholesky forward substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cholfwd8x8x1s,"Cholesky forward substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cholfwdmxnxps,"Cholesky forward substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cholbkw2x1s,"Cholesky back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cholbkw4x1s,"Cholesky back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cholbkw8x1s,"Cholesky back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cholbkwnxps,"Cholesky back substitution (streaming order complex 16-bit data)");


ANNOTATE_FUN(cholmxnn  ,"Cholesky decomposition (block order complex 16-bit data)");    
ANNOTATE_FUN(chol8x8n  ,"Cholesky decomposition (block order complex 16-bit data)");
ANNOTATE_FUN(chol16x16n,"Cholesky decomposition (block order complex 16-bit data)");
ANNOTATE_FUN(chol32x32n,"Cholesky decomposition (block order complex 16-bit data)");
ANNOTATE_FUN(cholfwdmxnxpn   ,"Cholesky forward substitution (block order complex 16-bit data)");
ANNOTATE_FUN(cholfwd8x8x1n   ,"Cholesky forward substitution (block order complex 16-bit data)");
ANNOTATE_FUN(cholfwd16x16x1n ,"Cholesky forward substitution (block order complex 16-bit data)");
ANNOTATE_FUN(cholfwd32x32x1n ,"Cholesky forward substitution (block order complex 16-bit data)");

ANNOTATE_FUN(cholbkwnxpn  ,"Cholesky back substitution (block order complex 16-bit data)");
ANNOTATE_FUN(cholbkw8x1n  ,"Cholesky back substitution (block order complex 16-bit data)");
ANNOTATE_FUN(cholbkw16x1n ,"Cholesky back substitution (block order complex 16-bit data)");
ANNOTATE_FUN(cholbkw32x1n ,"Cholesky back substitution (block order complex 16-bit data)");

ANNOTATE_FUN(chol2x2sf ,"Cholesky decomposition (streaming order complex floating-point data)");
ANNOTATE_FUN(chol3x3sf ,"Cholesky decomposition (streaming order complex floating-point data)");
ANNOTATE_FUN(chol4x4sf ,"Cholesky decomposition (streaming order complex floating-point data)");
ANNOTATE_FUN(cholmxnsf ,"Cholesky decomposition (streaming order complex floating-point data)");
ANNOTATE_FUN(rchol2x2sf,"Cholesky decomposition (streaming order real floating-point data)");
ANNOTATE_FUN(rchol3x3sf,"Cholesky decomposition (streaming order real floating-point data)");
ANNOTATE_FUN(rchol4x4sf,"Cholesky decomposition (streaming order real floating-point data)");
ANNOTATE_FUN(rcholmxnsf,"Cholesky decomposition (streaming order real floating-point data)");
ANNOTATE_FUN(cholfwd2x2x1sf,"Cholesky forward substitution (streaming order complex floating-point data)");
ANNOTATE_FUN(cholfwd3x3x1sf,"Cholesky forward substitution (streaming order complex floating-point data)");
ANNOTATE_FUN(cholfwd4x4x1sf,"Cholesky forward substitution (streaming order complex floating-point data)");
ANNOTATE_FUN(cholfwdmxnxpsf,"Cholesky forward substitution (streaming order complex floating-point data)");
ANNOTATE_FUN(rcholfwd2x2x1sf,"Cholesky forward substitution (streaming order real floating-point data)");
ANNOTATE_FUN(rcholfwd3x3x1sf,"Cholesky forward substitution (streaming order real floating-point data)");
ANNOTATE_FUN(rcholfwd4x4x1sf,"Cholesky forward substitution (streaming order real floating-point data)");
ANNOTATE_FUN(rcholfwdmxnxpsf,"Cholesky forward substitution (streaming order real floating-point data)");

ANNOTATE_FUN(cholbkw2x1sf,"Cholesky back substitution (streaming order complex floating-point data)");
ANNOTATE_FUN(cholbkw3x1sf,"Cholesky back substitution (streaming order complex floating-point data)");
ANNOTATE_FUN(cholbkw4x1sf,"Cholesky back substitution (streaming order complex floating-point data)");
ANNOTATE_FUN(cholbkwnxpsf,"Cholesky back substitution (streaming order complex floating-point data)");
ANNOTATE_FUN(rcholbkw2x1sf,"Cholesky back substitution (streaming order real floating-point data)");
ANNOTATE_FUN(rcholbkw3x1sf,"Cholesky back substitution (streaming order real floating-point data)");
ANNOTATE_FUN(rcholbkw4x1sf,"Cholesky back substitution (streaming order real floating-point data)");
ANNOTATE_FUN(rcholbkwnxpsf,"Cholesky back substitution (streaming order real floating-point data)");

ANNOTATE_FUN(cholmmse2x2x1sf,"Cholesky MMSE (streaming order complex floating-point data)");
ANNOTATE_FUN(cholmmse3x3x1sf,"Cholesky MMSE (streaming order complex floating-point data)");
ANNOTATE_FUN(cholmmse4x4x1sf,"Cholesky MMSE (streaming order complex floating-point data)");
ANNOTATE_FUN(cholmmsemxnxpsf,"Cholesky MMSE (streaming order complex floating-point data)");
ANNOTATE_FUN(rcholmmse2x2x1sf,"Cholesky MMSE(streaming order real floating-point data)");
ANNOTATE_FUN(rcholmmse3x3x1sf,"Cholesky MMSE(streaming order real floating-point data)");
ANNOTATE_FUN(rcholmmse4x4x1sf,"Cholesky MMSE(streaming order real floating-point data)");
ANNOTATE_FUN(rcholmmsemxnxpsf,"Cholesky MMSE(streaming order real floating-point data)");

ANNOTATE_FUN(chol8x8nf   ,"Cholesky decomposition (block order complex floating-point data)");
ANNOTATE_FUN(chol16x16nf ,"Cholesky decomposition (block order complex floating-point data)");
ANNOTATE_FUN(cholmxnnf   ,"Cholesky decomposition (block order complex floating-point data)");
ANNOTATE_FUN(rchol8x8nf  ,"Cholesky decomposition (block order real floating-point data)");
ANNOTATE_FUN(rchol16x16nf,"Cholesky decomposition (block order real floating-point data)");
ANNOTATE_FUN(rcholmxnnf  ,"Cholesky decomposition (block order real floating-point data)");
ANNOTATE_FUN(cholfwd8x8x1nf   ,"Cholesky forward substitution (block order complex floating-point data)");
ANNOTATE_FUN(cholfwd16x16x1nf ,"Cholesky forward substitution (block order complex floating-point data)");
ANNOTATE_FUN(cholfwdmxnxpnf   ,"Cholesky forward substitution (block order complex floating-point data)");
ANNOTATE_FUN(rcholfwd8x8x1nf  ,"Cholesky forward substitution (block order real floating-point data)");
ANNOTATE_FUN(rcholfwd16x16x1nf,"Cholesky forward substitution (block order real floating-point data)");
ANNOTATE_FUN(rcholfwdmxnxpnf  ,"Cholesky forward substitution (block order real floating-point data)");

ANNOTATE_FUN(cholbkw8x1nf  ,"Cholesky MMSE (block order complex floating-point data)");
ANNOTATE_FUN(cholbkw16x1nf ,"Cholesky MMSE (block order complex floating-point data)");
ANNOTATE_FUN(cholbkwnxpnf  ,"Cholesky MMSE (block order complex floating-point data)");
ANNOTATE_FUN(rcholbkw8x1nf ,"Cholesky back substitution (block order real floating-point data)");
ANNOTATE_FUN(rcholbkw16x1nf,"Cholesky back substitution (block order real floating-point data)");
ANNOTATE_FUN(rcholbkwnxpnf ,"Cholesky back substitution (block order real floating-point data)");

ANNOTATE_FUN(cholmmse8x8x1nf    ,"Cholesky MMSE (block order, complex floating-point data)");
ANNOTATE_FUN(cholmmse16x16x1nf  ,"Cholesky MMSE (block order, complex floating-point data)");
ANNOTATE_FUN(cholmmsemxnxpnf    ,"Cholesky MMSE (block order, complex floating-point data)");
ANNOTATE_FUN(rcholmmse8x8x1nf   ,"Cholesky MMSE (block order, real floating-point data)");
ANNOTATE_FUN(rcholmmse16x16x1nf ,"Cholesky MMSE (block order, real floating-point data)");
ANNOTATE_FUN(rcholmmsemxnxpnf   ,"Cholesky MMSE (block order, real floating-point data)");

ANNOTATE_FUN(bcholwxnn  ,"Cholesky decomposition (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bchol16x32n,"Cholesky decomposition (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bchol8x32n ,"Cholesky decomposition (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bchol8x16n ,"Cholesky decomposition (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bchol4x32n ,"Cholesky decomposition (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bchol4x16n ,"Cholesky decomposition (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bcholfwdwxnxpn  ,"Cholesky forward substitution (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bcholfwd16x32x1n,"Cholesky forward substitution (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bcholfwd8x32x1n ,"Cholesky forward substitution (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bcholfwd8x16x1n ,"Cholesky forward substitution (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bcholfwd4x32x1n ,"Cholesky forward substitution (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bcholfwd4x16x1n ,"Cholesky forward substitution (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bcholbkwwxnxpn  ,"Cholesky back substitution (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bcholbkw16x32x1n,"Cholesky back substitution (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bcholbkw8x32x1n ,"Cholesky back substitution (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bcholbkw8x16x1n ,"Cholesky back substitution (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bcholbkw4x32x1n ,"Cholesky back substitution (banded matrices complex 16-bit data)");
ANNOTATE_FUN(bcholbkw4x16x1n ,"Cholesky back substitution (banded matrices complex 16-bit data)");

ANNOTATE_FUN(cqr2x2s           ,"QR decomposition (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_r3x3s   ,"QR decomposition (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_r4x2s   ,"QR decomposition (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_r4x4s   ,"QR decomposition (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_r5x5s   ,"QR decomposition (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_r6x6s   ,"QR decomposition (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_r7x7s   ,"QR decomposition (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_r8x4s   ,"QR decomposition (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_r8x8s   ,"QR decomposition (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_r16x8s  ,"QR decomposition (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_r16x16s ,"QR decomposition (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_rmxns   ,"QR decomposition (streaming order complex 16-bit data)");
ANNOTATE_FUN(qr2x2s            ,"QR decomposition (streaming order real 16-bit data)");
ANNOTATE_FUN(qr_build_r4x4s    ,"QR decomposition (streaming order real 16-bit data)");
ANNOTATE_FUN(qr_build_r8x8s    ,"QR decomposition (streaming order real 16-bit data)");
ANNOTATE_FUN(qr_build_rmxns    ,"QR decomposition (streaming order real 16-bit data)");


ANNOTATE_FUN(cqr_calc_qb3x3x1s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb3x3x3s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb4x2x1s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb4x4x1s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb4x4x4s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb5x5x1s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb5x5x5s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb6x6x1s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb6x6x6s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb7x7x1s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb7x7x7s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb8x4x1s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb8x8x1s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb8x8x8s    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb16x8x1s   ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb16x16x1s  ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb16x16x16s ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qbmxnxps    ,"Rotation by Housholder vectors (streaming order complex 16-bit data)");
ANNOTATE_FUN(qr_calc_qb4x4x1s     ,"Rotation by Housholder vectors (streaming order real 16-bit data)");
ANNOTATE_FUN(qr_calc_qb4x4x4s     ,"Rotation by Housholder vectors (streaming order real 16-bit data)");
ANNOTATE_FUN(qr_calc_qb8x8x1s     ,"Rotation by Housholder vectors (streaming order real 16-bit data)");
ANNOTATE_FUN(qr_calc_qb8x8x8s     ,"Rotation by Housholder vectors (streaming order real 16-bit data)");
ANNOTATE_FUN(qr_calc_qbmxnxps     ,"Rotation by Housholder vectors (streaming order real 16-bit data)");

ANNOTATE_FUN(cqr_bkw2x1s  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw3x1s  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw3x3s  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw4x1s  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw4x4s  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw5x1s  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw5x5s  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw6x1s  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw6x6s  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw7x1s  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw7x7s  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw8x1s  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw8x8s  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw16x1s ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw16x16s,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkwnxps  ,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw2x2x1s,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw2x2x2s,"QR-based back substitution (streaming order complex 16-bit data)");
ANNOTATE_FUN(qr_bkw4x1s   ,"QR-based back substitution (streaming order real 16-bit data)");
ANNOTATE_FUN(qr_bkw4x4s   ,"QR-based back substitution (streaming order real 16-bit data)");
ANNOTATE_FUN(qr_bkw8x1s   ,"QR-based back substitution (streaming order real 16-bit data)");
ANNOTATE_FUN(qr_bkw8x8s   ,"QR-based back substitution (streaming order real 16-bit data)");
ANNOTATE_FUN(qr_bkwnxps   ,"QR-based back substitution (streaming order real 16-bit data)");
ANNOTATE_FUN(qr_bkw2x2x1s ,"QR-based back substitution (streaming order real 16-bit data)");
ANNOTATE_FUN(qr_bkw2x2x2s ,"QR-based back substitution (streaming order real 16-bit data)");

ANNOTATE_FUN(cqr_build_rmxnn     ,"QR decomposition (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_r8x8n     ,"QR decomposition (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_r16x16n   ,"QR decomposition (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_build_r32x32n   ,"QR decomposition (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qbmxnxpn   ,"Rotation by Housholder vectors (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb8x8x1n   ,"Rotation by Housholder vectors (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb8x8x8n   ,"Rotation by Housholder vectors (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb16x16x1n ,"Rotation by Housholder vectors (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb16x16x16n,"Rotation by Housholder vectors (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb32x32x1n ,"Rotation by Housholder vectors (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_calc_qb32x32x32n,"Rotation by Housholder vectors (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkwmxnxpn       ,"QR-based back substitution (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw8x8x1n       ,"QR-based back substitution (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw8x8x8n       ,"QR-based back substitution (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw16x16x1n     ,"QR-based back substitution (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw16x16x16n    ,"QR-based back substitution (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw32x32x1n     ,"QR-based back substitution (block order complex 16-bit data)");
ANNOTATE_FUN(cqr_bkw32x32x32n    ,"QR-based back substitution (block order complex 16-bit data)");

ANNOTATE_FUN(cqr_build_r2x2sf  ,"QR decomposition (streaming order complex floating-point data)");
ANNOTATE_FUN(cqr_build_r3x3sf  ,"QR decomposition (streaming order complex floating-point data)");
ANNOTATE_FUN(cqr_build_r4x4sf  ,"QR decomposition (streaming order complex floating-point data)");
ANNOTATE_FUN(cqr_build_rmxnsf  ,"QR decomposition (streaming order complex floating-point data)");
ANNOTATE_FUN(qr_build_r2x2sf   ,"QR decomposition (streaming order real floating-point data)");
ANNOTATE_FUN(qr_build_r3x3sf   ,"QR decomposition (streaming order real floating-point data)");
ANNOTATE_FUN(qr_build_r4x4sf   ,"QR decomposition (streaming order real floating-point data)");
ANNOTATE_FUN(qr_build_rmxnsf   ,"QR decomposition (streaming order real floating-point data)");

ANNOTATE_FUN(cqr_calc_qb2x2x1sf,"Rotation by Housholder vectors (streaming order complex floating-point data)");
ANNOTATE_FUN(cqr_calc_qb3x3x1sf,"Rotation by Housholder vectors (streaming order complex floating-point data)");
ANNOTATE_FUN(cqr_calc_qb4x4x1sf,"Rotation by Housholder vectors (streaming order complex floating-point data)");
ANNOTATE_FUN(cqr_calc_qbmxnxpsf,"Rotation by Housholder vectors (streaming order complex floating-point data)");
ANNOTATE_FUN(qr_calc_qb2x2x1sf ,"Rotation by Housholder vectors (streaming order real floating-point data)");
ANNOTATE_FUN(qr_calc_qb3x3x1sf ,"Rotation by Housholder vectors (streaming order real floating-point data)");
ANNOTATE_FUN(qr_calc_qb4x4x1sf ,"Rotation by Housholder vectors (streaming order real floating-point data)");
ANNOTATE_FUN(qr_calc_qbmxnxpsf ,"Rotation by Housholder vectors (streaming order real floating-point data)");

ANNOTATE_FUN(cqr_bkw2x1sf,"QR-based back substitution (streaming order complex floating-point data)");
ANNOTATE_FUN(cqr_bkw3x1sf,"QR-based back substitution (streaming order complex floating-point data)");
ANNOTATE_FUN(cqr_bkw4x1sf,"QR-based back substitution (streaming order complex floating-point data)");
ANNOTATE_FUN(cqr_bkwnxpsf,"QR-based back substitution (streaming order complex floating-point data)");
ANNOTATE_FUN(qr_bkw2x1sf ,"QR-based back substitution (streaming order real floating-point data)");
ANNOTATE_FUN(qr_bkw3x1sf ,"QR-based back substitution (streaming order real floating-point data)");
ANNOTATE_FUN(qr_bkw4x1sf ,"QR-based back substitution (streaming order real floating-point data)");
ANNOTATE_FUN(qr_bkwnxpsf ,"QR-based back substitution (streaming order real floating-point data)");

ANNOTATE_FUN(cqr_build_rmxnnf   ,"QR decomposition (block order complex floating-point data)");
ANNOTATE_FUN(cqr_build_r8x8nf   ,"QR decomposition (block order complex floating-point data)");
ANNOTATE_FUN(cqr_build_r16x16nf ,"QR decomposition (block order complex floating-point data)");
ANNOTATE_FUN(qr_build_rmxnnf    ,"QR decomposition (block order real floating-point data)");
ANNOTATE_FUN(qr_build_r8x8nf    ,"QR decomposition (block order real floating-point data)");
ANNOTATE_FUN(qr_build_r16x16nf  ,"QR decomposition (block order real floating-point data)");
ANNOTATE_FUN(cqr_calc_qbmxnxpnf   ,"Rotation by Housholder vectors (block order complex floating-point data)");
ANNOTATE_FUN(cqr_calc_qb8x8x1nf   ,"Rotation by Housholder vectors (block order complex floating-point data)");
ANNOTATE_FUN(cqr_calc_qb16x16x1nf ,"Rotation by Housholder vectors (block order complex floating-point data)");
ANNOTATE_FUN(qr_calc_qbmxnxpnf    ,"Rotation by Housholder vectors (block order real floating-point data)");
ANNOTATE_FUN(qr_calc_qb8x8x1nf    ,"Rotation by Housholder vectors (block order real floating-point data)");
ANNOTATE_FUN(qr_calc_qb16x16x1nf  ,"Rotation by Housholder vectors (block order real floating-point data)");
ANNOTATE_FUN(cqr_bkwmxnxpnf    ,"QR-based back substitution (block order complex floating-point data)");
ANNOTATE_FUN(cqr_bkw8x8x1nf    ,"QR-based back substitution (block order complex floating-point data)");
ANNOTATE_FUN(cqr_bkw16x16x1nf  ,"QR-based back substitution (block order complex floating-point data)");
ANNOTATE_FUN(qr_bkwmxnxpnf     ,"QR-based back substitution (block order real floating-point data)");
ANNOTATE_FUN(qr_bkw8x8x1nf     ,"QR-based back substitution (block order real floating-point data)");
ANNOTATE_FUN(qr_bkw16x16x1nf   ,"QR-based back substitution (block order real floating-point data)");

ANNOTATE_FUN( svdmxnnf  , "SVD (block order complex floating-point data)");
ANNOTATE_FUN( svd16x16nf, "SVD (block order complex floating-point data)");
ANNOTATE_FUN( svd8x8nf  , "SVD (block order complex floating-point data)");
ANNOTATE_FUN( svdmxnsf  , "SVD (streaming order complex floating-point data)");
ANNOTATE_FUN( svd4x4sf  , "SVD (streaming order complex floating-point data)");
ANNOTATE_FUN( svd3x3sf  , "SVD (streaming order complex floating-point data)");
ANNOTATE_FUN( svd2x2sf  , "SVD (streaming order complex floating-point data)");

ANNOTATE_FUN( rsvdmxnnf  , "SVD (block order real floating-point data)");
ANNOTATE_FUN( rsvd16x16nf, "SVD (block order real floating-point data)");
ANNOTATE_FUN( rsvd8x8nf  , "SVD (block order real floating-point data)");
ANNOTATE_FUN( rsvdmxnsf  , "SVD (streaming order real floating-point data)");
ANNOTATE_FUN( rsvd4x4sf  , "SVD (streaming order real floating-point data)");
ANNOTATE_FUN( rsvd3x3sf  , "SVD (streaming order real floating-point data)");
ANNOTATE_FUN( rsvd2x2sf  , "SVD (streaming order real floating-point data)");

ANNOTATE_FUN( eigennxnnf  , "Eigenvalues/eigenvectors (block order complex floating-point data)");
ANNOTATE_FUN( eigen16x16nf, "Eigenvalues/eigenvectors (block order complex floating-point data)");
ANNOTATE_FUN( eigen8x8nf  , "Eigenvalues/eigenvectors (block order complex floating-point data)");
ANNOTATE_FUN( eigennxnsf  , "Eigenvalues/eigenvectors (streaming order complex floating-point data)");
ANNOTATE_FUN( eigen4x4sf  , "Eigenvalues/eigenvectors (streaming order complex floating-point data)");
ANNOTATE_FUN( eigen3x3sf  , "Eigenvalues/eigenvectors (streaming order complex floating-point data)");
ANNOTATE_FUN( eigen2x2sf  , "Eigenvalues/eigenvectors (streaming order complex floating-point data)");
ANNOTATE_FUN(reigennxnnf  , "Eigenvalues/eigenvectors (block order real floating-point data)");
ANNOTATE_FUN(reigen16x16nf, "Eigenvalues/eigenvectors (block order real floating-point data)");
ANNOTATE_FUN(reigen8x8nf  , "Eigenvalues/eigenvectors (block order real floating-point data)");
ANNOTATE_FUN(reigennxnsf  , "Eigenvalues/eigenvectors (streaming order real floating-point data)");
ANNOTATE_FUN(reigen4x4sf  , "Eigenvalues/eigenvectors (streaming order real floating-point data)");
ANNOTATE_FUN(reigen3x3sf  , "Eigenvalues/eigenvectors (streaming order real floating-point data)");
ANNOTATE_FUN(reigen2x2sf  , "Eigenvalues/eigenvectors (streaming order real floating-point data)");
