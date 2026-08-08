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
 * Communications
 * Annotations
 */
#include "NatureDSP_types.h"
#include "NatureDSP_Baseband_fit.h"
#include "common.h"

ANNOTATE_FUN(qammap_lte_qpsk   ,"LTE QPSK mapper");
ANNOTATE_FUN(qammap_lte_qam16  ,"LTE QAM16 mapper");
ANNOTATE_FUN(qammap_lte_qam64  ,"LTE QAM64 mapper");
ANNOTATE_FUN(softdemap_lte_qpsk  ,"LTE QPSK soft demapper");
ANNOTATE_FUN(softdemap_lte_16qam ,"LTE QAM16 soft demapper");
ANNOTATE_FUN(softdemap_lte_64qam ,"LTE QAM64 soft demapper");
ANNOTATE_FUN(zcseq_gen   , "Zadov-Chu generator (given u,v)");
ANNOTATE_FUN(zcseq_gen2  , "Zadov-Chu generator (u-th root)");
ANNOTATE_FUN(lteprs_reset, "LTE PRS sequence generator reset");
ANNOTATE_FUN(lteprs_gen1 , "LTE PRS sequence (codeword length 1)");
ANNOTATE_FUN(lteprs_gen2 , "LTE PRS sequence (codeword length 2)");
ANNOTATE_FUN(lteprs_gen4 , "LTE PRS sequence (codeword length 4)");
ANNOTATE_FUN(lteprs_gen6 , "LTE PRS sequence (codeword length 6)");
ANNOTATE_FUN(lteprs_gen8 , "LTE PRS sequence (codeword length 8)");
ANNOTATE_FUN(lteprs_gen10, "LTE PRS sequence (codeword length 10)");
ANNOTATE_FUN(lteprs_gen16, "LTE PRS sequence (codeword length 16)");
ANNOTATE_FUN(lteprs_gen1_fast, "LTE PRS sequence (codeword length 1)");
ANNOTATE_FUN(lteprs_gen2_fast, "LTE PRS sequence (codeword length 2)");
ANNOTATE_FUN(lteprs_gen4_fast, "LTE PRS sequence (codeword length 4)");
ANNOTATE_FUN(lteprs_gen6_fast, "LTE PRS sequence (codeword length 6)");
ANNOTATE_FUN(lteprs_gen8_fast, "LTE PRS sequence (codeword length 8)");
ANNOTATE_FUN(lteprs_gen10_fast, "LTE PRS sequence (codeword length 10)");
ANNOTATE_FUN(lteprs_gen16_fast, "LTE PRS sequence (codeword length 16)");
ANNOTATE_FUN(lteprs_scramble1  , "Scrambling with LTE PRS sequence (codeword length 1)");
ANNOTATE_FUN(lteprs_scramble2  , "Scrambling with LTE PRS sequence (codeword length 2)");
ANNOTATE_FUN(lteprs_scramble4  , "Scrambling with LTE PRS sequence (codeword length 4)");
ANNOTATE_FUN(lteprs_scramble6  , "Scrambling with LTE PRS sequence (codeword length 6)");
ANNOTATE_FUN(lteprs_scramble8  , "Scrambling with LTE PRS sequence (codeword length 8)");
ANNOTATE_FUN(lteprs_scramble10 , "Scrambling with LTE PRS sequence (codeword length 10)");
ANNOTATE_FUN(lteprs_scramble16 , "Scrambling with LTE PRS sequence (codeword length 16)");
ANNOTATE_FUN(lteprs_scramble1_fast , "Scrambling with LTE PRS sequence (codeword length 1)");
ANNOTATE_FUN(lteprs_scramble2_fast , "Scrambling with LTE PRS sequence (codeword length 2)");
ANNOTATE_FUN(lteprs_scramble4_fast , "Scrambling with LTE PRS sequence (codeword length 4)");
ANNOTATE_FUN(lteprs_scramble6_fast , "Scrambling with LTE PRS sequence (codeword length 6)");
ANNOTATE_FUN(lteprs_scramble8_fast , "Scrambling with LTE PRS sequence (codeword length 8)");
ANNOTATE_FUN(lteprs_scramble10_fast, "Scrambling with LTE PRS sequence (codeword length 10)");
ANNOTATE_FUN(lteprs_scramble16_fast, "Scrambling with LTE PRS sequence (codeword length 16)");
ANNOTATE_FUN(convenc_process2,"Convolution encoder, rate 2");
ANNOTATE_FUN(convenc_process3,"Convolution encoder, rate 3");
ANNOTATE_FUN(convenc_process4,"Convolution encoder, rate 4");
ANNOTATE_FUN(crc_process, "CRC");
ANNOTATE_FUN(gcode      , "Universal Grey coding/QAM mapping");
ANNOTATE_FUN(sliceh_lte_qpsk  ,"LTE QPSK slicer");
ANNOTATE_FUN(sliceh_lte_16qam ,"LTE QAM16 slicer");
ANNOTATE_FUN(sliceh_lte_64qam ,"LTE QAM64 slicer");
ANNOTATE_FUN(bsegmnt1 ,"Bitstream segmentation, codeword length 1 bits/word");
ANNOTATE_FUN(bsegmnt2 ,"Bitstream segmentation, codeword length 2 bits/word");
ANNOTATE_FUN(bsegmnt3 ,"Bitstream segmentation, codeword length 3 bits/word");
ANNOTATE_FUN(bsegmnt4 ,"Bitstream segmentation, codeword length 4 bits/word");
ANNOTATE_FUN(bsegmnt6 ,"Bitstream segmentation, codeword length 6 bits/word");
ANNOTATE_FUN(bsegmnt8 ,"Bitstream segmentation, codeword length 8 bits/word");
ANNOTATE_FUN(bpack1 ,"Bitstream packing, codeword length 1 bits/word");
ANNOTATE_FUN(bpack2 ,"Bitstream packing, codeword length 2 bits/word");
ANNOTATE_FUN(bpack3 ,"Bitstream packing, codeword length 3 bits/word");
ANNOTATE_FUN(bpack4 ,"Bitstream packing, codeword length 4 bits/word");
ANNOTATE_FUN(bpack6 ,"Bitstream packing, codeword length 6 bits/word");
ANNOTATE_FUN(bpack8 ,"Bitstream packing, codeword length 8 bits/word");
ANNOTATE_FUN(sfd2x2 ,"Alamouti SFD, diversity 2");
ANNOTATE_FUN(sfd4x4 ,"Alamouti SFD, diversity 4");
