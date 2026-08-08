/* ------------------------------------------------------------------------ */
/* Copyright (c) 2016 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
/*          Copyright (C) 2009-2016 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */

/*
	NatureDSP_Baseband library. FFT part.
    Radix-2 forward/inverse FFT on real data, auto scaling
    Twiddle factor tables for real-to-complex spectrum conversion routines
	IntegrIT, 2006-2016
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Twiddle factor table declarations. */
#include "fft_tw.h"
/* Common utility declarations. */
#include "common.h"

const int16_t ALIGN(32) fft_tw_tab_rfft_512[2*256+16] = {
  // N = 512;
  // twd = reshape(1j*exp(-2*pi*1j*1/2*(0:N/4-1)/(N/2)),8,N/32).';
  // twd = [(1-twd)/2,(1+twd)/2].';
  // twd = reshape([real(twd(:))';imag(twd(:))'],4,N/4)';
  // twd = max(-2^15,min(2^15-1,floor(twd*2^15+0.5)));
  // twd = mod(2^16+twd,2^16);
  (int16_t)0x4000,(int16_t)0xc000,(int16_t)0x3f37,(int16_t)0xc001,
  (int16_t)0x3e6e,(int16_t)0xc005,(int16_t)0x3da5,(int16_t)0xc00b,
  (int16_t)0x3cdc,(int16_t)0xc014,(int16_t)0x3c13,(int16_t)0xc01f,
  (int16_t)0x3b4b,(int16_t)0xc02c,(int16_t)0x3a82,(int16_t)0xc03c,
  (int16_t)0x4000,(int16_t)0x4000,(int16_t)0x40c9,(int16_t)0x3fff,
  (int16_t)0x4192,(int16_t)0x3ffb,(int16_t)0x425b,(int16_t)0x3ff5,
  (int16_t)0x4324,(int16_t)0x3fec,(int16_t)0x43ed,(int16_t)0x3fe1,
  (int16_t)0x44b5,(int16_t)0x3fd4,(int16_t)0x457e,(int16_t)0x3fc4,
  (int16_t)0x39ba,(int16_t)0xc04f,(int16_t)0x38f2,(int16_t)0xc064,
  (int16_t)0x382a,(int16_t)0xc07b,(int16_t)0x3763,(int16_t)0xc095,
  (int16_t)0x369c,(int16_t)0xc0b1,(int16_t)0x35d5,(int16_t)0xc0d0,
  (int16_t)0x350f,(int16_t)0xc0f1,(int16_t)0x3449,(int16_t)0xc115,
  (int16_t)0x4646,(int16_t)0x3fb1,(int16_t)0x470e,(int16_t)0x3f9c,
  (int16_t)0x47d6,(int16_t)0x3f85,(int16_t)0x489d,(int16_t)0x3f6b,
  (int16_t)0x4964,(int16_t)0x3f4f,(int16_t)0x4a2b,(int16_t)0x3f30,
  (int16_t)0x4af1,(int16_t)0x3f0f,(int16_t)0x4bb7,(int16_t)0x3eeb,
  (int16_t)0x3384,(int16_t)0xc13b,(int16_t)0x32bf,(int16_t)0xc163,
  (int16_t)0x31fa,(int16_t)0xc18e,(int16_t)0x3136,(int16_t)0xc1bb,
  (int16_t)0x3073,(int16_t)0xc1eb,(int16_t)0x2fb0,(int16_t)0xc21d,
  (int16_t)0x2eee,(int16_t)0xc251,(int16_t)0x2e2d,(int16_t)0xc288,
  (int16_t)0x4c7c,(int16_t)0x3ec5,(int16_t)0x4d41,(int16_t)0x3e9d,
  (int16_t)0x4e06,(int16_t)0x3e72,(int16_t)0x4eca,(int16_t)0x3e45,
  (int16_t)0x4f8d,(int16_t)0x3e15,(int16_t)0x5050,(int16_t)0x3de3,
  (int16_t)0x5112,(int16_t)0x3daf,(int16_t)0x51d3,(int16_t)0x3d78,
  (int16_t)0x2d6c,(int16_t)0xc2c1,(int16_t)0x2cac,(int16_t)0xc2fd,
  (int16_t)0x2bed,(int16_t)0xc33b,(int16_t)0x2b2e,(int16_t)0xc37b,
  (int16_t)0x2a70,(int16_t)0xc3be,(int16_t)0x29b4,(int16_t)0xc403,
  (int16_t)0x28f7,(int16_t)0xc44a,(int16_t)0x283c,(int16_t)0xc493,
  (int16_t)0x5294,(int16_t)0x3d3f,(int16_t)0x5354,(int16_t)0x3d03,
  (int16_t)0x5413,(int16_t)0x3cc5,(int16_t)0x54d2,(int16_t)0x3c85,
  (int16_t)0x5590,(int16_t)0x3c42,(int16_t)0x564c,(int16_t)0x3bfd,
  (int16_t)0x5709,(int16_t)0x3bb6,(int16_t)0x57c4,(int16_t)0x3b6d,
  (int16_t)0x2782,(int16_t)0xc4df,(int16_t)0x26c9,(int16_t)0xc52d,
  (int16_t)0x2611,(int16_t)0xc57e,(int16_t)0x2559,(int16_t)0xc5d0,
  (int16_t)0x24a3,(int16_t)0xc625,(int16_t)0x23ee,(int16_t)0xc67c,
  (int16_t)0x233a,(int16_t)0xc6d5,(int16_t)0x2287,(int16_t)0xc731,
  (int16_t)0x587e,(int16_t)0x3b21,(int16_t)0x5937,(int16_t)0x3ad3,
  (int16_t)0x59ef,(int16_t)0x3a82,(int16_t)0x5aa7,(int16_t)0x3a30,
  (int16_t)0x5b5d,(int16_t)0x39db,(int16_t)0x5c12,(int16_t)0x3984,
  (int16_t)0x5cc6,(int16_t)0x392b,(int16_t)0x5d79,(int16_t)0x38cf,
  (int16_t)0x21d5,(int16_t)0xc78f,(int16_t)0x2124,(int16_t)0xc7ee,
  (int16_t)0x2074,(int16_t)0xc850,(int16_t)0x1fc6,(int16_t)0xc8b5,
  (int16_t)0x1f19,(int16_t)0xc91b,(int16_t)0x1e6d,(int16_t)0xc983,
  (int16_t)0x1dc3,(int16_t)0xc9ee,(int16_t)0x1d19,(int16_t)0xca5b,
  (int16_t)0x5e2b,(int16_t)0x3871,(int16_t)0x5edc,(int16_t)0x3812,
  (int16_t)0x5f8c,(int16_t)0x37b0,(int16_t)0x603a,(int16_t)0x374b,
  (int16_t)0x60e7,(int16_t)0x36e5,(int16_t)0x6193,(int16_t)0x367d,
  (int16_t)0x623d,(int16_t)0x3612,(int16_t)0x62e7,(int16_t)0x35a5,
  (int16_t)0x1c72,(int16_t)0xcac9,(int16_t)0x1bcb,(int16_t)0xcb3a,
  (int16_t)0x1b26,(int16_t)0xcbad,(int16_t)0x1a82,(int16_t)0xcc21,
  (int16_t)0x19e0,(int16_t)0xcc98,(int16_t)0x193f,(int16_t)0xcd11,
  (int16_t)0x18a0,(int16_t)0xcd8c,(int16_t)0x1802,(int16_t)0xce08,
  (int16_t)0x638e,(int16_t)0x3537,(int16_t)0x6435,(int16_t)0x34c6,
  (int16_t)0x64da,(int16_t)0x3453,(int16_t)0x657e,(int16_t)0x33df,
  (int16_t)0x6620,(int16_t)0x3368,(int16_t)0x66c1,(int16_t)0x32ef,
  (int16_t)0x6760,(int16_t)0x3274,(int16_t)0x67fe,(int16_t)0x31f8,
  (int16_t)0x1766,(int16_t)0xce87,(int16_t)0x16cb,(int16_t)0xcf07,
  (int16_t)0x1632,(int16_t)0xcf8a,(int16_t)0x159b,(int16_t)0xd00e,
  (int16_t)0x1505,(int16_t)0xd094,(int16_t)0x1471,(int16_t)0xd11c,
  (int16_t)0x13df,(int16_t)0xd1a6,(int16_t)0x134e,(int16_t)0xd231,
  (int16_t)0x689a,(int16_t)0x3179,(int16_t)0x6935,(int16_t)0x30f9,
  (int16_t)0x69ce,(int16_t)0x3076,(int16_t)0x6a65,(int16_t)0x2ff2,
  (int16_t)0x6afb,(int16_t)0x2f6c,(int16_t)0x6b8f,(int16_t)0x2ee4,
  (int16_t)0x6c21,(int16_t)0x2e5a,(int16_t)0x6cb2,(int16_t)0x2dcf,
  (int16_t)0x12bf,(int16_t)0xd2bf,(int16_t)0x1231,(int16_t)0xd34e,
  (int16_t)0x11a6,(int16_t)0xd3df,(int16_t)0x111c,(int16_t)0xd471,
  (int16_t)0x1094,(int16_t)0xd505,(int16_t)0x100e,(int16_t)0xd59b,
  (int16_t)0x0f8a,(int16_t)0xd632,(int16_t)0x0f07,(int16_t)0xd6cb,
  (int16_t)0x6d41,(int16_t)0x2d41,(int16_t)0x6dcf,(int16_t)0x2cb2,
  (int16_t)0x6e5a,(int16_t)0x2c21,(int16_t)0x6ee4,(int16_t)0x2b8f,
  (int16_t)0x6f6c,(int16_t)0x2afb,(int16_t)0x6ff2,(int16_t)0x2a65,
  (int16_t)0x7076,(int16_t)0x29ce,(int16_t)0x70f9,(int16_t)0x2935,
  (int16_t)0x0e87,(int16_t)0xd766,(int16_t)0x0e08,(int16_t)0xd802,
  (int16_t)0x0d8c,(int16_t)0xd8a0,(int16_t)0x0d11,(int16_t)0xd93f,
  (int16_t)0x0c98,(int16_t)0xd9e0,(int16_t)0x0c21,(int16_t)0xda82,
  (int16_t)0x0bad,(int16_t)0xdb26,(int16_t)0x0b3a,(int16_t)0xdbcb,
  (int16_t)0x7179,(int16_t)0x289a,(int16_t)0x71f8,(int16_t)0x27fe,
  (int16_t)0x7274,(int16_t)0x2760,(int16_t)0x72ef,(int16_t)0x26c1,
  (int16_t)0x7368,(int16_t)0x2620,(int16_t)0x73df,(int16_t)0x257e,
  (int16_t)0x7453,(int16_t)0x24da,(int16_t)0x74c6,(int16_t)0x2435,
  (int16_t)0x0ac9,(int16_t)0xdc72,(int16_t)0x0a5b,(int16_t)0xdd19,
  (int16_t)0x09ee,(int16_t)0xddc3,(int16_t)0x0983,(int16_t)0xde6d,
  (int16_t)0x091b,(int16_t)0xdf19,(int16_t)0x08b5,(int16_t)0xdfc6,
  (int16_t)0x0850,(int16_t)0xe074,(int16_t)0x07ee,(int16_t)0xe124,
  (int16_t)0x7537,(int16_t)0x238e,(int16_t)0x75a5,(int16_t)0x22e7,
  (int16_t)0x7612,(int16_t)0x223d,(int16_t)0x767d,(int16_t)0x2193,
  (int16_t)0x76e5,(int16_t)0x20e7,(int16_t)0x774b,(int16_t)0x203a,
  (int16_t)0x77b0,(int16_t)0x1f8c,(int16_t)0x7812,(int16_t)0x1edc,
  (int16_t)0x078f,(int16_t)0xe1d5,(int16_t)0x0731,(int16_t)0xe287,
  (int16_t)0x06d5,(int16_t)0xe33a,(int16_t)0x067c,(int16_t)0xe3ee,
  (int16_t)0x0625,(int16_t)0xe4a3,(int16_t)0x05d0,(int16_t)0xe559,
  (int16_t)0x057e,(int16_t)0xe611,(int16_t)0x052d,(int16_t)0xe6c9,
  (int16_t)0x7871,(int16_t)0x1e2b,(int16_t)0x78cf,(int16_t)0x1d79,
  (int16_t)0x792b,(int16_t)0x1cc6,(int16_t)0x7984,(int16_t)0x1c12,
  (int16_t)0x79db,(int16_t)0x1b5d,(int16_t)0x7a30,(int16_t)0x1aa7,
  (int16_t)0x7a82,(int16_t)0x19ef,(int16_t)0x7ad3,(int16_t)0x1937,
  (int16_t)0x04df,(int16_t)0xe782,(int16_t)0x0493,(int16_t)0xe83c,
  (int16_t)0x044a,(int16_t)0xe8f7,(int16_t)0x0403,(int16_t)0xe9b4,
  (int16_t)0x03be,(int16_t)0xea70,(int16_t)0x037b,(int16_t)0xeb2e,
  (int16_t)0x033b,(int16_t)0xebed,(int16_t)0x02fd,(int16_t)0xecac,
  (int16_t)0x7b21,(int16_t)0x187e,(int16_t)0x7b6d,(int16_t)0x17c4,
  (int16_t)0x7bb6,(int16_t)0x1709,(int16_t)0x7bfd,(int16_t)0x164c,
  (int16_t)0x7c42,(int16_t)0x1590,(int16_t)0x7c85,(int16_t)0x14d2,
  (int16_t)0x7cc5,(int16_t)0x1413,(int16_t)0x7d03,(int16_t)0x1354,
  (int16_t)0x02c1,(int16_t)0xed6c,(int16_t)0x0288,(int16_t)0xee2d,
  (int16_t)0x0251,(int16_t)0xeeee,(int16_t)0x021d,(int16_t)0xefb0,
  (int16_t)0x01eb,(int16_t)0xf073,(int16_t)0x01bb,(int16_t)0xf136,
  (int16_t)0x018e,(int16_t)0xf1fa,(int16_t)0x0163,(int16_t)0xf2bf,
  (int16_t)0x7d3f,(int16_t)0x1294,(int16_t)0x7d78,(int16_t)0x11d3,
  (int16_t)0x7daf,(int16_t)0x1112,(int16_t)0x7de3,(int16_t)0x1050,
  (int16_t)0x7e15,(int16_t)0x0f8d,(int16_t)0x7e45,(int16_t)0x0eca,
  (int16_t)0x7e72,(int16_t)0x0e06,(int16_t)0x7e9d,(int16_t)0x0d41,
  (int16_t)0x013b,(int16_t)0xf384,(int16_t)0x0115,(int16_t)0xf449,
  (int16_t)0x00f1,(int16_t)0xf50f,(int16_t)0x00d0,(int16_t)0xf5d5,
  (int16_t)0x00b1,(int16_t)0xf69c,(int16_t)0x0095,(int16_t)0xf763,
  (int16_t)0x007b,(int16_t)0xf82a,(int16_t)0x0064,(int16_t)0xf8f2,
  (int16_t)0x7ec5,(int16_t)0x0c7c,(int16_t)0x7eeb,(int16_t)0x0bb7,
  (int16_t)0x7f0f,(int16_t)0x0af1,(int16_t)0x7f30,(int16_t)0x0a2b,
  (int16_t)0x7f4f,(int16_t)0x0964,(int16_t)0x7f6b,(int16_t)0x089d,
  (int16_t)0x7f85,(int16_t)0x07d6,(int16_t)0x7f9c,(int16_t)0x070e,
  (int16_t)0x004f,(int16_t)0xf9ba,(int16_t)0x003c,(int16_t)0xfa82,
  (int16_t)0x002c,(int16_t)0xfb4b,(int16_t)0x001f,(int16_t)0xfc13,
  (int16_t)0x0014,(int16_t)0xfcdc,(int16_t)0x000b,(int16_t)0xfda5,
  (int16_t)0x0005,(int16_t)0xfe6e,(int16_t)0x0001,(int16_t)0xff37,
  (int16_t)0x7fb1,(int16_t)0x0646,(int16_t)0x7fc4,(int16_t)0x057e,
  (int16_t)0x7fd4,(int16_t)0x04b5,(int16_t)0x7fe1,(int16_t)0x03ed,
  (int16_t)0x7fec,(int16_t)0x0324,(int16_t)0x7ff5,(int16_t)0x025b,
  (int16_t)0x7ffb,(int16_t)0x0192,(int16_t)0x7fff,(int16_t)0x00c9,
  // Shift-and-reverse pattern for select.
  16,17,14,15,12,13,10,11,8,9,6,7,4,5,2,3
};
