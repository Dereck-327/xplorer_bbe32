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
 */

#ifndef __NATUREDSP_BASEBAND_COMM_H
#define __NATUREDSP_BASEBAND_COMM_H

#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
Communications
qammap      QAM Mapper: QPSK, QAM16, QAM64
zcseq_gen   Zadoff-Chu Sequence
softdemap   Soft Demapper for QPSK and QAM Constellations
convenc     Convolution Encoder
lteprs      LTE Pseudo-random Sequence Generation and Scrambling
crc         Cyclic redundancy check
gcode       Universal Grey Coding/QAM Mapping
sliceh      Hard demapper for QPSK and QAM constellations
bsegmnt     Bit segmentation, bitstream of 1, 2, 3, 4, 6, 8 bits to 16 bits.   
bpack       Bit packing, formatted bitstream to raw bitstream.   
sfd         Alamouti SFD
===========================================================================*/

/*-------------------------------------------------------------------------
  QAM mapper (QPSK, QAM16, QAM64).

  Translates a number of binary n-tuples (where n is 2 for QPSK, 4 for QAM16
  or 6 for QAM64) into an array of complex points.

  For detailed information on mapping bits onto constellations, please
  refer to
  [1] ETSI TS 136 211 V10.0.0 (2011-01); LTE; Evolved Universal Terrestrial
      Radio Access (E-UTRA); Physical channels and modulation (3GPP TS
      36.211 version 10.0.0 Release 10);
 
  Methods:
    qammap_lte_qpsk()    -  LTE QPSK  mapper
    qammap_lte_qam16()   -  LTE QAM16 mapper
    qammap_lte_qam64()   -  LTE QAM64 mapper
  Restrictions:
    All mapper variants require a multiple of 16 for the number of 
    points N, and 32 bytes alignment for input/output arrays d, r.
---------------------------------------------------------------------------*/

#define QAMMAP_QPSK_TYPE_A    0
#define QAMMAP_QPSK_TYPE_B    1
#define QAMMAP_QPSK_TYPE_C    2
#define QAMMAP_QPSK_LTE       QAMMAP_QPSK_TYPE_C

/*-------------------------------------------------------------------------
  LTE QPSK mapper

  Translate binary 2-tuples into QPSK points. QPSK LTE constellation: 
                          
                          [10]    [00]
                          -1+j     1+j
                     
                          -1-j     1-j
                          [11]    [01]
  Input:
    d[N]         Binary 2-tuples, each one stored at the least significant bits
                 of a 16-bit word (0..1). 14 MSBs of the word (2..15) must hold
                 zero
  Output:
    r[2*N]       QPSK points, Q5.10. In-phase (I) and quadrature (Q) phase 
                 components are interleaved with the I-component going first
                 (at even indices).
  Restrictions:
    N            Multiple of 16
    d[N],r[2*N]  Aligned on 32-byte boundary
---------------------------------------------------------------------------*/
void qammap_lte_qpsk ( int16_t * restrict r,
                 const int16_t *          d,
                 int N );


/*-------------------------------------------------------------------------
  Fast LTE QAM16 mapper

  Translate binary 4-tuples into LTE QAM16 points.

  Input:
    d[N]     Binary 4-tuples, each one stored at the least significant bits
             of a 16-bit word (0..3). 12 MSBs of the word (4..15) must hold
             zero
  Output:
   r[2*N]    QAM16 points, Q5.10. In-phase (I) and quadrature (Q) phase 
             components are interleaved with the I-component going first
             (at even indices).
  Restrictions:
    N        Multiple of 16
    d[],r[]  Aligned on 32-byte boundary
---------------------------------------------------------------------------*/
void qammap_lte_qam16 ( int16_t * restrict r,
                  const int16_t *          d,
                  int N );

/*-------------------------------------------------------------------------
  LTE QAM64 mapper

  Translate binary 6-tuples into LTE QAM64 points.

  Input:
    d[N]     Binary 6-tuples, each one stored at the least significant bits
             of a 16-bit word (0..5). 10 MSBs of the word (6..15) must hold
             zero
  Output:
   r[2*N]    QAM64 points, Q5.10. In-phase (I) and quadrature (Q) phase 
             components are interleaved with the I-component going first
             (at even indices).
  Restrictions:
    N        Multiple of 16
    d[],r[]  Aligned on 32-byte boundary
---------------------------------------------------------------------------*/
void qammap_lte_qam64 ( int16_t * restrict r,
                  const int16_t *          d,
                  int N );


/*-------------------------------------------------------------------------
Soft demapper

  Soft demapper algorithm converts equalized complex data into the sequence
  of log-likelihood ratios (LLRs). It supports QAM/QPSK constellations from
  LTE standards.

  For detailed information on mapping bits onto constellations, please refer
  to
  [1] ETSI TS 136 211 V10.0.0 (2011-01); LTE; Evolved Universal Terrestrial
      Radio Access (E-UTRA); Physical channels and modulation (3GPP TS
      36.211 version 10.0.0 Release 10);

  Soft decoding algorithm is run over a block of I/Q samples. For each input
  sample s, the demapper computes approximations of sub-optimal simplified
  log-likelihood ratios (LLR) for K bit positions, where K equals the base-2 
  logarithm of the constellation size. Let S(0,k,s) denote the constellation
  point with 0 at k-th bit position which is closest to s, while S(1,k,s) is
  the closest point with 1 at k-th bit position, k=0..K-1. Then the value of
  sub-optimal simplified LLR is defined as:

    LLR(s,k) = (|s-S(0,k,s)|^2 - |s-S(1,k,s)|^2),

          |...|^2 - squared absolute value of a complex number

  I/Q components are represented in Q10

  Methods:
    softdemap_lte_qpsk()       - LTE QPSK [1] 7.1.2
    softdemap_lte_16qam()      - LTE 16QAM [1] 7.1.3
    softdemap_lte_64qam()      - LTE 64QAM [1] 7.1.4
  Input:
     s[2*N]     Input complex samples. I and Q components are interleaved
                with the real part going first (at even indices). Fixed point
                position depends on the constellation size, as described
                above
  Output:
     llr[K*N]   LLR estimations, Q0. K equals the base-2 logarithm of the
                constellation size. llr[0] holds the LLR value for the left-most
                bit position of a K-tuple that matches s[0]+j*s[1].
  Restrictions:
    s[], llr[]  Must be aligned on 32-byte boundary
    N           Must be a multiple of 8
  Performance restrictions:
    None
---------------------------------------------------------------------------*/
/* LTE QPSK [2] 7.1.2 */
void softdemap_lte_qpsk ( int16_t * restrict llr,
                    const int16_t *          s,
                    int N );
/* LTE 16QAM [2] 7.1.3 */
void softdemap_lte_16qam ( int16_t * restrict llr,
                     const int16_t *          s,
                     int N );
/* LTE 64QAM [2] 7.1.4 */
void softdemap_lte_64qam ( int16_t * restrict llr,
                     const int16_t *          s,
                     int N );

/*-------------------------------------------------------------------------
Zadov-Chu sequence generator

Function generates a complex exponential with some properties. Implements
the base algorithm from LTE standard see para 5.5.1.1, 3GPP TS 36.211
V8.8.0 (2009-09)

Accuracy: 9 (2.7e-4)

Input:
  u       Group number, 0..29
  v       Base sequence number within the group, 0..1
  M       Size of sequence (36..1320 in steps of 12)

Output:
  r[2*M]  Output complex Zadoff-Chu sequence, Q15

Return value:
  none

Restrictions:
  r[]     Must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
void zcseq_gen ( int16_t * restrict r, int u, int v, int M );

/*-------------------------------------------------------------------------
  zcseq_gen2()

  Operation : This function generates u'th root Zadoff-Chu(ZC) sequence
              as defined in section 5.7.2 of 3GPP 36.211 V8.4.0(2008-09)
              xu(n) = e^(-j*pi*u*n*(n+1)/Nzc)
              where 0 <= n <= Nzc-1

  Accuracy: 44 (1.3e-3)

  Output:
    r    Pointer to output complex ZC sequence.
         Real and imaginary outputs are interleaved.
         Output is in Q15 format.
  Input:
    u    Root no. to be used for ZC sequence generation
         Range should be (1 <= u <= Nzc-1)
    Nzc  Length of ZC sequence (839 or 139)

  Restrictions:
    r[]  Must be aligned on 32-byte boundary
-------------------------------------------------------------------------*/
void zcseq_gen2 ( int16_t * restrict r, const uint16_t u, uint16_t Nzc ); 

/*-------------------------------------------------------------------------
LTE Pseudo-random Sequence Generation

Functions generate pseudo-random sequence defined in the LTE standard (para
7.2) with given length and codeword size. Supported sizes are 1, 2, 4, 6, 8
and 10 bits. Scrambler function also makes bitwise modulo 2 summation
between generating sequence and input codewords.
-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------
Reset LTE PRS generator

Function lteprs_reset() prepares shift registers and fills them with proper
values derived from initial values (1 for first register and cinit for second
one).

Algorithm:
  See para 7.2 of 3GPP TS 36.211 V8.8.0 (2009-09)

Output:
  r[2]   32-bit LFSR states

Input:
  cinit  Initialization value for second shift register
  M      Codeword size (1, 2, 4, 6, 8 10 and 16 bits)

Return value:
  none
-------------------------------------------------------------------------*/
void lteprs_reset ( uint32_t * r, uint32_t cinit, int16_t M );

/*-------------------------------------------------------------------------
Generate LTE PRS sequence 

Functions generate pseudo-random sequence defined in the LTE standard with
given length and codeword size. Supported sizes are 1, 2, 4, 6, 8 10 and 16 
bits.
Two versions of routines are available: regular versions (lteprs_gen1, 
lteprs_gen2, lteprs_gen4, lteprs_gen6, lteprs_gen8, lteprs_gen10, lteprs_gen16) 
work with arbitrary arguments, faster versions (lteprs_gen1_fast, 
lteprs_gen2_fast, lteprs_gen4_fast, lteprs_gen6_fast, lteprs_gen8_fast,
lteprs_gen10_fast, lteprs_gen16_fast) apply some restrictions.

Algorithm:
  See para 7.2 of 3GPP TS 36.211 V8.8.0 (2009-09)

Input/Output:
  r[2]   32-bit LFSR states

Output:
  c[N]   Generated sequence

Return value:
  none

Restrictions:
  For fast versions:
  N    Multiple of 16
  c[]  Aligned on 32-byte boundary
-------------------------------------------------------------------------*/
void lteprs_gen1  ( uint32_t * r, uint16_t * restrict c, int N );
void lteprs_gen2  ( uint32_t * r, uint16_t * restrict c, int N );
void lteprs_gen4  ( uint32_t * r, uint16_t * restrict c, int N );
void lteprs_gen6  ( uint32_t * r, uint16_t * restrict c, int N );
void lteprs_gen8  ( uint32_t * r, uint16_t * restrict c, int N );
void lteprs_gen10 ( uint32_t * r, uint16_t * restrict c, int N );
void lteprs_gen16 ( uint32_t * r, uint16_t * restrict c, int N );

void lteprs_gen1_fast  ( uint32_t * r, uint16_t * restrict c, int N );
void lteprs_gen2_fast  ( uint32_t * r, uint16_t * restrict c, int N );
void lteprs_gen4_fast  ( uint32_t * r, uint16_t * restrict c, int N );
void lteprs_gen6_fast  ( uint32_t * r, uint16_t * restrict c, int N );
void lteprs_gen8_fast  ( uint32_t * r, uint16_t * restrict c, int N );
void lteprs_gen10_fast ( uint32_t * r, uint16_t * restrict c, int N );
void lteprs_gen16_fast ( uint32_t * r, uint16_t * restrict c, int N );

/*-------------------------------------------------------------------------
Scramble by PRS

Functions generate pseudo-random sequence defined in the LTE standard with
given length and codeword size, when, make bitwise modulo 2 summation
between generated sequence and input codewords. Supported sizes are 1, 2,
4, 6, 8 10 and 16 bits. 
Two versions of routines are available: regular versions (lteprs_scramble1, 
lteprs_scramble2, lteprs_scramble4, lteprs_scramble6, lteprs_scramble8,
lteprs_scramble10,lteprs_scramble16) work with arbitrary arguments, faster versions 
(lteprs_scramble1_fast, lteprs_scramble2_fast, lteprs_scramble4_fast, 
lteprs_scramble6_fast, lteprs_scramble8_fast, lteprs_scramble10_fast,
lteprs_scramble16_fast) 
apply some restrictions.

Algorithm:
  See para 7.2 of 3GPP TS 36.211 V8.8.0 (2009-09)

Input/Output:
  r[2]   32-bit LFSR states

Output:
  c[N]   Output codewords

Input:
  b[N]   Input codewords

Return value:
  none

Restrictions:
  For fast versions: 
  N        Multiple of 16
  c[],b[]  Aligned by 32-byte boundary
-------------------------------------------------------------------------*/

void lteprs_scramble1  ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );
void lteprs_scramble2  ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );
void lteprs_scramble4  ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );
void lteprs_scramble6  ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );
void lteprs_scramble8  ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );
void lteprs_scramble10 ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );
void lteprs_scramble16 ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );

void lteprs_scramble1_fast  ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );
void lteprs_scramble2_fast  ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );
void lteprs_scramble4_fast  ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );
void lteprs_scramble6_fast  ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );
void lteprs_scramble8_fast  ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );
void lteprs_scramble10_fast ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );
void lteprs_scramble16_fast ( uint32_t * r, uint16_t * restrict c, const uint16_t * b, int N );

/*---------------------------------------------------------------------------
Convolution encoder
Description: Allows to encode data with arbitrary selected polynomials with constraint length 
from 3 to 9 and the code rate up to 4.
---------------------------------------------------------------------------*/

typedef void* convenc_handle_t;

/*---------------------------------------------------------------------------
Object allocation
Input:
none
Output:
none
Returns: size of memory in bytes to be allocated
---------------------------------------------------------------------------*/
size_t convenc_alloc();

/*---------------------------------------------------------------------------
Object initialization

Input:
objmem  Allocated memory block
K       Constraint length (3-9)
R       Code rate (2-4)
poly[R] Polynomials in octal form
Output:
None
Returns: handle to the object
---------------------------------------------------------------------------*/
convenc_handle_t convenc_init ( void * objmem, 
                                int K, int R,
                                const int16_t * restrict poly );

/*---------------------------------------------------------------------------
Encoding

Input:
s[N]            Unformatted bitstream (16 bits per word)
N               Number of input 16-bit words
Output:
e[N*R]          Unformatted bitstream
Returned value  Updated state of encoder
---------------------------------------------------------------------------*/
int16_t convenc_process2 ( convenc_handle_t handle, 
                           int16_t * restrict e, const int16_t * s, int N );
int16_t convenc_process3 ( convenc_handle_t handle, 
                           int16_t * restrict e, const int16_t * s, int N );
int16_t convenc_process4 ( convenc_handle_t handle, 
                           int16_t * restrict e, const int16_t * s, int N );

/*-------------------------------------------------------------------------
CRC 

Description: This function calculates 8-, 16-, 24- or 32-bit CRC checksum

Parameters:
Input:
objmem          Allocated memory block
bitstream[N]    Bitstream
order           Order of polynomial (8,16,24 or 32)
poly            Polynomial without higher degree, for example: 
                x^16 + x^12 + x^5 + 1 corresponds to 0x1021
reg             Original crc register value

Returned value  Updated crc register 

Restrictions:
bistream must be aligned on 2-byte boundary
-------------------------------------------------------------------------*/

typedef void* crc_handle_t ;
size_t crc_alloc ( int order, uint32_t poly ); // Returns: size of memory in bytes to be allocated
crc_handle_t crc_init ( void * objmem, int order, uint32_t poly ); // Returns: handle to the object
uint32_t crc_process ( crc_handle_t handle, uint32_t reg, const uint8_t * bitstream, int N ); // Returns: crc value

/*-------------------------------------------------------------------------
  Universal Grey coding/QAM mapping

  Map data n-tuples (n=1..8) onto a user-defined QAM<2^n> constellation.

  User provides a constellation by passing an array S of 16-bit complex numbers,
  where the correspondence between a QAM point s and its related data item m
  is established through the array index: s_I = S[2*m+0] (in-phase component), 
  s_Q = S[2*m+1] (quadrature phase component).

  Input:
  S[2*(2^n)]  Complex constellation
  d[N]        Input n-tuples

  Output:
  r[2*N]      Coded complex points

  Scratch:
  pScr        Scratch memory for storing temporary data. Allocated space must 
              be at least GCODE_SCRATCH() bytes

  Restrictions:
  N           Multiple of 16. 
  pScr,r,d    Aligned on 32-byte boundary.
---------------------------------------------------------------------------*/
void gcode ( void * pScr,
             int16_t * restrict r,
       const int16_t *          d,
       const int16_t *          S,
       int N );

#define GCODE_SCRATCH (0)

/*-------------------------------------------------------------------------
Hard demapper

Hard demapper algorithm converts equalized complex data into the sequence
decoded bits. It supports QAM/QPSK constellations from LTE standards.

For detailed information on mapping bits onto constellations, please refer
to
[1] ETSI TS 136 211 V10.0.0 (2011-01); LTE; Evolved Universal Terrestrial
    Radio Access (E-UTRA); Physical channels and modulation (3GPP TS
    36.211 version 10.0.0 Release 10);

Methods:
  sliceh_lte_qpsk()       - LTE QPSK [1] 7.1.2
  sliceh_lte_16qam()      - LTE 16QAM [1] 7.1.3
  sliceh_lte_64qam()      - LTE 64QAM [1] 7.1.4
Input:
   s[2*N]  Input complex samples. I and Q components are interleaved
           with the real part going first (at even indices). Fixed point
           position depends on the constellation size, as described
           above
Output:
   b[N]    Bitstream formatted as sequence of 16-bit words each containing 
           K-tuples
Restrictions:
  s[]      Aligned on 32-byte boundary
  b        Aligned on 32-byte boundary
  N        Multiple of 8
Performance restrictions:
  None
---------------------------------------------------------------------------*/
/* LTE QPSK [1] 7.1.2 */
void sliceh_lte_qpsk  ( int16_t * restrict b,
                  const int16_t *          s,
                  int N);
/* LTE 16QAM [1] 7.1.3 */
void sliceh_lte_16qam ( int16_t * restrict b,
                  const int16_t *          s,
                  int N );
/* LTE 64QAM [1] 7.1.4 */
void sliceh_lte_64qam ( int16_t * restrict b,
                  const int16_t *          s,
                  int N );

/*-------------------------------------------------------------------------
Bit segmentation

Convert stream of 1, 2, 4, 6, 8 bits to stream of 16-bit words

Input:
b[N*M/8]  Unformatted bitstream
N         Number of symbols
M         Number of bits per symbol
Output:
s[N]      Bitstream formatted as M bits per word

Restrictions:
M   1, 2, 3, 4, 6 or 8
N   Multiple of 16
s   Aligned on 32-byte boundary
b   Aligned on 32-byte boundary
-------------------------------------------------------------------------*/
void bsegmnt1 ( int16_t * restrict s, const uint8_t* b, int N );
void bsegmnt2 ( int16_t * restrict s, const uint8_t* b, int N );
void bsegmnt3 ( int16_t * restrict s, const uint8_t* b, int N );
void bsegmnt4 ( int16_t * restrict s, const uint8_t* b, int N );
void bsegmnt6 ( int16_t * restrict s, const uint8_t* b, int N );
void bsegmnt8 ( int16_t * restrict s, const uint8_t* b, int N );

/*-------------------------------------------------------------------------
Bit packing

Bitstream packing - convert stream of 16-bit words to stream of 1, 2, 3, 4, 6,
8 bits. 

Input:
s[N]      Bitstream formatted as M bits per word
N         Number of symbols
M         Number of bits per symbol
Output:
b[N*M/8]  Unformatted bitstream

Restrictions:
M         1, 2, 3, 4, 6 or 8
N         Multiple of 16
s[]       Aligned on 32-byte boundary
b[]       Aligned on 32-byte boundary
-------------------------------------------------------------------------*/
void bpack1 ( uint8_t * restrict b, const int16_t* s, int N );
void bpack2 ( uint8_t * restrict b, const int16_t* s, int N );
void bpack3 ( uint8_t * restrict b, const int16_t* s, int N );
void bpack4 ( uint8_t * restrict b, const int16_t* s, int N );
void bpack6 ( uint8_t * restrict b, const int16_t* s, int N );
void bpack8 ( uint8_t * restrict b, const int16_t* s, int N );

/*-------------------------------------------------------------------------
Alamouti SFD

Function calculates maximum likelihood decision by Alamouti formulas
for given Tx/Rx diversity. See 
A Simple Transmit Diversity Technique for Wireless Communications. Siavash 
M. Alamouti, IEEE JOURNAL ON SELECT AREAS IN COMMUNICATIONS, VOL. 16, NO. 
8, OCTOBER 1998
LTE - The UMTS Long Term Evolution: From Theory to Practice. Stefania Sesia,
Issam Toufik, Matthew Baker

Input:
r0,r1,r2,r3[2*N]  Input signals, complex vectors of N elements, 
                  fixed point presentation Qr
h[2*2][N][2]      Impulse response, sequence of N 2x2 complex matrices
                  stored in the streaming order, fixed point 
                  presentation Qh
N                 Size of vectors
rsh               Additional right shift amount for output data

Output:
s0,s1[2*N]        Output signals, complex vectors of N elements, 
                  fixed point presentation Qr+Qh-rsh

Restrictions:
Arrays must not overlap
Arrays must be aligned on 32-byte boundary
N must be a multiple of 8
-------------------------------------------------------------------------*/
void sfd2x2 ( int16_t * restrict s0,
              int16_t * restrict s1,
        const int16_t * restrict r0,
        const int16_t * restrict r1,
        const int16_t * restrict r2,
        const int16_t * restrict r3,
        const int16_t * restrict h,
        int N, int rsh );

/*-------------------------------------------------------------------------
SFD for combination of FSTD and SFBC technologies. In fact this function
calculates two independents sfd2x2.

See 
A Simple Transmit Diversity Technique for Wireless Communications. Siavash 
M. Alamouti, IEEE JOURNAL ON SELECT AREAS IN COMMUNICATIONS, VOL. 16, NO. 
8, OCTOBER 1998
LTE - The UMTS Long Term Evolution: From Theory to Practice. Stefania Sesia,
Issam Toufik, Matthew Baker

Input:
r0,r1,r2,r3,
r4,r5,r6,r7[2*N]  Input signals, complex vectors of N elements, 
                  fixed point presentation Qr
h[2*2*2][N][2]    Impulse response, sequence of N 2x2x2 complex matrices
                  stored in the streaming order, fixed point 
                  presentation Qh
N                 Size of vectors
rsh               Additional right shift amount for output data

Output:
s0,s1,s2,s3[2*N]  Output signals, complex vectors of N elements, 
                  fixed point presentation Qr+Qh-rsh

Restrictions:
Arrays must not overlap
Arrays must be aligned on 32-byte boundary
N must be a multiple of 8
-------------------------------------------------------------------------*/
void sfd4x4 ( int16_t * restrict s0, int16_t * restrict s1,
              int16_t * restrict s2, int16_t * restrict s3,
        const int16_t * restrict r0,
        const int16_t * restrict r1,
        const int16_t * restrict r2,
        const int16_t * restrict r3,
        const int16_t * restrict r4,
        const int16_t * restrict r5,
        const int16_t * restrict r6,
        const int16_t * restrict r7,
        const int16_t * restrict h,
        int N, int rsh );

#ifdef __cplusplus
};
#endif

#endif /* __NATUREDSP_BASEBAND_COMM_H */
