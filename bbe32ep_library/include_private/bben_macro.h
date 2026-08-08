#ifndef __XTFFT_BBEN_MACRO_H__
#define __XTFFT_BBEN_MACRO_H__
#include "common.h"
/*
 * Outermost macro definitions follow this naming scheme:
 *  <Rank>_[<Unrolling>_]<DFT_kernel>_[inner-tensor-I_]_<twiddles_buffer_placement>
 *
 *  For example: R3_U2S_DFT4_SKEW_T_XI4_iw_o
 *   This is rank 3 (R3) with outer loop (skew loop) unrolled by two (U2S),
 *   radix 4 (DFT4) with twiddles (SKEW) decimation in freqency (T stands for transposed),
 *   inner loop with 4 iterations (XI4), this loop iterates over vector elements therefore not
 *   explicitly present in the code,
 *   twiddle factors expected in the same DRAM as input (iw_o).
 * 
 */

/*
 * Radix 2
 */

#define DFT2_SKEW_EX(p_src, p_dst, is, os, ivs, ovs, tw1, iinc, oinc) \
{                                                                        \
  VT t0, t1, t2, t3;/*, t4, t5, t6;    */                \
  LD_INC(t0, p_src, is);                                                \
  LD_INC(t1, p_src, -is + ivs + iinc);                                        \
  /* t2 = t0 + t1 + A + B; A==0 && B==0 */                                \
  t2 = BBE_FFTADD4SABN_2XCQ15(t0, t1, 0, 0);                                \
  /* t3 = t0 - t1 + A - B; A==0 && B==0 */                                \
  t3 = BBE_FFTADD4SABN_2XCQ15(t0, t1, 2, 0);                                \
                                                                        \
  MUL(t3, tw1);                                                                \
                                                                        \
  RST_INC(t2, p_dst, os);                                                \
  RST_INC(t3, p_dst, -os + ovs + oinc);                                        \
}

#define DFT2_SKEW(p_src, p_dst, is, os, ivs, ovs, tw1) \
  DFT2_SKEW_EX(p_src, p_dst, is, os, ivs, ovs, tw1, 0, 0)

/* rank 2 skew DFT 2 */

#define R2_DFT2_SKEW_T_XI4_iw_o(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int count = (_count) / 4;                                                \
  int tail = (_count) % 4;                                                \
  int i, is = (_is), ivs = (_ivs), os = (_os),  ovs = (_ovs);                \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)(_p_in);                                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
                                                                        \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  VT tw1, tw2, tw3, tw4;                                                \
      xb_vecNx16 zero = BBE_ZERONX16();                                    \
  MVA(zero); MVB(zero);                            \
                                                                        \
  for (i=0; i<count; i++) {                                                \
    LD_IX4_TW_UNPACK_4(p_tw, tw1, tw2, tw3, tw4);                        \
    DFT2_SKEW(p_src, p_dst, is, os, ivs, ovs, tw1);                        \
    DFT2_SKEW(p_src, p_dst, is, os, ivs, ovs, tw2);                        \
    DFT2_SKEW(p_src, p_dst, is, os, ivs, ovs, tw3);                        \
    DFT2_SKEW(p_src, p_dst, is, os, ivs, ovs, tw4);                        \
  }                                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
  for (i=0; i<tail; i++) {                                                \
    LDA_IX4_TW_UNPACK_1(align, p_tw, tw1);                                \
    DFT2_SKEW(p_src, p_dst, is, os, ivs, ovs, tw1);                        \
  }

#define R2_DFT2_SKEW_T_XI4_i_ow(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int count = (_count) / 4;                                                \
  int tail = (_count) % 4;                                                \
  int i, is = (_is), ivs = (_ivs), os = (_os),  ovs = (_ovs);                \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)(_p_in);                                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
                                                                        \
                                                                        \
  VT tw1, tw2, tw3, tw4;                                                \
  VT zero = BBE_ZERON_2XCQ15();                                    \
  MVA(zero); MVB(zero);                            \
  __Pragma("ymemory(p_dst)");                        \
  __Pragma("ymemory(p_tw)");                        \
                                                                        \
  for (i=0; i<count; i++) {                                                \
    LD_IX4_TW_UNPACK_4(p_tw, tw1, tw2, tw3, tw4);                        \
    DFT2_SKEW(p_src, p_dst, is, os, ivs, ovs, tw1);                        \
    DFT2_SKEW(p_src, p_dst, is, os, ivs, ovs, tw2);                        \
    DFT2_SKEW(p_src, p_dst, is, os, ivs, ovs, tw3);                        \
    DFT2_SKEW(p_src, p_dst, is, os, ivs, ovs, tw4);                        \
  }                                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
  for (i=0; i<tail; i++) {                                                \
    LDA_IX4_TW_UNPACK_1(align, p_tw, tw1);                                \
    DFT2_SKEW(p_src, p_dst, is, os, ivs, ovs, tw1);                        \
  }

/*
 * Radix 3
 */

// r3_tw1 = sign * (0.0, 0.8660254037844385)

#define DFT3_TW_RE   0
#define DFT3_TW_IM   0.8660254037844385

/*
 * DFT3 macro assumes states A and C contain same value input as t0 and
 * state D is set zero.
 */

#define DFT3(t0, t1, t2, r3_tw, scaleid)                                \
  {                                                                        \
    VT  _r3_1, _r3_2;                                                        \
                                                                        \
    BBE_FFTAVGN_2XCQ15SB(t1, t2);                                        \
                                                                        \
    /* t0 = t1 + t2 + C + D; D==0; */                                        \
    t0 = BBE_FFTADD4SCDN_2XCQ15(t1, t2, 0, scaleid);                        \
                                                                        \
    MUL(t1, r3_tw);                                                        \
    MUL(t2, r3_tw);                                                        \
                                                                        \
    /* _r3_1 = t1 - t2 + A - B; */                                        \
    _r3_1 = BBE_FFTADD4SABN_2XCQ15(t1, t2, 2, scaleid);                        \
    /* _r3_2 = t2 - t1 + A - B; */                                        \
    _r3_2 = BBE_FFTADD4SABN_2XCQ15(t2, t1, 2, scaleid);                        \
    t1 = _r3_1;                                                                \
    t2 = _r3_2;                                                                \
  }


#define DFT3_SKEW_EX(p_src, p_dst, is, os, ivs, ovs, r3_tw, tw1, tw2, iinc, oinc, scaleid) \
{                                                                        \
  VT t0, t1, t2;                                                        \
  LD_INC(t2, p_src, -is);                                                \
  LD_INC(t1, p_src, -is);                                                \
  LD_C(p_src, 0);                                                        \
  LD_A(p_src, 2*is + ivs + iinc);                                        \
                                                                        \
  DFT3(t0, t1, t2, r3_tw, scaleid);                                        \
                                                                        \
  MUL(t1, tw1);                                                                \
  MUL(t2, tw2);                                                                \
                                                                        \
  RST_INC(t0, p_dst, os);                                                \
  RST_INC(t1, p_dst, os);                                                \
  RST_INC(t2, p_dst, -2*os + ovs + oinc);                                \
}


#define DFT3_SKEW(p_src, p_dst, is, os, ivs, ovs, r3_tw, tw1, tw2)        \
  DFT3_SKEW_EX(p_src, p_dst, is, os, ivs, ovs, r3_tw, tw1, tw2, 0, 0, 0)


/* rank 2 skew DFT 3 */

#define R2_DFT3_SKEW_T_XI4_iw_o(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int count = (_count) / 2;                                                \
  int tail = (_count) % 2;                                                \
  int i, is = (_is), ivs = (_ivs), os = (_os),  ovs = (_ovs);                \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)((char *)(_p_in) + 2*is);                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
  xb_vecNx16 zero = BBE_ZERONX16();                                \
                                                                        \
                                                                        \
  VT tw1, tw2, tw4, tw5;                                                \
  VT r3_tw;                                                                \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  RX3_TW(_sign, r3_tw);                                                        \
  BBE_MOVSDV(zero);                            \
                                                                        \
  for (i=0; i<count; i++) {                                                \
    LD_IX4_TW_UNPACK_4(p_tw, tw1, tw2, tw4, tw5);                        \
    DFT3_SKEW(p_src, p_dst, is, os, ivs, ovs, r3_tw, tw1, tw2);                \
    DFT3_SKEW(p_src, p_dst, is, os, ivs, ovs, r3_tw, tw4, tw5);                \
  }                                                                        \
  if (tail) {                                                                \
    align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
    LDA_IX4_TW_UNPACK_2(align, p_tw, tw1, tw2);                                \
    DFT3_SKEW(p_src, p_dst, is, os, ivs, ovs, r3_tw, tw1, tw2);                \
  }

#define R2_DFT3_SKEW_T_XI4_i_ow(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int count = (_count) / 2;                                                \
  int tail = (_count) % 2;                                                \
  int i, is = (_is), ivs = (_ivs), os = (_os),  ovs = (_ovs);                \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)((char *)(_p_in) + 2*is);                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
                                                                        \
  VT tw1, tw2, tw4, tw5;                        \
  VT r3_tw;                                \
    xb_vecNx16 zero = BBE_ZERONX16();                                 \
  __Pragma("ymemory(p_tw)");                        \
  __Pragma("ymemory(p_dst)");                        \
                                    \
  RX3_TW(_sign, r3_tw);                            \
  BBE_MOVSDV(zero);                            \
                                                                        \
  for (i=0; i<count; i++) {                                                \
    LD_IX4_TW_UNPACK_4(p_tw, tw1, tw2, tw4, tw5);                        \
    DFT3_SKEW(p_src, p_dst, is, os, ivs, ovs, r3_tw, tw1, tw2);                \
    DFT3_SKEW(p_src, p_dst, is, os, ivs, ovs, r3_tw, tw4, tw5);                \
  }                                                                        \
  if (tail) {                                                                \
    align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
    LDA_IX4_TW_UNPACK_2(align, p_tw, tw1, tw2);                                \
    DFT3_SKEW(p_src, p_dst, is, os, ivs, ovs, r3_tw, tw1, tw2);                \
  }


/* rank 3 skew DFT 3 */

#define R3_DFT3_SKEW_T_XI4_iw_o(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                                                \
  int count1 = (_count1);                                                \
  int count2 = (_count2);                                                \
  int is = (_is), ivs1 = (_ivs1), ivs2 = (_ivs2);                        \
  int os = (_os), ovs1 = (_ovs1), ovs2 = (_ovs2);                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)((char *)(_p_in) + 2*is);                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
                                                                          \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  VT r3_tw;                                                                \
                                                                        \
  RX3_TW(_sign, r3_tw);                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
                                                                        \
  BBE_MOVSDV(0);                                                        \
                                                                        \
  for (j=0; j<count2; j++) {                                                \
    VT tw1, tw2;                                                        \
    LDA_IX4_TW_UNPACK_2(align, p_tw, tw1, tw2);                                \
    for (i=0; i<count1; i++) {                                                \
      DFT3_SKEW(p_src, p_dst, is, os, ivs1, ovs1, r3_tw, tw1, tw2);        \
    }                                                                        \
    p_src += (-count1*ivs1 + ivs2)/VSIZE;                                \
    p_dst += (-count1*ovs1 + ovs2)/VSIZE;                                \
  }

#define R3_DFT3_SKEW_T_XI4_i_ow(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                                                \
  int count1 = (_count1);                                                \
  int count2 = (_count2);                                                \
  int is = (_is), ivs1 = (_ivs1), ivs2 = (_ivs2);                        \
  int os = (_os), ovs1 = (_ovs1), ovs2 = (_ovs2);                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)((char *)(_p_in) + 2*is);                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
  xb_vecNx16 zero = BBE_ZERONX16();                                     \
                                    \
  VT r3_tw;                                \
  __Pragma("ymemory(p_tw)");                        \
  __Pragma("ymemory(p_dst)");                        \
                                    \
  RX3_TW(_sign, r3_tw);                            \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                    \
                                    \
  BBE_MOVSDV(zero);                            \
                                                                        \
  for (j=0; j<count2; j++) {                                                \
    VT tw1, tw2;                                                        \
    LDA_IX4_TW_UNPACK_2(align, p_tw, tw1, tw2);                                \
    for (i=0; i<count1; i++) {                                                \
      DFT3_SKEW(p_src, p_dst, is, os, ivs1, ovs1, r3_tw, tw1, tw2);        \
    }                                                                        \
    p_src += (-count1*ivs1 + ivs2)/VSIZE;                                \
    p_dst += (-count1*ovs1 + ovs2)/VSIZE;                                \
  }

/* R3_DFT3_SKEW_T_XI4 with outer loop unrolled by 2 */

#define R3_U2S_DFT3_SKEW_T_XI4_iw_o(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                                                \
  int tail = (_count2) % 2;                                                \
  int count1 = (_count1);                                                \
  int count2 = (_count2) / 2;                                                \
  int is = (_is), ivs1 = (_ivs1), ivs2 = (_ivs2);                        \
  int os = (_os), ovs1 = (_ovs1), ovs2 = (_ovs2);                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src0 = (VT *)((char *)(_p_in) + 2*is);                                \
  VT * p_dst0 = (VT *)(_p_out);                                                \
  VT * p_src1 = (VT *)((char *)(_p_in) + 2*is + ivs2);                        \
  VT * p_dst1 = (VT *)((char *)(_p_out) + ovs2);                        \
  valign align;                                                                \
  VT r3_tw;                                \
  xb_vecNx16 zero = BBE_ZERONX16();                                     \
  __Pragma("ymemory(p_dst0)");                        \
  __Pragma("ymemory(p_dst1)");                        \
                                    \
                                    \
  RX3_TW(_sign, r3_tw);                            \
  BBE_MOVSDV(zero);                            \
                                                                        \
  for (j=0; j<count2; j++) {                                                \
    VT tw1, tw2, tw4, tw5;                                                \
    LD_IX4_TW_UNPACK_4(p_tw, tw1, tw2, tw4, tw5);                        \
    __Pragma("super_swp ii=8, unroll=1");                                \
    for (i=0; i<count1; i++) {                                                \
      DFT3_SKEW(p_src0, p_dst0, is, os, ivs1, ovs1, r3_tw, tw1, tw2);        \
      DFT3_SKEW(p_src1, p_dst1, is, os, ivs1, ovs1, r3_tw, tw4, tw5);        \
    }                                                                        \
    p_src0 += (-count1*ivs1 + 2*ivs2)/VSIZE;                                \
    p_dst0 += (-count1*ovs1 + 2*ovs2)/VSIZE;                                \
    p_src1 += (-count1*ivs1 + 2*ivs2)/VSIZE;                                \
    p_dst1 += (-count1*ovs1 + 2*ovs2)/VSIZE;                                \
  }                                                                        \
  if (tail) {                                                                \
    VT tw1, tw2;                                                        \
    align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
    LDA_IX4_TW_UNPACK_2(align, p_tw, tw1, tw2);                                \
    for (i=0; i<count1; i++) {                                                \
      DFT3_SKEW(p_src0, p_dst0, is, os, ivs1, ovs1, r3_tw, tw1, tw2);        \
    }                                                                        \
  }

/* R3_DFT3_SKEW_T_XI4 with outer loop unrolled by 2 */

#define R3_U2S_DFT3_SKEW_T_XI4_i_ow(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                                                \
  int tail = (_count2) % 2;                                                \
  int count1 = (_count1);                                                \
  int count2 = (_count2) / 2;                                                \
  int is = (_is), ivs1 = (_ivs1), ivs2 = (_ivs2);                        \
  int os = (_os), ovs1 = (_ovs1), ovs2 = (_ovs2);                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src0 = (VT *)((char *)(_p_in) + 2*is);                                \
  VT * p_dst0 = (VT *)(_p_out);                                                \
  VT * p_src1 = (VT *)((char *)(_p_in) + 2*is + ivs2);                        \
  VT * p_dst1 = (VT *)((char *)(_p_out) + ovs2);                        \
  valign align;                                                                \
  xb_vecNx16 zero = BBE_ZERONX16();                                \
                                    \
  VT r3_tw;                                \
  __Pragma("ymemory(p_tw)");                        \
  __Pragma("ymemory(p_dst0)");                        \
  __Pragma("ymemory(p_dst1)");                        \
                                    \
  RX3_TW(_sign, r3_tw);                            \
  BBE_MOVSDV(zero);                            \
                                                                        \
  for (j=0; j<count2; j++) {                                                \
    VT tw1, tw2, tw4, tw5;                                                \
    LD_IX4_TW_UNPACK_4(p_tw, tw1, tw2, tw4, tw5);                        \
    __Pragma("super_swp ii=8, unroll=1");                                \
    for (i=0; i<count1; i++) {                                                \
      DFT3_SKEW(p_src0, p_dst0, is, os, ivs1, ovs1, r3_tw, tw1, tw2);        \
      DFT3_SKEW(p_src1, p_dst1, is, os, ivs1, ovs1, r3_tw, tw4, tw5);        \
    }                                                                        \
    p_src0 += (-count1*ivs1 + 2*ivs2)/VSIZE;                                \
    p_dst0 += (-count1*ovs1 + 2*ovs2)/VSIZE;                                \
    p_src1 += (-count1*ivs1 + 2*ivs2)/VSIZE;                                \
    p_dst1 += (-count1*ovs1 + 2*ovs2)/VSIZE;                                \
  }                                                                        \
  if (tail) {                                                                \
    VT tw1, tw2;                                                        \
    align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
    LDA_IX4_TW_UNPACK_2(align, p_tw, tw1, tw2);                                \
    for (i=0; i<count1; i++) {                                                \
      DFT3_SKEW(p_src0, p_dst0, is, os, ivs1, ovs1, r3_tw, tw1, tw2);        \
    }                                                                        \
  }

/* special case of R3_DFT3_SKEW_T_XI4 with inner loop unrolled (_count1 == 4) */

#define R3_U4_DFT3_SKEW_T_XI4_iw_o(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                                                \
  int count2 = (_count2);                                                \
  int is = (_is), ivs1 = (_ivs1), ivs2 = (_ivs2);                        \
  int os = (_os), ovs1 = (_ovs1), ovs2 = (_ovs2);                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)((char *)(_p_in) + 2*is);                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
  VT r3_tw;                                                                \
                                                                        \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  RX3_TW(_sign, r3_tw);                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
                                                                        \
  BBE_MOVSDV(0);                                                        \
                                                                        \
  for (j=0; j<count2; j++) {                                                \
    VT tw1, tw2;                                                        \
    LDA_IX4_TW_UNPACK_2(align, p_tw, tw1, tw2);                                \
    DFT3_SKEW(p_src, p_dst, is, os, ivs1, ovs1, r3_tw, tw1, tw2);        \
    DFT3_SKEW(p_src, p_dst, is, os, ivs1, ovs1, r3_tw, tw1, tw2);        \
    DFT3_SKEW(p_src, p_dst, is, os, ivs1, ovs1, r3_tw, tw1, tw2);        \
    DFT3_SKEW_EX(p_src, p_dst, is, os, ivs1, ovs1, r3_tw, tw1, tw2, (-4*ivs1 + ivs2), (-4*ovs1 + ovs2), 0); \
  }


/* special case of R3_DFT3_SKEW_T_XI4 with inner loop unrolled (_count1 == 4) */

#define R3_U4_DFT3_SKEW_T_XI4_i_ow(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                                                \
  int count2 = (_count2);                                                \
  int is = (_is);                                                        \
  int ivs1 = (_ivs1);                                                        \
  int ivs2 = (_ivs2);                                                        \
  int os = (_os);                                                        \
  int ovs1 = (_ovs1);                                                        \
  int ovs2 = (_ovs2);                                                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)((char *)(_p_in) + 2*is);                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
  VT r3_tw;                                                                \
                                                                        \
  __Pragma("ymemory(p_tw)");                                                \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  RX3_TW(_sign, r3_tw);                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
                                                                        \
  BBE_MOVSDV(0);                                                        \
                                                                        \
  for (j=0; j<count2; j++) {                                                \
    VT tw1, tw2;                                                        \
    LDA_IX4_TW_UNPACK_2(align, p_tw, tw1, tw2);                                \
    DFT3_SKEW(p_src, p_dst, is, os, ivs1, ovs1, r3_tw, tw1, tw2);        \
    DFT3_SKEW(p_src, p_dst, is, os, ivs1, ovs1, r3_tw, tw1, tw2);        \
    DFT3_SKEW(p_src, p_dst, is, os, ivs1, ovs1, r3_tw, tw1, tw2);        \
    DFT3_SKEW_EX(p_src, p_dst, is, os, ivs1, ovs1, r3_tw, tw1, tw2, (-4*ivs1 + ivs2), (-4*ovs1 + ovs2), 0); \
  }

/*
 *  DFT(4)
 */
#ifndef IS_INV_FFT

#define DFT4(y0, y1, y2, y3, x0, x1, S1, S2, scaleid)            \
  y0 = BBE_FFTADD4S ## S1 ## S2 ## N_2XCQ15(x0, x1, 0, scaleid);    \
  y1 = BBE_FFTADD4S ## S1 ## S2 ## N_2XCQ15(x0, x1, 1, scaleid);    \
  y2 = BBE_FFTADD4S ## S1 ## S2 ## N_2XCQ15(x0, x1, 2, scaleid);    \
  y3 = BBE_FFTADD4S ## S1 ## S2 ## N_2XCQ15(x0, x1, 3, scaleid);

#else
#define DFT4(y0, y1, y2, y3, x0, x1, S1, S2, scaleid)            \
  y0 = BBE_FFTADD4S ## S1 ## S2 ## N_2XCQ15(x0, x1, 0, scaleid);    \
  y3 = BBE_FFTADD4S ## S1 ## S2 ## N_2XCQ15(x0, x1, 1, scaleid);    \
  y2 = BBE_FFTADD4S ## S1 ## S2 ## N_2XCQ15(x0, x1, 2, scaleid);    \
  y1 = BBE_FFTADD4S ## S1 ## S2 ## N_2XCQ15(x0, x1, 3, scaleid);

#endif //#ifndef IS_INV_FFT

#define DFT4_SKEW_T_EX(S1, S2, p_src, p_dst, is, os, ivs, ovs, tw1, tw2, tw3, iinc, oinc) \
{                                                                        \
  VT t0, t1, t2, t3, t4, t5;                                                \
  LD_INC(t0, p_src, is);                                                \
  LD_INC(t1, p_src, is);                                                \
  LD_ ## S1(p_src, is);                                                        \
  LD_ ## S2(p_src, -3*is + ivs + iinc);                                        \
                                                                          \
  DFT4(t2, t3, t4, t5, t0, t1, S1, S2, 0);                                \
                                                                          \
  MUL(t3, tw1);                                                                \
  MUL(t4, tw2);                                                                \
  MUL(t5, tw3);                                                                \
                                                                          \
  RST_INC(t2, p_dst, os);                                                \
  RST_INC(t3, p_dst, os);                                                \
  RST_INC(t4, p_dst, os);                                                \
  RST_INC(t5, p_dst, -3*os + ovs + oinc);                                \
}

#define DFT4_SKEW_T(S1, S2, p_src, p_dst, is, os, ivs, ovs, tw1, tw2, tw3) \
  DFT4_SKEW_T_EX(S1, S2, p_src, p_dst, is, os, ivs, ovs, tw1, tw2, tw3, 0, 0)


/* rank 2 skew DFT 4 */

#define R2_DFT4_SKEW_T_XI4_iw_o(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int tail = (_count) % 2;                                                \
  int count = (_count) / 2;                                                \
  int i, is = (_is), ivs = (_ivs), os = (_os), ovs = (_ovs);                \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)(_p_in);                                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
                                                                        \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
                                                                        \
  for (i=0; i<count; i++) {                                                \
    VT tw1, tw2, tw3;                                                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
    DFT4_SKEW_T(A, B, p_src, p_dst, is, os, ivs, ovs, tw1, tw2, tw3);        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
    DFT4_SKEW_T(C, D, p_src, p_dst, is, os, ivs, ovs, tw1, tw2, tw3);        \
  }                                                                        \
  if (tail) {                                                                \
    VT tw1, tw2, tw3;                                                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
    DFT4_SKEW_T(A, B, p_src, p_dst, is, os, ivs, ovs, tw1, tw2, tw3);        \
  }  


#define R2_DFT4_SKEW_T_XI4_i_ow(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int tail = (_count) % 2;                                                \
  int count = (_count) / 2;                                                \
  int i, is = (_is), ivs = (_ivs), os = (_os), ovs = (_ovs);                \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)(_p_in);                                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
                                                                        \
  __Pragma("ymemory(p_tw)");                                                \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
                                                                        \
  for (i=0; i<count; i++) {                                                \
    VT tw1, tw2, tw3;                                                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
    DFT4_SKEW_T(A, B, p_src, p_dst, is, os, ivs, ovs, tw1, tw2, tw3);        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
    DFT4_SKEW_T(C, D, p_src, p_dst, is, os, ivs, ovs, tw1, tw2, tw3);        \
  }                                                                        \
  if (tail) {                                                                \
    VT tw1, tw2, tw3;                                                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
    DFT4_SKEW_T(A, B, p_src, p_dst, is, os, ivs, ovs, tw1, tw2, tw3);        \
  }  


/* rank 3 skew DFT 4 */

#define R3_DFT4_SKEW_T_XI4_iw_o(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                                                \
  int tail = (_count1) % 2;                                                \
  int count1 = (_count1) / 2;                                                \
  int count2 = (_count2);                                                \
  int is = (_is), ivs1 = (_ivs1), ivs2 = (_ivs2);                        \
  int os = (_os), ovs1 = (_ovs1), ovs2 = (_ovs2);                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)(_p_in);                                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
                                                                        \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
                                                                        \
  if (tail == 0) {                                                        \
    for (j=0; j<count2; j++) {                                                \
      VT tw1, tw2, tw3;                                                        \
      LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
      for (i=0; i<count1; i++) {                                        \
        DFT4_SKEW_T(A, B, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
        DFT4_SKEW_T(C, D, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
      }                                                                        \
      p_src += (-2*count1*ivs1 + ivs2)/VSIZE;                                \
      p_dst += (-2*count1*ovs1 + ovs2)/VSIZE;                                \
    }                                                                        \
  } else {                                                                \
    for (j=0; j<count2; j++) {                                                \
      VT tw1, tw2, tw3;                                                        \
      LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
      for (i=0; i<count1; i++) {                                        \
        DFT4_SKEW_T(A, B, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
        DFT4_SKEW_T(C, D, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
      }                                                                        \
      DFT4_SKEW_T_EX(A, B, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3, (-(2*count1+1)*ivs1 + ivs2), (-(2*count1+1)*ovs1 + ovs2)); \
    }                                                                        \
  }

#define R3_DFT4_SKEW_T_XI4_i_ow(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                                                \
  int tail = (_count1) % 2;                                                \
  int count1 = (_count1) / 2;                                                \
  int count2 = (_count2);                                                \
  int is = (_is), ivs1 = (_ivs1), ivs2 = (_ivs2);                        \
  int os = (_os), ovs1 = (_ovs1), ovs2 = (_ovs2);                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)(_p_in);                                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
                                                                        \
  __Pragma("ymemory(p_tw)");                                                \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
                                                                        \
  if (tail == 0) {                                                        \
    for (j=0; j<count2; j++) {                                                \
      VT tw1, tw2, tw3;                                                        \
      LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
      for (i=0; i<count1; i++) {                                        \
        DFT4_SKEW_T(A, B, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
        DFT4_SKEW_T(C, D, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
      }                                                                        \
      p_src += (-2*count1*ivs1 + ivs2)/VSIZE;                                \
      p_dst += (-2*count1*ovs1 + ovs2)/VSIZE;                                \
    }                                                                        \
  } else {                                                                \
    for (j=0; j<count2; j++) {                                                \
      VT tw1, tw2, tw3;                                                        \
      LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
      for (i=0; i<count1; i++) {                                        \
        DFT4_SKEW_T(A, B, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
        DFT4_SKEW_T(C, D, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
      }                                                                        \
      DFT4_SKEW_T_EX(A, B, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3, (-(2*count1+1)*ivs1 + ivs2), (-(2*count1+1)*ovs1 + ovs2)); \
    }                                                                        \
  }

/* rank 3 skew DFT 4 (R3_DFT4_SKEW_T_XI4) with outer loop unrolled by two  */

#define R3_U2S_DFT4_SKEW_T_XI4_iw_o(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                                                \
  int tail = (_count2) % 2;                                                \
  int count1 = (_count1);                                                \
  int count2 = (_count2) / 2;                                                \
  int is = (_is), ivs1 = (_ivs1), ivs2 = (_ivs2);                        \
  int os = (_os), ovs1 = (_ovs1), ovs2 = (_ovs2);                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src0 = (VT *)(_p_in);                                                \
  VT * p_dst0 = (VT *)(_p_out);                                                \
  VT * p_src1 = (VT *)((char *)(_p_in) + ivs2);                                \
  VT * p_dst1 = (VT *)((char *)(_p_out) + ovs2);                        \
  valign align;                                                                \
                                                                        \
  __Pragma("ymemory(p_dst0)");                                                \
  __Pragma("ymemory(p_dst1)");                                                \
                                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
                                                                        \
  for (j=0; j<count2; j++) {                                                \
    VT tw1, tw2, tw3, tw5, tw6, tw7;                                        \
                                                                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw5, tw6, tw7);                        \
                                                                            \
    for (i=0; i<count1; i++) {                                                \
      DFT4_SKEW_T(A, B, p_src0, p_dst0, is, os, ivs1, ovs1, tw1, tw2, tw3); \
      DFT4_SKEW_T(C, D, p_src1, p_dst1, is, os, ivs1, ovs1, tw5, tw6, tw7); \
    }                                                                        \
    p_src0 += (-count1*ivs1 + 2*ivs2)/VSIZE;                                \
    p_dst0 += (-count1*ovs1 + 2*ovs2)/VSIZE;                                \
    p_src1 += (-count1*ivs1 + 2*ivs2)/VSIZE;                                \
    p_dst1 += (-count1*ovs1 + 2*ovs2)/VSIZE;                                \
  }                                                                        \
  if (tail) {                                                                \
    VT tw1, tw2, tw3;                                                        \
                                                                        \
    tail = count1 % 2;                                                        \
    count1 /= 2;                                                        \
                                                                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
                                                                            \
    for (i=0; i<count1; i++) {                                                \
      DFT4_SKEW_T(A, B, p_src0, p_dst0, is, os, ivs1, ovs1, tw1, tw2, tw3); \
      DFT4_SKEW_T(C, D, p_src0, p_dst0, is, os, ivs1, ovs1, tw1, tw2, tw3); \
    }                                                                        \
    if (tail) {                                                                \
      DFT4_SKEW_T(A, B, p_src0, p_dst0, is, os, ivs1, ovs1, tw1, tw2, tw3); \
    }                                                                        \
  }

#define R3_U2S_DFT4_SKEW_T_XI4_i_ow(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                                                \
  int tail = (_count2) % 2;                                                \
  int count1 = (_count1);                                                \
  int count2 = (_count2) / 2;                                                \
  int is = (_is), ivs1 = (_ivs1), ivs2 = (_ivs2);                        \
  int os = (_os), ovs1 = (_ovs1), ovs2 = (_ovs2);                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src0 = (VT *)(_p_in);                                                \
  VT * p_dst0 = (VT *)(_p_out);                                                \
  VT * p_src1 = (VT *)((char *)(_p_in) + ivs2);                                \
  VT * p_dst1 = (VT *)((char *)(_p_out) + ovs2);                        \
  valign align;                                                                \
                                                                        \
  __Pragma("ymemory(p_tw)");                                                \
  __Pragma("ymemory(p_dst0)");                                                \
  __Pragma("ymemory(p_dst1)");                                                \
                                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
                                                                        \
  for (j=0; j<count2; j++) {                                                \
    VT tw1, tw2, tw3, tw5, tw6, tw7;                                        \
                                                                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw5, tw6, tw7);                        \
                                                                            \
    for (i=0; i<count1; i++) {                                                \
      DFT4_SKEW_T(A, B, p_src0, p_dst0, is, os, ivs1, ovs1, tw1, tw2, tw3); \
      DFT4_SKEW_T(C, D, p_src1, p_dst1, is, os, ivs1, ovs1, tw5, tw6, tw7); \
    }                                                                        \
    p_src0 += (-count1*ivs1 + 2*ivs2)/VSIZE;                                \
    p_dst0 += (-count1*ovs1 + 2*ovs2)/VSIZE;                                \
    p_src1 += (-count1*ivs1 + 2*ivs2)/VSIZE;                                \
    p_dst1 += (-count1*ovs1 + 2*ovs2)/VSIZE;                                \
  }                                                                        \
  if (tail) {                                                                \
    VT tw1, tw2, tw3;                                                        \
                                                                        \
    tail = count1 % 2;                                                        \
    count1 /= 2;                                                        \
                                                                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
                                                                            \
    for (i=0; i<count1; i++) {                                                \
      DFT4_SKEW_T(A, B, p_src0, p_dst0, is, os, ivs1, ovs1, tw1, tw2, tw3); \
      DFT4_SKEW_T(C, D, p_src0, p_dst0, is, os, ivs1, ovs1, tw1, tw2, tw3); \
    }                                                                        \
    if (tail) {                                                                \
      DFT4_SKEW_T(A, B, p_src0, p_dst0, is, os, ivs1, ovs1, tw1, tw2, tw3); \
    }                                                                        \
  }

/* rank 3 skew DFT 4 (R3_DFT4_SKEW_T_XI4) with inner loop unrolled (_count1 == 4)  */

#define R3_U4_DFT4_SKEW_T_XI4_iw_o(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int /*i,*/ j;                                \
  int count2 = (_count2);                                                \
  int is = (_is);                                                        \
  int ivs1 = (_ivs1);                                                        \
  int ivs2 = (_ivs2);                                                        \
  int os = (_os);                                                        \
  int ovs1 = (_ovs1);                                                        \
  int ovs2 = (_ovs2);                                                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)(_p_in);                                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
                                                                        \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
                                                                        \
  for (j=0; j<count2; j++) {                                                \
    VT tw1, tw2, tw3;                                                        \
                                                                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
                                                                            \
    DFT4_SKEW_T(A, B, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
    DFT4_SKEW_T(C, D, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
    DFT4_SKEW_T(A, B, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
    DFT4_SKEW_T_EX(C, D, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3, (-4*ivs1 + ivs2), (-4*ovs1 + ovs2)); \
  }

#define R3_U4_DFT4_SKEW_T_XI4_i_ow(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                                                \
  int count2 = (_count2);                                                \
  int is = (_is);                                                        \
  int ivs1 = (_ivs1);                                                        \
  int ivs2 = (_ivs2);                                                        \
  int os = (_os);                                                        \
  int ovs1 = (_ovs1);                                                        \
  int ovs2 = (_ovs2);                                                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)(_p_in);                                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
                                                                        \
  __Pragma("ymemory(p_tw)");                                                \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
                                                                        \
  for (j=0; j<count2; j++) {                                                \
    VT  tw1, tw2, tw3;                                                        \
                                                                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
                                                                            \
    DFT4_SKEW_T(A, B, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
    DFT4_SKEW_T(C, D, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
    DFT4_SKEW_T(A, B, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3); \
    DFT4_SKEW_T_EX(C, D, p_src, p_dst, is, os, ivs1, ovs1, tw1, tw2, tw3, (-4*ivs1 + ivs2), (-4*ovs1 + ovs2)); \
  }

/*
 * Radix 5 
 */

// r5_tw1 = sign * (-0.29389262614623656, 0.47552825814757676)
// r5_tw2 = (0.55901699437494751, 0);

#define DFT5_TW1_RE (-0.29389262614623656)
#define DFT5_TW1_IM ( 0.47552825814757676)
#define DFT5_TW2_RE ( 0.55901699437494751)
#define DFT5_TW2_IM ( 0)

#define DFT5_SKEW_T_EX(S1, S2, S3, S4, p_src, p_dst, is, os, ivs, ovs, r5_tw1, r5_tw2, r5_tw3, tw1, tw2, tw3, tw4, iinc, oinc) \
{                                                                        \
  VT t0, t1, t2, t3, t4, t5, t6, t7;                                        \
  LD_INC(t0, p_src, is);                                                \
  LD_INC(t1, p_src, is);                                                \
  LD_INC(t2, p_src, is);                                                \
  LD_ ## S2 (p_src, is);                                                \
  LD_ ## S1 (p_src, -4*is + ivs + iinc);                                \
                                                                        \
  t0 = BBE_FFTSRAN_2XCQ15(t0);                                                \
                                                                        \
  /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */                \
  DFT4(t3, t4, t5, t6, t1, t2, S1, S2, 1);                                \
                                                                        \
  t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);                                        \
  {  xb_vecN_2xcq15 tmpc;\
    tmpc = BBE_SRAIN_2XCQ15(t3, 2);   \
    t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            \
    }                                    \
  MUL(t4, r5_tw1);                                                        \
  MUL(t5, r5_tw2);                                                        \
  MUL(t6, r5_tw3);                                                        \
                                                                        \
  MV ## S3 (t5);                                                        \
  MV ## S4 (t6);                                                        \
                                                                        \
  t3 = t4;                                                                \
                                                                        \
  DFT4(t4, t5, t6, t7, t2, t3, S3, S4, 0);                                \
                                                                        \
  MUL(t4, tw1);                                                                \
  MUL(t5, tw2);                                                                \
  MUL(t7, tw3);                                                                \
  MUL(t6, tw4);                                                                \
                                                                        \
  RST_INC(t1, p_dst, os);                                                \
  RST_INC(t4, p_dst, os);                                                \
  RST_INC(t5, p_dst, os);                                                \
  RST_INC(t7, p_dst, os);                                                \
  RST_INC(t6, p_dst, -4*os + ovs + oinc);                                \
}

#define DFT5_SKEW_T(p_src, p_dst, is, os, ivs, ovs, r5_tw1, r5_tw2, r5_tw3, tw1, tw2, tw3, tw4) \
  DFT5_SKEW_T_EX(A, B, C, D, p_src, p_dst, is, os, ivs, ovs, r5_tw1, r5_tw2, r5_tw3, tw1, tw2, tw3, tw4, 0, 0)


#define R2_DFT5_SKEW_T_XI4_iw_o(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int i, count = (_count), is = (_is), ivs = (_ivs), os = (_os), ovs = (_ovs); \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)(_p_in);                                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  VT r5_tw1, r5_tw2, r5_tw3;                                                \
                                                                        \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  RX5_TW(_sign, r5_tw1, r5_tw2, r5_tw3);                                \
                                                                        \
  for (i=0; i<count; i++) {                                                \
    VT tw1, tw2, tw3, tw4;                                                \
    LD_IX4_TW_UNPACK_4(p_tw, tw1, tw2, tw3, tw4);                        \
    DFT5_SKEW_T(p_src, p_dst, is, os, ivs, ovs, r5_tw1, r5_tw2, r5_tw3, tw1, tw2, tw3, tw4); \
  }


#define R2_DFT5_SKEW_T_XI4_i_ow(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int i, count = (_count), is = (_is), ivs = (_ivs), os = (_os), ovs = (_ovs); \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)(_p_in);                                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  VT r5_tw1, r5_tw2, r5_tw3;                                                \
                                                                        \
  __Pragma("ymemory(p_tw)");                                                \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  RX5_TW(_sign, r5_tw1, r5_tw2, r5_tw3);                                \
                                                                        \
  for (i=0; i<count; i++) {                                                \
    VT tw1, tw2, tw3, tw4;                                                \
    LD_IX4_TW_UNPACK_4(p_tw, tw1, tw2, tw3, tw4);                        \
    DFT5_SKEW_T(p_src, p_dst, is, os, ivs, ovs, r5_tw1, r5_tw2, r5_tw3, tw1, tw2, tw3, tw4); \
  }


#define R3_DFT5_SKEW_T_XI4_iw_o(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                                                \
  int count1 = (_count1);                                                \
  int count2 = (_count2);                                                \
  int is = (_is), ivs1 = (_ivs1), ivs2 = (_ivs2);                        \
  int os = (_os), ovs1 = (_ovs1), ovs2 = (_ovs2);                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)(_p_in);                                                \
  VT * p_dst = (VT *)(_p_out);                                                \
                                                                        \
  VT r5_tw1, r5_tw2, r5_tw3;                                                \
                                                                        \
  __Pragma("ymemory(p_dst)");                                                \
                                                                        \
  RX5_TW(_sign, r5_tw1, r5_tw2, r5_tw3);                                \
                                                                        \
  for (j=0; j<count2; j++) {                                                \
    VT tw1, tw2, tw3, tw4;                                                \
    LD_IX4_TW_UNPACK_4(p_tw, tw1, tw2, tw3, tw4);                        \
    for (i=0; i<count1; i++) {                                                \
      DFT5_SKEW_T(p_src, p_dst, is, os, ivs1, ovs1, r5_tw1, r5_tw2, r5_tw3, tw1, tw2, tw3, tw4); \
    }                                                                        \
    p_src += (-count1*ivs1 + ivs2)/VSIZE;                                \
    p_dst += (-count1*ovs1 + ovs2)/VSIZE;                                \
  }

#define R3_DFT5_SKEW_T_XI4_i_ow(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                                                \
  int count1 = (_count1);                                                \
  int count2 = (_count2);                                                \
  int is = (_is), ivs1 = (_ivs1), ivs2 = (_ivs2);                        \
  int os = (_os), ovs1 = (_ovs1), ovs2 = (_ovs2);                        \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)(_p_in);                                                \
  VT * p_dst = (VT *)(_p_out);                                                \
                                                                        \
                                    \
  VT r5_tw1, r5_tw2, r5_tw3;                        \
  __Pragma("ymemory(p_tw)");                        \
  __Pragma("ymemory(p_dst)");                        \
                                                                        \
  RX5_TW(_sign, r5_tw1, r5_tw2, r5_tw3);                                \
                                                                        \
  for (j=0; j<count2; j++) {                                                \
    VT tw1, tw2, tw3, tw4;                                                \
    LD_IX4_TW_UNPACK_4(p_tw, tw1, tw2, tw3, tw4);                        \
    for (i=0; i<count1; i++) {                                                \
      DFT5_SKEW_T(p_src, p_dst, is, os, ivs1, ovs1, r5_tw1, r5_tw2, r5_tw3, tw1, tw2, tw3, tw4); \
    }                                                                        \
    p_src += (-count1*ivs1 + ivs2)/VSIZE;                                \
    p_dst += (-count1*ovs1 + ovs2)/VSIZE;                                \
  }


/* rank 2 skew DFT 6 */

#define DFT6(t0, t1, t2, t3, t4, t5, r3_tw)                                \
  {                                                                        \
    VT _r60 = t0, _r61 = t3, _r62 = t4, _r63 = t1, _r64 = t2, _r65 = t5; \
                                                                        \
    BFLY(_r60, _r61);                                                        \
    BFLY(_r62, _r63);                                                        \
    BFLY(_r64, _r65);                                                        \
                                                                        \
    MVC(_r60);                                                                \
    MVA(_r60);                                                                \
    DFT3(_r60, _r62, _r64, r3_tw, 1);                                        \
                                                                        \
    MVC(_r61);                                                                \
    MVA(_r61);                                                                \
    DFT3(_r61, _r63, _r65, r3_tw, 1);                                        \
                                                                        \
    t0 = _r60; t1 = _r65; t2 = _r62; t3 = _r61; t4 = _r64; t5 = _r63;        \
  }


/* rank 2 skew DFT 6 */
/* xb_vecNx16 input/output */
#define __DFT6(t0, t1, t2, t3, t4, t5, r3_tw)                                \
  {                                                                        \
    VT _r60 = BBE_MOVN_2XCQ15_FROMNX16(t0), _r61 = BBE_MOVN_2XCQ15_FROMNX16(t3), _r62 = BBE_MOVN_2XCQ15_FROMNX16(t4), \
       _r63 = BBE_MOVN_2XCQ15_FROMNX16(t1), _r64 = BBE_MOVN_2XCQ15_FROMNX16(t2), _r65 = BBE_MOVN_2XCQ15_FROMNX16(t5); \
                                                                        \
    BFLY(_r60, _r61);                                                        \
    BFLY(_r62, _r63);                                                        \
    BFLY(_r64, _r65);                                                        \
                                                                        \
    MVC(_r60);                                                                \
    MVA(_r60);                                                                \
    DFT3(_r60, _r62, _r64, r3_tw, 1);                                        \
                                                                        \
    MVC(_r61);                                                                \
    MVA(_r61);                                                                \
    DFT3(_r61, _r63, _r65, r3_tw, 1);                                        \
                                                                        \
    t0 = BBE_MOVNX16_FROMN_2XCQ15(_r60); t1 = BBE_MOVNX16_FROMN_2XCQ15(_r65); t2 = BBE_MOVNX16_FROMN_2XCQ15(_r62); \
    t3 = BBE_MOVNX16_FROMN_2XCQ15(_r61); t4 = BBE_MOVNX16_FROMN_2XCQ15(_r64); t5 = BBE_MOVNX16_FROMN_2XCQ15(_r63);        \
  }

#define R2_DFT6_SKEW_T_XI4_iw_o(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int count = (_count);                                                        \
  int i, is = (_is), ivs = (_ivs), os = (_os), ovs = (_ovs);                \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src = (VT *)(_p_in);                                                \
  VT * p_dst = (VT *)(_p_out);                                                \
  valign align;                                                                \
  xb_vecNx16 zero = BBE_ZERONX16();                                    \
                                    \
  VT r3_tw;                                \
  __Pragma("ymemory(p_dst)");                        \
                                    \
  RX3_TW(_sign, r3_tw);                            \
  BBE_MOVSDV(zero);                            \
                                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
                                                                        \
  for (i=0; i<count; i++) {                                                \
    VT tw1, tw2, tw3, tw4, tw5;                                                \
    VT t0, t1, t2, t3, t4, t5;                                                \
    LD_INC(t0, p_src, is);                                                \
    LD_INC(t1, p_src, is);                                                \
    LD_INC(t2, p_src, is);                                                \
    LD_INC(t3, p_src, is);                                                \
    LD_INC(t4, p_src, is);                                                \
    LD_INC(t5, p_src, -5*is + ivs);                                        \
                                                                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
    LDA_IX4_TW_UNPACK_2(align, p_tw, tw4, tw5);                                \
                                                                        \
    DFT6(t0, t1, t2, t3, t4, t5, r3_tw);                                \
                                                                        \
    MUL(t1, tw1);                                                        \
    MUL(t2, tw2);                                                        \
    MUL(t3, tw3);                                                        \
    MUL(t4, tw4);                                                        \
    MUL(t5, tw5);                                                        \
                                                                        \
    RST_INC(t0, p_dst, os);                                                \
    RST_INC(t1, p_dst, os);                                                \
    RST_INC(t2, p_dst, os);                                                \
    RST_INC(t3, p_dst, os);                                                \
    RST_INC(t4, p_dst, os);                                                \
    RST_INC(t5, p_dst, -5*os + ovs);                                        \
   }                                                                        \

#define R2_DFT6_SKEW_T_XI4_i_ow(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int count = (_count);                                                        \
  int i, is = (_is), ivs = (_ivs), os = (_os), ovs = (_ovs);                \
  VT * p_tw = (VT *)(_p_tw);                                                \
  VT * p_src0 = (VT *)(_p_in);                                                \
  VT * p_src1 = (VT *)((char *)_p_in + is);                                \
  VT * p_dst0 = (VT *)(_p_out);                                                \
  VT * p_dst1 = (VT *)((char *)_p_out + os);                                \
  valign align;                                                                \
  VT r3_tw;                                \
  xb_vecNx16 zero = BBE_ZERONX16();                                    \
  __Pragma("ymemory(p_dst0)");                        \
  __Pragma("ymemory(p_dst1)");                        \
  __Pragma("ymemory(p_tw)");                        \
                                    \
                                    \
  RX3_TW(_sign, r3_tw);                            \
  BBE_MOVSDV(zero);                            \
                                                                        \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                                        \
                                                                        \
  for (i=0; i<count; i++) {                                                \
    VT tw1, tw2, tw3, tw4, tw5;                                                \
    VT t0, t1, t2, t3, t4, t5;                                                \
    LD_INC(t0, p_src0, 2*is);                                                \
    LD_INC(t1, p_src1, 2*is);                                                \
    LD_INC(t2, p_src0, 2*is);                                                \
    LD_INC(t3, p_src1, 2*is);                                                \
    LD_INC(t4, p_src0, -4*is + ivs);                                        \
    LD_INC(t5, p_src1, -4*is + ivs);                                        \
                                                                        \
    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);                        \
    LDA_IX4_TW_UNPACK_2(align, p_tw, tw4, tw5);                                \
                                                                        \
    DFT6(t0, t1, t2, t3, t4, t5, r3_tw);                                \
                                                                        \
    MUL(t1, tw1);                                                        \
    MUL(t2, tw2);                                                        \
    MUL(t3, tw3);                                                        \
    MUL(t4, tw4);                                                        \
    MUL(t5, tw5);                                                        \
                                                                        \
    RST_INC(t0, p_dst0, 2*os);                                                \
    RST_INC(t1, p_dst1, 2*os);                                                \
    RST_INC(t2, p_dst0, 2*os);                                                \
    RST_INC(t3, p_dst1, 2*os);                                                \
    RST_INC(t4, p_dst0, -4*os + ovs);                                        \
    RST_INC(t5, p_dst1, -4*os + ovs);                                        \
   }                                                                        \



/* rank 2 skew DFT 8, in input t4 expected in state A, t6 in state B,
 * t5 in state C and t7 in state D.
 */

#define DFT8(t0, t1, t2, t3, t4, t5, t6, t7, r8_e4, r8_e8)        \
  {                                                                \
    VT rt0, rt1, rt4, rt5;                                        \
                                                                    \
    rt0 = t0; rt1 = t2;                                                        \
    rt4 = t1; rt5 = t4;                                                        \
                                                                            \
    DFT4(t0, t1, t2, t3, rt0, rt1, A, B, 1);                                \
    DFT4(t4, t5, t6, t7, rt4, rt5, C, D, 1);                                \
                                                                        \
    MUL(t5, r8_e8);                                                        \
    MUL(t6, r8_e4);                                                        \
    MUL_CONJ(t7, r8_e8);                                                \
                                                                            \
    BFLY(t0, t4);                                                        \
    BFLY(t1, t5);                                                        \
    BFLY(t2, t6);                                                        \
    BFLY(t3, t7);                                                        \
  }



/*
 * mode register initialization;
 *
 *  N: radix size;
 *  sign: direction (-1 for forward; 1 for inverse);
 *  skew: non zero if it's a pass with twiddle factors multiplication;
 *  range: NSA returned by BBE_RRANGE();
 *  scaling: inout variable, incremented acording to shift amount set in BBE_MODE; 
 *
 */

#if 0

/* this version probably could work fast if we could keep all masks in the WVEC and
   extract as we need them instead of loading using l32r. This is a challenge though
   because compiler will try to eliminate common subexpressions and will end up
   with significantly suboptimal code */

#define RANGE_BEGIN(N, sign, skew, range, scaling)        \
  {                                                        \
    unsigned _r, _mode, _shift, _sign;                        \
                                                        \
    _sign = (sign) & 0x10;                                \
    switch (N) {                                        \
      /* radix 2  */                                                        \
    case 2: if (skew) {        /* shifting by 2 if range==0 to avoid overflow on twiddle multiplications */ \
        _mode = _shift = 0x00012;                                        \
      } else {        /* last pass, shifting only if range==0 */                \
        _mode = _shift = 0x00001;                                        \
      }                                                                        \
      break;                                                                \
      /* radix 3 is always skewed codelet */                                \
    case 3: _mode = _shift = 0x00122; break;                                \
    case 4: if (skew) {        /* taking into account possible overflow on twiddles, \
                           hence shifting by 3 if range is 0 */                \
        _mode = _shift = 0x00123;                                        \
      } else { /* last pass has no twiddles */                                \
        _mode = _shift = 0x00012;                                        \
      }                                                                        \
      break;                                                                \
      /* radix 5 uses fftadd on second stage so can shift only by 1 there */ \
    case 5: _mode = 0x00159; _shift = 0x00123; _sign = 0x10; break;        \
      /* radix 6 appears on last pass only */                \
    case 6: _mode = 0x00489; _shift = 0x00123; break;        \
      /* radix 8 appears on last pass only with fftadds at the end of calculations */ \
    case 8: _mode = 0x00159; _shift = 0x00123; break;        \
    default:                                                \
      _mode = 0; _shift = 0; break;                        \
    }                                                        \
     /*_r = ((range) <= 4 ? (range) : 4) << 2;        */        \
     _r = (range) << 2;                                                     \
    _shift = (_shift >> _r) & 0xF;                        \
    _mode = (_mode >> _r) & 0xF;                        \
    (scaling) += _shift;                                \
    BBE_FFTWMODE(_mode | _sign);                        \
    BBE_WRANGE(4);                                        \
  }

#else

#define RANGE_BEGIN(N, sign, skew, range, scaling)                        \
  {                                                                        \
    unsigned _sign = (sign) & 0x10;                                        \
    int _range = (range);                                                \
    int _shift, _mode;                                                        \
    switch (N) {                                                        \
      /* radix 2  */                                                        \
    case 2: if (skew) {        /* shifting by 2 if range==0 to avoid overflow on twiddle multiplications */ \
        _shift = 2 - _range;                                                \
      } else {        /* last pass, shifting only if range==0 */                \
        _shift = 1 - _range;                                                \
      }                                                                        \
      _shift = _shift < 0 ? 0 : _shift;                                        \
      _mode = _shift; /* mode mask is same as shift amount */                \
      break;                                                                \
      /* radix 3 is always skewed codelet */                                \
    case 3:                                                                \
      _shift = 3 - _range;                                                \
      _range = _range == 0 ? -1 : 0;                                        \
      _shift = (_shift < 0 ? 0 : _shift) + _range;                        \
      _mode = _shift; /* mode mask is p4 edit same as shift amount */        \
      break;                                                                \
    case 4: if (skew) {        /* taking into account possible overflow on twiddles, \
                           hence shifting by 3 if range is 0 */                \
        _shift = 3 - _range;                                                \
      } else { /* last pass has no twiddles */                                \
        _shift = 2 - _range;                                                \
      }                                                                        \
      _shift = _shift < 0 ? 0 : _shift;                                        \
      _mode = _shift; /* mode mask is same as shift amount */                \
      break;                                                                \
      /* radix 5 uses fftadd on second stage so can shift only by 1 there */ \
    case 5:                                                                \
      _shift = 3 - _range;                                                \
      _shift = _shift < 0 ? 0 : _shift;                                        \
      _mode = 9 - 4*_range;                                                \
      _mode = _mode < 0 ? 0 : _mode;                                        \
      _sign = 0x10;                                                        \
      break;                                                                \
      /* radix 6 appears on last pass only */                                \
    case 6:                                                                \
      _shift = 3 - _range;                                                \
      _shift = _shift < 0 ? 0 : _shift;                                        \
      _mode = (0x00489 >> (4*_range)) & 0xF;                                \
      break;                                                                \
      /* radix 8 appears on last pass only with fftadds at the end of calculations.  */ \
      /* Taking into account possible overflow on twiddles, between  butterfly radix 4 and  butterfly radix 2 */                                                                              \
    case 8:                                                                \
      _shift = 4 - _range;                                                \
      _shift = _shift < 0 ? 0 : _shift;                                        \
      _mode = 13 - 4*_range;                                                \
      _mode = _mode < 0 ? 0 : _mode;                                        \
      break;                                                                \
    default:                                                                \
      _mode = 0; _shift = 0; break;                                        \
    }                                                                        \
    (scaling) += _shift;                                                \
    BBE_FFTWMODE(_mode | _sign);                                        \
    BBE_WRANGE(4);                                                        \
  }
#endif /* RANGE_BEGIN */

#define RANGE_END(range) \
  (range) = BBE_RRANGE();


#define LDST_STRETCH_DFT4_SKEW_T_L64_16_with_lshift(S1, S2, p_src0, p_src1, p_src2, p_src3, a0, a1, a2, a3, p_dst, p_tw, part, os, ovs, shift33) \
{                                                                        \
  VT t0, t1, t2, t3, t4, t5, tw1, tw2, tw3;                                \
  xb_vecNx16 _t0, _t1, _t2, _t3;                        \
  BBE_LAVN_2XCQ15_XP(t0, a0, p_src0, part);                                \
  BBE_LAVN_2XCQ15_XP(t1, a1, p_src1, part);                                \
  BBE_LAVN_2XCQ15_XP(t2, a2, p_src2, part);                                \
  BBE_LAVN_2XCQ15_XP(t3, a3, p_src3, part);                                \
                                                                                                        \
   _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                               \
   _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                               \
   _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                               \
   _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                               \
  _t0 = BBE_SRANX16(_t0, shift33);                                   \
  _t1 = BBE_SRANX16(_t1, shift33);                                   \
  _t2 = BBE_SRANX16(_t2, shift33);                                   \
  _t3 = BBE_SRANX16(_t3, shift33);                                   \
   t0 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               \
   t1 = BBE_MOVN_2XCQ15_FROMNX16(_t1);                               \
   t2 = BBE_MOVN_2XCQ15_FROMNX16(_t2);                               \
   t3 = BBE_MOVN_2XCQ15_FROMNX16(_t3);                               \
  LD_INC_IMM(tw1, p_tw, VSIZE);                                                \
  LD_INC_IMM(tw2, p_tw, VSIZE);                                                \
  LD_INC_IMM(tw3, p_tw, VSIZE);                                                \
                                                                            \
  MV ## S1(t2);                                                                \
  MV ## S2(t3);                                                                \
                                                                        \
  __DFT4_SKEW_T(t2, t3, t4, t5, t0, t1, S1, S2, tw1, tw2, tw3)                \
                                                                        \
  INTLV(t2, t4);                                                        \
  INTLV(t3, t5);                                                        \
  INTLV(t2, t3);                                                        \
  INTLV(t4, t5);                                                        \
                                                                        \
  RST_INC(t2, p_dst, os);                                                \
  RST_INC(t3, p_dst, os);                                                \
  RST_INC(t4, p_dst, os);                                                \
  RST_INC(t5, p_dst, -3*os + ovs);                                        \
}


#define R1_DFT4_SKEW_T_L64_16_STRETCH_iw_o_with_lshift(_count, _is, _os, dummy, _ovs, _p_in, _p_out, _p_tw, _shift_) \
  int rest = ((_count) % VLEN) * ELM_SIZE;                                \
  int count = (_count) / VLEN;                                                \
  int tail = count % 2;                                                        \
  int i, is = (_is), os = (_os), ovs = (_ovs);                                \
  VT *  p_tw = (VT *)(_p_tw);                                                \
  VT *  p_src0 = (VT *)(_p_in);                                                \
  VT *  p_src1 = (VT *)((char *)(_p_in) + is);                                \
  VT *  p_src2 = (VT *)((char *)(_p_in) + 2*is);                        \
  VT *  p_src3 = (VT *)((char *)(_p_in) + 3*is);                        \
  VT *  p_dst = (VT *)(_p_out);                                                \
  valign a0, a1, a2, a3;                                                \
  vsaN shift22=_shift_;                                                                         \
  __Pragma("ymemory(p_dst)");                                                \
  /*shift22 = BBE_MOVVSA32(_shift_);  */                                                                         \
  a0 = BBE_LAN_2XCQ15_PP(p_src0);                                        \
  a1 = BBE_LAN_2XCQ15_PP(p_src1);                                        \
  a2 = BBE_LAN_2XCQ15_PP(p_src2);                                        \
  a3 = BBE_LAN_2XCQ15_PP(p_src3);                                        \
                                                                        \
  for (i=0; i<count/2; i++) {                                                \
    LDST_STRETCH_DFT4_SKEW_T_L64_16_with_lshift(A, B, p_src0, p_src1, p_src2, p_src3, a0, a1, a2, a3, p_dst, p_tw, VSIZE, os, ovs, shift22); \
    LDST_STRETCH_DFT4_SKEW_T_L64_16_with_lshift(C, D, p_src0, p_src1, p_src2, p_src3, a0, a1, a2, a3, p_dst, p_tw, VSIZE, os, ovs, shift22); \
  }                                                                        \
  if (tail) {                                                                \
    LDST_STRETCH_DFT4_SKEW_T_L64_16_with_lshift(A, B, p_src0, p_src1, p_src2, p_src3, a0, a1, a2, a3, p_dst, p_tw, VSIZE, os, ovs, shift22); \
    LDST_STRETCH_DFT4_SKEW_T_L64_16_with_lshift(C, D, p_src0, p_src1, p_src2, p_src3, a0, a1, a2, a3, p_dst, p_tw, rest, os, ovs, shift22); \
  } else {                                                                \
    LDST_STRETCH_DFT4_SKEW_T_L64_16_with_lshift(A, B, p_src0, p_src1, p_src2, p_src3, a0, a1, a2, a3, p_dst, p_tw, rest, os, ovs, shift22); \
  }



#define PRINT_RESULT(p, cnt) \
  {                                                \
    short *_p = (short *)p;                        \
    int i;                                        \
    printf ("result:\n");                        \
    for (i=0; i<(cnt)/2; i++) {                        \
      printf ("%d, %d,\n", _p[2*i], _p[2*i+1]);        \
    }                                                \
  }

#endif /* __XTFFT_BBEN_MACRO_H__ */
