#ifndef __XTFFT_BBE_16XCQ15_MACRO_H__
#define __XTFFT_BBE_16XCQ15_MACRO_H__

#include "common.h"

#define VLEN 16
#define ELM_SIZE 4

#define VSIZE (VLEN*ELM_SIZE)

#define VT_TW_ALIGN(tw_count) ((((tw_count)*sizeof(DT) + (VSIZE) - 1)/(VSIZE))*(VSIZE/sizeof(DT)))

#define VDIVIDES(a)  ((a) != 0 && (a) % VLEN == 0)

#define VT xb_vecN_2xcq15
#define WVT xb_vecN_2xcq9_30

#define LD(v, addr, offs)             v = BBE_LVN_2XCQ15_X(addr, offs)
#define LD_IMM(v, addr, imm_offs)     v = BBE_LVN_2XCQ15_I(addr, imm_offs)
#define _LD_IMM(addr, imm_offs)       BBE_LVN_2XCQ15_I(addr, imm_offs)

#define LD_INC(v, addr, inc)          BBE_LVN_2XCQ15_XP(v, addr, inc)
#define LD_INC_IMM(v, addr, imm_inc)  BBE_LVN_2XCQ15_IP(v, addr, imm_inc)

#define ST(v, addr, offs)             BBE_SVN_2XCQ15_X(v, addr, offs)
#define ST_IMM(v, addr, imm_offs)     BBE_SVN_2XCQ15_I(v, addr, imm_offs)

#define ST_INC(v, addr, inc)          BBE_SVN_2XCQ15_XP(v, addr, inc)
#define ST_INC_IMM(v, addr, imm_inc)  BBE_SVN_2XCQ15_IP(v, addr, imm_inc)

#define RST_INC(v, addr, inc)         BBE_SVRN_2XCQ15_XP(v, addr, inc)

#define REPL(a, b)   BBE_REPN_2XCQ15(a,b)

#ifndef IS_INV_FFT
#define MUL(a, b)       a = BBE_MULN_2XCQ15PACKQ(a,b)
#define MUL_CONJ(a, b)  a = BBE_MULN_2XCQ15JPACKQ(a,b)
#define _MUL(a, b)       a = BBE_MULNX16CPACKQ(a,b)
#define _MUL_CONJ(a, b)  a = BBE_MULNX16JPACKQ(a,b)

#else
#define MUL_CONJ(a, b)       a = BBE_MULN_2XCQ15PACKQ(a,b)
#define MUL(a, b)  a = BBE_MULN_2XCQ15JPACKQ(a,b)
#define _MUL(a, b)       a = BBE_MULNX16JPACKQ(a,b)
#define _MUL_CONJ(a, b)  a = BBE_MULNX16CPACKQ(a,b)
#endif

#define BFLY(t0, t1)                    \
  {                            \
    VT _vt0 = t0, _vt1 = t1;                \
    t0 = BBE_FFTADDSSRN_2XCQ15(_vt0, _vt1);        \
    t1 = BBE_FFTSUBSSRN_2XCQ15(_vt0, _vt1);        \
  }

#define INTLV(t0, t1) BBE_INTLV2XN_2XCQ15(t1, t0, t1, t0, BBE_INTLVC_2X2)

#define INTLV_SEL(t0, t1)              \
  {                          \
    VT _tis0 = BBE_SELN_2XCQ15I(t1, t0, BBE_SELI_INTERLEAVE_2_LO);  \
    VT _tis1 = BBE_SELN_2XCQ15I(t1, t0, BBE_SELI_INTERLEAVE_2_HI);  \
    t0 = _tis0;                      \
    t1 = _tis1;                      \
  }


#define INTLV_XI4(t0, t1) BBE_INTLV2XN_2XCQ15(t1, t0, t1, t0, BBE_INTLVC_4X4)

#define LD_A(addr, inc) BBE_LVA_XP(addr, inc)
#define LD_B(addr, inc) BBE_LVB_XP(addr, inc)
#define LD_C(addr, inc) BBE_LVC_XP(addr, inc)
#define LD_D(addr, inc) BBE_LVD_XP(addr, inc)

#define MVA(t) { xb_vecNx16 z = BBE_MOVNX16_FROMN_2XCQ15(t);  BBE_MOVSAV(z); }
#define MVB(t) { xb_vecNx16 z = BBE_MOVNX16_FROMN_2XCQ15(t);  BBE_MOVSBV(z); }
#define MVC(t) { xb_vecNx16 z = BBE_MOVNX16_FROMN_2XCQ15(t);  BBE_MOVSCV(z); }
#define MVD(t) { xb_vecNx16 z = BBE_MOVNX16_FROMN_2XCQ15(t);  BBE_MOVSDV(z); }


/* FIXME: there is an instruction to do this */
#define BROADCAST_INC(p, re, im) \
  { \
  DT __re = re; \
  DT __im = im; \
  int __i; \
  for (__i=0; __i<VLEN; __i++) {\
  *p++ = __re; \
  *p++ = __im; \
  } \
  }

#define VCOMPACT(_p_src, _p_dst, _blocks, _dst_count)    \
  {                            \
    valign align;                    \
    VT *p_src = (VT *)(_p_src);                \
    VT *p_dst = (VT *)(_p_dst);                \
    int i, j;                        \
    int block_count = (_dst_count)/VLEN;        \
    int slack = 4*((_dst_count) % VLEN);        \
    int blocks = (_blocks);                \
    /*VT t;        */                \
    __Pragma("ymemory(p_dst)");                \
                                \
    align =  BBE_ZALIGN();                \
    for (i=0; i<blocks; i++) {                        \
      for (j=0; j<block_count; j++) {                    \
    BBE_SAVN_2XCQ15_XP(*p_src++, align, p_dst, VSIZE);        \
      }                                    \
      BBE_SAVN_2XCQ15_XP(*p_src++, align, p_dst, slack);        \
    }                                    \
    BBE_SAN_2XCQ15POS_FC(align, p_dst);                    \
  }

#define VSTRETCH(_p_src, _p_dst, _blocks, _src_count)    \
  {                            \
    valign align;                    \
    VT *p_src = (VT * )(_p_src);            \
    VT *p_dst = (VT * )(_p_dst);            \
    int i, j;                        \
    int block_count = (_src_count)/VLEN;        \
    int slack = 4*((_src_count) % VLEN);        \
    int blocks = (_blocks);                \
    __Pragma("ymemory(p_dst)");                \
                                \
    align =  BBE_LAN_2XCQ15_PP(p_src);            \
                                \
    for (i=0; i<blocks; i++) {                        \
      VT temp;                                \
      for (j=0; j<block_count; j++) {                    \
    BBE_LAVN_2XCQ15_XP(temp, align, p_src, VSIZE);            \
    ST_INC_IMM(temp, p_dst, VSIZE);                    \
      }                                    \
      BBE_LAVN_2XCQ15_XP(temp, align, p_src, slack);            \
      ST_INC_IMM(temp, p_dst, VSIZE);                    \
    }                                    \
  }


#define LDA_IX4_TW_UNPACK_1(align, p_tw, tw1)                \
  { VT __tw;                                \
    BBE_LAVN_2XCQ15_XP(__tw, align, p_tw, 4*4);                \
    tw1 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_0);        \
  }

#define LDA_IX4_TW_UNPACK_2(align, p_tw, tw1, tw2)            \
  { VT __tw;                                \
    BBE_LAVN_2XCQ15_XP(__tw, align, p_tw, 4*4);                \
    tw1 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_0);        \
    tw2 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_1);        \
  }

/* xb_vecNx16 output */
#define __LDA_IX4_TW_UNPACK_2(align, p_tw, tw1, tw2)            \
  { xb_vecNx16 __tw;                                \
    BBE_LAVNX16_XP(__tw, align, p_tw, 4*4);                \
    tw1 = BBE_SHFLNX16I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_0);        \
    tw2 = BBE_SHFLNX16I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_1);        \
  }

/* xb_vecNx16 output */
#define __LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3)            \
  { xb_vecNx16 __tw;                                \
  BBE_LAVNX16_XP(__tw, align, p_tw, 6*4);                \
  tw1 = BBE_SHFLNX16I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_0);        \
  tw2 = BBE_SHFLNX16I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_1);        \
  tw3 = BBE_SHFLNX16I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_2);        \
  }

#define LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3)            \
  { VT __tw;                                \
  BBE_LAVN_2XCQ15_XP(__tw, align, p_tw, 6*4);                \
  tw1 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_0);        \
  tw2 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_1);        \
  tw3 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_2);        \
  }

#define LD_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3)            \
  { VT __tw;                                \
  BBE_LVN_2XCQ15_XP(__tw, p_tw, 2*BBE_SIMD_WIDTH);                \
  tw1 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_0);        \
  tw2 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_1);        \
  tw3 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_2);        \
  }


#define LDA_IX4_TW_UNPACK_4(align, p_tw, tw1, tw2, tw3, tw4)        \
  { VT __tw;                                \
  BBE_LAVN_2XCQ15_XP(__tw, align, p_tw, 16*4);                \
  tw1 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_0);        \
  tw2 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_1);        \
  tw3 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_2);        \
  tw4 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_3);        \
  }

#define LD_IX4_TW_UNPACK_4(p_tw, tw1, tw2, tw3, tw4)            \
  { VT __tw;                                \
    LD_INC_IMM(__tw, p_tw, VSIZE);                    \
    tw1 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_0);        \
    tw2 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_1);        \
    tw3 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_2);        \
    tw4 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_3);        \
  }

/*
 * common twiddles, replicated in to the whole vector
 */

#define __rx3_tw1_idx 0
#define __rx5_tw1_idx 1
#define __rx5_tw2_idx 2
#define __rx5_tw3_idx 3
#define __rx8_tw1_idx 4
#define __rx8_tw2_idx 5
#define __rx9_tw1_idx 8
#define __rx9_tw2_idx 9
#define __rx9_tw3_idx 10

#define RX_TWIDDLES(sign) \
  _LD_IMM( (VT*)( (sign) < 0 ? __xtfft_fwd_twiddles : __xtfft_inv_twiddles ), 0 )

#define RX3_TW(sign, tw)                             \
{                                          \
  VT t = RX_TWIDDLES(sign);                \
  tw = BBE_REPN_2XCQ15 (t, __rx3_tw1_idx); \
}

#define RX5_TW(sign, tw1, tw2, tw3)            \
  {                            \
    VT __tw = RX_TWIDDLES(sign);            \
    tw1 = BBE_REPN_2XCQ15 (__tw, __rx5_tw1_idx);    \
    tw2 = BBE_REPN_2XCQ15 (__tw, __rx5_tw2_idx);    \
    tw3 = BBE_REPN_2XCQ15 (__tw, __rx5_tw3_idx);    \
  }
  
#define RX8_TW(sign, tw1, tw2)                \
  {                            \
    VT __tw = RX_TWIDDLES(sign);            \
    tw1 = BBE_REPN_2XCQ15 (__tw, __rx8_tw1_idx);    \
    tw2 = BBE_REPN_2XCQ15 (__tw, __rx8_tw2_idx);    \
  }

#define RX9_TW(tw1, tw2, tw3)                                                                   \
{                                                                                               \
    VT __tw = BBE_LVN_2XCQ15_I((VT*)__xtfft_fwd_twiddles, (__rx9_tw1_idx>>3)*2*BBE_SIMD_WIDTH); \
    tw1 = BBE_REPN_2XCQ15 (__tw, __rx9_tw1_idx &(BBE_SIMD_WIDTH/2-1));                         \
    tw2 = BBE_REPN_2XCQ15 (__tw, __rx9_tw2_idx &(BBE_SIMD_WIDTH/2-1));                         \
    tw3 = BBE_REPN_2XCQ15 (__tw, __rx9_tw3_idx &(BBE_SIMD_WIDTH/2-1));                         \
}


#include "bben_macro.h"


#define __DFT4_SKEW_T(y0, y1, y2, y3, x0, x1, S1, S2, tw1, tw2, tw3) \
{                                    \
  DFT4(y0, y1, y2, y3, x0, x1, S1, S2, 0);                \
                                    \
  MUL(y1, tw1);                                \
  MUL(y2, tw2);                                \
  MUL(y3, tw3);                                \
}

#ifdef NO_INTLV_STORE
/* This version of the first pass R1_DFT4_SKEW_T_L64_16 used when no
   store with interleaving available. It is 1 cycle over load-store
   bound per radix 4.
   Using SELI instead of INTLV gives better loop schedule but if loop
   is unrolled INTLV is better */

#define LDST_DFT4_SKEW_T_L64_16(S1, S2, p_src, p_dst, p_tw, is, os, ivs, ovs) \
{                                    \
  VT t0, t1, t2, t3, t4, t5, tw1, tw2, tw3;                \
  LD_INC(t4, p_src, is);                        \
  LD_INC(t5, p_src, is);                        \
  LD_ ## S1(p_src, is);                            \
  LD_ ## S2(p_src, -3*is + ivs);                    \
                                    \
  LD_INC_IMM(tw1, p_tw, VSIZE);                        \
  LD_INC_IMM(tw2, p_tw, VSIZE);                        \
  LD_INC_IMM(tw3, p_tw, VSIZE);                        \
                                    \
  __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, S1, S2, tw1, tw2, tw3)        \
                                    \
  INTLV(t0, t2);                            \
  INTLV(t1, t3);                            \
  INTLV_SEL(t0, t1);                            \
  INTLV_SEL(t2, t3);                            \
                                    \
  RST_INC(t0, p_dst, os);                        \
  RST_INC(t1, p_dst, os);                        \
  RST_INC(t2, p_dst, os);                        \
  RST_INC(t3, p_dst, -3*os + ovs);                    \
}

#define R1_DFT4_SKEW_T_L64_16_iw_o(_count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int count = (_count) / (VLEN);                    \
  int i, is = (_is), ivs = (_ivs), os = (_os), ovs = (_ovs);        \
  VT * p_tw = (VT *)(_p_tw);                        \
  VT * p_src = (VT *)(_p_in);                        \
  VT * p_dst = (VT *)(_p_out);                        \
                                    \
  __Pragma("ymemory(p_dst)");                        \
                                    \
  for (i=0; i<count; i++) {                        \
    LDST_DFT4_SKEW_T_L64_16(A, B, p_src, p_dst, p_tw, is, os, ivs, ovs) \
  }                                    \

#define R1_DFT4_SKEW_T_L64_16_i_ow(_count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int count = (_count) / (VLEN);                    \
  int i, is = (_is), ivs = (_ivs), os = (_os), ovs = (_ovs);        \
  VT * p_tw = (VT *)(_p_tw);                        \
  VT * p_src = (VT *)(_p_in);                        \
  VT * p_dst = (VT *)(_p_out);                        \
                                    \
  __Pragma("ymemory(p_tw)");                        \
  __Pragma("ymemory(p_dst)");                        \
                                    \
  for (i=0; i<count; i++) {                        \
    LDST_DFT4_SKEW_T_L64_16(A, B, p_src, p_dst, p_tw, is, os, ivs, ovs) \
  }                                    \

#endif /* NO_INTLV_STORE */

#define R1_DFT4_SKEW_T_L64_16_iw_o(_count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int count = (_count) / VLEN;                                        \
  int i, is = (_is), ivs = (_ivs), os = (_os), ovs = (_ovs);        \
  VT * p_tw = (VT *)(_p_tw);                                        \
  VT * p_src = (VT *)(_p_in);                                       \
  xb_vecNx16 * p_dst = (xb_vecNx16 *)(_p_out);                      \
  valign uu0, uu1;                                                  \
  VT t0, t1, t2, t3, t4, t5, tw1, tw2, tw3;                         \
  xb_vecNx16 _t0, _t1, _t2, _t3;                                    \
  __Pragma("ymemory(p_dst)");                                       \
                                                                    \
  LD_INC(t4, p_src, is);                                            \
  LD_INC(t5, p_src, is);                                            \
  LD_A(p_src, is);                                                  \
  LD_B(p_src, -3*is + ivs);                                         \
                                                                    \
  LD_INC_IMM(tw1, p_tw, VSIZE);                                     \
  LD_INC_IMM(tw2, p_tw, VSIZE);                                     \
  LD_INC_IMM(tw3, p_tw, VSIZE);                                     \
                                                                    \
  __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);       \
                                                                    \
  INTLV(t0, t2);                                                    \
  INTLV(t1, t3);                                                    \
    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                             \
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                             \
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                             \
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                             \
                                                                    \
        uu0 = BBE_MOVUVR(_t0);                                      \
        uu1 = BBE_MOVUVR(_t2);                                      \
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*os, 0);              \
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  -os, 0);              \
  for (i=0; i<count-1; i++) {                                       \
    xb_vecNx16 _t0, _t1, _t2, _t3;                                  \
    LD_INC(t4, p_src, is);                                          \
    LD_INC(t5, p_src, is);                                          \
    LD_A(p_src, is);                                                \
    LD_B(p_src, -3*is + ivs);                                       \
                                                                    \
    LD_INC_IMM(tw1, p_tw, VSIZE);                                   \
    LD_INC_IMM(tw2, p_tw, VSIZE);                                   \
    LD_INC_IMM(tw3, p_tw, VSIZE);                                   \
                                                                    \
    __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);        \
                                                                    \
    INTLV(t0, t2);                                                  \
    INTLV(t1, t3);                                                  \
    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                             \
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                             \
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                             \
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                             \
                                                                    \
    BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*os);                        \
    BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, -3*os + ovs);              \
    BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*os, 0);                  \
    BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  -os, 0);                  \
  }                                                                 \
  _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                               \
  BBE_SALIGNVRNX16_XP(_t3, uu0, p_dst, 2*os);                       \
  BBE_SALIGNVRNX16_XP(_t3, uu1, p_dst, 2*os);


#define R1_DFT4_SKEW_T_L64_16_iw_o_with_lshift(_count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw, _shift) \
  int count = (_count) / VLEN;                                        \
  int i, is = (_is), ivs = (_ivs), os = (_os), ovs = (_ovs);        \
  VT * p_tw = (VT *)(_p_tw);                                        \
  xb_vecNx16 * p_src = (xb_vecNx16 *)(_p_in);                                       \
  xb_vecNx16 * p_dst = (xb_vecNx16 *)(_p_out);                      \
  valign uu0, uu1;                                                  \
  VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         \
  xb_vecNx16 _t0, _t1, _t2, _t3;                                    \
  /*vsaN lshift = BBE_MOVVSA32(_shift);       */         \
  __Pragma("ymemory(p_dst)");                                       \
                                                                    \
  BBE_LVNX16_XP(_t0, p_src, is);                                    \
  BBE_LVNX16_XP(_t1, p_src, is);                                    \
  BBE_LVNX16_XP(_t2, p_src, is);                                    \
  BBE_LVNX16_XP(_t3, p_src, -3*is + ivs);                           \
                                                                    \
  _t0 = BBE_SRANX16(_t0, _shift);                                   \
  _t1 = BBE_SRANX16(_t1, _shift);                                   \
  _t2 = BBE_SRANX16(_t2, _shift);                                   \
  _t3 = BBE_SRANX16(_t3, _shift);                                   \
                                                                    \
  t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               \
  t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);                               \
  BBE_MOVSAV(_t2);                                                  \
  BBE_MOVSBV(_t3);                                                  \
                                                                    \
  LD_INC_IMM(tw1, p_tw, VSIZE);                                     \
  LD_INC_IMM(tw2, p_tw, VSIZE);                                     \
  LD_INC_IMM(tw3, p_tw, VSIZE);                                     \
                                                                    \
  __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);       \
                                                                    \
  INTLV(t0, t2);                                                    \
  INTLV(t1, t3);                                                    \
    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                             \
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                             \
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                             \
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                             \
                                                                    \
        uu0 = BBE_MOVUVR(_t0);                                      \
        uu1 = BBE_MOVUVR(_t2);                                      \
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*os, 0);              \
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  -os, 0);              \
  for (i=0; i<count-1; i++) {                                       \
    xb_vecNx16 _t0, _t1, _t2, _t3;                                  \
      BBE_LVNX16_XP(_t0, p_src, is);                                    \
      BBE_LVNX16_XP(_t1, p_src, is);                                    \
      BBE_LVNX16_XP(_t2, p_src, is);                                    \
      BBE_LVNX16_XP(_t3, p_src, -3*is + ivs);                           \
                                                                        \
      _t0 = BBE_SRANX16(_t0, _shift);                                   \
      _t1 = BBE_SRANX16(_t1, _shift);                                   \
      _t2 = BBE_SRANX16(_t2, _shift);                                   \
      _t3 = BBE_SRANX16(_t3, _shift);                                   \
                                                                        \
      t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               \
      t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);                               \
      BBE_MOVSAV(_t2);                                                  \
      BBE_MOVSBV(_t3);                                                  \
                                                                    \
    LD_INC_IMM(tw1, p_tw, VSIZE);                                   \
    LD_INC_IMM(tw2, p_tw, VSIZE);                                   \
    LD_INC_IMM(tw3, p_tw, VSIZE);                                   \
                                                                    \
    __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);        \
                                                                    \
    INTLV(t0, t2);                                                  \
    INTLV(t1, t3);                                                  \
    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                             \
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                             \
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                             \
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                             \
                                                                    \
    BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*os);                        \
    BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, -3*os + ovs);              \
    BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*os, 0);                  \
    BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  -os, 0);                  \
  }                                                                 \
  _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                               \
  BBE_SALIGNVRNX16_XP(_t3, uu0, p_dst, 2*os);                       \
  BBE_SALIGNVRNX16_XP(_t3, uu1, p_dst, 2*os);


#define R1_DFT4_SKEW_T_L64_16_i_ow(_count, _is, _os, _ivs, _ovs, _p_in, _p_out, _p_tw) \
  int count = (_count) / VLEN;                        \
  int i, is = (_is), ivs = (_ivs), os = (_os), ovs = (_ovs);        \
  VT * p_tw = (VT *)(_p_tw);                        \
  VT * p_src = (VT *)(_p_in);                        \
  VT * p_dst = (VT *)(_p_out);                        \
  valign uu0, uu1;                            \
  xb_vecNx16 _t0, _t1, _t2, _t3;                                                             \
                                    \
  __Pragma("ymemory(p_dst)");                        \
  __Pragma("ymemory(p_tw)");                        \
                                    \
  VT t0, t1, t2, t3, t4, t5, tw1, tw2, tw3;                \
  LD_INC(t4, p_src, is);                        \
  LD_INC(t5, p_src, is);                        \
  LD_A(p_src, is);                            \
  LD_B(p_src, -3*is + ivs);                        \
                                    \
  LD_INC_IMM(tw1, p_tw, VSIZE);                        \
  LD_INC_IMM(tw2, p_tw, VSIZE);                        \
  LD_INC_IMM(tw3, p_tw, VSIZE);                        \
                                    \
  __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);        \
                                    \
  INTLV(t0, t2);                            \
  INTLV(t1, t3);                            \
    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                 \
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                 \
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                 \
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                 \
                                    \
  uu0 = BBE_MOVUVR((_t0));            \
  uu1 = BBE_MOVUVR((_t2));            \
  BBE_SVINTLARNX16_XP((_t1), uu0, p_dst, 2*os, 0); \
  BBE_SVINTLARNX16_XP((_t3), uu1, p_dst,  -os, 0); \
                                    \
  for (i=0; i<count-1; i++) {                        \
                                    \
    LD_INC(t4, p_src, is);                        \
    LD_INC(t5, p_src, is);                        \
    LD_A(p_src, is);                            \
    LD_B(p_src, -3*is + ivs);                        \
                                    \
    LD_INC_IMM(tw1, p_tw, VSIZE);                    \
    LD_INC_IMM(tw2, p_tw, VSIZE);                    \
    LD_INC_IMM(tw3, p_tw, VSIZE);                    \
                                    \
    __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);        \
                                    \
    INTLV(t0, t2);                            \
    INTLV(t1, t3);                            \
    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                 \
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                 \
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                 \
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                 \
                                    \
    BBE_SALIGNVRNX16_XP((_t0), uu0, p_dst, 2*os);    \
    BBE_SALIGNVRNX16_XP((_t2), uu1, p_dst, -3*os + ovs); \
    BBE_SVINTLARNX16_XP((_t1), uu0, p_dst, 2*os, 0); \
    BBE_SVINTLARNX16_XP((_t3), uu1, p_dst,  -os, 0); \
  }                                    \
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                 \
  BBE_SALIGNVRNX16_XP((_t3), uu0, p_dst, 2*os);    \
  BBE_SALIGNVRNX16_XP((_t3), uu1, p_dst, 2*os);


/* version with zero padding for DFTs */

#define LDST_STRETCH_DFT4_SKEW_T_L64_16(S1, S2, p_src0, p_src1, p_src2, p_src3, a0, a1, a2, a3, p_dst, p_tw, part, os, ovs) \
{                                    \
  VT t0, t1, t2, t3, t4, t5, tw1, tw2, tw3;                \
  BBE_LAVN_2XCQ15_XP(t0, a0, p_src0, part);                \
  BBE_LAVN_2XCQ15_XP(t1, a1, p_src1, part);                \
  BBE_LAVN_2XCQ15_XP(t2, a2, p_src2, part);                \
  BBE_LAVN_2XCQ15_XP(t3, a3, p_src3, part);                \
                                    \
  LD_INC_IMM(tw1, p_tw, VSIZE);                        \
  LD_INC_IMM(tw2, p_tw, VSIZE);                        \
  LD_INC_IMM(tw3, p_tw, VSIZE);                        \
                                        \
  MV ## S1(t2);                                \
  MV ## S2(t3);                                \
                                    \
  __DFT4_SKEW_T(t2, t3, t4, t5, t0, t1, S1, S2, tw1, tw2, tw3)        \
                                    \
  INTLV(t2, t4);                            \
  INTLV(t3, t5);                            \
  INTLV(t2, t3);                            \
  INTLV(t4, t5);                            \
                                    \
  RST_INC(t2, p_dst, os);                        \
  RST_INC(t3, p_dst, os);                        \
  RST_INC(t4, p_dst, os);                        \
  RST_INC(t5, p_dst, -3*os + ovs);                    \
}



#define R1_DFT4_SKEW_T_L64_16_STRETCH_iw_o(_count, _is, _os, dummy, _ovs, _p_in, _p_out, _p_tw) \
  int rest = ((_count) % VLEN) * ELM_SIZE;                \
  int count = (_count) / VLEN;                        \
  int tail = count % 2;                            \
  int i, is = (_is), os = (_os), ovs = (_ovs);                \
  VT *  p_tw = (VT *)(_p_tw);                        \
  VT *  p_src0 = (VT *)(_p_in);                        \
  VT *  p_src1 = (VT *)((char *)(_p_in) + is);                \
  VT *  p_src2 = (VT *)((char *)(_p_in) + 2*is);            \
  VT *  p_src3 = (VT *)((char *)(_p_in) + 3*is);            \
  VT *  p_dst = (VT *)(_p_out);                        \
  valign a0, a1, a2, a3;                        \
                                    \
  __Pragma("ymemory(p_dst)");                        \
                                    \
  a0 = BBE_LAN_2XCQ15_PP(p_src0);                    \
  a1 = BBE_LAN_2XCQ15_PP(p_src1);                    \
  a2 = BBE_LAN_2XCQ15_PP(p_src2);                    \
  a3 = BBE_LAN_2XCQ15_PP(p_src3);                    \
                                    \
  for (i=0; i<count/2; i++) {                        \
    LDST_STRETCH_DFT4_SKEW_T_L64_16(A, B, p_src0, p_src1, p_src2, p_src3, a0, a1, a2, a3, p_dst, p_tw, VSIZE, os, ovs); \
    LDST_STRETCH_DFT4_SKEW_T_L64_16(C, D, p_src0, p_src1, p_src2, p_src3, a0, a1, a2, a3, p_dst, p_tw, VSIZE, os, ovs); \
  }                                    \
  if (tail) {                                \
    LDST_STRETCH_DFT4_SKEW_T_L64_16(A, B, p_src0, p_src1, p_src2, p_src3, a0, a1, a2, a3, p_dst, p_tw, VSIZE, os, ovs); \
    LDST_STRETCH_DFT4_SKEW_T_L64_16(C, D, p_src0, p_src1, p_src2, p_src3, a0, a1, a2, a3, p_dst, p_tw, rest, os, ovs); \
  } else {                                \
    LDST_STRETCH_DFT4_SKEW_T_L64_16(A, B, p_src0, p_src1, p_src2, p_src3, a0, a1, a2, a3, p_dst, p_tw, rest, os, ovs); \
  }




#define DFT4_L16_4_xI4(t0, t1, t2, t3, S1, S2) \
{                                    \
  VT _t0, _t1;                                \
                                    \
  INTLV_XI4(t0, t2);                            \
  INTLV_XI4(t1, t3);                            \
  INTLV_XI4(t0, t1);                            \
  INTLV_XI4(t2, t3);                            \
                                    \
  _t0 = t0;                                \
  _t1 = t1;                                \
  MV ## S1 (t2);                            \
  MV ## S2 (t3);                            \
                                    \
  DFT4(t0, t1, t2, t3, _t0, _t1, S1, S2, 0);                \
}

#define LDST_DFT4_L16_4_xI4(S1, S2, p_src, p_dst, is, os, ivs, ovs)    \
{                                    \
  VT t0, t1, t2, t3;                            \
  LD_INC(t0, p_src, is);                        \
  LD_INC(t1, p_src, is);                        \
  LD_INC(t2, p_src, is);                        \
  LD_INC(t3, p_src, -3*is + ivs);                    \
                                    \
  DFT4_L16_4_xI4(t0, t1, t2, t3, S1, S2);                \
                                    \
  RST_INC(t0, p_dst, os);                        \
  RST_INC(t1, p_dst, os);                        \
  RST_INC(t2, p_dst, os);                        \
  RST_INC(t3, p_dst, -3*os + ovs);                    \
}

#define LDST_COMPACT_DFT4_L16_4_xI4(S1, S2, p_src, p_dst0, p_dst1, p_dst2, p_dst3, a0, a1, a2, a3, is, ivs, part) \
{                                    \
  VT t0, t1, t2, t3;                            \
  LD_INC(t0, p_src, is);                        \
  LD_INC(t1, p_src, is);                        \
  LD_INC(t2, p_src, is);                        \
  LD_INC(t3, p_src, -3*is + ivs);                    \
                                    \
  DFT4_L16_4_xI4(t0, t1, t2, t3, S1, S2);                \
                                    \
  BBE_SAVRN_2XCQ15_XP(t0, a0, p_dst0, part);                \
  BBE_SAVRN_2XCQ15_XP(t1, a1, p_dst1, part);                \
  BBE_SAVRN_2XCQ15_XP(t2, a2, p_dst2, part);                \
  BBE_SAVRN_2XCQ15_XP(t3, a3, p_dst3, part);                \
}

#define R2_DFT4_L_i_o(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out) \
  int tail = ((_count) / VLEN) % 2;                    \
  int count = (_count) / (2 * VLEN);                    \
  int i, is = (_is), ivs = (_ivs), os = (_os), ovs = (_ovs);        \
  VT * p_src = (VT *)(_p_in);                        \
  VT * p_dst = (VT *)(_p_out);                        \
                                    \
  __Pragma("ymemory(p_dst)");                        \
                                    \
  for (i=0; i<count; i++) {                        \
    LDST_DFT4_L16_4_xI4(A, B, p_src, p_dst, is, os, ivs, ovs);        \
    LDST_DFT4_L16_4_xI4(C, D, p_src, p_dst, is, os, ivs, ovs);        \
  }                                    \
  if (tail) {                                \
    LDST_DFT4_L16_4_xI4(A, B, p_src, p_dst, is, os, ivs, ovs);        \
  }

#define R2_DFT4_L_COMPACT_i_o(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out) \
  int rest = ((_count) % VLEN) * ELM_SIZE;                \
  int count = (_count) / VLEN;                        \
  int tail = count % 2;                            \
  int i, is = (_is), ivs = (_ivs), os = (_os);                \
  VT * p_src = (VT *)(_p_in);                        \
  VT * p_dst0 = (VT *)(_p_out);                        \
  VT * p_dst1 = (VT *)((char *)(_p_out) + os);                \
  VT * p_dst2 = (VT *)((char *)(_p_out) + 2*os);            \
  VT * p_dst3 = (VT *)((char *)(_p_out) + 3*os);            \
  valign a0, a1, a2, a3;                        \
                                    \
  __Pragma("ymemory(p_dst0)");                        \
  __Pragma("ymemory(p_dst1)");                        \
  __Pragma("ymemory(p_dst2)");                        \
  __Pragma("ymemory(p_dst3)");                        \
                                    \
  a0 =  BBE_ZALIGN();                            \
  a1 =  BBE_ZALIGN();                            \
  a2 =  BBE_ZALIGN();                            \
  a3 =  BBE_ZALIGN();                            \
                                    \
  for (i=0; i<count/2; i++) {                        \
    LDST_COMPACT_DFT4_L16_4_xI4(A, B, p_src, p_dst0, p_dst1, p_dst2, p_dst3, a0, a1, a2, a3, is, ivs, VSIZE); \
    LDST_COMPACT_DFT4_L16_4_xI4(C, D, p_src, p_dst0, p_dst1, p_dst2, p_dst3, a0, a1, a2, a3, is, ivs, VSIZE); \
  }                                    \
  if (tail) {                                \
    LDST_COMPACT_DFT4_L16_4_xI4(A, B, p_src, p_dst0, p_dst1, p_dst2, p_dst3, a0, a1, a2, a3, is, ivs, VSIZE); \
    LDST_COMPACT_DFT4_L16_4_xI4(C, D, p_src, p_dst0, p_dst1, p_dst2, p_dst3, a0, a1, a2, a3, is, ivs, rest); \
  } else {                                \
    LDST_COMPACT_DFT4_L16_4_xI4(A, B, p_src, p_dst0, p_dst1, p_dst2, p_dst3, a0, a1, a2, a3, is, ivs, rest); \
  }                                    \
  BBE_SAN_2XCQ15POS_FC(a0, p_dst0);                    \
  BBE_SAN_2XCQ15POS_FC(a1, p_dst1);                    \
  BBE_SAN_2XCQ15POS_FC(a2, p_dst2);                    \
  BBE_SAN_2XCQ15POS_FC(a3, p_dst3);



#define R2_DFT8_L_i_o(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out) \
  int i;                                \
  int count = (_count) / VLEN;                        \
  int is = (_is);                            \
  int ivs = (_ivs);                            \
  int os = (_os);                            \
  int ovs = (_ovs);                            \
  VT * p_src0 = (VT *)((char *)_p_in + 3*is);                \
  VT * p_src1 = (VT *)((char *)_p_in);                    \
  VT * p_dst0 = (VT *)(_p_out);                        \
  VT * p_dst1 = (VT *)((char *)_p_out + 4*os);                \
  VT e4, e8;                                \
                                    \
  __Pragma("ymemory(p_dst0)");                        \
  __Pragma("ymemory(p_dst1)");                        \
                                    \
  RX8_TW(_sign, e4, e8);                        \
                                    \
  for (i=0; i<count; i++) {                        \
    VT t0, t1, t2, t3, t4, t5, t6, t7, rt0, rt1, rt4, rt5;        \
                                    \
    LD_INC(t3, p_src0, 4*is);                        \
    LD_INC(t7, p_src0,-5*is);                        \
    LD_INC(t2, p_src0, 4*is);                        \
    LD_INC(t6, p_src0,-3*is + ivs);                    \
    LD_INC(t0, p_src1, 4*is);                        \
    LD_INC(t4, p_src1,-3*is);                        \
    LD_INC(t1, p_src1, 4*is);                        \
    LD_INC(t5, p_src1,-5*is + ivs);                    \
                                    \
    INTLV_XI4(t0, t4);                            \
    INTLV_XI4(t1, t5);                            \
    INTLV_XI4(t2, t6);                            \
    INTLV_XI4(t3, t7);                            \
                                    \
    INTLV_XI4(t0, t2);                            \
    INTLV_XI4(t4, t6);                            \
    INTLV_XI4(t1, t3);                            \
    INTLV_XI4(t5, t7);                            \
                                    \
    rt0 = t0; rt1 = t4; MVA(t1); MVB(t5);                \
    rt4 = t2; rt5 = t6; MVC(t3); MVD(t7);                \
                                    \
    DFT4(t0, t1, t2, t3, rt0, rt1, A, B, 1);                \
    DFT4(t4, t5, t6, t7, rt4, rt5, C, D, 1);                \
                                    \
    MUL(t5, e8);                            \
    MUL(t6, e4);                            \
    MUL_CONJ(t7, e8);                            \
                                        \
    BFLY(t0, t4);                            \
    BFLY(t1, t5);                            \
    BFLY(t2, t6);                            \
    BFLY(t3, t7);                            \
                                        \
    RST_INC(t0, p_dst0, os);                        \
    RST_INC(t1, p_dst0, os);                        \
    RST_INC(t2, p_dst0, os);                        \
    RST_INC(t7, p_dst0, -3*os + ovs);                    \
    RST_INC(t4, p_dst1, os);                        \
    RST_INC(t5, p_dst1, os);                        \
    RST_INC(t6, p_dst1, os);                        \
    RST_INC(t3, p_dst1, -3*os + ovs);                    \
  }

/*
 * last pass radix 4 with half of L(16,4) x I(4) moved to the previos pass
 * where it's implemented using store with interleaving.
 */

#define DFT4_L16_4_split_xI4(t0, t1, t2, t3, S1, S2) \
{                                    \
  VT _t0, _t1, _t2, _t3;                        \
                                    \
  _t0 = BBE_SELN_2XCQ15I(t2, t0, BBE_SELI_EXTRACT_LO_HALVES);        \
  _t1 = BBE_SELN_2XCQ15I(t2, t0, BBE_SELI_EXTRACT_HI_HALVES);        \
  _t2 = BBE_SELN_2XCQ15I(t3, t1, BBE_SELI_EXTRACT_LO_HALVES);        \
  _t3 = BBE_SELN_2XCQ15I(t3, t1, BBE_SELI_EXTRACT_HI_HALVES);        \
                                    \
  MV ## S1 (_t2);                            \
  MV ## S2 (_t3);                            \
                                    \
  DFT4(t0, t1, t2, t3, _t0, _t1, S1, S2, 0);                \
}

#define LDST_DFT4_L16_4_split_xI4(S1, S2, p_src, p_dst, is, os, ivs, ovs)    \
{                                    \
  VT t0, t1, t2, t3;                            \
  LD_INC(t0, p_src, is);                        \
  LD_INC(t1, p_src, is);                        \
  LD_INC(t2, p_src, is);                        \
  LD_INC(t3, p_src, -3*is + ivs);                    \
                                    \
  DFT4_L16_4_split_xI4(t0, t1, t2, t3, S1, S2);                \
                                    \
  RST_INC(t0, p_dst, os);                        \
  RST_INC(t1, p_dst, os);                        \
  RST_INC(t2, p_dst, os);                        \
  RST_INC(t3, p_dst, -3*os + ovs);                    \
}

#define R2_DFT4_L_split_i_o(_sign, _count, _is, _os, _ivs, _ovs, _p_in, _p_out) \
  int tail = ((_count) / VLEN) % 2;                    \
  int count = (_count) / (2 * VLEN);                    \
  int i, is = (_is), ivs = (_ivs), os = (_os), ovs = (_ovs);        \
  VT * p_src = (VT *)(_p_in);                        \
  VT * p_dst = (VT *)(_p_out);                        \
                                    \
  __Pragma("ymemory(p_dst)");                        \
                                    \
  for (i=0; i<count; i++) {                        \
    LDST_DFT4_L16_4_split_xI4(A, B, p_src, p_dst, is, os, ivs, ovs);    \
    LDST_DFT4_L16_4_split_xI4(C, D, p_src, p_dst, is, os, ivs, ovs);    \
  }                                    \
  if (tail) {                                \
    LDST_DFT4_L16_4_split_xI4(A, B, p_src, p_dst, is, os, ivs, ovs);    \
  }

/* Second to last pass with part of L(16,4) implemented using stores with interleaving.
   Assuming next pass (the last one is R2_DFT4_L_split_i_o), this one has same formula as
   R3_DFT4_SKEW_T_XI4 with _count2 == 1 && (_count1 % 2) == 0
*/

#define R3_DFT4_SKEW_T_XI4_split_iw_o(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                \
  int count1 = (_count1) / 2;                        \
  int is = (_is), ivs1 = (_ivs1), ivs2 = (_ivs2);            \
  int os = (_os), ovs1 = (_ovs1), ovs2 = (_ovs2);            \
  VT * p_tw = (VT *)(_p_tw);                        \
  VT * p_src = (VT *)((char *)_p_in + 3*is);                \
  VT * p_dst = (VT *)((char *)_p_out);                    \
  valign align, uu0, uu1;                        \
                                    \
  __Pragma("ymemory(p_dst)");                        \
                                    \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                    \
                                    \
  VT x0, x1, x2, x3, y0, y1, y2, y3, t4, t5, tw1, tw2, tw3;        \
  LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);            \
                                    \
  LD_B(p_src, -is);                            \
  LD_A(p_src, -is);                            \
  LD_INC(t5, p_src, -is);                        \
  LD_INC(t4, p_src, 3*is + ivs1);                    \
                                    \
  __DFT4_SKEW_T(x0, x1, x2, x3, t4, t5, A, B, tw1, tw2, tw3);        \
                                    \
  LD_D(p_src, -is);                            \
  LD_C(p_src, -is);                            \
  LD_INC(t5, p_src, -is);                        \
  LD_INC(t4, p_src, 3*is + ivs1);                    \
                                    \
  __DFT4_SKEW_T(y0, y1, y2, y3, t4, t5, C, D, tw1, tw2, tw3);        \
                                    \
  uu0 = BBE_MOVUVR(BBE_MOVNX16_FROMN_2XCQ15(x0));            \
  uu1 = BBE_MOVUVR(BBE_MOVNX16_FROMN_2XCQ15(x1));            \
  BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y0), uu0, p_dst,  os, 1); \
                                    \
  for (i=0; i<count1-1; i++) {                        \
    BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y1), uu1, p_dst, -os + ovs1, 1); \
    BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(x2), uu0, p_dst,  os);    \
    BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(x3), uu1, p_dst,  os - ovs1); \
    BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y2), uu0, p_dst,  os, 1); \
    BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y3), uu1, p_dst, -os + ovs1, 1); \
                                    \
    LD_B(p_src, -is);                            \
    LD_A(p_src, -is);                            \
    LD_INC(t5, p_src, -is);                        \
    LD_INC(t4, p_src, 3*is + ivs1);                    \
                                    \
    __DFT4_SKEW_T(x0, x1, x2, x3, t4, t5, A, B, tw1, tw2, tw3);        \
                                    \
    LD_D(p_src, -is);                            \
    LD_C(p_src, -is);                            \
    LD_INC(t5, p_src, -is);                        \
    LD_INC(t4, p_src, 3*is + ivs1);                    \
                                    \
    __DFT4_SKEW_T(y0, y1, y2, y3, t4, t5, C, D, tw1, tw2, tw3);        \
                                    \
    BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(x0), uu0, p_dst,  os);    \
    BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(x1), uu1, p_dst, -3*os + ovs1); \
    BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y0), uu0, p_dst,  os, 1); \
  }                                    \
  BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y1), uu1, p_dst, -os + ovs1, 1); \
  BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(x2), uu0, p_dst,  os);    \
  BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(x3), uu1, p_dst,  os - ovs1); \
  BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y2), uu0, p_dst,  os, 1); \
  BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y3), uu1, p_dst, -os + ovs1, 1); \
  BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y3), uu0, p_dst, os);    \
  BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y3), uu1, p_dst, os);    \

#define R3_DFT4_SKEW_T_XI4_split_i_ow(_sign, _count1, _count2, _is, _os, _ivs1, _ovs1, _ivs2, _ovs2, _p_in, _p_out, _p_tw) \
  int i, j;                                \
  int count1 = (_count1) / 2;                        \
  int is = (_is), ivs1 = (_ivs1), ivs2 = (_ivs2);            \
  int os = (_os), ovs1 = (_ovs1), ovs2 = (_ovs2);            \
  VT * p_tw = (VT *)(_p_tw);                        \
  VT * p_src = (VT *)((char *)_p_in + 3*is);                \
  VT * p_dst = (VT *)((char *)_p_out);                    \
  valign align, uu0, uu1;                        \
                                    \
  __Pragma("ymemory(p_tw)");                        \
  __Pragma("ymemory(p_dst)");                        \
                                    \
  align =  BBE_LAN_2XCQ15_PP(p_tw);                    \
                                    \
  VT x0, x1, x2, x3, y0, y1, y2, y3, t4, t5, tw1, tw2, tw3;        \
  LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3);            \
                                    \
  LD_B(p_src, -is);                            \
  LD_A(p_src, -is);                            \
  LD_INC(t5, p_src, -is);                        \
  LD_INC(t4, p_src, 3*is + ivs1);                    \
                                    \
  __DFT4_SKEW_T(x0, x1, x2, x3, t4, t5, A, B, tw1, tw2, tw3);        \
                                    \
  LD_D(p_src, -is);                            \
  LD_C(p_src, -is);                            \
  LD_INC(t5, p_src, -is);                        \
  LD_INC(t4, p_src, 3*is + ivs1);                    \
                                    \
  __DFT4_SKEW_T(y0, y1, y2, y3, t4, t5, C, D, tw1, tw2, tw3);        \
                                    \
  uu0 = BBE_MOVUVR(BBE_MOVNX16_FROMN_2XCQ15(x0));            \
  uu1 = BBE_MOVUVR(BBE_MOVNX16_FROMN_2XCQ15(x1));            \
  BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y0), uu0, p_dst,  os, 1); \
                                    \
  for (i=0; i<count1-1; i++) {                        \
    BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y1), uu1, p_dst, -os + ovs1, 1); \
    BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(x2), uu0, p_dst,  os);    \
    BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(x3), uu1, p_dst,  os - ovs1); \
    BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y2), uu0, p_dst,  os, 1); \
    BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y3), uu1, p_dst, -os + ovs1, 1); \
                                    \
    LD_B(p_src, -is);                            \
    LD_A(p_src, -is);                            \
    LD_INC(t5, p_src, -is);                        \
    LD_INC(t4, p_src, 3*is + ivs1);                    \
                                    \
    __DFT4_SKEW_T(x0, x1, x2, x3, t4, t5, A, B, tw1, tw2, tw3);        \
                                    \
    LD_D(p_src, -is);                            \
    LD_C(p_src, -is);                            \
    LD_INC(t5, p_src, -is);                        \
    LD_INC(t4, p_src, 3*is + ivs1);                    \
                                    \
    __DFT4_SKEW_T(y0, y1, y2, y3, t4, t5, C, D, tw1, tw2, tw3);        \
                                    \
    BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(x0), uu0, p_dst,  os);    \
    BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(x1), uu1, p_dst, -3*os + ovs1); \
    BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y0), uu0, p_dst,  os, 1); \
  }                                    \
  BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y1), uu1, p_dst, -os + ovs1, 1); \
  BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(x2), uu0, p_dst,  os);    \
  BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(x3), uu1, p_dst,  os - ovs1); \
  BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y2), uu0, p_dst,  os, 1); \
  BBE_SVINTLARNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y3), uu1, p_dst, -os + ovs1, 1); \
  BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y3), uu0, p_dst, os);    \
  BBE_SALIGNVRNX16_XP(BBE_MOVNX16_FROMN_2XCQ15(y3), uu1, p_dst, os);    \


#endif /* __XTFFT_BBE_16XCQ15_MACRO_H__ */
