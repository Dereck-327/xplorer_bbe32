#ifndef _CLANGD_TIE_SHIM_H_
#define _CLANGD_TIE_SHIM_H_
typedef int immediate;
/* clangd-only: xt-clang has these as frontend builtins. Not compiled. */
typedef long long              _TIE_xt_bbe32_valign;
typedef short  __attribute__((ext_vector_type(16))) _TIE_xt_bbe32_vboolN;
typedef short  __attribute__((ext_vector_type(8)))  _TIE_xt_bbe32_vboolN_2;
typedef int   _TIE_xt_bbe32_xb_c16, _TIE_xt_bbe32_xb_cq11_20, _TIE_xt_bbe32_xb_cq15;
typedef int   _TIE_xt_bbe32_xb_cq19_20, _TIE_xt_bbe32_xb_cq1_30, _TIE_xt_bbe32_xb_cq5_10;
typedef long long _TIE_xt_bbe32_xb_c32, _TIE_xt_bbe32_xb_c40, _TIE_xt_bbe32_xb_cq9_30;
typedef short _TIE_xt_bbe32_xb_int16; typedef unsigned short _TIE_xt_bbe32_xb_int16U;
typedef short _TIE_xt_bbe32_xb_q15;
typedef int   _TIE_xt_bbe32_xb_int32, _TIE_xt_bbe32_xb_q11_20, _TIE_xt_bbe32_xb_q19_20;
typedef int   _TIE_xt_bbe32_xb_q1_30, _TIE_xt_bbe32_xb_q5_10;
typedef long long _TIE_xt_bbe32_xb_int40, _TIE_xt_bbe32_xb_q9_30;
#define _TIE_V(n,base,name) typedef base __attribute__((ext_vector_type(n))) name;
_TIE_V(16,short,_TIE_xt_bbe32_xb_vecNx16) typedef _TIE_xt_bbe32_xb_vecNx16 _TIE_xt_bbe32_xb_vecNx16U;
_TIE_V(16,short,_TIE_xt_bbe32_xb_vecNxq15)
_TIE_V(16,int,_TIE_xt_bbe32_xb_vecNx32) typedef _TIE_xt_bbe32_xb_vecNx32 _TIE_xt_bbe32_xb_vecNx32U;
_TIE_V(16,int,_TIE_xt_bbe32_xb_vecNxq11_20) _TIE_V(16,int,_TIE_xt_bbe32_xb_vecNxq19_20)
_TIE_V(16,int,_TIE_xt_bbe32_xb_vecNxq1_30) _TIE_V(16,int,_TIE_xt_bbe32_xb_vecNxq5_10)
_TIE_V(16,long long,_TIE_xt_bbe32_xb_vecNx40) _TIE_V(16,long long,_TIE_xt_bbe32_xb_vecNxq9_30)
_TIE_V(16,int,_TIE_xt_bbe32_xb_vecNxc16) _TIE_V(16,int,_TIE_xt_bbe32_xb_vecNxcq15)
_TIE_V(16,long long,_TIE_xt_bbe32_xb_vecNxc32) _TIE_V(16,long long,_TIE_xt_bbe32_xb_vecNxc40)
_TIE_V(16,int,_TIE_xt_bbe32_xb_vecNxcq11_20) _TIE_V(16,int,_TIE_xt_bbe32_xb_vecNxcq19_20)
_TIE_V(16,int,_TIE_xt_bbe32_xb_vecNxcq1_30) _TIE_V(16,int,_TIE_xt_bbe32_xb_vecNxcq5_10)
_TIE_V(16,long long,_TIE_xt_bbe32_xb_vecNxcq9_30)
_TIE_V(8,int,_TIE_xt_bbe32_xb_vecN_2xc16) _TIE_V(8,int,_TIE_xt_bbe32_xb_vecN_2xcq15)
_TIE_V(8,long long,_TIE_xt_bbe32_xb_vecN_2xc32) _TIE_V(8,long long,_TIE_xt_bbe32_xb_vecN_2xc40)
_TIE_V(8,int,_TIE_xt_bbe32_xb_vecN_2xcq11_20) _TIE_V(8,int,_TIE_xt_bbe32_xb_vecN_2xcq19_20)
_TIE_V(8,int,_TIE_xt_bbe32_xb_vecN_2xcq1_30) _TIE_V(8,int,_TIE_xt_bbe32_xb_vecN_2xcq5_10)
_TIE_V(8,long long,_TIE_xt_bbe32_xb_vecN_2xcq9_30)
/* vsa and vsel are shift-amount/select registers; xt-clang converts them
   to/from vectors implicitly, so alias them rather than use a scalar. */
typedef _TIE_xt_bbe32_xb_vecNx16 _TIE_xt_bbe32_vsaN, _TIE_xt_bbe32_vsaNC, _TIE_xt_bbe32_vsaN_2C;
typedef _TIE_xt_bbe32_xb_vecNx16 _TIE_xt_bbe32_vsaS, _TIE_xt_bbe32_vsaSC;
typedef _TIE_xt_bbe32_xb_vecNx16 _TIE_xt_bbe32_vselN, _TIE_xt_bbe32_vselN_2;
#undef _TIE_V
typedef int _TIE_xtbool, _TIE_xtbool2, _TIE_xtbool4, _TIE_xtbool8, _TIE_xtbool16;

/* Private headers (common.h, *_common.h) assume the .c already pulled in
   NatureDSP_types.h for inline_/float32_t. Pull it in so opening a header
   standalone resolves too. */
#include "NatureDSP_types.h"
#endif
