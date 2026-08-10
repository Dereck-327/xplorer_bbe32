/**
    @file       bbe_type.h
    @brief      处理链上下文, 处理链各阶段统一函数签名`int32_t Fn(bbe_CtxType *const aCtx)`
                新增参数时流水线管道不需要修改

    @version    0.1.0
    @date       2026-08-10
    @author     hjk

*/

#ifndef BBE_TYPE_H__
#define BBE_TYPE_H__

#include "NatureDSP_types.h"
#include "simulate/param.h"
#include <stdint.h>


/*==========================[ 状态 ]=========================*/

typedef uint8_t StateType;
#define STATE_UNINIT		(0U)	/* 未 Init */
#define STATE_READY		(1U)	/* 可接受 SubmitFrame */

/*==========================[ 输入 / 中间 / 输出 ]=========================*/

/* 单个 chirp 的输入。bexp 是块指数, 定点谱的定标随帧变化。 */
typedef struct sBbe_FrameType
{
    uint16_t mag2[MAG2_SIZE]    COMP_ALIGN(VEC_ALIGN) ;	/* U(16,15) 功率谱 */
    uint8_t  dataBexp;		    /* U(8,0) FFT 输出的块指数 */
    uint8_t  magsqBexp;		    /* U(8,0) 幅度谱 取模平方引入的额外指数 */
} Bbe_FrameType;


/* 中间缓冲区, 带 _PAD 后缀的缓冲两端各有 RA_PADDING 个哨兵样点,
 * 有效数据从 [RA_PADDING] 开始。*/
typedef struct sBbe_WorkType
{
    uint16_t mag2[CHIRP_NUM][MAG2_SIZE]         COMP_ALIGN(VEC_ALIGN) ;	/* U(16,15) 各 chirp 暂存 */
    int16_t mag2Bexp[CHIRP_NUM]                COMP_ALIGN(VEC_ALIGN) ;	/* S(16,0)  各 chirp 块指数 */
    int16_t acc[MAG2_SIZE_PAD]                 COMP_ALIGN(VEC_ALIGN) ;	/* S(16,14) 非相干累积结果 */
    int16_t accEq[MAG2_SIZE_PAD]               COMP_ALIGN(VEC_ALIGN) ;	/* S(16,14) 去底噪 + 遮蔽后 */
    int16_t accEqMf[MAG2_SIZE_PAD]             COMP_ALIGN(VEC_ALIGN) ;	/* S(16,14) 匹配滤波后 */
    uint16_t peakIdx[PEAKS_NUM]                 COMP_ALIGN(VEC_ALIGN) ;	/* U(16,0)  谱内峰位 */
    int16_t  peakLeftW[PEAKS_NUM]              COMP_ALIGN(VEC_ALIGN) ;	/* S(16,0)  左半高宽 */
    int16_t  peakRightW[PEAKS_NUM]             COMP_ALIGN(VEC_ALIGN) ;	/* S(16,0)  右半高宽 */
    int16_t  peakFrac[PEAKS_NUM]                COMP_ALIGN(VEC_ALIGN) ;	/* S(16,12) 亚 bin 偏移 */
    int16_t accEqBexp;		                                            /* S(16,0) 累积结果的块指数 */
} Bbe_WorkType;

typedef struct sBbe_ResultType
{
    int16_t peaksIdx[PEAKS_NUM]		COMP_ALIGN(VEC_ALIGN) ;	/* S(16,0)  相对中心的峰位 */
    int16_t peaksFrac[PEAKS_NUM]	COMP_ALIGN(VEC_ALIGN) ;	/* S(16,12) 亚 bin 偏移 */
    uint16_t peaksFwhm[PEAKS_NUM]	COMP_ALIGN(VEC_ALIGN) ;	/* U(16,0)  半高全宽 */
    uint16_t peaksSnr[PEAKS_NUM]	COMP_ALIGN(VEC_ALIGN) ;	/* U(16,0)  信噪比 */
} Bbe_ResultType;

/*==========================[ context ]=========================*/

typedef struct sBbe_CtxType
{
    ParamType_t param;
    Bbe_FrameType frame;        /* 当前 chirp 输入 */
    Bbe_WorkType work;          /* 中间缓冲 */
    Bbe_ResultType result;      /* 最近一帧的结果 */
    uint8_t chirpIdx;           /* chirp计数 */
    StateType state;
} Bbe_CtxType;

/* 上下文按向量对齐放置, 静态实例用这个宏声明:
 *   BBE_CTX_DEFINE(static, g_bbeCtx);
 * 结构里的缓冲已各自对齐, 这里再对齐整体是为了保证首字段偏移不被打破。 */
#define BBE_CTX_DEFINE(scope, name) \
    scope Bbe_CtxType name COMP_ALIGN(VEC_ALIGN)

#endif /* BBE_TYPE_H__ */