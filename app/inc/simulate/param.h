/**
    @file       param.h
    @brief      simulate 运行时可调参数 描述符表
                添加参数 e.g.
                1. ParamType 里加字段
                param.c ParamDesc[] 里加一行

    @version    0.1.0
    @date       2026-08-10
    @author     hjk

*/

#ifndef SIMULATE_PARAM_H__
#define SIMULATE_PARAM_H__

#include "cfg.h"
#include "utils/compiler.h"
#include "utils/err.h"

/* 亚 bin 插值的取宽策略 */
typedef uint8_t CompModeType;
#define COMP_MODE_AUTO	(0U)	/* 用实测峰宽 */
#define COMP_MODE_MANUAL	(1U)	/* 用 compWidth 固定宽度 */
#define COMP_MODE_HYBRID	(2U)	/* 实测宽 > compWidth 时用实测, 否则用固定 */
#define COMP_MODE_MAX	(COMP_MODE_HYBRID)

/* 描述符表里一行的数据类型 */
typedef uint8_t ParamKindType;
#define PARAM_KIND_U8	(0U)
#define PARAM_KIND_S16	(1U)
#define PARAM_KIND_U16	(2U)
#define PARAM_KIND_S32	(3U)


/*==========================[ 参数结构 ]=========================*/

typedef struct ParamType
{
    uint16_t fftEq[MAG2_SIZE]               COMP_ALIGN(VEC_ALIGN) ;     /* U(16,15) 底噪 */
    uint16_t blankTable[MAX_BLANK_SIZE]     COMP_ALIGN(VEC_ALIGN) ;     /* U(16,0) 各区间屏蔽中心 */
    int16_t  mfCoeff[MF_COEFF_SIZE]         COMP_ALIGN(VEC_ALIGN) ;    /* S(16,15) 匹配滤波半边系数 */

    int16_t  powerShift;		/* S(16,0)  SNR 定标, 饱和到 RA_SNR_SHIFT_MAX */
    int16_t  compWidth;			/* S(16,0)  MANUAL/HYBRID 模式的固定半宽 */
    int16_t  scaleComp;			/* S(16,12) 亚 bin 插值增益 */
    uint16_t blankRange;		/* U(16,0)  遮蔽中心两侧各置零多少点 */
    uint8_t  blankValidSize;	/* U(8,0)   blankTable 实际有效条数 */
    CompModeType compMode;	    /* COMP_MODE_* */
    uint8_t  iqReverse;			/* U(8,0)   IQ是否反转 非零则输出索引与小数部分取反 */
} ParamType_t;


/*==========================[ 描述符表 ]=========================*/

/* 表里一行 = 一个参数的全部元信息。min/max 用 int32_t 装得下所有 kind。 */
typedef struct ParamDescType
{
	const char			*name;		/* 外部键名, 与 params.txt 里的 key 一致 */
	uint16_t			 offset;	/* offsetof(ParamType, 字段) */
	ParamKindType	     kind;
	uint16_t			 count;		/* 元素个数, 1 表示标量 */
	int32_t				 min;		/* 合法闭区间下界 */
	int32_t				 max;		/* 合法闭区间上界 */
} ParamDescType_t;

typedef ErrorType (ParamSrcFnType)(void *const aSrcCtx, const char **const aKey, char **const aVal);

/*==========================[ 接口 ]=========================*/

/** @brief 按描述符表从任意来源填充参数结构
    @param[out] aParam   目标参数结构
    @param[in]  aSrcFn   数据源回调
    @param[in]  aSrcCtx  传给回调的上下文 (如 FILE*)
    @retval #ERR_OK            全部项已加载且逐项通过范围校验
    @retval #ERR_INVAL_PARAMS  aParam/aSrcFn 为空, 或某项越界
    @retval #ERR_DATA_INTEG    数据源读取失败
    @post 未在数据源里出现的字段保持调用前的值, 调用者应先清零 */
ErrorType ParamLoad(ParamType_t *const aParam,
                          ParamSrcFnType aSrcFn,
                          void *const aSrcCtx);

/** @brief 逐项范围校验 + 跨字段不变式校验
    @param[in] aParam 待校验参数
    @retval #ERR_OK            全部合法
    @retval #ERR_INVAL_PARAMS  某项越界或不变式不成立 */
ErrorType ParamCheck(const ParamType_t *const aParam);

#endif /* SIMULATE_PARAM_H__ */
