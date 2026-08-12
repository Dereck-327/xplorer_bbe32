/**
    @file       debug_trace.h
    @brief      debug print info (array)

    @version    0.1.0
    @date       2026-08-10
    @author     hjk

*/

#ifndef DEBUG_TRACE_H__
#define DEBUG_TRACE_H__

#include "utils/err.h"
#include "simulate/cfg.h"

#ifndef TRACE_MASK
	#define TRACE_MASK	TRACE_DEBUG
#endif

#define TRACE_ENABLED	(0U != TRACE_MASK)

/*==========================[ 数组类型 ]=========================*/

typedef uint8_t ArrayKindType;
#define ARRAY_KIND_U8		(0U)
#define ARRAY_KIND_S16		(1U)
#define ARRAY_KIND_U16		(2U)
#define ARRAY_KIND_S32		(3U)
#define ARRAY_KIND_CPLX16	(4U)

/*==========================[ 接口 ]=========================*/

/** @brief 打印一段数组 (仅 TRACE_ENABLED 时有实现)
    @param[in] aArr   数组首地址
    @param[in] aKind  元素类型
    @param[in] aCount 元素个数
    @param[in] aLabel 打印前缀
    @retval #ERR_OK            打印完成
    @retval #ERR_INVAL_PARAMS  aArr 为空或 aKind 未知 */
ErrorType DumpArray(const void *const aArr,
                          const ArrayKindType aKind,
                          const uint16_t aCount,
                          const char *const aLabel);


#if TRACE_ENABLED
	#define TRACE(aArr, aKind, aCount, aLabel)						\
				(void) DumpArray((aArr), (aKind), (aCount), (aLabel));	\

#else
	#define TRACE(aArr, aKind, aCount, aLabel)	do { } while (0)
#endif

#endif /* DEBUG_TRACE_H__ */