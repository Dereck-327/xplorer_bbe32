/**
    @file       stage_emit.c
    @brief      索引归心 + IQ Reverse
                谱内下标减去 RA_IDX_CORR 变成
                以中心为零点的带符号索引; iqReverse 为真时取反 (IQ 接反的板子靠这个纠回来)。
    @version    0.1.0
    @date       2026-08-11
    @author     hjk

*/

#include "utils/err.h"
#include "utils/bbe_type.h"
#include "utils/debug_trace.h"
#include "simulate/dsp_ops.h"

//==========================[ Function Implementations ]=========================
ErrorType StageEmit(Bbe_CtxType *const aCtx)
{
	ErrorType ret = ERR_OK;

	if (NULL == aCtx)
	{
		ret = ERR_INVAL_PARAMS;
	}
	else
	{
		DSP_CalcPeaksIdxFrac(aCtx->work.peakIdx,
		                     aCtx->work.peakFrac,
		                     PEAKS_NUM,
		                     IDX_CORR,
		                     aCtx->param.iqReverse,
		                     aCtx->result.peaksIdx,
		                     aCtx->result.peaksFrac);

		TRACE(aCtx->result.peaksIdx,
		         ARRAY_KIND_S16, PEAKS_NUM, "peaks_idx:\n");
		TRACE(aCtx->result.peaksFrac,
		         ARRAY_KIND_S16, PEAKS_NUM, "peaks_frac:\n");
	}

	return ret;
}