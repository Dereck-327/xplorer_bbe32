/**
    @file       stage_refine.c
    @brief      按 compmode 做亚bin的插值

    @version    0.1.0
    @date       2026-08-11
    @author     hjk

*/

#include "simulate/cfg.h"
#include "utils/err.h"
#include "utils/debug_trace.h"
#include "utils/bbe_type.h"
#include "simulate/dsp_ops.h"

//==========================[ Static Function Prototypes ]=========================
static void width_select(const Bbe_CtxType *const aCtx,
                         const uint16_t aIdx,
                         int16_t *const aLeftW,
                         int16_t *const aRightW);

//==========================[ Static Function Implementations ]=========================
/* AUTO   : 用实测半高宽
 * MANUAL : 一律用 compWidth
 * HYBRID : 实测宽度大于 compWidth 时才信实测, 否则退回 compWidth
 *          (窄峰的实测宽度噪声大, 插值会跑偏) */
static void width_select(const Bbe_CtxType *const aCtx,
                         const uint16_t aIdx,
                         int16_t *const aLeftW,
                         int16_t *const aRightW)
{
	uint8_t useMeasured;

	if (COMP_MODE_AUTO == aCtx->param.compMode)
	{
		useMeasured = 1U;
	}
	else if (COMP_MODE_HYBRID == aCtx->param.compMode)
	{
		useMeasured = (aCtx->work.peakLeftW[aIdx] > aCtx->param.compWidth) ? 1U : 0U;
	}
	else
	{
		useMeasured = 0U;		/* MANUAL */
	}

	if (0U != useMeasured)
	{
		*aLeftW = aCtx->work.peakLeftW[aIdx];
		*aRightW = aCtx->work.peakRightW[aIdx];
	}
	else
	{
		*aLeftW = aCtx->param.compWidth;
		*aRightW = aCtx->param.compWidth;
	}
}


//==========================[ Function Implementations ]=========================

ErrorType StageRefine(Bbe_CtxType *const aCtx)
{
    ErrorType ret = ERR_OK;
	uint16_t i;
	int16_t leftW;
	int16_t rightW;

	if (NULL == aCtx)
	{
		ret = ERR_INVAL_PARAMS;
	}
	else
	{
		for (i = 0U; i < PEAKS_NUM; ++i)
		{
			width_select(aCtx, i, &leftW, &rightW);

			/* 用回灌的 accEq (未被遮蔽的滤波结果), 峰位加 RA_PADDING 换成缓冲下标 */
			DSP_PeakIdxComp(aCtx->work.accEq,
			                aCtx->work.peakIdx[i] + PADDING,
			                leftW,
			                rightW,
			                aCtx->param.scaleComp,
			                &aCtx->work.peakFrac[i]);
		}

		TRACE(aCtx->work.peakLeftW,
		         ARRAY_KIND_S16, PEAKS_NUM, "peaks_left_width:\n");
		TRACE(aCtx->work.peakRightW,
		         ARRAY_KIND_S16, PEAKS_NUM, "peaks_right_width:\n");
		TRACE(aCtx->work.peakFrac,
		         ARRAY_KIND_S16, PEAKS_NUM, "peaks_frac_buffer:\n");
	}

    return ret;
}