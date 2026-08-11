/**
    @file       stage_peak.c
    @brief      搜峰 + 峰宽 + FWHM + SNR

    @version    0.1.0
    @date       2026-08-10
    @author     hjk

*/

#include "simulate/cfg.h"
#include "utils/bbe_type.h"
#include "utils/debug_trace.h"
#include "simulate/dsp_ops.h"

ErrorType StagePeak(Bbe_CtxType *const aCtx)
{
	ErrorType ret = ERR_OK;

	if (NULL == aCtx)
	{
		ret = ERR_INVAL_PARAMS;
	}
	else
	{
		/* 传的是 frame.dataBexp: SNR 的定标用的是当前 chirp 的块指数,
		 * 而不是累积后的 work.accEqBexp。原 second_stage.c 亦然。 */
		DSP_CalcPeaksInfo(aCtx->work.accEqMf,
		                  MAG2_SIZE,
		                  aCtx->work.peakLeftW,
		                  aCtx->work.peakRightW,
		                  PEAKS_NUM,
		                  aCtx->work.peakIdx,
		                  aCtx->frame.dataBexp,
		                  aCtx->param.powerShift,
		                  aCtx->result.peaksFwhm,
		                  aCtx->result.peaksSnr);

		TRACE(TRACE_DEBUG, aCtx->work.peakIdx,
		         ARRAY_KIND_U16, PEAKS_NUM, "peaks_idx_buffer:\n");
		TRACE(TRACE_DEBUG, aCtx->result.peaksFwhm,
		         ARRAY_KIND_U16, PEAKS_NUM, "peaks_fwhm:\n");
		TRACE(TRACE_DEBUG, aCtx->result.peaksSnr,
		         ARRAY_KIND_U16, PEAKS_NUM, "peaks_snr:\n");
	}

	return ret;
}