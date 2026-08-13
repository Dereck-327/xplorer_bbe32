/**
    @file       stage_integrate.c
    @brief      非相干累积 + fft_eq

    @version    0.1.0
    @date       2026-08-10
    @author     hjk

*/


#include "simulate/pipeline.h"
#include "utils/debug_trace.h"
#include "simulate/cfg.h"
#include "simulate/dsp_ops.h"

//==========================[ Function Implementations ]=========================

ErrorType StageIntegrate(Bbe_CtxType *const aCtx)
{
    ErrorType ret = ERR_OK;
	if (NULL == aCtx)
	{
		ret = ERR_INVAL_PARAMS;
	}
	else
	{
		/* 四路 mag2 按最大块指数对齐后相加 -> acc S(16,14) */
		DSP_NciFourFft(&aCtx->work.mag2[0][0],
		               aCtx->work.mag2Bexp,
		               &aCtx->work.acc[PADDING],
		               &aCtx->work.accEqBexp,
		               MAG2_SIZE);

		TRACE(&aCtx->work.acc[PADDING],
		         ARRAY_KIND_S16, MAG2_SIZE, "mag2_acc:\n");
		TRACE(&aCtx->work.accEqBexp,
		         ARRAY_KIND_S16, 1U, "mag2_acc_eq_bexp:\n");

		/* 均衡: (acc S(16,14) * fftEq U(16,15)) >> 15 -> accEq S(16,14) */
		DSP_MulS16U16(&aCtx->work.acc[PADDING],
		              aCtx->param.fftEq,
		              &aCtx->work.accEq[PADDING],
		              EQ_SHIFT,
		              MAG2_SIZE);

		TRACE(&aCtx->work.accEq[PADDING],
		         ARRAY_KIND_S16, MAG2_SIZE, "mag2_acc_eq:\n");
	}


    return ret;
}