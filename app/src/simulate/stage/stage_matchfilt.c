/**
    @file       stage_matchfilt.c
    @brief      15 tap 匹配滤波
                *   1. pad_reset() 的尾部清零打在 [RA_MAG2_SIZE, +RA_PADDING) 上, 而有效数据占
                *      [RA_PADDING, RA_PADDING + RA_MAG2_SIZE)。也就是清掉了谱最后 16 个真实 bin,
                *      而真正的尾部 padding [RA_MAG2_SIZE_PAD - RA_PADDING, RA_MAG2_SIZE_PAD)
                *      始终没被清, conv 读到的是残留值。少加一个 RA_PADDING 偏移所致。
                *   2. pad_reset() 的头部清零调 DSP_ZeroPad(buf, 0, RA_PADDING), 而 DSP_ZeroPad 内部
                *      从 &input[input_samples - last_samples] 即 &buf[-16] 起读写一个向量 —— 负下标。
                *      读回来原样写回, 效果无害, 但确实越界。
                *   3. DSP_ConvS16Vw15 取 samples = RA_MAG2_SIZE_PAD 时按 258 个向量步进, 输入最远读到
                *      下标 4150、输出最远写到 4143, 两者都超出 4128 长的缓冲。输出溢出落在紧随其后的
                *      peakIdx/peakLeftW 等数组上, 而它们在 RA_StagePeak 里会被整体重写, 所以现在看不出
                *      问题 —— 字段顺序一旦调整就会炸。给 conv 传 RA_MAG2_SIZE 才是它真正该处理的长度。

    @version    0.1.0
    @date       2026-08-10
    @author     hjk

*/

#include "simulate/pipeline.h"
#include "simulate/dsp_ops.h"
#include "utils/debug_trace.h"

//==========================[ Static Function Prototypes ]=========================
static void pad_reset(int16_t *const aBuf);

//==========================[ Static Function Implementations ]=========================
/* 原 reset_acc_buffers()。区间照抄, 见文件头 ponytail 第 1、2 条。 */
static void pad_reset(int16_t *const aBuf)
{
	DSP_ZeroPad(aBuf, 0U, PADDING);
	DSP_ZeroPad(aBuf, MAG2_SIZE, MAG2_SIZE + PADDING);
}

//==========================[ Function Implementations ]=========================
ErrorType StageMatchFilt(Bbe_CtxType *const aCtx)
{
	ErrorType ret = ERR_OK;

	if (NULL == aCtx)
	{
		ret = ERR_INVAL_PARAMS;
	}
	else
	{
		pad_reset(aCtx->work.accEq);
		pad_reset(aCtx->work.accEqMf);

		DSP_ConvS16Vw15(aCtx->work.accEq,
		                &aCtx->work.accEqMf[PADDING],
		                MAG2_SIZE_PAD,
		                aCtx->param.mfCoeff);

		TRACE(&aCtx->work.accEqMf[PADDING],
		         ARRAY_KIND_S16, MAG2_SIZE, "mag2_acc_eq_mf:\n");

		/* 滤波结果回灌 accEq: 后面 RA_StagePeak 的遮蔽会就地破坏 accEqMf,
		 * 而 RA_StageRefine 做亚 bin 插值时需要一份没被遮蔽的滤波结果。 */
		DSP_CopyConvRes(aCtx->work.accEqMf,
		                aCtx->work.accEq,
		                MAG2_SIZE_PAD);
	}

	return ret;
}