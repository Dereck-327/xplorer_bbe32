/**
    @file       stage_prepare.c
    @brief      拷贝频谱 + 算本 chirp 的块指数
                每个 chirp 都跑

    @version    0.1.0
    @date       2026-08-10
    @author     hjk

*/

#include "simulate/pipeline.h"
#include "simulate/dsp_ops.h"

ErrorType StagePrepare(Bbe_CtxType *const aCtx)
{
    ErrorType ret = ERR_OK;

    if (NULL == aCtx)
    {
        ret = ERR_INVAL_PARAMS;
    }
    else if (aCtx->chirpIdx >= (uint8_t) CHIRP_NUM)
    {
        ret = ERR_INVAL_STATE;		/* 越界写 work.mag2[] 的唯一入口, 挡在这里 */
    }
    else
    {
        DSP_PrepareData(aCtx->frame.mag2,
                        aCtx->work.mag2[aCtx->chirpIdx],
                        MAG2_SIZE);

        /* 定点谱的定标随帧变化: 32 位满量程减去 FFT 与取模平方各自的块指数。
         * 乘 2 是因为取模平方让指数翻倍。 */
        aCtx->work.mag2Bexp[aCtx->chirpIdx] = (int16_t) (32U -
            ((uint16_t) aCtx->frame.dataBexp << 1U) -
            ((uint16_t) aCtx->frame.magsqBexp << 1U));
    }

    return ret;
}