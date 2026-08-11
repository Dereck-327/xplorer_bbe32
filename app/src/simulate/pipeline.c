#include "simulate/pipeline.h"
#include "simulate/param.h"



//==========================[ Macros ]=========================
#define STAGE_NUM	((uint16_t) (sizeof(StageTbl) / sizeof(StageTbl[0])))


//==========================[ Static Variables Definition ]=========================
/* 链路即数据: 增删阶段 = 增删一行, 顺序在一处可读。
 * 下标必须与 pipeline.h 的 STAGE_* 一致 (也是 TRACE_MASK 的位号)。 */
static const StageDescType_t StageTbl[] =
{
	{ "prepare",	StagePrepare,	WHEN_EVERY_CHIRP	},
	{ "integrate",	StageIntegrate,	WHEN_LAST_CHIRP	},
	{ "matchfilt",	StageMatchFilt,	WHEN_LAST_CHIRP	},
	{ "peak",		StagePeak,		WHEN_LAST_CHIRP	},
	{ "refine",		StageRefine,	WHEN_LAST_CHIRP	},
	{ "emit",		StageEmit,		WHEN_LAST_CHIRP	},
};

COMP_STATIC_ASSERT(STAGE_NUM == PIPELINE_STAGE_NUM, stage_tbl_matches_indices);

//==========================[ Function Implementations ]=========================

ErrorType Pipeline_Init(Bbe_CtxType *const aCtx, const ParamType_t *const aParam)
{
    ErrorType ret = ERR_OK;

	if ((NULL == aCtx) || (NULL == aParam))
	{
		ret = ERR_INVAL_PARAMS;
	}
	else
	{
		ret = ParamCheck(aParam);		/* 非法参数在这里会失败, 不进链路 */

		if (ERR_OK == ret)
		{
			Mem_Set(aCtx, 0U, (uint32_t) sizeof(Bbe_CtxType));
			Mem_Copy(&aCtx->param, aParam, (uint32_t) sizeof(ParamType_t));

			aCtx->chirpIdx = 0U;
			aCtx->state = STATE_READY;
		}
	}
    return ret;
}

ErrorType Pipeline_SubmitFrame(Bbe_CtxType *const aCtx,
                               const uint16_t *const aMag2,
                               const uint8_t aDataBexp,
                               const uint8_t aMagsqBexp)
{
    ErrorType ret = ERR_OK;

	if ((NULL == aCtx) || (NULL == aMag2))
	{
		ret = ERR_INVAL_PARAMS;
	}
	else if (STATE_READY != aCtx->state)
	{
		ret = ERR_UNINIT;
	}
	else if (0U != (((uint32_t) (const void *) aMag2) % VEC_ALIGN))
	{
		ret = ERR_INVAL_PARAMS;			/* 未对齐的输入会让向量加载崩掉 */
	}
	else
	{
		Mem_Copy(aCtx->frame.mag2, aMag2, (uint32_t) sizeof(aCtx->frame.mag2));
		aCtx->frame.dataBexp = aDataBexp;
		aCtx->frame.magsqBexp = aMagsqBexp;
	}
    return ret;
}

ErrorType Pipeline_ProcessFrame(Bbe_CtxType *const aCtx)
{
    ErrorType ret = ERR_OK;
    const StageDescType_t *stage;
    uint8_t isLast = 0U;
    uint16_t i;
    if (NULL == aCtx)
    {
        ret = ERR_INVAL_PARAMS;
    } else if (STATE_READY != aCtx->state)
    {
        ret = ERR_UNINIT;
    }
    else
    {
        isLast = (aCtx->chirpIdx == (uint8_t) LAST_CHIRP) ? 1U : 0U;
        for (i = 0U; i < STAGE_NUM; ++i) {
            stage = &StageTbl[i];
            if ((WHEN_LAST_CHIRP == stage->when) && (0U == isLast))
            {
                continue;				/* 没攒满就跳过只在末尾chirp跑的阶段 */
            }

            ret = stage->fn(aCtx);
            if (ERR_OK != ret)
            {
                break;					/* 遇错即停, 错误码上抛 */
            }

        }

		if (ERR_OK == ret)
		{
			if (0U != isLast)
			{
				aCtx->chirpIdx = 0U;	/* 一帧完成, 开始下一帧的累积 */
			}
			else
			{
				++aCtx->chirpIdx;
				ret = ERR_BUSY;			/* 还要继续 submit */
			}
		}

    }
    return ret;
}

ErrorType Pipeline_GetResult(const Bbe_CtxType *const aCtx, Bbe_ResultType *const aOut)
{
	ErrorType ret = ERR_OK;

	if ((NULL == aCtx) || (NULL == aOut))
	{
		ret = ERR_INVAL_PARAMS;
	}
	else if (STATE_READY != aCtx->state)
	{
		ret = ERR_UNINIT;
	}
	else
	{
		Mem_Copy(aOut, &aCtx->result, (uint32_t) sizeof(Bbe_ResultType));
	}

	return ret;
}