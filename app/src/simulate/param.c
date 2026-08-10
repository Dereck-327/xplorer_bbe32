#include "simulate/param.h"
#include "utils/debug_trace.h"
#include "utils/compiler.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#if TRACE_ENABLED
    #include <stdio.h>
#endif


//==========================[ Macros ]=========================
#define OF(field)		((uint16_t) offsetof(ParamType_t, field))
#define DESC_NUM		((uint16_t) (sizeof(ParamDesc) / sizeof(ParamDesc[0])))

/* 各类型的取值范围 */

#define LIM_S16_MIN		(-32768)
#define LIM_S16_MAX		(32767)
#define LIM_U16_MAX		(65535)
#define LIM_S32_MIN		((int32_t) 0x80000000)
#define LIM_S32_MAX		((int32_t) 0x7FFFFFFF)

//==========================[ Static Variables Definition ]=========================
/*  name                 offset               kind                count               min          max
 *  ------------------------------------------------------------------------------------------------------------- */

static const ParamDescType_t ParamDesc[] =
{
	{ "fft_eq",           OF(fftEq),           PARAM_KIND_U16, MAG2_SIZE,           0,           LIM_U16_MAX          },
	{ "blank_table",      OF(blankTable),      PARAM_KIND_U16, MAX_BLANK_SIZE,      0,           MAG2_SIZE - 1U    },
	{ "mf_coeff",         OF(mfCoeff),         PARAM_KIND_S16, MF_COEFF_SIZE,       LIM_S16_MIN, LIM_S16_MAX          },
	{ "power_shift",      OF(powerShift),      PARAM_KIND_S16, 1U,                  0,           SNR_SHIFT_MAX     },
	{ "comp_width",       OF(compWidth),       PARAM_KIND_S16, 1U,                  0,           BBE_SIMD_WIDTH       },
	{ "scale_comp",       OF(scaleComp),       PARAM_KIND_S16, 1U,                  LIM_S16_MIN, LIM_S16_MAX          },
	{ "blank_range",      OF(blankRange),      PARAM_KIND_U16, 1U,                  0,           MAG2_SIZE / 2U    },
	{ "blank_valid_size", OF(blankValidSize),  PARAM_KIND_U8,  1U,                  0,           MAX_BLANK_SIZE    },
	{ "comp_mode",        OF(compMode),        PARAM_KIND_U8,  1U,                  0,           COMP_MODE_MAX     },
	{ "iq_reverse",       OF(iqReverse),       PARAM_KIND_U8,  1U,                  0,           1                    },
};

/* 表必须覆盖结构里每个字段, 漏一行就是"参数加载了但没校验"。
 * 字段数变了这里会编译不过, 提醒去补表。 */
COMP_STATIC_ASSERT(DESC_NUM == 10U, param_desc_covers_all_fields);

//==========================[ Static Function Prototypes ]=========================
static uint8_t elem_size(const ParamKindType aKind);
static int32_t elem_get(const void *const aBase, const ParamKindType aKind, const uint16_t aIdx);
static void elem_set(void *const aBase, const ParamKindType aKind, const uint16_t aIdx, const int32_t aVal);
static ErrorType desc_find(const char *const aKey, const ParamDescType_t **const aOut);
static ErrorType desc_store(const ParamDescType_t *const aDesc, void *const aField, char *const aVal);

//==========================[ Static Function Implementations ]=========================
static uint8_t elem_size(const ParamKindType aKind)
{
	uint8_t size;

	if (PARAM_KIND_U8 == aKind)
	{
		size = (uint8_t) sizeof(uint8_t);
	}
	else if (PARAM_KIND_S32 == aKind)
	{
		size = (uint8_t) sizeof(int32_t);
	}
	else
	{
		size = (uint8_t) sizeof(uint16_t);	/* S16 与 U16 同宽 */
	}

	return size;
}

static int32_t elem_get(const void *const aBase, const ParamKindType aKind, const uint16_t aIdx)
{
	int32_t val;

	if (PARAM_KIND_U8 == aKind)
	{
		val = (int32_t) ((const uint8_t *) aBase)[aIdx];
	}
	else if (PARAM_KIND_S16 == aKind)
	{
		val = (int32_t) ((const int16_t *) aBase)[aIdx];
	}
	else if (PARAM_KIND_U16 == aKind)
	{
		val = (int32_t) ((const uint16_t *) aBase)[aIdx];
	}
	else
	{
		val = ((const int32_t *) aBase)[aIdx];
	}

	return val;
}

static void elem_set(void *const aBase, const ParamKindType aKind, const uint16_t aIdx, const int32_t aVal)
{
	if (PARAM_KIND_U8 == aKind)
	{
		((uint8_t *) aBase)[aIdx] = (uint8_t) aVal;
	}
	else if (PARAM_KIND_S16 == aKind)
	{
		((int16_t *) aBase)[aIdx] = (int16_t) aVal;
	}
	else if (PARAM_KIND_U16 == aKind)
	{
		((uint16_t *) aBase)[aIdx] = (uint16_t) aVal;
	}
	else
	{
		((int32_t *) aBase)[aIdx] = aVal;
	}
}

static ErrorType desc_find(const char *const aKey, const ParamDescType_t **const aOut)
{
	ErrorType ret = ERR_NOT_FOUND;
	uint16_t i;

	for (i = 0U; i < DESC_NUM; ++i)
	{
		if (0 == strcmp(aKey, ParamDesc[i].name))
		{
			*aOut = &ParamDesc[i];
			ret = ERR_OK;
			break;
		}
	}

	return ret;
}

static ErrorType desc_store(const ParamDescType_t *const aDesc, void *const aField, char *const aVal)
{
	ErrorType ret = ERR_OK;
	char *tok;
	uint16_t n = 0U;
	int32_t v;

	tok = strtok(aVal, ",");
	while ((NULL != tok) && (n < aDesc->count))
	{
		v = (int32_t) strtol(tok, NULL, 10);

		if ((v < aDesc->min) || (v > aDesc->max))
		{
			ret = ERR_INVAL_PARAMS;		/* 越界的值一律不写入 */
			break;
		}

		elem_set(aField, aDesc->kind, n, v);
		++n;
		tok = strtok(NULL, ",");
	}

	if ((ERR_OK == ret) && (n != aDesc->count))
	{
		ret = ERR_DATA_INTEG;			/* 给的个数不对 */
	}

	return ret;
}

//==========================[ Function Implementations ]=========================

ErrorType ParamLoad(ParamType_t *const aParam,
                          ParamSrcFnType aSrcFn,
                          void *const aSrcCtx)
{
	ErrorType ret = ERR_OK;
	ErrorType srcRet;
	const ParamDescType_t *desc = NULL;
	const char *key = NULL;
	char *val = NULL;
	uint16_t seen = 0U;

	if ((NULL == aParam) || (NULL == aSrcFn))
	{
		ret = ERR_INVAL_PARAMS;
	}
	else
	{
		/* 先清零: 数据源没给的字段处于确定状态, 而不是上一帧的残留 */
		Mem_Set(aParam, 0U, (uint32_t) sizeof(ParamType_t));

		for (;;)
		{
			srcRet = aSrcFn(aSrcCtx, &key, &val);

			if (ERR_NOT_FOUND == srcRet)
			{
				break;					/* 数据源读完, 正常结束 */
			}
			else if (ERR_OK != srcRet)
			{
				ret = ERR_DATA_INTEG;
				break;
			}
			else if (ERR_OK != desc_find(key, &desc))
			{
				continue;				/* 表里没有的键忽略, 便于参数文件向前兼容 */
			}
			else
			{
				ret = desc_store(desc, (void *) ((uint8_t *) aParam + desc->offset), val);

				if (ERR_OK != ret)
				{
					break;
				}

				++seen;
			}
		}

		if ((ERR_OK == ret) && (seen != DESC_NUM))
		{
            printf("ParamLoad: seen=%d, DESC_NUM=%d\n", seen, DESC_NUM);
			ret = ERR_DATA_INTEG;		/* 有参数没在数据源里出现 */
		}

		if (ERR_OK == ret)
		{
			ret = ParamCheck(aParam);
		}
	}

	return ret;
}

ErrorType ParamCheck(const ParamType_t *const aParam)
{
	ErrorType ret = ERR_OK;
	const ParamDescType_t *desc;
	const void *field;
	uint16_t i;
	uint16_t n;
	int32_t v;

	if (NULL == aParam)
	{
		ret = ERR_INVAL_PARAMS;
	}
	else
	{
		/* 逐项范围校验 */
		for (i = 0U; (i < DESC_NUM) && (ERR_OK == ret); ++i)
		{
			desc = &ParamDesc[i];
			field = (const void *) ((const uint8_t *) aParam + desc->offset);

			for (n = 0U; n < desc->count; ++n)
			{
				v = elem_get(field, desc->kind, n);

				if ((v < desc->min) || (v > desc->max))
				{
					ret = ERR_INVAL_PARAMS;
					break;
				}
			}
		}

		/* 遮蔽半宽必须落在 padding 能兜住的范围内, 否则置零会写到缓冲之外 */
		if (ERR_OK == ret)
		{
			if (aParam->blankRange >= MAG2_SIZE)
			{
				ret = ERR_INVAL_PARAMS;
			}
		}
	}

	return ret;
}
