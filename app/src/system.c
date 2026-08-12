#include "system.h"
#include "utils/bbe_type.h"
#include "simulate/pipeline.h"

#include <stdio.h>
#include <time.h>

#define DATA_DIR        "data"
//==========================[ Data Structures ]=========================
typedef struct sTxtSrcType
{
	FILE	*fp;
	char	 key[TXT_KEY_MAX];
	char	 val[TXT_VAL_MAX];
} TxtSrcType;


//==========================[ Static Variables Definition ]=========================
BBE_CTX_DEFINE(static, g_ctx);

static ParamType_t g_param COMP_ALIGN(VEC_ALIGN);
static uint16_t g_mag2[MAG2_SIZE] COMP_ALIGN(VEC_ALIGN);
static TxtSrcType g_src;
static char g_path[PATH_MAX_LEN];

//==========================[ Static Function Prototypes ]=========================
static ErrorType txt_next(void *const aSrcCtx, const char **const aKey, char **const aVal);

static ErrorType frame_load(const char *const aPath,
                               uint16_t *const aMag2,
                               uint8_t *const aDataBexp,
                               uint8_t *const aMagsqBexp);

static uint16_t file_count(const char *const aDir);

//==========================[ Static Function Implementations ]=========================

static ErrorType txt_next(void *const aSrcCtx, const char **const aKey, char **const aVal)
{
	ErrorType ret = ERR_NOT_FOUND;
	TxtSrcType *src = (TxtSrcType *) aSrcCtx;
	char *eq;
	size_t len;

	while (NULL != fgets(src->val, TXT_VAL_MAX, src->fp))
	{
		eq = strchr(src->val, '=');

		if (NULL == eq)
		{
			continue;					/* 空行或注释, 跳过 */
		}

		len = (size_t) (eq - src->val);

		if (len >= TXT_KEY_MAX)
		{
            printf("txt_next: key too long: %zu\n", len);
			ret = ERR_DATA_INTEG;		/* 键名过长 */
			break;
		}

		Mem_Copy(src->key, src->val, (uint32_t) len);
		src->key[len] = '\0';

		*aKey = src->key;
		*aVal = eq + 1;
		ret = ERR_OK;
		break;
	}

	return ret;
}

/* 一个 ra_data_N.txt 就是一个 chirp: 谱 + 两个块指数 */
static ErrorType frame_load(const char *const aPath,
                               uint16_t *const aMag2,
                               uint8_t *const aDataBexp,
                               uint8_t *const aMagsqBexp)
{
	ErrorType ret = ERR_OK;
	ErrorType srcRet;
	const char *key = NULL;
	char *val = NULL;
	char *tok;
	uint16_t n;
	uint8_t seen = 0U;

	g_src.fp = fopen(aPath, "r");

	if (NULL == g_src.fp)
	{
		ret = ERR_NOT_FOUND;
	}
	else
	{
		for (;;)
		{
			srcRet = txt_next(&g_src, &key, &val);

			if (ERR_NOT_FOUND == srcRet)
			{
				break;
			}
			else if (ERR_OK != srcRet)
			{
				ret = ERR_DATA_INTEG;
				break;
			}
			else if (0 == strcmp(key, "ra_data"))
			{
				n = 0U;
				tok = strtok(val, ",");

				while ((NULL != tok) && (n < MAG2_SIZE))
				{
					aMag2[n] = (uint16_t) strtoul(tok, NULL, 10);
					++n;
					tok = strtok(NULL, ",");
				}

				if (MAG2_SIZE != n)
				{
					ret = ERR_DATA_INTEG;
					break;
				}

				++seen;
			}
			else if (0 == strcmp(key, "ra_data_bexp"))
			{
				*aDataBexp = (uint8_t) strtoul(val, NULL, 10);
				++seen;
			}
			else if (0 == strcmp(key, "ra_data_magsq_bexp"))
			{
				*aMagsqBexp = (uint8_t) strtoul(val, NULL, 10);
				++seen;
			}
			else
			{
				continue;				/* 未知键忽略 */
			}
		}

		(void) fclose(g_src.fp);
		g_src.fp = NULL;

		if ((ERR_OK == ret) && (FRAME_KEY_NUM != seen))
		{
			ret = ERR_DATA_INTEG;		/* 三样缺一不可 */
		}
	}

	return ret;
}

static uint16_t file_count(const char *const aDir)
{
	uint16_t n = 0U;
	FILE *fp;

	for (;;)
	{
		(void) snprintf(g_path, sizeof(g_path), "%s/ra_data_%u.txt", aDir, (unsigned) n);

		fp = fopen(g_path, "r");

		if (NULL == fp)
		{
			break;
		}

		(void) fclose(fp);
		++n;
	}

	return n;
}



ErrorType system_init(void)
{
    ErrorType ret = ERR_OK;
#if TRACE_ENABLED
	if (ERR_OK != ParamSelfTest())
	{
		printf("RA_ParamSelfTest FAILED\n");
		return ERR_DATA_INTEG;
	}
	if (ERR_OK != PipelineSelfTest())
	{
		printf("PipelineSelfTest FAILED\n");
		return ERR_DATA_INTEG;
	}
#endif
    do
    {
        (void) snprintf(g_path, sizeof(g_path), "%s/params.txt", DATA_DIR);
        g_src.fp = fopen(g_path, "r");

        if (NULL == g_src.fp)
        {
            ret = ERR_NOT_FOUND;
            break;
        }

        ret = ParamLoad(&g_param, txt_next, &g_src);
        (void) fclose(g_src.fp);
        g_src.fp = NULL;

        if (ERR_OK != ret)
        {
            printf("ParamLoad failed: %d\n", (int) ret);
            break;
        }

        (void) ParamDump(&g_param, "System Parameters");
        ret = Pipeline_Init(&g_ctx, &g_param);
        if (ret != ERR_OK) {
            printf("Pipeline_Init failed: %d\n", (int) ret);
            break;
        }

    } while (0);

    return ret;
}


void system_run(void)
{
    ErrorType ret = ERR_OK;
    Bbe_ResultType out;
    clock_t t0, t1;
    double sec;
    uint16_t i = 0U;
    uint8_t dataBexp = 0U;
    uint8_t magsqBexp = 0U;
    uint16_t frames = 0U;
    uint16_t total = file_count(DATA_DIR);
    printf("Total chirps to process: %u\n", (unsigned) total);
    t0 = clock();
    while (i < total)
    {
        (void) snprintf(g_path, sizeof(g_path), "%s/ra_data_%u.txt", DATA_DIR, (unsigned) i);
        ret = frame_load(g_path, g_mag2, &dataBexp, &magsqBexp);
        if (ret != ERR_OK) {
            printf("frame_load failed: %d\n", (int) ret);
            break;
        }
        ret = Pipeline_SubmitFrame(&g_ctx, g_mag2, dataBexp, magsqBexp);
        if (ret != ERR_OK) {
            printf("Pipeline_SubmitFrame failed: %d\n", (int) ret);
            break;
        }
        ret = Pipeline_ProcessFrame(&g_ctx);
        if (ret == ERR_OK) {
			(void) Pipeline_GetResult(&g_ctx, &out);
			++frames;
        } else if (ret != ERR_BUSY) {
            printf("Pipeline_ProcessFrame failed: %d\n", (int) ret);
            break;
        } else {} /* chirp 还在积累 没到最后一个*/
        ++i;

    }
    t1 = clock();
    sec = (double)(t1 - t0) / CLOCKS_PER_SEC;
    printf("========================================\n");
    printf("Processed %u chirps in %.3f seconds\n", (unsigned) frames, sec);
	if (0U != i)
	{
		printf("Average time per chirp: %.6f seconds\n", sec / (double) i);
	}
    printf("========================================\n");
}

void system_deinit(void)
{
    printf("system over\n");
}
