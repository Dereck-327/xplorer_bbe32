#include "simulate/iq_missmatch/iq_missmatch.h"
#include "bbe/iq_correct_mag2.h"
#include "utils/bbe_type.h"
#include "utils/compiler.h"
#include "system.h"			/* TXT_VAL_MAX */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils/perf_stat.h"

/* 交织存储 [re, im, re, im, ...], 布局同 complex_fract16, 每路 2*MAG2_SIZE 个 int16 */
// static int16_t g_i_f[2U * MAG2_SIZE] COMP_ALIGN(VEC_ALIGN);
// static int16_t g_q_f[2U * MAG2_SIZE] COMP_ALIGN(VEC_ALIGN);
static complex_fract16 g_i_f[MAG2_SIZE] COMP_ALIGN(VEC_ALIGN);
static complex_fract16 g_q_f[MAG2_SIZE] COMP_ALIGN(VEC_ALIGN);

static char    g_line[TXT_VAL_MAX];

/* 逗号分隔整数 -> aDst, 返回解析到的个数 */
static uint32_t parse_i16_csv(char *aVal, complex_fract16 *const aDst, const uint32_t aMax)
{
	uint32_t n = 0U;
	char *tok = strtok(aVal, ",");

	while ((NULL != tok) && (n < aMax))
	{
		aDst[n].s.re = (int16_t) strtol(tok, NULL, 10);
        tok = strtok(NULL, ",");
        if (NULL == tok) { break; }  // 缺少匹配的叙述
		aDst[n].s.im = (int16_t) strtol(tok, NULL, 10);
		tok = strtok(NULL, ",");
		++n;
	}

	return n;
}

/*
 * IQ missmatch data file
 * i_f(交织 2*MAG2_SIZE int16): re, im, re, im, ...
 * q_f(交织 2*MAG2_SIZE int16): re, im, re, im, ...
 * fft_bexp : real value = value / (2 ** 15) * (2 ** bexp)

 * aRho 默认值 0.2325524544260392
 * aGainCorr 默认值 1.0897237313241515
*/
ErrorType iqMissmatch_frame_load(const char *const aPath,
                                 const float aRho,
                                 const float aGainCorr,
                                 Bbe_FrameType *const aFrame)
{
	FILE *fp = NULL;
	char *eq;
	uint8_t seen = 0U;			/* bit0 i_f, bit1 q_f, bit2 fft_bexp */
	int    fftBexp = 0;
	int16_t rhoQ15, gainQ14;

	if ((NULL == aPath) || (NULL == aFrame))
	{
		return ERR_INVAL_PARAMS;
	}

	fp = fopen(aPath, "r");
	if (NULL == fp)
	{
		return ERR_NOT_FOUND;
	}
    PERF_BEGIN(chirp_read);
	while (NULL != fgets(g_line, (int) sizeof(g_line), fp))
	{
		eq = strchr(g_line, '=');
		if (NULL == eq)
		{
			continue;
		}
		*eq = '\0';

		if (0 == strcmp(g_line, "i_f"))
		{
			if (MAG2_SIZE == parse_i16_csv(eq + 1, g_i_f, MAG2_SIZE))
			{
				seen |= 0x01U;
			}
		}
		else if (0 == strcmp(g_line, "q_f"))
		{
			if (MAG2_SIZE == parse_i16_csv(eq + 1, g_q_f, MAG2_SIZE))
			{
				seen |= 0x02U;
			}
		}
		else if (0 == strcmp(g_line, "fft_bexp"))
		{
			fftBexp = (int) strtol(eq + 1, NULL, 10);
			seen |= 0x04U;
		}
		else
		{
			/* 未知键忽略 */
		}
	}

	(void) fclose(fp);
    PERF_END(chirp_read);

	if (0x07U != seen)
	{
        printf("seen error :: 0X%0x\n", seen);
		return ERR_DATA_INTEG;		/* i_f / q_f / fft_bexp 三样缺一不可 */
	}

	/* 系数转定点: rho Q15, gainCorr Q14 (见 iq_correct_mag2.h) */
	rhoQ15  = (int16_t) lroundf(aRho * (float) (1 << IQ_RHO_Q));
	gainQ14 = (int16_t) lroundf(aGainCorr * (float) (1 << IQ_GAIN_Q));

    PERF_BEGIN(iq_correct_mag2);
	/* 全定点校正 + 取模平方 -> U(16,15) 功率谱
	 * g_i_f/g_q_f 是 int16_t[], 但布局与 complex_fract16[] 等价 (交织 [re,im,...]) */
	iq_correct_mag2((const complex_fract16 *) g_i_f,
	                (const complex_fract16 *) g_q_f,
	                rhoQ15, gainQ14,
	                aFrame->mag2, &aFrame->magsqBexp, (uint16_t) MAG2_SIZE);

	aFrame->dataBexp = (uint8_t) fftBexp;
    PERF_END(iq_correct_mag2);
	return ERR_OK;
}