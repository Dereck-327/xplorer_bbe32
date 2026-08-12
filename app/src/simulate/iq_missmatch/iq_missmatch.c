#include "simulate/iq_missmatch/iq_missmatch.h"
#include "utils/bbe_type.h"
#include "utils/compiler.h"
#include "system.h"			/* TXT_VAL_MAX */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* 交织存储 [re, im, re, im, ...], 布局同 complex_fract16, 每路 2*MAG2_SIZE 个 int16 */
static int16_t g_i_f[2U * MAG2_SIZE] COMP_ALIGN(VEC_ALIGN);
static int16_t g_q_f[2U * MAG2_SIZE] COMP_ALIGN(VEC_ALIGN);
static char    g_line[TXT_VAL_MAX];

/* 逗号分隔整数 -> aDst, 返回解析到的个数 */
static uint32_t parse_i16_csv(char *aVal, int16_t *const aDst, const uint32_t aMax)
{
	uint32_t n = 0U;
	char *tok = strtok(aVal, ",");

	while ((NULL != tok) && (n < aMax))
	{
		aDst[n] = (int16_t) strtol(tok, NULL, 10);
		++n;
		tok = strtok(NULL, ",");
	}

	return n;
}

/* 校正一个 bin: Q_corr = gainCorr*(Q - rho*I); sig = I + j*Q_corr; 返回 |sig|^2
 * j*Q_corr 让实虚部交换: sig_re = I_re - Qc_im, sig_im = I_im + Qc_re */
static float bin_mag2(const uint32_t k, const float aRho, const float aGainCorr)
{
	float iRe = (float) g_i_f[2U * k];
	float iIm = (float) g_i_f[2U * k + 1U];
	float qRe = (float) g_q_f[2U * k];
	float qIm = (float) g_q_f[2U * k + 1U];

	float qcRe = aGainCorr * (qRe - (aRho * iRe));
	float qcIm = aGainCorr * (qIm - (aRho * iIm));

	float sRe = iRe - qcIm;
	float sIm = iIm + qcRe;

	return (sRe * sRe) + (sIm * sIm);
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
	float  peakP = 0.0F;
	int    magsqBexp;
	float  invScale;
	uint32_t k;

	if ((NULL == aPath) || (NULL == aFrame))
	{
		return ERR_INVAL_PARAMS;
	}

	fp = fopen(aPath, "r");
	if (NULL == fp)
	{
		return ERR_NOT_FOUND;
	}

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
			if ((2U * MAG2_SIZE) == parse_i16_csv(eq + 1, g_i_f, 2U * MAG2_SIZE))
			{
				seen |= 0x01U;
			}
		}
		else if (0 == strcmp(g_line, "q_f"))
		{
			if ((2U * MAG2_SIZE) == parse_i16_csv(eq + 1, g_q_f, 2U * MAG2_SIZE))
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

	if (0x07U != seen)
	{
		return ERR_DATA_INTEG;		/* i_f / q_f / fft_bexp 三样缺一不可 */
	}

	/* 第一趟: 找校正后功率峰, 定块指数 */
	for (k = 0U; k < MAG2_SIZE; ++k)
	{
		float p = bin_mag2(k, aRho, aGainCorr);
		if (p > peakP)
		{
			peakP = p;
		}
	}

	/* U(16,15): mag2 = round(P / 2^(15 + 2*magsqBexp)) 需落进 uint16。
	 * 解 2^(15+2*b) >= peakP/65535 -> b = ceil((log2(peakP) - 31)/2)。dataBexp 用 fft_bexp。 */
	magsqBexp = 0;
	if (peakP > 0.0F)
	{
		magsqBexp = (int) ceilf((log2f(peakP) - 31.0F) / 2.0F);
		if (magsqBexp < 0)
		{
			magsqBexp = 0;
		}
	}

	/* 第二趟: 量化写出 */
	invScale = ldexpf(1.0F, -(15 + (2 * magsqBexp)));	/* 1 / 2^(15+2*magsqBexp) */
	for (k = 0U; k < MAG2_SIZE; ++k)
	{
		float m = bin_mag2(k, aRho, aGainCorr) * invScale;
		long  v = lroundf(m);

		if (v < 0L)      { v = 0L; }
		if (v > 65535L)  { v = 65535L; }
		aFrame->mag2[k] = (uint16_t) v;
	}

	aFrame->dataBexp  = (uint8_t) fftBexp;
	aFrame->magsqBexp = (uint8_t) magsqBexp;

	return ERR_OK;
}