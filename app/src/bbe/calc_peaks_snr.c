/*
 * calc_peaks_snr.c
 *
 *  Created on: 2025.11.11
 *      Author: qirui
 */

//==========================[ Headers ]=========================
#include "calc_peaks_snr.h"

//==========================[ Macros ]=========================
#define LIMIT_POWER_SHIFT(sh) ((sh) > 8 ? 8 : (sh))
#define SATURATE_U16(val) ((val) > UINT16_MAX ? UINT16_MAX : (uint16_t)(val))

//==========================[ Function Implementations ]=========================
void calc_peaks_snr(
	const int16_t *restrict input,
	const uint16_t samples,
	const int16_t *restrict peaks_val,
	uint16_t *restrict peaks_snr,
	const uint8_t peaks_num,
	const uint8_t ra_data_bexp,
	const int16_t power_shift)
{
	int16_t zeros_samples, zeros_total = 0, power_shift_correct;
	uint16_t idx, cycles, divisor, peaks_ns;
	uint32_t snr_calc;

	xb_vecNx16 input_16v, ones_16v, zeros_16v, mask_16v;
	xb_vecNx40 input_40v, sum_40v;
	vboolN mask_vboolN;

    const xb_vecNx16 *restrict p_input = (const xb_vecNx16 *restrict) input;

	cycles = samples / BBE_SIMD_WIDTH;
	ones_16v = BBE_MOVVA16(1);
	zeros_16v = BBE_ZERONX16();
	sum_40v = BBE_ZERONX40();
	for (idx = 0; idx < cycles; idx++)
	{
		BBE_LVNX16_IP(input_16v, p_input, 2 * BBE_SIMD_WIDTH);
		mask_vboolN = BBE_GENX16(input_16v, ones_16v);
		mask_16v = BBE_MOVNX16T(zeros_16v, ones_16v, mask_vboolN);
		zeros_samples = BBE_RADDNX16(mask_16v);
		zeros_total += zeros_samples;
		input_40v = BBE_UNPKSNX16(input_16v);
		sum_40v = BBE_ADDNX40(sum_40v, input_40v);
	}
	sum_40v = BBE_RADDNX40(sum_40v);

	if (!ra_data_bexp)
	{
	    for (idx = 0; idx < peaks_num; idx++)
	        peaks_snr[idx] = UINT16_MAX;
	    return;
	}
	power_shift_correct = LIMIT_POWER_SHIFT(power_shift);
	divisor = samples - zeros_total;
	peaks_ns = (divisor == 0) ? 1 : (uint32_t)(((int64_t)sum_40v << 4) / divisor);
	for (idx = 0; idx < peaks_num; idx++)
	{
	    snr_calc = ((uint32_t)peaks_val[idx] << (power_shift_correct + 4)) / (peaks_ns + (1U << 4));
	    peaks_snr[idx] = SATURATE_U16(snr_calc);
	}
}
