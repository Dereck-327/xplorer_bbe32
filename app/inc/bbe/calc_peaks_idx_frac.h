/*
 * calc_peaks_idx_frac.h
 *
 *  Created on: 2026.1.28
 *      Author: qirui
 */

#ifndef CALC_PEAKS_IDX_FRAC_H_
#define CALC_PEAKS_IDX_FRAC_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"

//==========================[ Function Prototypes ]=========================
void calc_peaks_idx_frac(
	const uint16_t *restrict peaks_idx_buffer,
	const int16_t *restrict peaks_frac_buffer,
	const uint8_t peaks_num,
	const uint16_t idx_corr,
	const uint8_t iq_reverse,
	int16_t *restrict peaks_idx,
	int16_t *restrict peaks_frac);

#endif /* CALC_PEAKS_IDX_FRAC_H_ */
