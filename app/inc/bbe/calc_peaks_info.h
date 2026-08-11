/*
 * calc_peaks_info.h
 *
 *  Created on: 2025.11.7
 *      Author: qirui
 */

#ifndef CALC_PEAKS_INFO_H_
#define CALC_PEAKS_INFO_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"

//==========================[ Function Prototypes ]=========================
void calc_peaks_info(
	int16_t *restrict input,
	const uint16_t samples,
	int16_t *restrict peaks_left_width,
	int16_t *restrict peaks_right_width,
	const uint8_t peaks_num,
	uint16_t *restrict peaks_idx_buffer,
	const uint8_t ra_data_bexp,
	const int16_t power_shift,
	uint16_t *restrict peaks_fwhm,
	uint16_t *restrict peaks_snr);

#endif /* CALC_PEAKS_INFO_H_ */
