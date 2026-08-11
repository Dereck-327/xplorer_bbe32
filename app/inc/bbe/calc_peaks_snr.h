/*
 * calc_peaks_snr.h
 *
 *  Created on: 2025.11.11
 *      Author: qirui
 */

#ifndef CALC_PEAKS_SNR_H_
#define CALC_PEAKS_SNR_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"

//==========================[ Function Prototypes ]=========================
void calc_peaks_snr(
	const int16_t *restrict input,
	const uint16_t samples,
	const int16_t *restrict peaks_val,
	uint16_t *restrict peaks_snr,
	const uint8_t peaks_num,
	const uint8_t ra_data_bexp,
	const int16_t power_shift);

#endif /* CALC_PEAKS_SNR_H_ */
