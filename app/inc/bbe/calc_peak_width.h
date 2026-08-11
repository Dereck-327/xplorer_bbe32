/*
 * calc_peak_width.h
 *
 *  Created on: 2025.11.7
 *      Author: qirui
 */

#ifndef CALC_PEAK_WIDTH_H_
#define CALC_PEAK_WIDTH_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"

//==========================[ Function Prototypes ]=========================
void calc_peak_width(
	int16_t *restrict input,
	const uint16_t samples,
	const int16_t peak_val,
	const uint16_t peak_idx,
	int16_t *restrict left_width,
	int16_t *restrict right_width);

#endif /* CALC_PEAK_WIDTH_H_ */
