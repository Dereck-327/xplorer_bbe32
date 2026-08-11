/*
 * peak_idx_comp.h
 *
 *  Created on: 2025.11.13
 *      Author: qirui
 */

#ifndef PEAK_IDX_COMP_H_
#define PEAK_IDX_COMP_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"

//==========================[ Function Prototypes ]=========================
void peak_idx_comp(
	const int16_t *restrict input,
	const uint16_t peak_idx,
	const int16_t peak_left_width,
	const int16_t peak_right_width,
	const int16_t scale_comp,
	int16_t *restrict peak_frac);

#endif /* PEAK_IDX_COMP_H_ */
