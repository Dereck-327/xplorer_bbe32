/*
 * blank_neighbours.h
 *
 *  Created on: 2025.11.7
 *      Author: qirui
 */

#ifndef BLANK_NEIGHBOURS_H_
#define BLANK_NEIGHBOURS_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"

//==========================[ Function Prototypes ]=========================
void blank_neighbours(
	int16_t *restrict input,
	const uint16_t peak_idx,
	const int16_t peak_left_width,
	const int16_t peakright_width);

#endif /* BLANK_NEIGHBOURS_H_ */
