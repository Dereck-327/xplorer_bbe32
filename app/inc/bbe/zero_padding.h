/*
 * zero_padding.h
 *
 *  Created on: 2025.11.4
 *      Author: qirui
 */

#ifndef ZERO_PADDING_H_
#define ZERO_PADDING_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"

//==========================[ Function Prototypes ]=========================
void zero_padding(
	int16_t *restrict input,
	const uint16_t input_samples,
	const uint16_t output_samples);

#endif /* ZERO_PADDING_H_ */
