/*
 * CONV_S16_S16_VW_SIZE15_H_.h
 *
 *  Created on: 2025.11.3
 *      Author: qirui
 */

#ifndef CONV_S16_S16_VW_SIZE15_H_
#define CONV_S16_S16_VW_SIZE15_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"

//==========================[ Function Prototypes ]=========================
void conv_S16_S16_vw_size15(
	const int16_t *restrict input,
    int16_t *restrict output,
    const uint16_t samples,
    const int16_t *restrict mf_coeff);

#endif /* CONV_S16_S16_VW_SIZE15_H_ */
