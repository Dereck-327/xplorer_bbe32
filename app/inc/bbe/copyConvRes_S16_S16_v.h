/*
 * copyConvRes_S16_S16_v.h
 *
 *  Created on: 2026.2.27
 *      Author: qirui
 */

#ifndef COPYCONVRES_S16_S16_V_H_
#define COPYCONVRES_S16_S16_V_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"

//==========================[ Function Prototypes ]=========================
void copyConvRes_S16_S16_v(
	const int16_t *restrict input,
	int16_t *restrict output,
	const uint16_t samples);

#endif /* COPYCONVRES_S16_S16_V_H_ */
