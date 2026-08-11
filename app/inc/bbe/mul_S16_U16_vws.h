/*
 * mul_S16_U16_vws.h
 *
 *  Created on: 2026.3.21
 *      Author: qirui
 */

#ifndef MUL_S16_U16_VWS_H_
#define MUL_S16_U16_VWS_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"

//==========================[ Function Prototypes ]=========================
void mul_S16_U16_vws(
	const int16_t *restrict x,
	const uint16_t *restrict y,
	int16_t *restrict z,
	const int16_t rsh,
	const uint16_t samples);

#endif /* MUL_S16_U16_VWS_H_ */
