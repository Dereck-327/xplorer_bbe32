/*
 * copyConvRes_S16_S16_v.c
 *
 *  Created on: 2026.2.27
 *      Author: qirui
 */

//==========================[ Headers ]=========================
#include "copyConvRes_S16_S16_v.h"

//==========================[ Function Implementations ]=========================
void copyConvRes_S16_S16_v(
	const int16_t *restrict input,
	int16_t *restrict output,
	const uint16_t samples)
{
	#pragma aligned(input, 32)
    #pragma aligned(output, 32)

	uint16_t i;
	xb_vecNx16 input_16v;

	const xb_vecNx16 *restrict p_input = (xb_vecNx16 *restrict) input;
	xb_vecNx16 *restrict p_output = (xb_vecNx16 *restrict) output;

	for (i = 0; i < (samples / BBE_SIMD_WIDTH); i++)
	{
		// Load input
		BBE_LVNX16_IP(input_16v, p_input, 2 * BBE_SIMD_WIDTH);

		// Store input
		BBE_SVNX16_IP(input_16v, p_output, 2 * BBE_SIMD_WIDTH);
	}
}
