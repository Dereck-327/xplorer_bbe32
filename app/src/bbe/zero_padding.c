/*
 * zero_padding.c
 *
 *  Created on: 2025.11.4
 *      Author: qirui
 */

//==========================[ Headers ]=========================
#include "zero_padding.h"

//==========================[ Function Implementations ]=========================
void zero_padding(
	int16_t *restrict input,
	const uint16_t input_samples,
	const uint16_t output_samples)
{
	#pragma aligned(input, 32)

	uint16_t last_samples, zero_padding_len, i;

	xb_vecNx16 zeros_16v, val_16v, aux_16v;

	vboolN comp_vboolN;

	valign load_valign, store_valign;

	xb_vecNx16 *restrict p_input;

	vselN select_vsel;

	zeros_16v = BBE_ZERONX16();

	zero_padding_len = output_samples - input_samples;

	// Calculate select_vsel to get last samples of input
	last_samples = 16 - (zero_padding_len % BBE_SIMD_WIDTH);

	aux_16v = BBE_SEQNX16();

	val_16v = BBE_MOVVA16(last_samples);

	comp_vboolN = BBE_GTNX16(val_16v, aux_16v);

	val_16v = BBE_MOVVA16(16);

	select_vsel = BBE_MOVNX16T(aux_16v, val_16v, comp_vboolN);

	// Load last input vector
	p_input = (xb_vecNx16 *restrict) &input[input_samples - last_samples];

	load_valign = BBE_LANX16_PP(p_input);

	BBE_LANX16_IP(val_16v, load_valign, p_input);

	// Select vector to write
	val_16v = BBE_SELNX16(zeros_16v, val_16v, select_vsel);

	// Reset position
	p_input = (xb_vecNx16 *restrict) &input[input_samples - last_samples];

	// Store alignment
	store_valign = BBE_ZALIGN();

	// Store
	BBE_SAVNX16_XP(val_16v, store_valign, p_input, 2 * BBE_SIMD_WIDTH);

	for (i = 0; i < (zero_padding_len / BBE_SIMD_WIDTH); i++)
	{
		// Store zeros
		BBE_SAVNX16_XP(zeros_16v, store_valign, p_input, 2 * BBE_SIMD_WIDTH);
	}

	// Flush store alignment
	BBE_SAPOS_FP(store_valign, p_input);
}
