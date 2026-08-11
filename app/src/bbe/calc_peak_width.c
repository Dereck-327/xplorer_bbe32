/*
 * calc_peak_width.c
 *
 *  Created on: 2025.11.7
 *      Author: qirui
 */

//==========================[ Headers ]=========================
#include "calc_peak_width.h"

//==========================[ Function Implementations ]=========================
void calc_peak_width(
	int16_t *restrict input,
	const uint16_t samples,
	const int16_t peak_val,
	const uint16_t peak_idx,
	int16_t *restrict peak_left_width,
	int16_t *restrict peak_right_width)
{
	#pragma aligned(input, 32)

	int16_t peak_min, peak_half;
	uint16_t peak_left_idx, peak_right_idx;

	xb_vecNx16 peak_left_16v, peak_right_16v, peak_min_16v, peak_half_16v, zeros_16v, seq_16v, ones_16v, left_seq_16v, right_seq_16v;
	valign peak_left_valign, peak_right_valign;
	vboolN peak_left_vboolN, peak_right_vboolN;

	peak_left_idx = peak_idx - (uint16_t) BBE_SIMD_WIDTH;
	peak_right_idx = peak_idx + 1;

	xb_vecNx16 *restrict p_peak_left = (xb_vecNx16 *restrict) &(input[peak_left_idx]);
	xb_vecNx16 *restrict p_peak_right = (xb_vecNx16 *restrict) &(input[peak_right_idx]);

	peak_left_valign = BBE_LANX16_PP(p_peak_left);
	peak_right_valign = BBE_LANX16_PP(p_peak_right);

	BBE_LANX16_IP(peak_left_16v, peak_left_valign, p_peak_left);
	BBE_LANX16_IP(peak_right_16v, peak_right_valign, p_peak_right);

	peak_min_16v = BBE_MINNX16(peak_left_16v, peak_right_16v);
	peak_min = BBE_RMINNX16(peak_min_16v);
	peak_half = (peak_val >> 1) + (peak_min >> 1);
	peak_half_16v = BBE_MOVVA16(peak_half);

	peak_left_vboolN = BBE_GENX16(peak_left_16v, peak_half_16v);
	peak_right_vboolN = BBE_GENX16(peak_right_16v, peak_half_16v);

	zeros_16v = BBE_ZERONX16();
	seq_16v = BBE_SEQNX16();
	ones_16v = BBE_MOVVA16(1);
	right_seq_16v = BBE_ADDNX16(seq_16v, ones_16v);
	left_seq_16v = BBE_SHFLNX16I(right_seq_16v, BBE_SHFLI_REVERSE_1);

	peak_left_16v = BBE_MOVNX16T(left_seq_16v, zeros_16v, peak_left_vboolN);
	peak_right_16v = BBE_MOVNX16T(right_seq_16v, zeros_16v, peak_right_vboolN);

	*peak_left_width = BBE_RMAXNX16(peak_left_16v);
	*peak_right_width = BBE_RMAXNX16(peak_right_16v);
}
