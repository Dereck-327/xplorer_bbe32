/*
 * peak_idx_comp.c
 *
 *  Created on: 2025.11.13
 *      Author: qirui
 */

//==========================[ Headers ]=========================
#include "peak_idx_comp.h"

//==========================[ Macros ]=========================
#define OUTPUT_SHIFT	(((uint16_t) 19) - ((uint16_t) 12))

//==========================[ Function Implementations ]=========================
void peak_idx_comp(
	const int16_t *restrict input,
	const uint16_t peak_idx,
	const int16_t peak_left_width,
	const int16_t peak_right_width,
	const int16_t scale_comp,
	int16_t *restrict peak_frac)
{
	#pragma aligned(input, 32)

	int64_t peak_num_sum_64, peak_den_sum_64, peak_div_64;
	xb_int40 peak_left_sum_40, peak_right_sum_40, peak_div_40, round_40;
	xb_vecNx16 peak_left_16v, peak_right_16v, seq_16v, peak_left_comp_16v, peak_right_comp_16v, peak_left_compare_16v, peak_right_compare_16v,
	seq_offset_16v, peak_left_seq_16v, zero_16v, peak_left_coeff_16v, peak_right_coeff_16v;
	xb_vecNx40 peak_left_40v, peak_right_40v, peak_left_mul_40v, peak_right_mul_40v;
	valign peak_left_valign, peak_right_valign;
	vboolN peak_left_vboolN, peak_right_vboolN;
	vsaS shift_vsaS;

	xb_vecNx16 *restrict p_peak_left = (xb_vecNx16 *restrict) &input[peak_idx - peak_left_width];
	xb_vecNx16 *restrict p_peak_right = (xb_vecNx16 *restrict) &input[peak_idx + 1];

	peak_left_valign = BBE_LANX16_PP(p_peak_left);
	peak_right_valign = BBE_LANX16_PP(p_peak_right);

	BBE_LAVNX16_XP(peak_left_16v, peak_left_valign, p_peak_left, 2 * peak_left_width);
	BBE_LAVNX16_XP(peak_right_16v, peak_right_valign, p_peak_right, 2 * peak_right_width);

	peak_left_40v = BBE_UNPKSNX16(peak_left_16v);
	peak_right_40v = BBE_UNPKSNX16(peak_right_16v);

	// S(16,14)
	peak_left_sum_40 = BBE_RADDNX40(peak_left_40v);
	peak_right_sum_40 = BBE_RADDNX40(peak_right_40v);

	peak_den_sum_64 = (int64_t) peak_left_sum_40 + (int64_t) peak_right_sum_40 + (int64_t) input[peak_idx];

	seq_16v = BBE_SEQNX16();

	peak_left_comp_16v = int16_rtor_xb_vecNxc16(peak_left_width);
	peak_right_comp_16v = int16_rtor_xb_vecNxc16(peak_right_width);

	peak_left_seq_16v = BBE_SUBNX16(seq_16v, peak_left_comp_16v);
	peak_left_compare_16v = int16_rtor_xb_vecNxc16(-1);
	seq_offset_16v = int16_rtor_xb_vecNxc16(1);
	peak_right_compare_16v = BBE_ADDNX16(seq_16v, seq_offset_16v);

	peak_left_vboolN = BBE_LENX16(peak_left_seq_16v, peak_left_compare_16v);
	peak_right_vboolN = BBE_GENX16(peak_right_comp_16v, peak_right_compare_16v);

	zero_16v = BBE_ZERONX16();

	peak_left_coeff_16v = BBE_MOVNX16T(peak_left_seq_16v, zero_16v, peak_left_vboolN);
	peak_right_coeff_16v = BBE_MOVNX16T(peak_right_compare_16v, zero_16v, peak_right_vboolN);

	// S(16,0) * S(16,14) = S(32,14)
	peak_left_mul_40v = BBE_MULNX16(peak_left_coeff_16v, peak_left_16v);
	peak_right_mul_40v = BBE_MULNX16(peak_right_coeff_16v, peak_right_16v);

	peak_left_sum_40 = BBE_RADDNX40(peak_left_mul_40v);
	peak_right_sum_40 = BBE_RADDNX40(peak_right_mul_40v);

	peak_left_sum_40 = BBE_OPERATOR_SLLI40(peak_left_sum_40, 15);
	peak_right_sum_40 = BBE_OPERATOR_SLLI40(peak_right_sum_40, 15);

	peak_num_sum_64 = (int64_t) peak_left_sum_40 + (int64_t) peak_right_sum_40;

	if (peak_num_sum_64 != 0)
	{
		// S(32,15) * S(16,12) = S(48,27)
		peak_div_64 = peak_num_sum_64 / peak_den_sum_64 * scale_comp;
		// S(48,27) -> S(40,19)
		peak_div_64 = ((uint64_t) peak_div_64 >> 8U);

		peak_div_40 = (xb_int40) peak_div_64;

		shift_vsaS = int16_rtor_vsaS(OUTPUT_SHIFT);

		round_40 = (xb_int40) (((uint32_t) 1) << ((uint16_t) (OUTPUT_SHIFT - 1U)));

		peak_div_40 = BBE_OPERATOR_ADD40(peak_div_40, round_40);

		*peak_frac = BBE_PACKV40(peak_div_40, shift_vsaS);
	}
}
