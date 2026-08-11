/*
 * calc_peaks_idx_frac.c
 *
 *  Created on: 2026.1.28
 *      Author: qirui
 */

//==========================[ Headers ]=========================
#include "calc_peaks_idx_frac.h"

//==========================[ Function Implementations ]=========================
void calc_peaks_idx_frac(
	const uint16_t *restrict peaks_idx_buffer,
	const int16_t *restrict peaks_frac_buffer,
	const uint8_t peaks_num,
	const uint16_t idx_corr,
	const uint8_t iq_reverse,
	int16_t *restrict peaks_idx,
	int16_t *restrict peaks_frac)
{
	xb_vecNx16 peaks_idx_16v, peaks_frac_16v, idx_corr_16v;
	valign idx_in_align, frac_in_align, idx_out_align, frac_out_align;

	// Input pointer
	xb_vecNx16 *restrict p_input_idx = (xb_vecNx16 *restrict) peaks_idx_buffer;
	xb_vecNx16 *restrict p_input_frac = (xb_vecNx16 *restrict) peaks_frac_buffer;
	// Output pointer
	xb_vecNx16U *restrict p_output_idx = (xb_vecNx16U *restrict) peaks_idx;
	xb_vecNx16U *restrict p_output_frac = (xb_vecNx16U *restrict) peaks_frac;

	// Idx
	idx_in_align = BBE_LANX16_PP(p_input_idx);
	BBE_LAVNX16_XP(peaks_idx_16v, idx_in_align, p_input_idx, peaks_num * 2U);
	idx_corr_16v = BBE_MOVVA16(idx_corr);
	// Substract
	peaks_idx_16v = BBE_SUBNX16(peaks_idx_16v, idx_corr_16v);

	// Frac
	frac_in_align = BBE_LANX16_PP(p_input_frac);
	BBE_LAVNX16_XP(peaks_frac_16v, frac_in_align, p_input_frac, peaks_num * 2U);

	// Negative
	if (iq_reverse)
	{
		peaks_idx_16v = BBE_NEGNX16(peaks_idx_16v);
		peaks_frac_16v = BBE_NEGNX16(peaks_frac_16v);
	}

	// Idx
	idx_out_align = BBE_ZALIGN();
	BBE_SAVNX16U_XP(peaks_idx_16v, idx_out_align, p_output_idx, peaks_num * 2U);
	BBE_SAVNX16UPOS_FP(idx_out_align, p_output_idx);

	// Frac
	frac_out_align = BBE_ZALIGN();
	BBE_SAVNX16U_XP(peaks_frac_16v, frac_out_align, p_output_frac, peaks_num * 2U);
	BBE_SAVNX16UPOS_FP(frac_out_align, p_output_frac);
}
