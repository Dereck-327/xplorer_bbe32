/*
 * blank_neighbours_v.c
 *
 *  Created on: 2025.11.7
 *      Author: qirui
 */

//==========================[ Headers ]=========================
#include "blank_neighbours.h"

//==========================[ Function Implementations ]=========================
void blank_neighbours(
	int16_t *restrict input,
	const uint16_t peak_idx,
	const int16_t peak_left_width,
	const int16_t peak_right_width)
{
    #pragma aligned(input, 32)

    xb_vecNx16 zeros_16v;
    valign out_align;

    xb_vecNx16 *restrict p_peak_left = (xb_vecNx16 *restrict) &(input[peak_idx - peak_left_width]);
    xb_vecNx16 *restrict p_peak_right = (xb_vecNx16 *restrict) &(input[peak_idx]);

    zeros_16v = BBE_MOVVA16(0);

    out_align = BBE_ZALIGN();

    BBE_SAVNX16_XP(zeros_16v, out_align, p_peak_left, peak_left_width << 1);
    BBE_SANX16POS_FP(out_align, p_peak_left);

    BBE_SAVNX16_XP(zeros_16v, out_align, p_peak_right, (peak_right_width << 1) + 2);
    BBE_SANX16POS_FP(out_align, p_peak_right);
}
