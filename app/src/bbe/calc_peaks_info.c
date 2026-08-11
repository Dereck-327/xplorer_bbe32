/*
 * calc_peaks_info.c
 *
 *  Created on: 2025.11.7
 *      Author: qirui
 */

//==========================[ Headers ]=========================
#include "calc_peaks_info.h"

#include "peak_search.h"
#include "calc_peak_width.h"
#include "blank_neighbours.h"
#include "calc_peaks_snr.h"

//==========================[ Macros ]=========================
#define PADDING			(16U)
#define EXPANSION_SIZE	(8U)

//==========================[ Function Implementations ]=========================
void calc_peaks_info(
	int16_t *restrict input,
	const uint16_t samples,
	int16_t *restrict peaks_left_width,
	int16_t *restrict peaks_right_width,
	const uint8_t peaks_num,
	uint16_t *restrict peaks_idx_buffer,
	const uint8_t ra_data_bexp,
	const int16_t power_shift,
	uint16_t *restrict peaks_fwhm,
	uint16_t *restrict peaks_snr)
{
	#pragma aligned(input, 32)

	int16_t peaks_val[peaks_num];
	uint16_t idx, peak_idx_offset;

    for (idx = 0; idx < peaks_num; idx++)
    {
    	peak_search(&(peaks_val[idx]), &(peaks_idx_buffer[idx]), &(input[PADDING]), samples);

    	peak_idx_offset = peaks_idx_buffer[idx] + (uint16_t) PADDING;

    	calc_peak_width(input, samples, peaks_val[idx], peak_idx_offset, &(peaks_left_width[idx]), &(peaks_right_width[idx]));

    	blank_neighbours(input, peak_idx_offset, peaks_left_width[idx] + EXPANSION_SIZE, peaks_right_width[idx] + EXPANSION_SIZE);

    	peaks_fwhm[idx] = (uint16_t) (peaks_left_width[idx] + peaks_right_width[idx] + 1);

    	peaks_left_width[idx] = peaks_right_width[idx] = (peaks_left_width[idx] > peaks_right_width[idx] ? peaks_left_width[idx] : peaks_right_width[idx]);
    }

    calc_peaks_snr(&(input[PADDING]), samples, peaks_val, peaks_snr, peaks_num, ra_data_bexp, power_shift);
}
