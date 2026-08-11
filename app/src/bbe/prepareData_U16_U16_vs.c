/*
 * prepareData_U16_U16_vs.c
 *
 *  Created on: 2026.3.21
 *      Author: qirui
 */

//==========================[ Headers ]=========================
#include "prepareData_U16_U16_vs.h"

//==========================[ Function Implementations ]=========================
void prepareData_U16_U16_vs(
    const uint16_t *restrict input,
    uint16_t *restrict output,
    const uint16_t samples)
{
    #pragma aligned(input, 32)
    #pragma aligned(output, 32)

    uint16_t i;
    xb_vecNx16U input_16v;

    const xb_vecNx16U *restrict p_input = (xb_vecNx16U *restrict) input;
    xb_vecNx16U *restrict p_output = (xb_vecNx16U *restrict) output;

    // Nature asserts
	NASSERT_ALIGN(input, (2 * BBE_SIMD_WIDTH));
	NASSERT_ALIGN(output, (2 * BBE_SIMD_WIDTH));

    for (i = 0; i < (samples / BBE_SIMD_WIDTH); i++)
    {
        // Load input
        BBE_LVNX16U_IP(input_16v, p_input, 2 * BBE_SIMD_WIDTH);
        // Store output
        BBE_SVNX16U_IP(input_16v, p_output, 2 * BBE_SIMD_WIDTH);
    }
}
