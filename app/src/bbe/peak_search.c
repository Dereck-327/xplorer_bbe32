/*
 * peak_search.c
 *
 *  Created on: 2025.11.7
 *      Author: qirui
 */

//==========================[ Headers ]=========================
#include "peak_search.h"

//==========================[ Function Implementations ]=========================
void peak_search(
	int16_t *restrict m,
	uint16_t *restrict idx,
	const int16_t *restrict x,
	const uint16_t N)
{
    #pragma aligned(x, 32)

    uint16_t i;
    uint32_t max_value;
    xb_vecNx16 value1_16v;
    xb_vecNx16 value2_16v;
    vboolN max_boolN;

    // Input pointer
    const xb_vecNx16 *restrict p_value = (const xb_vecNx16 *restrict) x;

    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT((N % BBE_SIMD_WIDTH) == 0);

    // Init BBE_MAX register to min signed value
    // Init BBE_MAXIDX register to zero
    BBE_SETDUALMAX(MIN_INT16);

    // Return if input vector is null
    if (N == 0U)
    {
        *m = 0;
        *idx = 0;
    }
    else
    {
        // Loop input vector reading 32 samples on each iteration
        for (i = 0; i < ((N / BBE_SIMD_WIDTH) / 2); i++)
        {
            // Load 32 samples from input vector
            BBE_LVNX16_IP(value1_16v, p_value, 2 * BBE_SIMD_WIDTH);
            BBE_LVNX16_IP(value2_16v, p_value, 2 * BBE_SIMD_WIDTH);
            // Single signed peak search
            BBE_DMAXNX16(value1_16v, value2_16v);
        }

        // Aggregate peak and index to a single vboolN vector
        max_boolN = BBE_GTMAXNX16();
        BBE_MOVDUALMAXT(max_boolN);

        // Extract max value and its idx
        BBE_RBDUALMAXR(max_value, max_boolN);

        *m = (uint16_t) max_value;
        *idx = BBE_SELMAXIDX(max_boolN, 0);
    }
}
