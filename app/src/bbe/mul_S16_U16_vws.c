/*
 * mul_S16_U16_vws.c
 *
 *  Created on: 2026.3.21
 *      Author: qirui
 */

//==========================[ Headers ]=========================
#include "mul_S16_U16_vws.h"

//==========================[ Function Implementations ]=========================
void mul_S16_U16_vws(
	const int16_t *restrict x,
	const uint16_t *restrict y,
	int16_t *restrict z,
	const int16_t rsh,
	const uint16_t samples)
{
    #pragma aligned(x, 32)
    #pragma aligned(y, 32)
    #pragma aligned(z, 32)

    uint16_t i;

    xb_vecNx16U x_16v, y_16v;
    xb_vecNx16 out_16v;
    xb_vecNx40 mul_40v;
    vsaN shift_vsaN;

    const xb_vecNx16 *restrict p_x = (xb_vecNx16 *restrict) x;
    const xb_vecNx16 *restrict p_y = (xb_vecNx16 *restrict) y;

    xb_vecNx16 *restrict p_z = (xb_vecNx16 *restrict) z;

    // Nature asserts
    NASSERT_ALIGN(x, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(y, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(z, (2 * BBE_SIMD_WIDTH));

    // Load shift
    shift_vsaN = BBE_MOVVSA32(rsh);

    for (i = 0; i < (samples / BBE_SIMD_WIDTH); i++)
    {
        // Load X vector
        BBE_LVNX16U_IP(x_16v, p_x, 2 * BBE_SIMD_WIDTH);

        // Load Y vector
        BBE_LVNX16U_IP(y_16v, p_y, 2 * BBE_SIMD_WIDTH);

        // Unsigned multiplication
        mul_40v = BBE_MULUUNX16(x_16v, y_16v);

        // Shift
        mul_40v = BBE_SRANX40(mul_40v, shift_vsaN);

        // Pack
        out_16v = BBE_PACKLNX40(mul_40v);

        // Store output
        BBE_SVNX16_IP(out_16v, p_z, 2 * BBE_SIMD_WIDTH);
    }
}
