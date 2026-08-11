/*
 * nci_U16_S16_withMaxBexpScaling_vsm__fourFFTs.c
 *
 *  Created on: 2026.3.21
 *      Author: qirui
 */

//==========================[ Headers ]=========================
#include "nci_U16_S16_withMaxBexpScaling_vsm__fourFFTs.h"

//==========================[ Macros ]=========================
#define MAG2_BEXP_SHIFT		(1)
#define LIST_SIZE_FOUR_FFTS	(4U)

//==========================[ Function Implementations ]=========================
void nci_U16_S16_withMaxBexpScaling_vsm__fourFFTs(
    const uint16_t *restrict mag2,
	const int16_t *restrict mag2_bexp,
    int16_t *restrict mag2_acc,
	int16_t *restrict mag2_acc_bexp,
    const uint16_t samples)
{
    #pragma aligned(mag2, 32)
    #pragma aligned(mag2_bexp, 32)
    #pragma aligned(mag2_acc, 32)

    const xb_vecNx16 *restrict p_mag1 = (const xb_vecNx16 *) (&(mag2[0]));
    const xb_vecNx16 *restrict p_mag2 = (const xb_vecNx16 *) (&(mag2[samples]));
    const xb_vecNx16 *restrict p_mag3 = (const xb_vecNx16 *) (&(mag2[samples * 2U]));
    const xb_vecNx16 *restrict p_mag4 = (const xb_vecNx16 *) (&(mag2[samples * 3U]));
    const xb_vecNx16 *restrict p_mag2_bexp = (const xb_vecNx16 *) mag2_bexp;

    xb_vecNx16 *restrict p_acc_mant = (xb_vecNx16 *) mag2_acc;
    // Current value registers
    xb_vecNx16U mag2_value_16Uv[LIST_SIZE_FOUR_FFTS];
    xb_vecNx16 mag2_value_16v[LIST_SIZE_FOUR_FFTS], bexp_16v, zero_16v = {0};
    // Shift registers for BExps
    vsaN shift_vsaN[LIST_SIZE_FOUR_FFTS];
    // Scalar BExp values
    xb_int16 max_bexp = 0, mag2_bexps[LIST_SIZE_FOUR_FFTS];
    uint16_t i, max_selector = 0x000F; // Boolean vector for max selection
    vboolN max_boolN;

    // max_bexp calculation
    // Load the bexps values
    bexp_16v = BBE_LVNX16_I(p_mag2_bexp, 0);
    // Load the max_selector to a boolN vector
    max_boolN = BBE_LBN_I((const vboolN *) &max_selector, 0);
    // Move only the values needed for max calculation
    bexp_16v = BBE_MOVNX16T(bexp_16v, zero_16v, max_boolN);
    // Get the maximum
    max_bexp = BBE_RMAXNX16(bexp_16v);
    // Generate mag2_acc_bexp adding by 2
    // max_bexp = BBE_OPERATOR_ADD16(max_bexp, 2);

    // Calculate bexp
    mag2_bexps[0] = BBE_OPERATOR_SUB16(BBE_OPERATOR_SUB16(mag2_bexp[0], max_bexp), MAG2_BEXP_SHIFT);
    mag2_bexps[1] = BBE_OPERATOR_SUB16(BBE_OPERATOR_SUB16(mag2_bexp[1], max_bexp), MAG2_BEXP_SHIFT);
    mag2_bexps[2] = BBE_OPERATOR_SUB16(BBE_OPERATOR_SUB16(mag2_bexp[2], max_bexp), MAG2_BEXP_SHIFT);
    mag2_bexps[3] = BBE_OPERATOR_SUB16(BBE_OPERATOR_SUB16(mag2_bexp[3], max_bexp), MAG2_BEXP_SHIFT);

    // Load shift registers with each bexp value
    shift_vsaN[0] = BBE_MOVVA16(mag2_bexps[0]);
    shift_vsaN[1] = BBE_MOVVA16(mag2_bexps[1]);
    shift_vsaN[2] = BBE_MOVVA16(mag2_bexps[2]);
    shift_vsaN[3] = BBE_MOVVA16(mag2_bexps[3]);

    // mag2_acc calculation
    for (i = 0; i < (samples / BBE_SIMD_WIDTH); i++)
    {
        // Load newest and previous values
        BBE_LVNX16U_IP(mag2_value_16Uv[0], p_mag1, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16U_IP(mag2_value_16Uv[1], p_mag2, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16U_IP(mag2_value_16Uv[2], p_mag3, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16U_IP(mag2_value_16Uv[3], p_mag4, 2 * BBE_SIMD_WIDTH);

        // Mag2 shift operand: make place for sign bit
        mag2_value_16v[0] = BBE_SLLNX16(mag2_value_16Uv[0], shift_vsaN[0]);
        mag2_value_16v[1] = BBE_SLLNX16(mag2_value_16Uv[1], shift_vsaN[1]);
        mag2_value_16v[2] = BBE_SLLNX16(mag2_value_16Uv[2], shift_vsaN[2]);
        mag2_value_16v[3] = BBE_SLLNX16(mag2_value_16Uv[3], shift_vsaN[3]);

        // Add the four mag2 values - with saturation
        mag2_value_16v[0] = BBE_ADDSNX16(mag2_value_16v[0], mag2_value_16v[1]);
        mag2_value_16v[2] = BBE_ADDSNX16(mag2_value_16v[2], mag2_value_16v[3]);
        mag2_value_16v[0] = BBE_ADDSNX16(mag2_value_16v[0], mag2_value_16v[2]);

        // Store result
        BBE_SVNX16U_IP(mag2_value_16v[0], p_acc_mant, 2 * BBE_SIMD_WIDTH);
    }

    *mag2_acc_bexp = max_bexp;
}
