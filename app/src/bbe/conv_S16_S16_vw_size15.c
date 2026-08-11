/*
 * conv_S16_S16_vw_size15.c
 *
 *  Created on: 2025.11.3
 *      Author: qirui
 */

//==========================[ Headers ]=========================
#include "conv_S16_S16_vw_size15.h"

//==========================[ Function Implementations ]=========================
void conv_S16_S16_vw_size15(
	const int16_t *restrict input,
    int16_t *restrict output,
    const uint16_t samples,
    const int16_t *restrict mf_coeff)
{
    #pragma aligned (input, 32)
    #pragma aligned (output, 32)

    uint16_t i;

    xb_vecNx16 in0_16v, in1_16v, in2_16v, in3_16v, in4_16v, in5_16v, in6_16v, in7_16v, in8_16v, in9_16v, in10_16v, in11_16v, in12_16v, in13_16v, in14_16v,
	a0_coeff_16v, a1_coeff_16v, a2_coeff_16v, a3_coeff_16v, a4_coeff_16v, a5_coeff_16v, a6_coeff_16v, a7_coeff_16v, add_16v;
    xb_vecNx40 add_40v;
    valign in0_align, in1_align, in2_align, in3_align, in4_align, in5_align, in6_align, in7_align, in8_align, in9_align, in10_align, in11_align, in12_align, in13_align, in14_align;
    vsaN shift_vsaN;

    // Init input pointers (15 samples)
    xb_vecNx16 *restrict p_in0  = (xb_vecNx16 *restrict) &input[9];
    xb_vecNx16 *restrict p_in1  = (xb_vecNx16 *restrict) &input[10];
    xb_vecNx16 *restrict p_in2  = (xb_vecNx16 *restrict) &input[11];
	xb_vecNx16 *restrict p_in3  = (xb_vecNx16 *restrict) &input[12];
	xb_vecNx16 *restrict p_in4  = (xb_vecNx16 *restrict) &input[13];
    xb_vecNx16 *restrict p_in5  = (xb_vecNx16 *restrict) &input[14];
	xb_vecNx16 *restrict p_in6  = (xb_vecNx16 *restrict) &input[15];
	xb_vecNx16 *restrict p_in7  = (xb_vecNx16 *restrict) &input[16];
	xb_vecNx16 *restrict p_in8  = (xb_vecNx16 *restrict) &input[17];
    xb_vecNx16 *restrict p_in9  = (xb_vecNx16 *restrict) &input[18];
	xb_vecNx16 *restrict p_in10 = (xb_vecNx16 *restrict) &input[19];
	xb_vecNx16 *restrict p_in11 = (xb_vecNx16 *restrict) &input[20];
	xb_vecNx16 *restrict p_in12 = (xb_vecNx16 *restrict) &input[21];
	xb_vecNx16 *restrict p_in13 = (xb_vecNx16 *restrict) &input[22];
	xb_vecNx16 *restrict p_in14 = (xb_vecNx16 *restrict) &input[23];

    // Init output pointer
    xb_vecNx16 *restrict p_out = (xb_vecNx16 *restrict) output;

    // Load multiplier coefficient
    a0_coeff_16v = BBE_MOVVA16(mf_coeff[0]);
    a1_coeff_16v = BBE_MOVVA16(mf_coeff[1]);
    a2_coeff_16v = BBE_MOVVA16(mf_coeff[2]);
    a3_coeff_16v = BBE_MOVVA16(mf_coeff[3]);
    a4_coeff_16v = BBE_MOVVA16(mf_coeff[4]);
    a5_coeff_16v = BBE_MOVVA16(mf_coeff[5]);
    a6_coeff_16v = BBE_MOVVA16(mf_coeff[6]);
    a7_coeff_16v = BBE_MOVVA16(mf_coeff[7]);

    // Align left, center and right pointers
    in0_align  = BBE_LANX16_PP(p_in0);
    in1_align  = BBE_LANX16_PP(p_in1);
    in2_align  = BBE_LANX16_PP(p_in2);
    in3_align  = BBE_LANX16_PP(p_in3);
    in4_align  = BBE_LANX16_PP(p_in4);
    in5_align  = BBE_LANX16_PP(p_in5);
    in6_align  = BBE_LANX16_PP(p_in6);
    in7_align  = BBE_LANX16_PP(p_in7);
    in8_align  = BBE_LANX16_PP(p_in8);
    in9_align  = BBE_LANX16_PP(p_in9);
    in10_align = BBE_LANX16_PP(p_in10);
    in11_align = BBE_LANX16_PP(p_in11);
    in12_align = BBE_LANX16_PP(p_in12);
    in13_align = BBE_LANX16_PP(p_in13);
    in14_align = BBE_LANX16_PP(p_in14);

    // Load shift to convert from S(32,29) to S(16,14)
	shift_vsaN = BBE_MOVVA16(15);

    for (i = 0; i < (samples / BBE_SIMD_WIDTH); i++)
    {
    	BBE_LANX16_IP(in0_16v,  in0_align,  p_in0);
        BBE_LANX16_IP(in1_16v,  in1_align,  p_in1);
        BBE_LANX16_IP(in2_16v,  in2_align,  p_in2);
        BBE_LANX16_IP(in3_16v,  in3_align,  p_in3);
        BBE_LANX16_IP(in4_16v,  in4_align,  p_in4);
        BBE_LANX16_IP(in5_16v,  in5_align,  p_in5);
        BBE_LANX16_IP(in6_16v,  in6_align,  p_in6);
        BBE_LANX16_IP(in7_16v,  in7_align,  p_in7);
        BBE_LANX16_IP(in8_16v,  in8_align,  p_in8);
        BBE_LANX16_IP(in9_16v,  in9_align,  p_in9);
        BBE_LANX16_IP(in10_16v, in10_align, p_in10);
        BBE_LANX16_IP(in11_16v, in11_align, p_in11);
        BBE_LANX16_IP(in12_16v, in12_align, p_in12);
        BBE_LANX16_IP(in13_16v, in13_align, p_in13);
        BBE_LANX16_IP(in14_16v, in14_align, p_in14);

        // Multiply center sample
        add_40v = BBE_MULNX16(in7_16v, a0_coeff_16v);
        // Multiply left samples and accumulate with previous result
        BBE_MULANX16(add_40v, in0_16v,  a7_coeff_16v);
        BBE_MULANX16(add_40v, in1_16v,  a6_coeff_16v);
        BBE_MULANX16(add_40v, in2_16v,  a5_coeff_16v);
        BBE_MULANX16(add_40v, in3_16v,  a4_coeff_16v);
        BBE_MULANX16(add_40v, in4_16v,  a3_coeff_16v);
        BBE_MULANX16(add_40v, in5_16v,  a2_coeff_16v);
        BBE_MULANX16(add_40v, in6_16v,  a1_coeff_16v);
        // Multiply right samples and accumulate with previous result
        BBE_MULANX16(add_40v, in8_16v,  a1_coeff_16v);
        BBE_MULANX16(add_40v, in9_16v,  a2_coeff_16v);
        BBE_MULANX16(add_40v, in10_16v, a3_coeff_16v);
        BBE_MULANX16(add_40v, in11_16v, a4_coeff_16v);
        BBE_MULANX16(add_40v, in12_16v, a5_coeff_16v);
        BBE_MULANX16(add_40v, in13_16v, a6_coeff_16v);
        BBE_MULANX16(add_40v, in14_16v, a7_coeff_16v);

        // Pack wide vector and convert to S(16,14) result
        add_16v = BBE_PACKVNX40(add_40v, shift_vsaN);

        // Store in output vector
        BBE_SVNX16_IP(add_16v, p_out, 2 * BBE_SIMD_WIDTH);
    }
}
