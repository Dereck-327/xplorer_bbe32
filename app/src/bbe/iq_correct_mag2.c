/*
 * iq_correct_mag2.c
 *
 * IQ 失配盲补偿
 *
 * 算法 (每个 bin):
 *   d  = Q - rho*I                    rho 为实标量 Q15     % 去相关
 *   Qc = gain * d                     gain 为实标量 Q14    % 增益平衡
 *   sig = I + j*Qc                    j 旋转: (re,im) -> (-im,re)
 *   mag2 = |sig|^2                    re^2 + im^2
 *
 * 两趟全向量:
 *   第一趟: 校正 + BBE_MAGINX16C + PRESCAN_SHIFT 压缩 + BBE_RMAXNX16 规约 -> bexp
 *   第二趟: 同样校正 + BBE_MAGINX16C + out_shift 量化 -> U(16,15) 输出
 *
 * 参考: bbe32ep_library/complex/vcabs_bbe32.c        (BBE_MAGINX16C)
 *       app/src/bbe/nci_...fourFFTs.c               (BBE_RMAXNX16)
 *       app/src/bbe/mul_S16_U16_vws.c               (MULNX16+SRANX40+PACKLNX40)
 */

//==========================[ Headers ]=========================
#include "bbe/iq_correct_mag2.h"
#include "utils/perf_stat.h"

/* 第一趟把 40-bit mag (最大 2*32767^2 ≈ 2^31) 右移到 int16 范围 (≤ 2^14)。
 * 恢复公式: bexp = max(0, (bits_16 - 14) / 2)
 *           其中 bits_16 = floor(log2(peak_16)) + 1                       */
#define IQ_PRESCAN_SHIFT    (17)

//==========================[ 向量化 IQ 校正 ]=========================
/* 把两趟相同的 12 步校正内联为宏, 维护方便 */
#define IQ_CORRECT_BODY(i0_, i1_, q0_, q1_, p_i_, p_q_)                    \
    BBE_LVNX16_IP(i0_, p_i_, 2 * BBE_SIMD_WIDTH);                          \
    BBE_LVNX16_IP(i1_, p_i_, 2 * BBE_SIMD_WIDTH);                          \
    BBE_LVNX16_IP(q0_, p_q_, 2 * BBE_SIMD_WIDTH);                          \
    BBE_LVNX16_IP(q1_, p_q_, 2 * BBE_SIMD_WIDTH);                          \
    BBE_DSELNX16I(iRe, iIm, i0_, i1_, BBE_DSELI_DEINTERLEAVE_2);           \
    BBE_DSELNX16I(qRe, qIm, q0_, q1_, BBE_DSELI_DEINTERLEAVE_2);           \
    mul40 = BBE_MULNX16(iRe, rho_sa);                                       \
    mul40 = BBE_SRANX40(mul40, shift_rho);                                  \
    dRe   = BBE_SUBNX16(qRe, BBE_PACKLNX40(mul40));                        \
    mul40 = BBE_MULNX16(iIm, rho_sa);                                       \
    mul40 = BBE_SRANX40(mul40, shift_rho);                                  \
    dIm   = BBE_SUBNX16(qIm, BBE_PACKLNX40(mul40));                        \
    mul40 = BBE_MULNX16(dRe, gain_sa);                                      \
    mul40 = BBE_SRANX40(mul40, shift_gain);                                 \
    qcRe  = BBE_PACKLNX40(mul40);                                           \
    mul40 = BBE_MULNX16(dIm, gain_sa);                                      \
    mul40 = BBE_SRANX40(mul40, shift_gain);                                 \
    qcIm  = BBE_PACKLNX40(mul40);                                           \
    sRe   = BBE_SUBNX16(iRe, qcIm);                                        \
    sIm   = BBE_ADDNX16(iIm, qcRe);                                        \
    BBE_DSELNX16I(s0, s1, sRe, sIm, BBE_DSELI_INTERLEAVE_2);              \
    mag40 = BBE_MAGINX16C(s0, s1)

//==========================[ Function Implementations ]=========================
void iq_correct_mag2(
    const complex_fract16 *restrict i_f,
    const complex_fract16 *restrict q_f,
    const int16_t rho_q15,
    const int16_t gain_q14,
    uint16_t *restrict mag2,
    uint8_t *restrict magsq_bexp,
    const uint16_t samples)
{
    #pragma aligned(i_f, 32)
    #pragma aligned(q_f, 32)
    #pragma aligned(mag2, 32)

    /* complex_fract16* -> xb_vecNx16*  */
    const xb_vecNx16 *restrict p_i;
    const xb_vecNx16 *restrict p_q;
    xb_vecNx16 *restrict p_out;

    xb_vecNx16 i0, i1, q0, q1;
    xb_vecNx16 iRe, iIm, qRe, qIm;
    xb_vecNx16 dRe, dIm, qcRe, qcIm;
    xb_vecNx16 sRe, sIm, s0, s1;
    xb_vecNx40 mul40, mag40;
    xb_vecNx16 mag16;
    xb_int16   vec_max;
    uint16_t   i_blk;
    int        bexp, shift, bits_16;

    const vsaN rho_sa     = BBE_MOVVA16(rho_q15);
    const vsaN gain_sa    = BBE_MOVVA16(gain_q14);
    const vsaN shift_rho  = BBE_MOVVA16((int16_t) IQ_RHO_Q);
    const vsaN shift_gain = BBE_MOVVA16((int16_t) IQ_GAIN_Q);
    const vsaN pre_shift  = BBE_MOVVA16((int16_t) IQ_PRESCAN_SHIFT);

    NASSERT_ALIGN(i_f,  (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(q_f,  (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(mag2, (2 * BBE_SIMD_WIDTH));
    NASSERT(samples % BBE_SIMD_WIDTH == 0);

    /* ==== 第一趟: 全向量化找峰 ==================================
     * mag40 最大约 2*32767^2 ≈ 2^31 (31 bit);
     * 右移 PRESCAN_SHIFT=17 后落进 [0, 2^14] ⊂ int16 范围;
     * BBE_RMAXNX16 每块规约, 跨块维护标量最大值。              */
    PERF_BEGIN(iq_scan_peak);

    vec_max = 0;
    p_i = (const xb_vecNx16 *) i_f;
    p_q = (const xb_vecNx16 *) q_f;

    for (i_blk = 0; i_blk < (samples / BBE_SIMD_WIDTH); ++i_blk)
    {
        IQ_CORRECT_BODY(i0, i1, q0, q1, p_i, p_q);    /* 校正 + mag40 */

        mag40 = BBE_SRANX40(mag40, pre_shift);          /* >> 17 压到 int16 */
        mag16 = BBE_PACKLNX40(mag40);                   /* saturate pack */

        {
            xb_int16 blk_max = BBE_RMAXNX16(mag16);    /* 块内规约 */
            if (BBE_OPERATOR_GT16(blk_max, vec_max))
            {
                vec_max = blk_max;
            }
        }
    }

    /* 从压缩后的 peak_16 恢复 bexp:
     * peak_true ≈ peak_16 << 17  →  bits_true = bits_16 + 17
     * bexp = max(0, ceil((bits_true - 32) / 2))
     *      = max(0, (bits_16 + 17 - 32 + 1) / 2)    [ceil: +1 再整除]
     *      = max(0, (bits_16 - 14) / 2)              */
    bexp = 0;
    bits_16 = 0;
    {
        int16_t t = (int16_t) vec_max;
        while (t > 0) { ++bits_16; t >>= 1; }
    }
    bexp = (bits_16 - 14) / 2;     /* integer floor division peak压缩过了 14 = 31-17 */
    if (bexp < 0) { bexp = 0; }
    *magsq_bexp = (uint8_t) bexp;

    PERF_END(iq_scan_peak);

    /* ==== 第二趟: 全向量化量化写出 ============================== */
    PERF_BEGIN(iq_quantize);

    shift = 15 + (2 * bexp);
    p_i   = (const xb_vecNx16 *) i_f;
    p_q   = (const xb_vecNx16 *) q_f;
    p_out = (xb_vecNx16 *) mag2;

    {
        const vsaN out_shift = BBE_MOVVA16((int16_t) shift);

        for (i_blk = 0; i_blk < (samples / BBE_SIMD_WIDTH); ++i_blk)
        {
            IQ_CORRECT_BODY(i0, i1, q0, q1, p_i, p_q);   /* 校正 + mag40 */

            mag40 = BBE_SRSNX40(mag40, out_shift);         /* 四舍五入右移 */
            mag16 = BBE_PACKLNX40(mag40);                  /* saturate -> int16 */

            BBE_SVNX16_IP(mag16, p_out, 2 * BBE_SIMD_WIDTH);
        }
    }

    PERF_END(iq_quantize);
}
