/*
 * iq_correct_mag2.c
 *
 * IQ 失配盲补偿 (全定点 BBE SIMD)。使用 complex_fract16 类型。
 *
 * 算法 (每个 bin):
 *   d  = Q - rho*I                    rho 为实标量 Q15
 *   Qc = gain * d                     gain 为实标量 Q14
 *   sig = I + j*Qc                    j 旋转: (re,im) -> (-im,re)
 *   mag2 = |sig|^2                    re^2 + im^2
 *
 * BBE 向量化:
 *   - complex_fract16* 直接强转为 xb_vecNx16* (见 vcabs_bbe32.c:89)
 *   - BBE_DSELNX16I 去交织/重交织
 *   - BBE_MULNX16 实标量乘法 (16x16->40 bit)
 *   - BBE_MAGINX16C 直接计算 |z|^2
 *
 * 参考: bbe32ep_library/complex/vcabs_bbe32.c
 *       app/src/bbe/nci_U16_S16_withMaxBexpScaling_vsm__fourFFTs.c
 */

//==========================[ Headers ]=========================
#include "bbe/iq_correct_mag2.h"
#include "utils/perf_stat.h"

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

    /* complex_fract16* 转为 xb_vecNx16* (标准 BBE 模式, 见 vcabs_bbe32.c:89) */
    const xb_vecNx16 *restrict p_i = (const xb_vecNx16 *) i_f;
    const xb_vecNx16 *restrict p_q = (const xb_vecNx16 *) q_f;
    xb_vecNx16 *restrict p_out = (xb_vecNx16 *) mag2;

    xb_vecNx16 i0, i1, q0, q1;      /* 交织载入 (每次 BBE_SIMD_WIDTH 复数) */
    xb_vecNx16 iRe, iIm, qRe, qIm;  /* 去交织 */
    xb_vecNx16 dRe, dIm;            /* Q - rho*I */
    xb_vecNx16 qcRe, qcIm;          /* Qc = gain*d */
    xb_vecNx16 sRe, sIm;            /* sig = I + j*Qc */
    xb_vecNx16 s0, s1;              /* 重交织给 MAGINX16C */
    xb_vecNx40 mul40;               /* 乘法中间结果 */
    xb_vecNx40 mag40;               /* |sig|^2 (40 bit) */
    xb_vecNx16 mag16;
    uint16_t i_blk;
    int bexp, shift;

    const vsaN rho_sa  = BBE_MOVVA16(rho_q15);
    const vsaN gain_sa = BBE_MOVVA16(gain_q14);
    const vsaN shift_rho  = BBE_MOVVA16((int16_t) IQ_RHO_Q);   /* 15 */
    const vsaN shift_gain = BBE_MOVVA16((int16_t) IQ_GAIN_Q);  /* 14 */

    NASSERT_ALIGN(i_f, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(q_f, (2 * BBE_SIMD_WIDTH));
    NASSERT_ALIGN(mag2, (2 * BBE_SIMD_WIDTH));
    NASSERT(samples % BBE_SIMD_WIDTH == 0);


    double sec;
    /* ---- 第一趟: 标量找峰 (向量峰值规约复杂, 暂用标量) ---- */
    {
        int64_t peakP = 0;
        uint16_t k;
        /* complex_fract16 直接访问 .s.re / .s.im */
        for (k = 0; k < samples; ++k)
        {
            int32_t iRe_s = (int32_t) i_f[k].s.re;
            int32_t iIm_s = (int32_t) i_f[k].s.im;
            int32_t qRe_s = (int32_t) q_f[k].s.re;
            int32_t qIm_s = (int32_t) q_f[k].s.im;

            int32_t dRe_s = qRe_s - (int32_t)(((int64_t) rho_q15 * iRe_s) >> IQ_RHO_Q);
            int32_t dIm_s = qIm_s - (int32_t)(((int64_t) rho_q15 * iIm_s) >> IQ_RHO_Q);
            int32_t qcRe_s = (int32_t)(((int64_t) gain_q14 * dRe_s) >> IQ_GAIN_Q);
            int32_t qcIm_s = (int32_t)(((int64_t) gain_q14 * dIm_s) >> IQ_GAIN_Q);

            int32_t sRe_s = iRe_s - qcIm_s;
            int32_t sIm_s = iIm_s + qcRe_s;

            int64_t p = (int64_t) sRe_s * sRe_s + (int64_t) sIm_s * sIm_s;
            if (p > peakP) { peakP = p; }
        }

        /* U(16,15): shift = 15 + 2*bexp, 使峰落进 uint16 */
        bexp = 0;
        if (peakP > 0)
        {
            int bits = 0;
            int64_t t = peakP;
            while (t > 0) { ++bits; t >>= 1; }
            bexp = ((bits - 32) + 1) / 2;
            if (bexp < 0) { bexp = 0; }
        }
        *magsq_bexp = (uint8_t) bexp;
    }

    /* ---- 第二趟: BBE 向量化校正 + 量化 ---- */
    shift = 15 + (2 * bexp);
    p_i = (const xb_vecNx16 *) i_f;
    p_q = (const xb_vecNx16 *) q_f;
    p_out = (xb_vecNx16 *) mag2;

    const vsaN out_shift = BBE_MOVVA16((int16_t) shift);

    for (i_blk = 0; i_blk < (samples / BBE_SIMD_WIDTH); ++i_blk)
    {
        /* 载入交织 complex_fract16: 每次 BBE_SIMD_WIDTH 复数 = 2*BBE_SIMD_WIDTH int16 */
        BBE_LVNX16_IP(i0, p_i, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(i1, p_i, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(q0, p_q, 2 * BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(q1, p_q, 2 * BBE_SIMD_WIDTH);

        /* 去交织: [re,im,re,im,...] -> [re0,re1,...], [im0,im1,...] */
        BBE_DSELNX16I(iRe, iIm, i0, i1, BBE_DSELI_DEINTERLEAVE_2);
        BBE_DSELNX16I(qRe, qIm, q0, q1, BBE_DSELI_DEINTERLEAVE_2);

        /* d = Q - rho*I (Q15 实标量乘) */
        mul40 = BBE_MULNX16(iRe, rho_sa);       /* rho * I_re -> 40 bit */
        mul40 = BBE_SRANX40(mul40, shift_rho);  /* >> 15 */
        dRe = BBE_PACKLNX40(mul40);             /* 回到 int16 */
        dRe = BBE_SUBNX16(qRe, dRe);            /* Q_re - rho*I_re */

        mul40 = BBE_MULNX16(iIm, rho_sa);
        mul40 = BBE_SRANX40(mul40, shift_rho);
        dIm = BBE_PACKLNX40(mul40);
        dIm = BBE_SUBNX16(qIm, dIm);            /* Q_im - rho*I_im */

        /* Qc = gain*d (Q14 实标量乘) */
        mul40 = BBE_MULNX16(dRe, gain_sa);
        mul40 = BBE_SRANX40(mul40, shift_gain); /* >> 14 */
        qcRe = BBE_PACKLNX40(mul40);

        mul40 = BBE_MULNX16(dIm, gain_sa);
        mul40 = BBE_SRANX40(mul40, shift_gain);
        qcIm = BBE_PACKLNX40(mul40);

        /* sig = I + j*Qc,  j*(re,im) = (-im,re) */
        sRe = BBE_SUBNX16(iRe, qcIm);           /* I_re - Qc_im */
        sIm = BBE_ADDNX16(iIm, qcRe);           /* I_im + Qc_re */

        /* 重交织成 [re,im,re,im,...] 给 BBE_MAGINX16C */
        BBE_DSELNX16I(s0, s1, sRe, sIm, BBE_DSELI_INTERLEAVE_2);

        /* |sig|^2 = sig_re^2 + sig_im^2 (BBE intrinsic 直接算) */
        mag40 = BBE_MAGINX16C(s0, s1);

        /* 右移舍入到 U(16,15) */
        mag40 = BBE_SRSNX40(mag40, out_shift);  /* 四舍五入右移 */
        mag16 = BBE_PACKLNX40(mag40);           /* saturate 到 int16 */

        /* 存为 uint16 (BBE_SVNX16 写 int16, 位模式与 uint16 等价) */
        BBE_SVNX16_IP(mag16, p_out, 2 * BBE_SIMD_WIDTH);
    }
}
