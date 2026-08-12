/*
 * iq_correct_mag2.c
 *
 * IQ 失配盲补偿, 全定点。rho/gain 为实标量 (Q15/Q14)。
 *
 * 每个 bin:
 *   Qc_re = gain*(Q_re - rho*I_re)          [实标量 x 值]
 *   Qc_im = gain*(Q_im - rho*I_im)
 *   sig_re = I_re - Qc_im                    [sig = I + j*Qc, j 使实虚交换+取负]
 *   sig_im = I_im + Qc_re
 *   mag2   = sig_re^2 + sig_im^2
 * 然后按峰值定块指数, 右移舍入成 U(16,15)。
 *
 * 本实现为标量定点 (已用主机验证, 见 test)。向量化映射 (Xplorer 内做):
 *   - 交织 [re,im] 载入 -> BBE_DSELNX16I(...DEINTERLEAVE_2) 去交织成 re/im 两路
 *   - (x<<15 - rho*x)   -> BBE_MULNX16 (16x16->40) + 移位相减
 *   - *gain >> (15+14)  -> BBE_MULNX16 + BBE_PACKQNX40
 *   - 平方与求和         -> BBE_MULNX16 + BBE_ADDNX40 (仿 nci kernel 的 ADDSNX16)
 *   - 峰值               -> BBE_RMAXNX16 / 40 位规约
 *   - 右移舍入写出       -> BBE_PACKQNX40 + BBE_SVNX16U_IP
 */

//==========================[ Headers ]=========================
#include "bbe/iq_correct_mag2.h"

//==========================[ Static Helpers ]=========================
/* 单 bin 校正后功率 |sig|^2 (纯定点, 中间用 32/64 位防溢出) */
static int64_t bin_power(const int16_t *restrict i_f,
                         const int16_t *restrict q_f,
                         const uint16_t k,
                         const int16_t rho_q15,
                         const int16_t gain_q14)
{
    int32_t iRe = (int32_t) i_f[2U * k];
    int32_t iIm = (int32_t) i_f[2U * k + 1U];
    int32_t qRe = (int32_t) q_f[2U * k];
    int32_t qIm = (int32_t) q_f[2U * k + 1U];

    /* Q - rho*I  在 Q15 域: (Q<<15 - rho*I) */
    int32_t dRe = (qRe << IQ_RHO_Q) - ((int32_t) rho_q15 * iRe);
    int32_t dIm = (qIm << IQ_RHO_Q) - ((int32_t) rho_q15 * iIm);

    /* *gain(Q14) 再退回整数域: >> (15+14) */
    int32_t qcRe = (int32_t) (((int64_t) gain_q14 * dRe) >> (IQ_RHO_Q + IQ_GAIN_Q));
    int32_t qcIm = (int32_t) (((int64_t) gain_q14 * dIm) >> (IQ_RHO_Q + IQ_GAIN_Q));

    int32_t sRe = iRe - qcIm;
    int32_t sIm = iIm + qcRe;

    return ((int64_t) sRe * sRe) + ((int64_t) sIm * sIm);
}

//==========================[ Function Implementations ]=========================
void iq_correct_mag2(
    const int16_t *restrict i_f,
    const int16_t *restrict q_f,
    const int16_t rho_q15,
    const int16_t gain_q14,
    uint16_t *restrict mag2,
    uint8_t *restrict magsq_bexp,
    const uint16_t samples)
{
    int64_t peakP = 0;
    int      bexp = 0;      /* magsq 块指数, 有效右移 = 15 + 2*bexp */
    int      shift;
    uint16_t k;

    /* 第一趟: 求校正后功率峰 */
    for (k = 0U; k < samples; ++k)
    {
        int64_t p = bin_power(i_f, q_f, k, rho_q15, gain_q14);
        if (p > peakP)
        {
            peakP = p;
        }
    }

    /* U(16,15): mag2 = round(P >> (15+2*bexp)) 需落进 uint16。
     * 需 15+2*bexp >= floor(log2(P))+1-16 -> bexp = ceil((bits-1-31)/2), bits=位宽 */
    if (peakP > 0)
    {
        int bits = 0;
        int64_t t = peakP;
        while (t > 0) { ++bits; t >>= 1; }         /* bits = floor(log2(P))+1 */
        bexp = ((bits - 32) + 1) / 2;              /* ceil((bits-32)/2), bits-32 可能为负 */
        if (bexp < 0) { bexp = 0; }
    }
    *magsq_bexp = (uint8_t) bexp;

    /* 第二趟: 右移舍入量化成 U(16,15) */
    shift = 15 + (2 * bexp);
    for (k = 0U; k < samples; ++k)
    {
        int64_t p = bin_power(i_f, q_f, k, rho_q15, gain_q14);
        int64_t v = (p + ((int64_t) 1 << (shift - 1))) >> shift;   /* 四舍五入 */

        if (v < 0)        { v = 0; }
        if (v > 65535)    { v = 65535; }
        mag2[k] = (uint16_t) v;
    }
}
