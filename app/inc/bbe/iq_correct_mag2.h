/*
 * iq_correct_mag2.h
 *
 * IQ 失配盲补偿 (定点/BBE 向量): 对交织的 I_F/Q_F 频谱做
 *   Q_corr = gainCorr*(Q - rho*I);  sig = I + j*Q_corr;  mag2 = |sig|^2
 * rho/gainCorr 为实标量常数, 分别按 Q15/Q14 定点。
 */

#ifndef IQ_CORRECT_MAG2_H_
#define IQ_CORRECT_MAG2_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"
#include "simulate/cfg.h"

//==========================[ 定点标度 ]=========================
#define IQ_RHO_Q        (15)    /* rho    < 1  -> Q15 */
#define IQ_GAIN_Q       (14)    /* gain ~ 1.09 -> Q14 (范围 +-2) */

//==========================[ Function Prototypes ]=========================
/* @param[in]  i_f        I 路 FFT, 交织 [re,im,...], 2*samples 个 int16
 * @param[in]  q_f        Q 路 FFT, 同上
 * @param[in]  rho_q15    去相关系数, Q15
 * @param[in]  gain_q14   增益平衡系数, Q14
 * @param[out] mag2       功率谱 U(16,15), samples 个
 * @param[out] magsq_bexp 取模平方引入的块指数; 有效右移 = 15 + 2*magsq_bexp
 *                        (与 stage_prepare 的 mag2Bexp = 32 - 2*dataBexp - 2*magsqBexp 对应)
 */
void iq_correct_mag2(
    const complex_fract16 *restrict i_f,
    const complex_fract16 *restrict q_f,
    const int16_t rho_q15,
    const int16_t gain_q14,
    uint16_t *restrict mag2,
    uint8_t *restrict magsq_bexp,
    const uint16_t samples);

#endif /* IQ_CORRECT_MAG2_H_ */
