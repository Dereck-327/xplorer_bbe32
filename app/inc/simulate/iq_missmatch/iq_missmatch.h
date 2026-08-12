/**
    @file       iq_missmatch.h
    @brief      IQ 相位失配盲补偿: 读 i_f/q_f 定点谱, 用 rho/gainCorr 校正, 输出功率谱

    @version    0.2.0
    @date       2026-08-12
    @author     hjk

*/


#ifndef SIMULATE_IQ_MISSMATCH__
#define SIMULATE_IQ_MISSMATCH__

#include "utils/bbe_type.h"
#include "utils/err.h"


/** @brief 读入一个 chirp 的 IQ 定点频谱, 做 Gram-Schmidt 失配校正, 生成功率谱帧

    文件格式 (ra_data_N.txt):
      i_f=re,im,re,im,...   I 路 FFT, 交织 int16, 共 2*MAG2_SIZE 个
      q_f=re,im,re,im,...   Q 路 FFT, 同上
      fft_bexp=<n>          两路共用块指数, 真值 = 尾数/2^15 * 2^fft_bexp

    校正 (与时域盲校准等价):
      Q_corr = gainCorr * (Q_F - rho * I_F)
      sig    = I_F + j * Q_corr
      mag2   = |sig|^2   -> 输出到 aFrame

    @param[in]   aPath      chirp 文件路径
    @param[in]   aRho       去相关系数 (I 泄漏到 Q 的投影)
    @param[in]   aGainCorr  增益平衡系数
    @param[out]  aFrame     输出帧: mag2(U16.15) + dataBexp + magsqBexp, 可直接 SubmitFrame
    @retval #ERR_OK
    @retval #ERR_INVAL_PARAMS  aPath/aFrame 为空
    @retval #ERR_NOT_FOUND     文件打不开
    @retval #ERR_DATA_INTEG    字段缺失或点数不符 */
ErrorType iqMissmatch_frame_load(const char *const aPath,
                                 const float aRho,
                                 const float aGainCorr,
                                 Bbe_FrameType *const aFrame);

#endif /* end of SIMULATE_IQ_MISSMATCH__ */
