/**
    @file       cfg.h
    @brief      simulate config

    @version    0.1.0
    @date       2026-08-10
    @author     hjk

*/

#ifndef SIMULATE_CFG_H__
#define SIMULATE_CFG_H__

#include "NatureDSP_types.h"
#include "common.h"			/* BBE_SIMD_WIDTH */

#define MAG2_SIZE        (4096U)		/* 每 chirp 的功率谱点数 */
#define CHIRP_NUM        (4U)		/* 非相干累积的 chirp 数 */
#define FIRST_CHIRP      (0U)
#define LAST_CHIRP       (CHIRP_NUM - 1U)
#define PEAKS_NUM        (4U)		/* 每帧输出的峰数 */

/* 累积缓冲两端各留 RA_PADDING 个哨兵样点: 搜峰后算宽度时会向左右各读一整个
 * 向量宽, 有 padding 才不越界。 */
#define PADDING          (16U)
#define MAG2_SIZE_PAD    (MAG2_SIZE + 2U * PADDING)

#define MAX_BLANK_SIZE   (31U)						/* blank_table 条目上限 */
#define MF_COEFF_SIZE    (8U)						/* 15 tap 对称 FIR 的半边 */

/* 向量对其字节数 */
#define VEC_ALIGN        (2U * BBE_SIMD_WIDTH)

/*==========================[  ]=========================*/

/* SNR: power_shift 的饱和上限, 超过按上限算 */
#define SNR_SHIFT_MAX	(8)

#endif /* SIMULATE_CFG_H__ */
