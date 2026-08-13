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

#include "utils/compiler.h"

/*==========================[ log ]=========================*/

#define TRACE_DEBUG (0)
#define IS_CALCU_IQ_BEXP (0)   // 是否需要计算iq相位失配块指数计算
// 耗时统计
#ifndef PERF_STAT_ENABLED
#define PERF_STAT_ENABLED   1
#endif
/* 时钟源: 0=clock(), 1=Xtensa cycle counter  */
#ifndef PERF_STAT_USE_CYCLES
#define PERF_STAT_USE_CYCLES  1
#endif

/*==========================[ input data ]=========================*/

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

/*==========================[ DSP Stage ]=========================*/

/* 非相干累积: 四路 mag2 相加前先各留 1 bit 进位余量 */
#define NCI_BEXP_SHIFT	(1)

/* 均衡: mag2 U(16,15) * fft_eq U(16,15) -> S(16,14) */
#define EQ_SHIFT			(15)

/* 匹配滤波: 累加器 S(32,29) -> 输出 S(16,14) */
#define MF_SHIFT			(15)

/* 峰宽遮蔽: 置零区间在左右宽度之外再各扩这么多点, 防止同一目标被判成两个峰 */
#define BLANK_EXPANSION	(8U)

/* 亚 bin 插值: 内部 S(40,19), scale_comp 是 S(16,12), 输出取 Q12 */
#define FRAC_Q_IN		(19)
#define FRAC_Q_OUT		(12)
#define FRAC_SHIFT		(RA_FRAC_Q_IN - RA_FRAC_Q_OUT)

/* SNR: power_shift 的饱和上限, 超过按上限算 */
#define SNR_SHIFT_MAX	(8)

/* 输出索引以谱中心为零点 */
#define IDX_CORR			(2048U)

/*==========================[ DSP Stage assert ]=========================*/

/* 算子按整向量步进, 谱长必须是向量宽的整数倍 */
COMP_STATIC_ASSERT((MAG2_SIZE % BBE_SIMD_WIDTH) == 0U, mag2_size_vec_multiple);

/* DSP_PeakSearch 每轮处理两个向量 (BBE_DMAXNX16), 需要 2 倍向量宽整除 */
COMP_STATIC_ASSERT((MAG2_SIZE % (2U * BBE_SIMD_WIDTH)) == 0U, mag2_size_dual_vec);

/* 算宽度时向左读一整个向量, padding 不足即越界读到别的缓冲 */
COMP_STATIC_ASSERT(PADDING >= BBE_SIMD_WIDTH, padding_covers_vec);

/* 中心点索引必须落在谱内 */
COMP_STATIC_ASSERT(IDX_CORR < MAG2_SIZE, idx_corr_in_range);

#endif /* SIMULATE_CFG_H__ */
