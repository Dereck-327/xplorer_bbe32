/**
    @file       pipeline.h
    @brief      处理链流水线

    @version    0.1.0
    @date       2026-08-10
    @author     hjk

*/

#ifndef SIMULATE_PIPELINE_H__
#define SIMULATE_PIPELINE_H__

#include "utils/bbe_type.h"
#include "utils/err.h"
#include "utils/bbe_type.h"
/*==========================[ 阶段表 ]=========================*/
typedef ErrorType (*StageFnType)(Bbe_CtxType *const aCtx);

/* 阶段的执行时机 */
typedef uint8_t StageWhenType;
#define WHEN_EVERY_CHIRP	(0U)	/* 每个 chirp 都跑 */
#define WHEN_LAST_CHIRP		(1U)	/* 只在累积满时跑 */

typedef struct sStageDescType
{
    const char			*name;	/* 阶段名, 同时是 trace 标签 */
    StageFnType		 fn;
    StageWhenType	 when;
} StageDescType_t;

/* 阶段在表中的下标, 同时是 RA_TRACE_MASK 的位号。
 * 增删 改这里 + 改 pipeline.c 里的 StageTbl[]。 */
#define STAGE_PREPARE	(0U)	/* 拷谱 + 算本 chirp 块指数 */
#define STAGE_INTEGRATE	(1U)	/* 非相干累积 + fft_eq 均衡去底噪 */
#define STAGE_MATCHFILT	(2U)	/* 15 tap 匹配滤波 */
#define STAGE_PEAK		(3U)	/* 搜峰 + 峰宽 + FWHM + SNR */
#define STAGE_REFINE	(4U)	/* 按 compMode 做亚 bin 插值 */
#define STAGE_EMIT		(5U)	/* 索引归心 + iqReverse */
#define PIPELINE_STAGE_NUM		(6U)

/*==========================[ 由 bbe/*.c 实现 ]=========================*/
ErrorType StagePrepare(Bbe_CtxType *const aCtx);
ErrorType StageIntegrate(Bbe_CtxType *const aCtx);
ErrorType StageMatchFilt(Bbe_CtxType *const aCtx);
ErrorType StagePeak(Bbe_CtxType *const aCtx);
ErrorType StageRefine(Bbe_CtxType *const aCtx);
ErrorType StageEmit(Bbe_CtxType *const aCtx);


/*==========================[ 对外接口 ]=========================*/

/** @brief 初始化上下文并载入参数
    @behavior Sync, Re-entrant (不同 aCtx 之间无共享状态)
    @param[out] aCtx   待初始化的上下文
    @param[in]  aParam 参数, 内部会做 ParamCheck 后拷入
    @retval #ERR_OK            初始化完成, 状态置 RA_STATE_READY
    @retval #ERR_INVAL_PARAMS  aCtx/aParam 为空, 或参数校验不通过
    @post 中间缓冲与 chirp 计数清零 */
ErrorType Pipeline_Init(Bbe_CtxType *const aCtx, const ParamType_t *const aParam);

/** @brief 提交一个 chirp 的谱数据
    @param[inout] aCtx       上下文
    @param[in]    aMag2      MAG2_SIZE 点功率谱, 需按 VEC_ALIGN 对齐
    @param[in]    aDataBexp  FFT 块指数
    @param[in]    aMagsqBexp 幅度谱, 取模平方的额外指数
    @retval #ERR_OK            已收下
    @retval #ERR_INVAL_PARAMS  aCtx/aMag2 为空, 或 aMag2 未对齐
    @retval #ERR_UNINIT        未 RA_Init */
ErrorType Pipeline_SubmitFrame(Bbe_CtxType *const aCtx,
                               const uint16_t *const aMag2,
                               const uint8_t aDataBexp,
                               const uint8_t aMagsqBexp);

/** @brief 推进一步: 遍历阶段表, 按 when 过滤, 遇错即停
    @param[inout] aCtx 上下文
    @retval #ERR_OK            本帧已完成, 结果可用 Pipeline_GetResult 取
    @retval #ERR_BUSY          chirp 还没攒满, 需要继续 Pipeline_SubmitFrame
    @retval #ERR_UNINIT        未 Pipeline_Init
    @retval #ERR_INVAL_PARAMS  aCtx 为空
    @post 返回 ERR_OK 时 chirp 计数归零, 开始下一帧的累积 */
ErrorType Pipeline_ProcessFrame(Bbe_CtxType *const aCtx);

/** @brief 取出最近一帧的结果
    @param[in]  aCtx 上下文
    @param[out] aOut 结果拷出目标
    @retval #ERR_OK            已拷出
    @retval #ERR_INVAL_PARAMS  aCtx/aOut 为空
    @retval #ERR_UNINIT        未 Pipeline_Init */
ErrorType Pipeline_GetResult(const Bbe_CtxType *const aCtx, Bbe_ResultType *const aOut);

/** @brief 自检: 前 CHIRP_NUM-1 次推进必须 ERR_BUSY, 第 CHIRP_NUM 次 ERR_OK;
    未初始化就推进必须 ERR_UNINIT。
    @retval #ERR_OK            自检通过
    @retval #ERR_DATA_INTEG    状态机行为与预期不符 */
ErrorType Pipeline_SelfTest(void);


#endif /* SIMULATE_PIPELINE_H__ */
