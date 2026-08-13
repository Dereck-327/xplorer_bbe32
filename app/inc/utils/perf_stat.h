/**
    @file       perf_stat.h
    @brief      零侵入式性能统计模块

    使用方式:
    ```c
    void my_function(void) {
        PERF_BEGIN(my_function);     // 自动计时开始
        // ... 函数体 ...
        PERF_END(my_function);       // 自动计时结束
    }

    // 或使用作用域自动计时 (gcc/clang)
    void another_function(void) {
        PERF_SCOPE(another_function);
        // ... 退出时自动停止计时 ...
    }
    ```

    特性:
    - 编译时可完全禁用 (零开销)
    - 自动收集 min/max/avg/count
    - 支持嵌套调用
    - 可用高精度时钟 (Xtensa cycle counter)

    @version    0.1.0
    @date       2026-08-12
    @author     hjk
*/

#ifndef UTILS_PERF_STAT_H_
#define UTILS_PERF_STAT_H_

#include <stdint.h>
#include <stdbool.h>
#include "simulate/cfg.h"

//==========================[ 配置 ]=========================
/* 编译时开关: 0=禁用 (零开销), 1=启用 */
#ifndef PERF_STAT_ENABLED
#define PERF_STAT_ENABLED   1
#endif

/* 最大统计点数量 */
#ifndef PERF_STAT_MAX_POINTS
#define PERF_STAT_MAX_POINTS  32
#endif

/* 时钟源: 0=clock(), 1=Xtensa cycle counter  */
#ifndef PERF_STAT_USE_CYCLES
#define PERF_STAT_USE_CYCLES  0
#endif

//==========================[ 数据类型 ]=========================
typedef uint64_t perf_time_t;

typedef struct {
    const char *name;           /* 统计点名称  */
    perf_time_t total;          /* 累计时间 */
    perf_time_t min;            /* 最小时间 */
    perf_time_t max;            /* 最大时间 */
    uint32_t    count;          /* 调用次数 */
    perf_time_t start;          /* 当前开始时间 (用于嵌套) */
    uint8_t     depth;          /* 嵌套深度 */
} PerfStatPoint;

//==========================[ API ]=========================
#if PERF_STAT_ENABLED

/* 初始化 */
void PerfStat_Init(void);

/* 开始计时 (用 PERF_BEGIN 宏) */
void PerfStat_Begin(const char *name);

/* 结束计时 (用 PERF_END 宏) */
void PerfStat_End(const char *name);

/* 打印所有统计信息 */
void PerfStat_Report(void);

/* 重置所有统计 */
void PerfStat_Reset(void);

/* 获取单个统计点 (不存在则返回 NULL ) */
const PerfStatPoint* PerfStat_Get(const char *name);

//==========================[ 宏接口 ]=========================
/* 显式开始/结束 */
#define PERF_BEGIN(tag)  PerfStat_Begin(#tag)
#define PERF_END(tag)    PerfStat_End(#tag)

/* 作用域自动计时 (gcc/clang with cleanup attribute) */
#if defined(__GNUC__) || defined(__clang__)
    typedef struct {
        const char *name;
    } _PerfScopeGuard;

    static inline void _perf_scope_cleanup(_PerfScopeGuard *guard) {
        if (guard->name) PerfStat_End(guard->name);
    }

    #define PERF_SCOPE(tag) \
        _PerfScopeGuard _perf_guard_##tag __attribute__((cleanup(_perf_scope_cleanup))) = {#tag}; \
        PerfStat_Begin(#tag)
#else
    /* 不支持 cleanup 的编译器: 退化为手动 BEGIN/END */
    #define PERF_SCOPE(tag)  PERF_BEGIN(tag); /* 需手动 PERF_END */
#endif

#else  /* PERF_STAT_ENABLED == 0 */

/* 禁用时: 所有宏编译为空 */
#define PerfStat_Init()       ((void)0)
#define PerfStat_Begin(name)  ((void)0)
#define PerfStat_End(name)    ((void)0)
#define PerfStat_Report()     ((void)0)
#define PerfStat_Reset()      ((void)0)
#define PerfStat_Get(name)    (NULL)

#define PERF_BEGIN(tag)       ((void)0)
#define PERF_END(tag)         ((void)0)
#define PERF_SCOPE(tag)       ((void)0)

#endif /* PERF_STAT_ENABLED */

#endif /* UTILS_PERF_STAT_H_ */
