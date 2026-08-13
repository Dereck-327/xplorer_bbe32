/**
    @file       perf_stat.c
    @brief      性能统计模块实现
*/

#include "utils/perf_stat.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

#if PERF_STAT_ENABLED

//==========================[ static fun ]=========================
static PerfStatPoint g_points[PERF_STAT_MAX_POINTS];
static uint32_t g_point_count = 0;

//==========================[ clock ]=========================
#if PERF_STAT_USE_CYCLES
#define XCHAL_CLOCK_FREQ_MHZ (800U)  // 800MHz
    /* Xtensa cycle counter */
    #include <xtensa/hal.h>
    static inline perf_time_t perf_get_time(void) {
        return xthal_get_ccount();
    }
    static inline uint64_t perf_time_to_us(perf_time_t cycles) {
        /* 假设 CPU 频率, 需根据实际调整 */
        return cycles / (XCHAL_CLOCK_FREQ_MHZ);
    }
#else
    /* 标准 clock() */
    static inline perf_time_t perf_get_time(void) {
        return (perf_time_t) clock();
    }
    static inline uint64_t perf_time_to_us(perf_time_t ticks) {
        return (ticks * 1000000ULL) / CLOCKS_PER_SEC;
    }
#endif

//==========================[ internal ]=========================
/* 查找或创建统计点 */
static PerfStatPoint* perf_find_or_create(const char *name)
{
    uint32_t i;

    for (i = 0; i < g_point_count; ++i)
    {
        if (0 == strcmp(g_points[i].name, name))
        {
            return &g_points[i];
        }
    }

    if (g_point_count >= PERF_STAT_MAX_POINTS)
    {
        return NULL;  /* 溢出 */
    }

    PerfStatPoint *pt = &g_points[g_point_count++];
    memset(pt, 0, sizeof(PerfStatPoint));
    pt->name = name;
    pt->min = (perf_time_t) -1;  /* UINT64_MAX */
    return pt;
}

//==========================[ API ]=========================
void PerfStat_Init(void)
{
    memset(g_points, 0, sizeof(g_points));
    g_point_count = 0;
}

void PerfStat_Begin(const char *name)
{
    PerfStatPoint *pt = perf_find_or_create(name);
    if (NULL == pt)
    {
        return;  /* 溢出, 静默忽略 */
    }

    if (0 == pt->depth)
    {
        /* 顶层调用: 记录开始时间 */
        pt->start = perf_get_time();
    }

    ++pt->depth;  /* 嵌套深度 +1 */
}

void PerfStat_End(const char *name)
{
    perf_time_t end = perf_get_time();  /* 尽早采样 */
    PerfStatPoint *pt = perf_find_or_create(name);

    if (NULL == pt || 0 == pt->depth)
    {
        return;  /* 不匹配的 End 调用 */
    }

    --pt->depth;

    if (0 == pt->depth)
    {
        /* 顶层调用结束: 更新统计 */
        perf_time_t elapsed = end - pt->start;

        pt->total += elapsed;
        ++pt->count;

        if (elapsed < pt->min)
        {
            pt->min = elapsed;
        }
        if (elapsed > pt->max)
        {
            pt->max = elapsed;
        }
    }
}

void PerfStat_Report(void)
{
    uint32_t i;

    printf("\n========== Performance Statistics ==========\n");
    printf("%-24s %8s %10s %10s %10s\n",
           "Name", "Count", "Avg(us)", "Min(us)", "Max(us)");
    printf("------------------------------------------------------------\n");

    for (i = 0; i < g_point_count; ++i)
    {
        const PerfStatPoint *pt = &g_points[i];
        if (0 == pt->count)
        {
            continue;  /* 跳过未使用的 */
        }

        uint64_t avg_us = perf_time_to_us(pt->total / pt->count);
        uint64_t min_us = perf_time_to_us(pt->min);
        uint64_t max_us = perf_time_to_us(pt->max);

        printf("%-24s %8u %10llu %10llu %10llu\n",
               pt->name, pt->count,
               (unsigned long long) avg_us,
               (unsigned long long) min_us,
               (unsigned long long) max_us);
    }

    printf("============================================\n\n");
}

void PerfStat_Reset(void)
{
    uint32_t i;
    for (i = 0; i < g_point_count; ++i)
    {
        g_points[i].total = 0;
        g_points[i].min = (perf_time_t) -1;
        g_points[i].max = 0;
        g_points[i].count = 0;
        g_points[i].depth = 0;
    }
}

const PerfStatPoint* PerfStat_Get(const char *name)
{
    uint32_t i;
    for (i = 0; i < g_point_count; ++i)
    {
        if (0 == strcmp(g_points[i].name, name))
        {
            return &g_points[i];
        }
    }
    return NULL;
}

#endif /* PERF_STAT_ENABLED */
