/**
    @file       system.h
    @brief      system enter

    @version    0.1.0
    @date       2026-08-10
    @author     hjk

*/

#ifndef SYSTEM_H__
#define SYSTEM_H__


#include "simulate/cfg.h"
#include "simulate/param.h"
#include "utils/err.h"

#define PATH_MAX_LEN	(512)
#define TXT_KEY_MAX		(64)

/* fft_eq 一行是 4096 个 "65535, ", 留足余量 */
#define TXT_VAL_MAX		(65536)
#define FRAME_KEY_NUM	(3U)	/* ra_data + 两个 bexp */


typedef struct sTxtSrcType;


ErrorType system_init(void);
void system_run(void);
void system_deinit(void);

#endif /* SYSTEM_H__ */