/**
    @file       compiler.h
    @brief      编译器抽象

    @version    0.1.0
    @date       2026-08-10
    @author     hjk

*/

#ifndef BBE_COMPILER_H__
#define BBE_COMPILER_H__

#include <stdint.h>
#include <string.h>

#if defined(__XCC__) || defined(__GNUC__)
	#define COMP_ALIGN(a)		__attribute__((aligned(a)))
	#define COMP_SECTION(s)		__attribute__((section(s)))
	#define COMP_INLINE			static inline
	#define COMP_UNUSED			__attribute__((unused))
#else
    #error "unsupported compiler: add COMP_* definitions for it here"
#endif

/* 编译期断言 */
#define COMP_STATIC_ASSERT(cond, tag) \
	typedef char comp_static_assert_##tag[(cond) ? 1 : -1] COMP_UNUSED

COMP_INLINE void Mem_Copy(void *const aDst, const void *const aSrc, const uint32_t aSize)
{
    (void) memcpy(aDst, aSrc, (size_t) aSize);
}

COMP_INLINE void Mem_Set(void *const aDst, const uint8_t aVal, const uint32_t aSize)
{
    (void) memset(aDst, (int) aVal, (size_t) aSize);
}

#endif /* BBE_COMPILER_H__ */
