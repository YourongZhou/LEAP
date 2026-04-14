/*
 * shift.h
 *
 * SSE / AVX 位移辅助函数。项目中的 lane 对齐、窗口提取和 mask 重定位
 * 都依赖这些跨字边界位移。
 */

#ifndef __SHIFT_H_
#define __SHIFT_H_

#include <stdint.h>
#include <x86intrin.h>
#include "leap_compat.h"

// 把 128 bit 向“逻辑右侧”移动若干 bit，用于 lane 对齐。
__m128i shift_right_sse(__m128i vec, int shift_num);
// 把 128 bit 向“逻辑左侧”移动若干 bit。
__m128i shift_left_sse(__m128i vec, int shift_num);
// 把 256 bit 向“逻辑右侧”移动若干 bit。
__m256i shift_right_avx(__m256i vec, int shift_num);
// 把 256 bit 向“逻辑左侧”移动若干 bit。
__m256i shift_left_avx(__m256i vec, int shift_num);

#endif /* __SHIFT_H_ */
