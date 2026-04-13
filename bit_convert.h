/*
 * bit_convert.h
 *
 * 将字符序列转换成项目内部使用的 bit-plane 表示。
 * 对 DNA 字符串而言，最终编码约定为：
 * A = 00, C = 01, G = 10, T = 11
 */

#ifndef BIT_CONVERT_H_
#define BIT_CONVERT_H_

#include <stdint.h>

#ifndef __aligned__
	#define __aligned__ __attribute__((aligned(32)))
#endif

// 把序列编码为紧凑的 2 bit 串表示。
void c_convert2bit(char *str, int length, uint8_t *bits);

// 把最多 128 bp 的序列编码为两个 SSE bit-plane。
void sse_convert2bit(char *str, uint8_t *bits0, uint8_t *bits1);

// 把最多 256 bp 的序列编码为两个 AVX bit-plane。
void avx_convert2bit(char *str, uint8_t *bits0, uint8_t *bits1);

#endif /* BIT_CONVERT_H_ */
