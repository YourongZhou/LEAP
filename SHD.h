/*
 * SHD.h
 *
 * Shifted Hamming Distance 预过滤器。它在正式进入 lane-based 搜索前，
 * 先用位并行的方式快速估计该对序列是否有机会落在目标阈值之内。
 */

#ifndef VECTOR_FILTER_H_
#define VECTOR_FILTER_H_

#ifndef __aligned__
#define __aligned__ __attribute__((aligned(16)))
#endif

#include <stdint.h>
#include <x86intrin.h>

// 直接基于 read / ref bit-plane 的 SSE 过滤接口。
int bit_vec_filter_sse(__m128i read_XMM0, __m128i read_XMM1,
		__m128i ref_XMM0, __m128i ref_XMM1, int length, int max_error);

// 直接基于 read / ref bit-plane 的 AVX 过滤接口。
int bit_vec_filter_avx(__m256i read_YMM0, __m256i read_YMM1,
		__m256i ref_YMM0, __m256i ref_YMM1, int length, int max_error);

// 直接基于预先算好的 lane mismatch mask 进行 AVX 过滤。
int bit_vec_filter_avx(__m256i *xor_masks, int length, int max_error);

#endif /* VECTOR_FILTER_H_ */
