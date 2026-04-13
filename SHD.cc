/*
 * SHD.cc
 *
 * Shifted Hamming Distance 预过滤实现。核心思路是对多个移位视图求交，
 * 用位并行方式快速判断这对序列是否可能落在目标阈值内。
 */

#include "print.h"
#include "SHD.h"
#include <stdio.h>
#include <string.h>
#include "popcount.h"
#include "bit_convert.h"
#include "mask.h"
#include "shift.h"

using namespace std;

/*
 * 由于这里的 bit 布局遵循 x86 小端序约定，源码中的“左 / 右移”语义
 * 需要结合最终的逻辑方向来理解。
 */

// 对 SSE 版本的 mismatch 向量做“填洞”，避免短零段被错误当作真实匹配。
void flip_false_zero(__m128i& vec) {

//	printf("vec: \t\t");
//	print128_bit(vec);

	// shuffle_epi8 会把高位当作控制位，因此这里先清掉高位。
	char temp_char = (char) 0x7f;
	__m128i boundary= _mm_set1_epi8(0x7f);
//	printf("MASK_7F: \t");
//	print128_bit(*boundary);

	__m128i shift = _mm_and_si128(boundary, vec);
	
//	printf("After and: \t");
//	print128_bit(shift);

	__m128i *mask = (__m128i *) MASK_SRS;

	shift = _mm_shuffle_epi8(*mask, shift);
	vec = _mm_or_si128(vec, shift);
	
//	printf("Last cases %d: \t", 0);
//	print128_bit(vec);

	int i;
	for (i = 1; i < 4; i++) {
		shift = _mm_srli_epi16(vec, i);
		shift = _mm_and_si128(boundary, shift);
//		printf("shift %d: \t", i);
//		print128_bit(shift);
		shift = _mm_shuffle_epi8(*mask, shift);
//		printf("shuffle %d: \t", i);
//		print128_bit(shift);
		shift = _mm_slli_epi16(shift, i);
		vec = _mm_or_si128(vec, shift);
//		printf("Last cases %d: \t", i);
//		print128_bit(vec);
	}

	// 处理跨 nibble / 跨字节边界时的传播情况。
	__m128i shifted_vec = shift_right_sse(vec, 4);
//	printf("shifted_vec: \t");
//	print128_bit(shifted_vec);

	shift = _mm_and_si128(boundary, shifted_vec);
//	printf("After and: \t");
//	print128_bit(shift);
	
	shift = _mm_shuffle_epi8(*mask, shift);
	shifted_vec = _mm_or_si128(shifted_vec, shift);
//	printf("Cross cases %d: \t", 0);
//	print128_bit(shifted_vec);

	for (i = 1; i < 4; i++) {
		shift = _mm_srli_epi16(shifted_vec, i);
		shift = _mm_and_si128(boundary, shift);
		shift = _mm_shuffle_epi8(*mask, shift);
		shift = _mm_slli_epi16(shift, i);
		shifted_vec = _mm_or_si128(shifted_vec, shift);
//		printf("Cross cases %d: \t", i);
//		print128_bit(shifted_vec);
	}

	shifted_vec = shift_left_sse(shifted_vec, 4);
	vec = _mm_or_si128(shifted_vec, vec);
//	printf("Final case: \t");
//	print128_bit(vec);
}

// AVX 版本的“填洞”处理，逻辑与 SSE 版本一致。
void flip_false_zero(__m256i& vec) {

//	printf("vec: \t\t");
//	print128_bit(vec);

	// shuffle_epi8 会把高位当作控制位，因此这里先清掉高位。
	__m256i boundary= _mm256_set1_epi8(0x7f);
//	printf("MASK_7F: \t");
//	print128_bit(*boundary);

	__m256i shift = _mm256_and_si256(boundary, vec);
	
//	printf("After and: \t");
//	print128_bit(shift);

	__m256i *mask = (__m256i *) MASK_SRS;

	shift = _mm256_shuffle_epi8(*mask, shift);
	vec = _mm256_or_si256(vec, shift);
	
//	printf("Last cases %d: \t", 0);
//	print128_bit(vec);

	int i;
	for (i = 1; i < 4; i++) {
		shift = _mm256_srli_epi16(vec, i);
		shift = _mm256_and_si256(boundary, shift);
//		printf("shift %d: \t", i);
//		print128_bit(shift);
		shift = _mm256_shuffle_epi8(*mask, shift);
//		printf("shuffle %d: \t", i);
//		print128_bit(shift);
		shift = _mm256_slli_epi16(shift, i);
		vec = _mm256_or_si256(vec, shift);
//		printf("Last cases %d: \t", i);
//		print128_bit(vec);
	}

	__m256i shifted_vec = shift_right_avx(vec, 4);
//	printf("shifted_vec: \t");
//	print128_bit(shifted_vec);

	shift = _mm256_and_si256(boundary, shifted_vec);
//	printf("After and: \t");
//	print128_bit(shift);
	
	shift = _mm256_shuffle_epi8(*mask, shift);
	shifted_vec = _mm256_or_si256(shifted_vec, shift);
//	printf("Cross cases %d: \t", 0);
//	print128_bit(shifted_vec);

	for (i = 1; i < 4; i++) {
		shift = _mm256_srli_epi16(shifted_vec, i);
		shift = _mm256_and_si256(boundary, shift);
		shift = _mm256_shuffle_epi8(*mask, shift);
		shift = _mm256_slli_epi16(shift, i);
		shifted_vec = _mm256_or_si256(shifted_vec, shift);
//		printf("Cross cases %d: \t", i);
//		print128_bit(shifted_vec);
	}

	shifted_vec = shift_left_avx(shifted_vec, 4);
	vec = _mm256_or_si256(shifted_vec, vec);
//	printf("Final case: \t");
//	print128_bit(vec);
}

// SSE 版本 SHD：直接从 read / ref bit-plane 构造多个移位视图并做过滤。
int bit_vec_filter_sse(__m128i read_XMM0, __m128i read_XMM1,
		__m128i ref_XMM0, __m128i ref_XMM1, int length, int max_error) {
	
	__m128i mask;
	if (length >= SSE_BIT_LENGTH)
		mask = _mm_set1_epi8(0xff);
	else
		mask = _mm_load_si128( (__m128i *) (MASK_SSE_END + (length *
										SSE_BYTE_LENGTH)));
	// 清除长度边界之外的无效 bit。
	read_XMM0 = _mm_and_si128(read_XMM0, mask);
	read_XMM1 = _mm_and_si128(read_XMM1, mask);
	ref_XMM0 = _mm_and_si128(ref_XMM0, mask);
	ref_XMM1 = _mm_and_si128(ref_XMM1, mask);

	int total_difference = 0;

	// 迭代所有允许的移位视图。
	int j;

	__m128i shift_XMM;
	__m128i diff_XMM;
	__m128i temp_diff_XMM;
	__m128i temp_shift_XMM;
	__m128i temp_mask;

	diff_XMM = _mm_xor_si128(read_XMM0, ref_XMM0);
	temp_diff_XMM = _mm_xor_si128(read_XMM1, ref_XMM1);
	diff_XMM = _mm_or_si128(diff_XMM, temp_diff_XMM);

	flip_false_zero(diff_XMM);

//	printf("diff_XMM: \t");
//	print128_bit(diff_XMM);

	for (j = 1; j <= max_error; j++) {
		temp_mask = _mm_load_si128( (__m128i *) (MASK_SSE_BEG + (j - 1) *
								SSE_BYTE_LENGTH));
		temp_mask = _mm_and_si128(temp_mask, mask);
		
		// 先对 read 侧做移位，再与 reference 比较。
		shift_XMM = shift_right_sse(read_XMM0, j);
		temp_diff_XMM = _mm_xor_si128(shift_XMM, ref_XMM0);
		shift_XMM = shift_right_sse(read_XMM1, j);
		temp_shift_XMM = _mm_xor_si128(shift_XMM, ref_XMM1);
		temp_diff_XMM = _mm_or_si128(temp_shift_XMM, temp_diff_XMM);
		temp_diff_XMM = _mm_and_si128(temp_diff_XMM, temp_mask);
//		printf("Before flip: \t");
//		print128_bit(temp_diff_XMM);
		flip_false_zero(temp_diff_XMM);
//		printf("After flip: \t");
//		print128_bit(temp_diff_XMM);
		diff_XMM = _mm_and_si128(diff_XMM, temp_diff_XMM);

//		printf("diff_XMM: \t");
//		print128_bit(diff_XMM);
		// 再对 reference 侧做移位，与原始 read 比较。
		shift_XMM = shift_right_sse(ref_XMM0, j);
		temp_diff_XMM = _mm_xor_si128(shift_XMM, read_XMM0);
		shift_XMM = shift_right_sse(ref_XMM1, j);
		temp_shift_XMM = _mm_xor_si128(shift_XMM, read_XMM1);
		temp_diff_XMM = _mm_or_si128(temp_shift_XMM, temp_diff_XMM);
		temp_diff_XMM = _mm_and_si128(temp_diff_XMM, temp_mask);
//		printf("Before flip: \t");
//		print128_bit(temp_diff_XMM);
		flip_false_zero(temp_diff_XMM);
//		printf("After flip: \t");
//		print128_bit(temp_diff_XMM);
		diff_XMM = _mm_and_si128(diff_XMM, temp_diff_XMM);
		
//		printf("diff_XMM: \t");
//		print128_bit(diff_XMM);
	}

	total_difference = popcount_SHD_sse(diff_XMM);

//	printf("total_difference: %d\n", total_difference);

	if (total_difference > (max_error) )
		return 0;
	else
		return 1;
}

// AVX 版本 SHD：直接从 read / ref bit-plane 构造多个移位视图并做过滤。
int bit_vec_filter_avx(__m256i read_YMM0, __m256i read_YMM1,
		__m256i ref_YMM0, __m256i ref_YMM1, int length, int max_error) {
	
	__m256i mask;
	if (length >= AVX_BIT_LENGTH)
		mask = _mm256_set1_epi8(0xff);
	else
		mask = _mm256_load_si256( (__m256i *) (MASK_AVX_END + (length *
										AVX_BYTE_LENGTH)));

	// 清除长度边界之外的无效 bit。
	read_YMM0 = _mm256_and_si256(read_YMM0, mask);
	read_YMM1 = _mm256_and_si256(read_YMM1, mask);
	ref_YMM0 = _mm256_and_si256(ref_YMM0, mask);
	ref_YMM1 = _mm256_and_si256(ref_YMM1, mask);

	int total_difference = 0;

	// 迭代所有允许的移位视图。
	int j;

	__m256i shift_YMM;
	__m256i diff_YMM;
	__m256i temp_diff_YMM;
	__m256i temp_shift_YMM;
	__m256i temp_mask;

	diff_YMM = _mm256_xor_si256(read_YMM0, ref_YMM0);
	temp_diff_YMM = _mm256_xor_si256(read_YMM1, ref_YMM1);
	diff_YMM = _mm256_or_si256(diff_YMM, temp_diff_YMM);

	flip_false_zero(diff_YMM);

#ifdef debug
	printf("diff_YMM: \t");
	print256_bit(diff_YMM);
#endif

	for (j = 1; j <= max_error; j++) {
		temp_mask = _mm256_load_si256( (__m256i *) (MASK_AVX_BEG + (j - 1) *
								AVX_BYTE_LENGTH));
		temp_mask = _mm256_and_si256(temp_mask, mask);
		
		// 先对 read 侧做移位，再与 reference 比较。
		shift_YMM = shift_right_avx(read_YMM0, j);
		temp_diff_YMM = _mm256_xor_si256(shift_YMM, ref_YMM0);
		shift_YMM = shift_right_avx(read_YMM1, j);
		temp_shift_YMM = _mm256_xor_si256(shift_YMM, ref_YMM1);
		temp_diff_YMM = _mm256_or_si256(temp_shift_YMM, temp_diff_YMM);
		temp_diff_YMM = _mm256_and_si256(temp_diff_YMM, temp_mask);
//		printf("Before flip: \t");
//		print128_bit(temp_diff_YMM);
		flip_false_zero(temp_diff_YMM);
//		printf("After flip: \t");
//		print128_bit(temp_diff_YMM);
		diff_YMM = _mm256_and_si256(diff_YMM, temp_diff_YMM);

//		printf("diff_YMM: \t");
//		print128_bit(diff_YMM);

		// 再对 reference 侧做移位，与原始 read 比较。
		shift_YMM = shift_right_avx(ref_YMM0, j);
		temp_diff_YMM = _mm256_xor_si256(shift_YMM, read_YMM0);
		shift_YMM = shift_right_avx(ref_YMM1, j);
		temp_shift_YMM = _mm256_xor_si256(shift_YMM, read_YMM1);
		temp_diff_YMM = _mm256_or_si256(temp_shift_YMM, temp_diff_YMM);
		temp_diff_YMM = _mm256_and_si256(temp_diff_YMM, temp_mask);
//		printf("Before flip: \t");
//		print128_bit(temp_diff_YMM);
		flip_false_zero(temp_diff_YMM);
#ifdef debug
		printf("After flip: \t");
		print256_bit(temp_diff_YMM);
#endif
		diff_YMM = _mm256_and_si256(diff_YMM, temp_diff_YMM);
		
#ifdef debug
		printf("diff_YMM: \t");
		print256_bit(diff_YMM);
#endif
	}

	total_difference = popcount_SHD_avx(diff_YMM);

#ifdef debug
	printf("total_difference: %d\n", total_difference);
#endif

	if (total_difference > (max_error) )
		return 0;
	else
		return 1;
}

// AVX 版本 SHD：直接复用预先计算好的 lane mismatch mask 做过滤。
int bit_vec_filter_avx(__m256i *xor_masks, int length, int max_error) {
	
	__m256i mask;
	if (length >= AVX_BIT_LENGTH)
		mask = _mm256_set1_epi8(0xff);
	else
		mask = _mm256_load_si256( (__m256i *) (MASK_AVX_END + (length *
										AVX_BYTE_LENGTH)));

	int total_difference = 0;

	__m256i diff_YMM;
	__m256i temp_mask;
	__m256i temp_diff;

	diff_YMM = _mm256_set1_epi8(0xff);

#ifdef debug
	printf("diff_YMM: \t");
	print256_bit(diff_YMM);
#endif

	for (int j = 0; j <= 2 * max_error; j++) {
		// 复用预先算好的 lane mask，并套上当前位置允许的边界掩码。
		int error = abs(j - max_error);
		temp_mask = _mm256_load_si256( (__m256i *) (MASK_AVX_BEG + (error - 1) *
								AVX_BYTE_LENGTH));
		temp_mask = _mm256_and_si256(temp_mask, mask);
		temp_diff = _mm256_and_si256(xor_masks[j], temp_mask);
		flip_false_zero(temp_mask);
//		printf("After flip: \t");
//		print128_bit(temp_diff_YMM);
		diff_YMM = _mm256_and_si256(diff_YMM, temp_diff);

#ifdef debug
		printf("diff_YMM: \t");
		print256_bit(diff_YMM);
#endif
	}

	total_difference = popcount_SHD_avx(diff_YMM);

#ifdef debug
	printf("total_difference: %d\n", total_difference);
#endif

	if (total_difference > (max_error) )
		return 0;
	else
		return 1;
}
