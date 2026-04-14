/*
 * SIMD_ED.h
 *
 * 项目的主实现。它在 lane-based farthest-reaching 搜索框架上加入
 * AVX 位并行匹配扩展，并可选接入 SHD 预过滤。
 */

#include <iostream>
#include <x86intrin.h>
#include <cstring>
#include <string>
#include <cstdlib>
#include "print.h"
#include "leap_compat.h"
#include "shift.h"
#include "bit_convert.h"

#ifndef _MAX_LENGTH_ 
#define _MAX_LENGTH_ 256
#endif

#ifndef _AFFINE_DEF_POS_ 
#define _AFFINE_DEF_POS_ 10000000
#endif

using namespace std;

#ifndef __ED_INFO_H_
#define __ED_INFO_H_

enum ED_TYPE {MISMATCH, A_INS, B_INS};

struct ED_INFO {
	ED_TYPE type;
	int id_length;	
};

#endif

#ifndef __SIMD_ED_H_
#define __SIMD_ED_H_

// 编辑距离模式：局部、全局，以及两个半自由边界模式。
enum ED_modes {ED_LOCAL, ED_GLOBAL, ED_SEMI_FREE_BEGIN, ED_SEMI_FREE_END};
enum OP_modes {SSE, AVX};

class SIMD_ED {
public:
	// 构造 / 释放内部缓存。
	SIMD_ED();
	~SIMD_ED();

	// 初始化普通 Levenshtein 模式。
	void init_levenshtein(int ED_threshold, ED_modes mode = ED_LOCAL, bool SHD_enable = true);
	// 初始化 affine gap 模式。
	void init_affine(int gap_threshold, int AF_threshold, ED_modes mode, int ms_penalty, int gap_open_penalty, int gap_ext_penalty, bool SHD_enable = false, int SHD_threshold = 10);
	// 在给定 lane 上统计从 start_pos 开始的连续匹配长度。
	int count_ID_length_avx(int lane_idx, int start_pos);

	// 把字符形式的输入转换成双 bit-plane 输出。
	void convert_reads(char *read, char *ref, int length, uint8_t *A0, uint8_t *A1, uint8_t *B0, uint8_t *B1);

	// 直接加载字符形式输入。
	void load_reads(char *read, char *ref, int length);
	// 直接加载已经编码好的 bit-plane 输入。
	void load_reads(uint8_t *A0, uint8_t *A1, uint8_t *B0, uint8_t *B1, int length);
	void load_reads(__m256i A0, __m256i A1, __m256i B0, __m256i B1, int length);
	
	// 分别加载 read 与 reference，适合批处理复用。
	void load_ref(__m256i B0, __m256i B1);
	void load_read(__m256i A0, __m256i A1, int length);
	
	// 为所有 lane 预先计算 mismatch mask。
	void calculate_masks();

	// 根据当前模式重置状态。
	void reset();	
	// 根据当前模式执行搜索。
	void run();
	// 返回当前任务是否通过阈值约束。
	bool check_pass();
	// 根据当前模式执行回溯。
	void backtrack();
	// 返回最终编辑距离或 affine 代价。
	int get_ED();
	// 输出回溯后的 CIGAR 风格字符串。
	string get_CIGAR();
private:
	// Levenshtein 模式内部实现。
    void reset_levenshtein();
	void run_levenshtein();
	void backtrack_levenshtein();

	// affine gap 模式内部实现。
	void run_affine();
    void reset_affine();
	void backtrack_affine();

    // 通用运行参数。
	int ED_t;
	__m256i *hamming_masks;
	ED_modes mode;
	bool SHD_enable;

    // affine gap 模式参数。
    int SHD_threshold;
    bool affine_mode;
    int gap_threshold;
    int af_threshold;
    int ms_penalty;
    int gap_open_penalty;
    int gap_ext_penalty;
	int **I_pos;
	int **D_pos;
	int ED_count;

	// lane 状态表。
	int *cur_ED;
	int **start;
	int **end;

	// 回溯与最终结果。
	bool ED_pass;
	int final_lane_idx;
	int final_ED;
	ED_INFO *ED_info;
	// 仅在全局 / 半自由模式下用于补齐边界收敛代价。
	int converge_ED;
	int converge_final_lane;

	// lane 布局参数。
	int mid_lane;
	int total_lanes;
	
	// 当前输入缓存。
	int buffer_length;
	
	LEAP_ALIGNAS(32) char A[_MAX_LENGTH_];
	LEAP_ALIGNAS(32) char B[_MAX_LENGTH_];

	LEAP_ALIGNAS(32) uint8_t A_bit0_t[_MAX_LENGTH_ / 4];
	LEAP_ALIGNAS(32) uint8_t A_bit1_t[_MAX_LENGTH_ / 4];
	LEAP_ALIGNAS(32) uint8_t B_bit0_t[_MAX_LENGTH_ / 4];
	LEAP_ALIGNAS(32) uint8_t B_bit1_t[_MAX_LENGTH_ / 4];

};

#endif
