/*
 * LV.h
 *
 * 标量版 Levenshtein 阈值搜索。该实现采用 Landau-Vishkin 风格的
 * “按对角线维护最远可达位置”的做法，作为 SIMD 版本的基线实现。
 */

#include <iostream>
#include <x86intrin.h>
#include <emmintrin.h>
#include <immintrin.h>
#include <cstring>
#include <string>
#include <cstdlib>

#ifndef _MAX_LENGTH_
#define _MAX_LENGTH_ 256
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

#ifndef __LV_H_
#define __LV_H_

class LV {
public:
	// 构造 / 释放内部状态缓存。
	LV();
	~LV();

	// 根据编辑距离阈值分配 lane 状态表。
	void init(int ED_threshold);

	// 加载待比较的 read 与 reference。
	void load_reads(char *read, char *ref, int length);

	// 重置一次匹配任务的运行状态。
	void reset();	
	// 在当前阈值下执行 Landau-Vishkin 搜索。
	void run();
	// 返回当前 read / reference 是否在阈值内通过。
	bool check_pass();
	// 根据记录下来的状态恢复编辑路径。
	void backtrack();
	// 返回最终编辑距离。
	int get_ED();
	// 把回溯结果转换为简化的 CIGAR 风格字符串。
	string get_CIGAR();
private:
	// 从给定 lane 和起点开始，向前统计连续匹配长度。
	int count_ID_length_sse(int lane_idx, int start_pos);

	int ED_t;

	// 每条 lane 在不同编辑距离预算下的状态。
	int *cur_ED;
	int **start;
	int **end;

	// 回溯所需的最终状态。
	bool ED_pass;
	int final_lane_idx;
	int final_ED;
	ED_INFO *ED_info;

	// lane 布局参数。
	int mid_lane;
	int total_lanes;
	
	// 当前输入缓存。
	int buffer_length;
	
	char A[_MAX_LENGTH_];
	char B[_MAX_LENGTH_];
};

#endif
