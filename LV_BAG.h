/*
 * LV_BAG.h
 *
 * 标量版 affine gap 基线实现。其状态组织方式与 SIMD_ED 的 affine 模式
 * 一致，适合做功能对照与性能对比。
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

enum ED_modes {ED_LOCAL, ED_GLOBAL, ED_SEMI_FREE_BEGIN, ED_SEMI_FREE_END};

class LV {
public:
	// 构造 / 释放 affine gap 搜索所需的状态表。
	LV();
	~LV();

	// 初始化 affine gap 参数与内部状态表。
	void init(int gap_threshold, int af_threshold, ED_modes mode, int ms_penalty, int gap_open_penalty, int gap_ext_penalty);

	// 加载待比较的 read 与 reference。
	void load_reads(char *read, char *ref, int length);

	// 重置一次运行的中间状态。
	void reset();	
	// 在当前阈值和打分参数下执行 affine gap 搜索。
	void run();
	// 返回当前任务是否通过阈值约束。
	bool check_pass();
	// 根据记录的状态恢复编辑路径。
	void backtrack();
	// 返回最终得分 / 距离。
	int get_ED();
	// 输出回溯后的 CIGAR 风格字符串。
	string get_CIGAR();
private:
	// 从给定 lane 与起点开始向前延伸连续匹配段。
	int count_ID_length(int lane_idx, int start_pos);

	// affine gap 状态表。
	int **start;
	int **end;
    int gap_threshold;
    int af_threshold;
    int ms_penalty;
    int gap_open_penalty;
    int gap_ext_penalty;
	int **I_pos;
	int **D_pos;
	int ED_count;

	// 回溯所需的最终状态。
	ED_modes mode;
	bool ED_pass;
	int final_lane_idx;
	int final_ED;
	ED_INFO *ED_info;
	int converge_ED;
	int converge_final_lane;

	// lane 布局参数。
	int mid_lane;
	int total_lanes;
	
	// 当前输入缓存。
	int buffer_length;
	
	char A[_MAX_LENGTH_];
	char B[_MAX_LENGTH_];
};

#endif
