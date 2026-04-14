/*
 * LV_BAG.cc
 *
 * 标量版 affine gap 搜索实现。其状态转移与 SIMD_ED 的 affine 模式
 * 基本对应，用于做功能验证和性能基线。
 */

#include "LV_BAG.h"
#include <cstdio>
#include <cassert>

// 沿当前 lane 从 start_pos 开始继续向前，直到遇到第一个不匹配字符。
int __attribute__((optimize("O0"))) LV::count_ID_length(int lane_idx, int start_pos) {
	int A_idx_offset = 0;
	int B_idx_offset = 0;

	if (lane_idx < mid_lane)
		A_idx_offset = mid_lane - lane_idx;
	else if (lane_idx > mid_lane) 
		B_idx_offset = lane_idx - mid_lane;

	while (start_pos < buffer_length && A[start_pos - A_idx_offset] == B[start_pos - B_idx_offset]) {
		start_pos++;
	}

	return start_pos;
}

// 仅初始化指针与默认状态，真正分配在 init() 中完成。
LV::LV() {
	start = NULL;
	end = NULL;

    I_pos = NULL;
    D_pos = NULL;

	mid_lane = 0;
	total_lanes = 0;
}

// 释放 affine gap 状态表和回溯缓存。
LV::~LV() {
	if (total_lanes != 0) {
		for (int i = 0; i < total_lanes; i++) {
			delete [] I_pos[i];
			delete [] D_pos[i];
		}

		delete [] I_pos;
        delete [] D_pos;

		for (int i = 0; i < total_lanes; i++) {
			delete [] start[i];
			delete [] end[i];
		}

		delete [] start;
		delete [] end;

		total_lanes = 0;

		start = NULL;
		end = NULL;

		I_pos = NULL;
		D_pos = NULL;
	}
}

// 根据 affine gap 参数与模式初始化 lane 状态表。
void LV::init(int gap_threshold, int af_threshold, ED_modes mode, int ms_penalty, int gap_open_penalty, int gap_ext_penalty) {
	if (total_lanes != 0)
		this->~LV();

    this->ms_penalty = ms_penalty;
    this->gap_open_penalty = gap_open_penalty;
    this->gap_ext_penalty = gap_ext_penalty;

 	this->gap_threshold = gap_threshold;
    this->af_threshold = af_threshold;

	this->mode = mode;
 
 	total_lanes = 2 * gap_threshold + 3;
	mid_lane = gap_threshold + 1;
	ED_info = new ED_INFO[af_threshold + 1];

	start = new int* [total_lanes];
	end = new int* [total_lanes];
	I_pos = new int* [total_lanes];
	D_pos = new int* [total_lanes];


	for (int i = 0; i < total_lanes; i++) {
		start[i] = new int [af_threshold + 1]();
		end[i] = new int [af_threshold + 1]();
        I_pos[i] = new int [af_threshold + 1]();
        D_pos[i] = new int [af_threshold + 1]();
	}

	for (int i = 0; i < total_lanes; i++) {
		for (int e = 0; e <= af_threshold; e++) {
			I_pos[i][e] = -2;
			D_pos[i][e] = -2;
			start[i][e] = -2;
			end[i][e] = -2;
		}
		int distance = abs(i - mid_lane);
		if (distance == 0 || mode == ED_LOCAL || mode == ED_SEMI_FREE_BEGIN)
			start[i][0] = distance;
	}
}

// 拷贝当前 read / reference 到内部缓冲区。
void LV::load_reads(char *read, char *ref, int length) {
	buffer_length = length;
	
	if (length > _MAX_LENGTH_)
		length = _MAX_LENGTH_;

	strncpy(A, read, length);
	strncpy(B, ref, length);

	//cout << "buffer_length: " << buffer_length << endl;
}

// 为一次新的 affine gap 搜索重置最终状态。
void LV::reset() {
	ED_pass = false;
	converge_ED = 1000000;
}

// affine gap 核心推进过程：同时维护匹配、插入链和删除链的最远可达位置。
void LV::run() {
	int top_offset = 0;
    int bot_offset = 0;

	for (int l = 1; l < total_lanes - 1; l++) {
		if (start[l][0] >= 0) {
	
			end[l][0] = count_ID_length(l, start[l][0]);;

#ifdef debug
			cout << "end[" << l << "][0]: " << end[l][0] << endl;
#endif
	
			if (end[l][0] == buffer_length) {
				final_lane_idx = l;
				final_ED = 0;
				ED_pass = true;
				return;
			}
		}
	}
	
	for (int e = 1; e <= af_threshold; e++) {

		for (int l = 1; l < total_lanes - 1; l++) {

			// top_offset 表示路径向下走时需要补的位移。
            if (l >= mid_lane)
                top_offset = 1;
			else
				top_offset = 0;	

			// bot_offset 表示路径向上走时需要补的位移。
            if (l <= mid_lane)
                bot_offset = 1;
			else
				bot_offset = 0;

			// 默认假设 gap_open_penalty 大于 gap_ext_penalty。
			if (e >= gap_open_penalty && end[l-1][e-gap_open_penalty] >= 0 && end[l-1][e-gap_open_penalty] > I_pos[l-1][e-gap_ext_penalty]) {
				I_pos[l][e] = end[l-1][e-gap_open_penalty] + top_offset;
#ifdef debug	
				cout << "Update I[" << l << "][" << e << "] from open from e[" << l-1 << "][" << e-gap_open_penalty << "]" << end[l-1][e-gap_open_penalty] << endl;
#endif
			}
			else if (e >= gap_ext_penalty && I_pos[l-1][e-gap_ext_penalty] >= 0) {
#ifdef debug	
				cout << "Update I[" << l << "][" << e << "] from ext from I[" << l-1 << "][" << e-gap_ext_penalty << "]" << I_pos[l-1][e-gap_ext_penalty] << endl;
#endif
				I_pos[l][e] = I_pos[l-1][e-gap_ext_penalty] + top_offset;
			}

			if (e >= gap_open_penalty && end[l+1][e-gap_open_penalty] >= 0 && end[l+1][e-gap_open_penalty] > D_pos[l+1][e-gap_ext_penalty])
				D_pos[l][e] = end[l+1][e-gap_open_penalty] + bot_offset;
			else if (e >= gap_ext_penalty && D_pos[l+1][e-gap_ext_penalty] >= 0)
				D_pos[l][e] = D_pos[l+1][e-gap_ext_penalty] + bot_offset;

			start[l][e] = -2;

			if (e >= ms_penalty && end[l][e-ms_penalty] >= 0) {
				start[l][e] = end[l][e-ms_penalty] + 1;
#ifdef debug	
				cout << "coming from end[" << l << "][" << e-ms_penalty << "]:" << end[l][e-ms_penalty] << endl;
#endif
			}

			if (I_pos[l][e] > start[l][e]) {
				start[l][e] = I_pos[l][e];
#ifdef debug	
				cout << "coming from I[" << l << "][" << e << "]:" << I_pos[l][e] << endl;
#endif
			}

			if (D_pos[l][e] > start[l][e]) {
				start[l][e] = D_pos[l][e];
#ifdef debug	
				cout << "coming from D[" << l << "][" << e << "]:" << D_pos[l][e] << endl;
#endif
			}

#ifdef debug	
				cout << "***start[" << l << "][" << e << "]:" << start[l][e] << endl;
#endif

			if (start[l][e] >= 0) {
				end[l][e] = count_ID_length(l, start[l][e]);

#ifdef debug	
				cout << "e: " << e << " l: " << l << endl;
				cout << "start[" << l << "][" << e << "]: " << start[l][e];
				cout << "   end[" << l << "][" << e << "]: " << end[l][e] << endl;
#endif

				if (end[l][e] == buffer_length) {
					if (mode == ED_GLOBAL || mode == ED_SEMI_FREE_BEGIN) {
						int lane_diff = abs(mid_lane - l);
						int temp_converge_ED = e;
						if (lane_diff != 0)
							temp_converge_ED += gap_open_penalty + (lane_diff - 1)	* gap_ext_penalty;
						if (temp_converge_ED <= af_threshold && temp_converge_ED < converge_ED) {
							final_lane_idx = l;
							final_ED = e;
							ED_pass = true;
							converge_ED = temp_converge_ED;
						}
					}
					else {
						final_lane_idx = l;
						final_ED = e;
						ED_pass = true;
					}
				}
			}
		}

		if (ED_pass)
			break;
	}
}

// 返回最近一次 run() 的通过结果。
bool LV::check_pass() {
	return ED_pass;
}

// 根据 affine gap 状态表逆向恢复编辑路径。
void LV::backtrack() {

	ED_count = 0;

	if (mode == ED_GLOBAL || mode == ED_SEMI_FREE_BEGIN) {
		for (int e = 0; e <abs(mid_lane - final_lane_idx); e++) {
			ED_info[ED_count].id_length = 0;
			if (final_lane_idx > mid_lane)
				ED_info[ED_count].type = B_INS;
			else
				ED_info[ED_count].type = A_INS;
			ED_count++;
		}
	}

	int lane_idx = final_lane_idx;
	int ED_probe = final_ED;

	int top_offset = 0;
	int bot_offset = 0;


	while (ED_probe != 0) {

#ifdef debug
		cout << "end[" << lane_idx << "][" << ED_probe  << "]: " << end[lane_idx][ED_probe];
		cout << "    start[" << lane_idx << "][" << ED_probe << "]: " << start[lane_idx][ED_probe] << endl;
#endif

		int match_count = end[lane_idx][ED_probe] - start[lane_idx][ED_probe];
		ED_info[ED_count].id_length = match_count;

		if (start[lane_idx][ED_probe] == I_pos[lane_idx][ED_probe])	{

			if (lane_idx >= mid_lane)
				top_offset = 1;
			else
				top_offset = 0;

			while (I_pos[lane_idx - 1][ED_probe-gap_ext_penalty] + top_offset == I_pos[lane_idx][ED_probe]) {
				ED_info[ED_count].type = A_INS;
				ED_count++;
				// 为下一次编辑动作预留一个新的输出槽位。
				ED_info[ED_count].id_length = 0;

				lane_idx--;
				ED_probe -= gap_ext_penalty;

				if (lane_idx >= mid_lane)
					top_offset = 1;
				else
					top_offset = 0;

			}
			// 如果不能继续延伸，则当前步骤一定对应一次 gap open。
			assert(end[lane_idx-1][ED_probe-gap_open_penalty] + top_offset == I_pos[lane_idx][ED_probe]);
			ED_info[ED_count].type = A_INS;
			ED_count++;

			lane_idx--;
			ED_probe -= gap_open_penalty;

		}
		else if (start[lane_idx][ED_probe] == D_pos[lane_idx][ED_probe]) {

			if (lane_idx <= mid_lane)
				bot_offset = 1;
			else
				bot_offset = 0;

			while (D_pos[lane_idx+1][ED_probe-gap_ext_penalty] + bot_offset == D_pos[lane_idx][ED_probe]) {
				ED_info[ED_count].type = B_INS;
				ED_count++;
				// 为下一次编辑动作预留一个新的输出槽位。
				ED_info[ED_count].id_length = 0;

				lane_idx++;
				ED_probe -= gap_ext_penalty;

				if (lane_idx <= mid_lane)
					bot_offset = 1;
				else
					bot_offset = 0;

			}
			// 如果不能继续延伸，则当前步骤一定对应一次 gap open。
			assert(end[lane_idx+1][ED_probe-gap_open_penalty] + bot_offset == D_pos[lane_idx][ED_probe]);
			ED_info[ED_count].type = B_INS;
			ED_count++;

			lane_idx++;
			ED_probe -= gap_open_penalty;
		}
		else {
			assert(start[lane_idx][ED_probe] == end[lane_idx][ED_probe - ms_penalty] + 1);
			ED_info[ED_count].type = MISMATCH;
			ED_count++;
			ED_probe -= ms_penalty;
		}
	}

	int match_count = end[lane_idx][ED_probe] - start[lane_idx][ED_probe];
	ED_info[ED_probe].id_length = match_count;
}

// 返回最近一次搜索得到的总代价。
int LV::get_ED() {
	return final_ED;
}

// 将 affine gap 回溯结果转换成简化的 CIGAR 风格字符串。
string LV::get_CIGAR() {
	char buffer[32];
	string CIGAR;
	//CIGAR = to_string(ED_info[0].id_length);
	sprintf(buffer, "%d", ED_info[0].id_length);
	CIGAR = string(buffer);
	for (int i = 1; i <= final_ED; i++) {
		switch (ED_info[i].type) {
		case MISMATCH:
			CIGAR += 'M';
			break;
		case A_INS:
			CIGAR += 'I';
			break;
		case B_INS:
			CIGAR += 'D';
			break;
		}

		sprintf(buffer, "%d", ED_info[0].id_length);
		CIGAR += string(buffer);
		//CIGAR += to_string(ED_info[i].id_length);
	}

	return CIGAR;
}
