/*
 * SIMD_ED.cc
 *
 * 主算法实现。该文件把 lane-based farthest-reaching 搜索与 AVX
 * 位并行匹配扩展结合起来，并提供普通 Levenshtein 与 affine gap 两种模式。
 */

#include <cstdio>
#include "SIMD_ED.h"
#include "SHD.h"
#include <cassert>

// 在指定 lane 上，从 start_pos 开始统计连续匹配长度。
// 这里通过对 mismatch mask 做位移，然后利用 tzcnt 找到第一个不匹配位。
int SIMD_ED::count_ID_length_avx(int lane_idx, int start_pos) {
	__m256i shifted_mask = shift_left_avx(hamming_masks[lane_idx], start_pos);
	
#ifdef debug	
	cout << "start_pos: " << start_pos << " ";
	print256_bit(shifted_mask);
#endif

	//unsigned long *byte_cast = (unsigned long*) &shifted_mask;
	uint64_t byte_cast [_MAX_LENGTH_/64] __aligned__;
	_mm256_store_si256((__m256i*) byte_cast, shifted_mask); 


	int length_result = 0;
	
	for (int i = 0; i <= (buffer_length - start_pos - 1) / (8 * sizeof(uint64_t) ); i++) {

#ifdef NO_BIT_VEC
		int id_length = 0;
		while ((id_length < 8 * sizeof(uint64_t) ) && ((byte_cast[i] & ((uint64_t)1 << id_length)) == 0) )
			id_length++;
/*
		cout << "id_length: " << id_length << endl;
		cout << "byte_cast[i]: " << byte_cast[i] << endl;
		int temp = byte_cast[i] & ((uint64_t)1 << id_length);
		cout << "temp: " << temp << endl;
*/
#else
		int id_length = _tzcnt_u64(byte_cast[i]);
#endif

		if (id_length == 8 * sizeof(uint64_t) && byte_cast[i] == 0) {
			//cout << "A" << endl;
			id_length = 8 * sizeof(uint64_t);
			length_result += id_length;
		}
		else {
			//cout << "B, byte_cast[" << i <<"]:" << byte_cast[i] << endl;
			length_result += id_length;
			break;
		}
	}

#ifdef debug	
	cout << "length result: " << length_result << endl;
#endif

	if (length_result < buffer_length - start_pos)
		return length_result;
	else
		return buffer_length - start_pos;
}

// 仅初始化默认状态，真正的缓冲区分配在 init_*() 中完成。
SIMD_ED::SIMD_ED() {
	ED_t = 0;

	hamming_masks = NULL;

	cur_ED = NULL;
	start = NULL;
	end = NULL;
	SHD_enable = false;

    affine_mode = false;
    I_pos = NULL;
    D_pos = NULL;

	mid_lane = 0;
	total_lanes = 0;
}

// 释放当前模式下分配的所有状态表与中间缓冲区。
SIMD_ED::~SIMD_ED() {
	if (total_lanes != 0) {
        if (affine_mode) {

		    for (int i = 0; i < total_lanes; i++) {
                delete [] I_pos[i];
                delete [] D_pos[i];
            }

            delete [] I_pos;
            delete [] D_pos;
        }
		else
			delete [] cur_ED;

		delete [] hamming_masks;

		for (int i = 0; i < total_lanes; i++) {
			delete [] start[i];
			delete [] end[i];
		}

		delete [] start;
		delete [] end;

		total_lanes = 0;

		hamming_masks = NULL;

		cur_ED = NULL;
		start = NULL;
		end = NULL;
		SHD_enable = false;

		affine_mode = false;
		I_pos = NULL;
		D_pos = NULL;
	}
}

// 把字符序列转换成两个 bit-plane，便于后续批量比较。
void SIMD_ED::convert_reads(char *read, char *ref, int length, uint8_t *A0, uint8_t *A1, uint8_t *B0, uint8_t *B1) {
	if (length > _MAX_LENGTH_)
		length = _MAX_LENGTH_;

	//cout << "length: " << length << " ref: " << ref << endl;;

	strncpy(A, read, length);
	sse_convert2bit(A, A_bit0_t, A_bit1_t);
	strncpy(B, ref, length);
	sse_convert2bit(B, B_bit0_t, B_bit1_t);

	memcpy(A0, A_bit0_t, (length - 1) / 8 + 1);
	memcpy(A1, A_bit1_t, (length - 1) / 8 + 1);
	memcpy(B0, B_bit0_t, (length - 1) / 8 + 1);
	memcpy(B1, B_bit1_t, (length - 1) / 8 + 1);
}

// 直接加载字符形式输入，并在内部完成 AVX bit-plane 编码。
void SIMD_ED::load_reads(char *read, char *ref, int length) {
	buffer_length = length;
	
	if (length > _MAX_LENGTH_)
		length = _MAX_LENGTH_;

	strncpy(A, read, length);
	avx_convert2bit(A, A_bit0_t, A_bit1_t);
	strncpy(B, ref, length);
	avx_convert2bit(B, B_bit0_t, B_bit1_t);

	//cout << "A: " << A  << endl;
	//cout << "B: " << B  << endl;
}

// 加载已经编码好的 bit-plane 输入。
void SIMD_ED::load_reads(uint8_t *A0, uint8_t *A1, uint8_t *B0, uint8_t *B1, int length) {
	buffer_length = length;
	memcpy(A_bit0_t, A0, (length - 1) / 8 + 1);
	memcpy(A_bit1_t, A1, (length - 1) / 8 + 1);
	memcpy(B_bit0_t, B0, (length - 1) / 8 + 1);
	memcpy(B_bit1_t, B1, (length - 1) / 8 + 1);
}

// 直接从 YMM 寄存器形式加载输入。
void SIMD_ED::load_reads(__m256i A0, __m256i A1, __m256i B0, __m256i B1, int length) {
	buffer_length = length;
	_mm256_store_si256((__m256i*) A_bit0_t, A0);
	_mm256_store_si256((__m256i*) A_bit1_t, A1);
	_mm256_store_si256((__m256i*) B_bit0_t, B0);
	_mm256_store_si256((__m256i*) B_bit1_t, B1);
}

// 只更新 read 侧，适合 reference 复用的场景。
void SIMD_ED::load_read(__m256i A0, __m256i A1, int length) {
	buffer_length = length;
	_mm256_store_si256((__m256i*) A_bit0_t, A0);
	_mm256_store_si256((__m256i*) A_bit1_t, A1);
}

// 只更新 reference 侧，适合批量比较多个 read。
void SIMD_ED::load_ref(__m256i B0, __m256i B1) {
	_mm256_store_si256((__m256i*) B_bit0_t, B0);
	_mm256_store_si256((__m256i*) B_bit1_t, B1);
}

// 为每条 lane 预计算 mismatch mask，后续搜索只需围绕 mask 推进。
void SIMD_ED::calculate_masks() {
	__m256i *A0 = (__m256i*) A_bit0_t;
	__m256i *A1 = (__m256i*) A_bit1_t;
	__m256i *B0 = (__m256i*) B_bit0_t;
	__m256i *B1 = (__m256i*) B_bit1_t;

	for (int i = 1; i < total_lanes - 1; i++) {
		__m256i shifted_A0 = *A0;
		__m256i shifted_A1 = *A1;
		__m256i shifted_B0 = *B0;
		__m256i shifted_B1 = *B1;

		int shift_amount = abs(i - mid_lane);

		if (i < mid_lane) {
			shifted_B0 = shift_right_avx(shifted_B0, shift_amount);
			shifted_B1 = shift_right_avx(shifted_B1, shift_amount);
		}
		else if (i > mid_lane) {
			shifted_A0 = shift_right_avx(shifted_A0, shift_amount);
			shifted_A1 = shift_right_avx(shifted_A1, shift_amount);
		}

		__m256i mask_bit0 = _mm256_xor_si256(shifted_A0, shifted_B0);
		__m256i mask_bit1 = _mm256_xor_si256(shifted_A1, shifted_B1);

		hamming_masks[i] = _mm256_or_si256(mask_bit0, mask_bit1);

		//cout << "hamming_masks[" << i << "]: ";
		//print128_bit(hamming_masks[i]);
		//cout << endl;
	}
}

// 初始化普通 Levenshtein 模式需要的 lane 状态表。
void SIMD_ED::init_levenshtein(int ED_threshold, ED_modes mode, bool SHD_enable) {
    // 如果之前处于 affine 模式，先释放对应状态表。
    this->~SIMD_ED();

	this->SHD_enable = SHD_enable;
	this->affine_mode = false;

	ED_t = ED_threshold;

	this->mode = mode;
	total_lanes = 2 * ED_t + 3;
	mid_lane = ED_t + 1;

	hamming_masks = new __m256i [total_lanes];

	cur_ED = new int[total_lanes];
	ED_info = new ED_INFO[ED_t + 1];

	start = new int* [total_lanes];
	end = new int* [total_lanes];

	for (int i = 0; i < total_lanes; i++) {
		start[i] = new int [ED_t + 1]();
		end[i] = new int [ED_t + 1]();
	}

	for (int i = 0; i < total_lanes; i++) {
		for (int j = 0; j < ED_t; j++) {
			start[i][j] = -2;
			end[i][j] = -2;
		}
	}

	for (int i = 1; i < total_lanes - 1; i++) {
		int ED = abs(i - mid_lane);
		if (mode == ED_GLOBAL || mode == ED_SEMI_FREE_END)
			start[i][ED] = ED;
		else
			start[i][0] = ED;
	}
}

// 为一次新的 Levenshtein 搜索重置当前编辑层信息。
void SIMD_ED::reset_levenshtein() {
	ED_pass = false;
	for (int i = 1; i < total_lanes - 1; i++) {
		if (mode == ED_GLOBAL || mode == ED_SEMI_FREE_END) {
			int ED = abs(i - mid_lane);
			cur_ED[i] = ED;
		}
		else {
			cur_ED[i] = 0;
		}
	}
}

// 普通 Levenshtein 模式：先可选跑 SHD，再按编辑距离层逐步扩展所有 lane。
void SIMD_ED::run_levenshtein() {
	if (SHD_enable && !bit_vec_filter_avx(hamming_masks+1, buffer_length, ED_t) ) {
		ED_pass = false;
		return;
	}

	int length;

	for (int l = 1; l < total_lanes - 1; l++) {
		if (cur_ED[l] == 0) {
			length = count_ID_length_avx(l, start[l][0]);
	
#ifdef debug
			cout << "length result: " << length << " buffer_length: " << buffer_length << endl;
#endif
	
			end[l][0] = length + start[l][0];
	
			if (end[l][0] == buffer_length) {
				final_lane_idx = l;
				final_ED = 0;
				ED_pass = true;
				return;
			}
			
			cur_ED[l]++;
		}
	}
	
	for (int e = 1; e <= ED_t; e++) {
		for (int l = 1; l < total_lanes - 1; l++) {
			if (cur_ED[l] == e) {
				
#ifdef debug	
				cout << "e: " << e << " l: " << l << endl;
#endif

				int top_offset = 0;
				int bot_offset = 0;

				if (l >= mid_lane)
					top_offset = 1;
				if (l <= mid_lane)
					bot_offset = 1;

				// 在三种可能转移中选出最靠前的起点。
				int max_start = end[l][e-1] + 1;
				if (end[l-1][e-1] + top_offset > max_start)
					max_start = end[l-1][e-1] + top_offset;
				if (end[l+1][e-1] + bot_offset > max_start)
					max_start = end[l+1][e-1] + bot_offset;

				start[l][e] = max_start;

				// 从该起点继续向前延伸连续匹配段。
				length = count_ID_length_avx(l, start[l][e]);

				end[l][e] = max_start + length;

#ifdef debug	
				cout << "start[" << l << "][" << e << "]: " << start[l][e];
				cout << "   end[" << l << "][" << e << "]: " << end[l][e] << endl;
#endif

				if (end[l][e] == buffer_length) {
					final_lane_idx = l;
					final_ED = e;
					ED_pass = true;
					
					break;
				}

				cur_ED[l]++;
			}
		}

		if (ED_pass)
			break;
	}

	if (mode == ED_GLOBAL || mode == ED_SEMI_FREE_BEGIN) {
		converge_ED = final_ED + abs(final_lane_idx - mid_lane);
		ED_pass = converge_ED <= ED_t;
	}
}

// 根据 Levenshtein 状态表逆向恢复编辑序列。
void SIMD_ED::backtrack_levenshtein() {

#ifdef debug	
				cout << "in backtrack!" << endl;
#endif

	ED_count = 0;

	if (mode == ED_GLOBAL || mode == ED_SEMI_FREE_BEGIN) {
		for (int e = converge_ED; e > final_ED; e--) {
#ifdef debug	
			cout << "e: " << e << endl;
#endif
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

	while (ED_probe != 0) {

#ifdef debug
		cout << "end[" << lane_idx << "][" << ED_probe  << "]: " << end[lane_idx][ED_probe];
		cout << "    start[" << lane_idx << "][" << ED_probe << "]: " << start[lane_idx][ED_probe] << endl;
#endif

		int match_count = end[lane_idx][ED_probe] - start[lane_idx][ED_probe];
		ED_info[ED_count].id_length = match_count;

		int top_offset = 0;
		int bot_offset = 0;

		if (lane_idx >= mid_lane)
			top_offset = 1;
		if (lane_idx <= mid_lane)
			bot_offset = 1;

		if (start[lane_idx][ED_probe] == (end[lane_idx][ED_probe - 1] + 1) ) {
#ifdef debug
			cout << "M" << endl;
#endif
			ED_info[ED_count].type = MISMATCH;
		}
		else if (start[lane_idx][ED_probe] == end[lane_idx - 1][ED_probe - 1] + top_offset) {
#ifdef debug
			cout << "I" << endl;
#endif
			lane_idx = lane_idx - 1;
			ED_info[ED_count].type = A_INS;
		}
		else if (start[lane_idx][ED_probe] == end[lane_idx + 1][ED_probe - 1] + bot_offset) {
#ifdef debug
			cout << "D" << endl;
#endif
			lane_idx = lane_idx + 1;
			ED_info[ED_count].type = B_INS;
		}
		else
			cerr << "Error! No lane!!" << endl;
		
		ED_probe--;
		ED_count++;
	}

	int match_count = end[lane_idx][ED_probe] - start[lane_idx][ED_probe];
	ED_info[ED_count].id_length = match_count;
#ifdef debug
		cout << "end[" << lane_idx << "][" << ED_probe  << "]: " << end[lane_idx][ED_probe];
		cout << "    start[" << lane_idx << "][" << ED_probe << "]: " << start[lane_idx][ED_probe] << endl;
	cout << "mid_lane: " << mid_lane << endl;
#endif
}

// 初始化 affine gap 模式需要的状态表与参数。
void SIMD_ED::init_affine(int gap_threshold, int af_threshold, ED_modes mode, int ms_penalty, int gap_open_penalty, int gap_ext_penalty, bool SHD_enable, int SHD_threshold) {
    // 如果之前处于普通模式，先释放对应状态表。
    this->~SIMD_ED();
    affine_mode = true;
    this->ms_penalty = ms_penalty;
    this->gap_open_penalty = gap_open_penalty;
    this->gap_ext_penalty = gap_ext_penalty;

 	this->gap_threshold = gap_threshold;
    this->af_threshold = af_threshold;

    this->SHD_enable = SHD_enable;
    this->SHD_threshold = SHD_threshold;

	this->mode = mode;
 
 	total_lanes = 2 * gap_threshold + 3;
	mid_lane = gap_threshold + 1;
	hamming_masks = new __m256i [total_lanes];
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

// 为一次新的 affine gap 搜索重置最终状态。
void SIMD_ED::reset_affine() {
	ED_pass = false;
	converge_ED = 1000000;
}

// affine gap 模式：同时推进匹配、插入链和删除链的最远可达位置。
void SIMD_ED::run_affine() {
	if (SHD_enable && !bit_vec_filter_avx(hamming_masks+1, buffer_length, SHD_threshold) ) {
		ED_pass = false;
		return;
	}

	int length;
	int top_offset = 0;
    int bot_offset = 0;

	for (int l = 1; l < total_lanes - 1; l++) {
		if (start[l][0] >= 0) {
			int lane_diff = abs(l - mid_lane);
			length = count_ID_length_avx(l, lane_diff);
	
#ifdef debug
			cout << "length result: " << length << " buffer_length: " << buffer_length << endl;
#endif
	
			end[l][0] = length + start[l][0];
	
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
				length = count_ID_length_avx(l, start[l][e]);
				end[l][e] = start[l][e] + length;

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

// 根据 affine gap 状态表逆向恢复 gap 打开 / 延伸与 mismatch 的组合路径。
void SIMD_ED::backtrack_affine() {

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

// 对外统一的 reset() 分发接口。
void SIMD_ED::reset() {
	if (affine_mode)
		reset_affine();
	else
		reset_levenshtein();
}

// 对外统一的 run() 分发接口。
void SIMD_ED::run() {
	if (affine_mode)
		run_affine();
	else
		run_levenshtein();
}

// 对外统一的 backtrack() 分发接口。
void SIMD_ED::backtrack() {
	if (affine_mode)
		backtrack_affine();
	else
		backtrack_levenshtein();
}

// 返回最近一次搜索的通过结果。
bool SIMD_ED::check_pass() {
	return ED_pass;
}

// 返回最近一次搜索得到的编辑距离或 affine 代价。
int SIMD_ED::get_ED() {
	if (mode == ED_GLOBAL || mode == ED_SEMI_FREE_BEGIN)
		return converge_ED;
	else
		return final_ED;
}

// 将回溯结果转换成简化的 CIGAR 风格字符串。
string SIMD_ED::get_CIGAR() {
	//char buffer[32];
	string CIGAR;
	CIGAR = to_string(ED_info[ED_count].id_length);
	//sprintf(buffer, "%d", ED_info[0].id_length);
	//CIGAR = string(buffer);
	for (int i = ED_count - 1; i >= 0; i--) {
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

		//sprintf(buffer, "%d", ED_info[0].id_length);
		//CIGAR += string(buffer);
		CIGAR += to_string(ED_info[i].id_length);
	}

	return CIGAR;
}
