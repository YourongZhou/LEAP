/*
 * RefDB.h
 *
 * 参考序列的预编码存储与查询接口。RefDB 会把 reference 转换成与 SIMD_ED
 * 相同的双 bit-plane 布局，以便查询时直接返回位对齐的 AVX 窗口。
 */

#ifndef REFDB_H_
#define REFDB_H_

#include <stdint.h>
#include <string>
#include <fstream>
#include <vector>
#include "bit_convert.h"
#include "shift.h"
#include "mask.h"

using namespace std;

struct chromoMeta{
	chromoMeta() {
		pos = 0;		
		loaded = false;
		length = 0;
		bit0 = NULL;
		bit1 = NULL;
	};

	~chromoMeta() {
		unload();
	};

	void unload() {
		if (loaded) {
			delete [] bit0;
			delete [] bit1;
			loaded = false;
			bit0 = NULL;
			bit1 = NULL;
		}
	}

	streampos pos;
	bool loaded;
	uint64_t length;
	uint8_t *bit0;
	uint8_t *bit1;
};

class RefDB {
public:
	RefDB();
	~RefDB();

	// 进入生成模式，准备写入新的 reference 数据库。
	void init_generate();
	// 向数据库中追加一条染色体或参考串。
	void add_chromo(char *chromo_string, uint64_t length);
	// 完成编码并写出最终数据库文件。
	void finish_and_store(string db_name);

	// 进入加载模式并读取数据库元信息。
	void init_load(string db_name);
	// 卸载所有已加载的染色体数据。
	void unload_all();
	// 按编号加载一条染色体。
	bool load_chromo(int chromo_num);
	// 按编号卸载一条染色体。
	void unload_chromo(int chromo_num);

	// 返回数据库中的染色体数量。
	int get_total_chromo_num();
	// 查询指定染色体上的一个窗口，并返回位对齐后的双 bit-plane。
	bool query(int chromo_num, int chromo_pos, int query_length, __m256i& bit0, __m256i& bit1);

	// 返回当前染色体长度。
	uint32_t get_chromo_length();
private:
	// 加载与查询阶段使用的元信息。
	chromoMeta * chromo_array;
	int chromo_total;

	// 生成阶段使用的临时缓存。
	uint64_t length_gen;
	uint8_t *bit0_gen;
	uint8_t *bit1_gen;
	vector<uint64_t> length_array;
	vector<streampos> pos_array;

	// 生成阶段的临时输出文件。
	ofstream temp_file;

	// 加载阶段使用的数据库文件句柄。
	ifstream db_file;
};

#endif /* REFDB_H_ */
