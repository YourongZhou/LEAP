# LEAP / ED 模块

这个仓库是一个面向短 DNA 序列的快速近似匹配实验代码库。核心目标是在给定错误阈值下，
快速判断两条序列是否足够接近，并在需要时回溯出编辑路径。

仓库中的主要实现包括：

- `LV`：标量版 Landau-Vishkin / Levenshtein 阈值搜索
- `LV_BAG`：标量版 affine gap 基线实现
- `SIMD_ED`：带 AVX 加速、可选 SHD 预过滤的主实现
- `RefDB`：把参考序列预编码后落盘，查询时直接返回位对齐的 SIMD 窗口
- 若干驱动程序和基准程序，例如 `testLV`、`vectorED`、`vectorSHD_ED`、`testNW`
- 两份 Bowtie2 集成实验目录：`bowtie2_augment/` 与 `bowtie2_skipED_augment/`
- vendored 的 Needleman-Wunsch 代码：`needleman_wunsch-0.3.5/`

## 总体结构

这个项目最适合按下面的流水线理解：

1. 输入编码
   `bit_convert.c/.h` 把 `A/C/G/T` 字符串转换成两个 bit-plane，后续阶段就可以通过
   XOR、位移和 popcount 一次处理很多个碱基。

2. lane 对齐与 mask 构造
   `shift.c/.h` 提供 SSE / AVX 的跨字边界位移能力。
   `SIMD_ED::calculate_masks()` 使用这些工具为每条对角 lane 预先构造 mismatch mask。

3. 可选的 SHD 预过滤
   `SHD.cc/.h` 实现 Shifted Hamming Distance 过滤器。
   如果一对序列明显不可能落在目标阈值内，就在进入主算法前直接剔除。

4. 核心阈值搜索
   `LV.cc/.h`、`LV_BAG.cc/.h` 和 `SIMD_ED.cc/.h` 都采用“按对角线维护最远可达前沿”
   的思路，而不是显式填完整个 DP 矩阵。

5. 可选回溯
   `backtrack()` 根据记录下来的状态恢复编辑动作，
   `get_CIGAR()` 进一步把编辑动作转换成紧凑的 CIGAR 风格字符串。

6. 可选集成层
   `RefDB.cc/.h` 提供参考库编码与查询能力；
   `bowtie2_augment/` 和 `bowtie2_skipED_augment/` 展示了该思路在 Bowtie2 中的集成实验。

## 关键数据结构

### 1. 基于对角线的搜索前沿

Levenshtein 和 affine gap 两套实现都围绕对角线 lane 展开：

- `mid_lane` 表示主对角线
- 偏离 `mid_lane` 的 lane 表示净插入 / 删除数量
- `start[l][e]` 表示 lane `l` 在误差或得分预算 `e` 下，从哪里开始继续扩展
- `end[l][e]` 表示从该状态能匹配到的最远位置

这就是整个项目最核心的优化点：它追踪的是“搜索前沿”，而不是完整 DP 表。

### 2. bit-plane 编码

每个碱基按 2 bit 编码：

- `A = 00`
- `C = 01`
- `G = 10`
- `T = 11`

仓库将两位拆到两个独立的 bit-plane 中保存。这样一来，某条 lane 的 mismatch mask
就可以通过下面三步得到：

1. 先把 read / reference 按该 lane 对齐
2. 对两个 bit-plane 分别做 XOR
3. 再把 XOR 结果 OR 起来

最终 mask 中的 `0` 表示匹配，`1` 表示不匹配。

## 主要文件

### 核心算法

- `SIMD_ED.h`, `SIMD_ED.cc`
- `LV.h`, `LV.cc`
- `LV_BAG.h`, `LV_BAG.cc`

### SIMD 支撑与过滤

- `bit_convert.h`, `bit_convert.c`
- `shift.h`, `shift.c`
- `SHD.h`, `SHD.cc`
- `mask.h`, `mask.c`
- `popcount.h`, `popcount.c`

### 参考库

- `RefDB.h`, `RefDB.cc`
- `RefDBMain.cc`

### 驱动 / 基准程序

- `testLV.cc`
- `testLV_BAG.cc`
- `vectorED.cc`
- `vectorSHD_ED.cc`
- `testNW.cc`

## 构建

项目使用顶层 `Makefile` 构建。仓库现在提供两份 conda 环境文件：

- `environment.yml`：手写的最小依赖说明，适合维护和后续调整
- `environment.lock.yml`：从当前已验证可工作的 `leap` 环境导出的锁定版本，适合直接复现

如果你是在 `linux-64 / x86_64` 上直接复现当前工作环境，优先用锁定版本：

```bash
conda env create -f environment.lock.yml
conda activate leap
```

如果你只是想按最小依赖重新求解环境，或者你准备自己调整包版本，再用：

```bash
conda env create -f environment.yml
conda activate leap
```

激活后，Makefile 会优先使用环境内的编译器工具链。

```bash
make smoke
make testLV
make testLV_BAG
make vectorSHD_ED
make vectorED
make testRefDB
```

说明：

- 默认构建目标不再依赖 `g++-5`
- Makefile 支持通过 `CC`、`CXX`、`OPT_FLAGS`、`ISA_FLAGS` 和 `USE_PARASAIL` 覆盖默认配置
- SIMD 编译选项默认假定机器支持 BMI、AVX2 和 SSE4.2
- `testNW` 需要显式启用 Parasail：`make testNW USE_PARASAIL=1`
- `environment.lock.yml` 是基于当前已验证的 `linux-64` conda 环境导出，不保证跨平台可直接复用
- Bowtie2 目录属于单独的集成实验，不走顶层默认构建流程

## 输入格式与运行方式

大部分驱动程序都从 `stdin` 读取输入，并采用“read 一行、reference 一行”的格式。
读取到 `end_of_file` 时停止。

示例输入：

```text
ACGTACGT
ACGTTCGT
ACGT
ACGT
end_of_file
```

仓库另外提供了一个可复用的基础输入集：

```bash
testdata/smoke_input.txt
```

## 当前状态

截至 2026-04-15，顶层核心构建已经迁移到现代 conda 工具链，当前可确认的状态如下：

- `testLV`：可构建，可运行
- `testLV_BAG`：可构建，可运行
- `vectorED`：可构建，可运行
- `testRefDB`：可构建，可运行
- `testNW`：在 `USE_PARASAIL=1` 时可构建，可运行
- `vectorSHD_ED`：可构建、可运行，但 Levenshtein 路径仍有正确性偏差，不能当成已验证结果

目前已复现的一个明确差异是，在同一组 smoke 输入上：

```bash
./testLV 1 < testdata/smoke_input.txt
# passNum: 5

./vectorSHD_ED 1 1 1 < testdata/smoke_input.txt
# passNum: 3
```

这说明当前待排查的问题在 SIMD Levenshtein 路径，而不是顶层构建本身。

## 最小 Smoke Case

如果只是想确认“环境、构建、最基本运行”是否正常，建议先跑这一组：

```bash
conda activate leap
make smoke

./testLV 1 < testdata/smoke_input.txt
./vectorED 2 0 < testdata/smoke_input.txt
./testRefDB
```

如果环境里装了 Parasail，再补一条：

```bash
make testNW USE_PARASAIL=1
./testNW 2 1 < testdata/smoke_input.txt
```

如果你要复现当前已知问题，再跑：

```bash
./vectorSHD_ED 1 1 1 < testdata/smoke_input.txt
```

### 标量 Levenshtein 基线

```bash
./testLV 2 < input.txt
```

### 标量 affine gap 基线

```bash
./testLV_BAG 2 < input.txt
```

### 带 SHD 的 SIMD 阈值搜索

```bash
./vectorSHD_ED 2 1 1 < input.txt
```

注意：这个入口当前只能作为“构建和运行是否通”的检查，不能作为已经完成验证的 Levenshtein 基准。

参数含义：

- `argv[1]`：错误阈值
- `argv[2]`：是否显式开启 / 关闭 SHD，`1` / `0`，可选
- `argv[3]`：使用 Levenshtein (`1`) 还是 affine gap (`0`)，可选

### SIMD affine 与 Needleman-Wunsch 对比

```bash
./vectorED 2 0 < input.txt
./vectorED 2 1 < input.txt
```

这里第二个参数用于选择：

- `0`：仓库内的 SIMD affine 实现
- `1`：Needleman-Wunsch 基线路径

## RefDB 使用流程

`RefDB` 的目标是预编码参考序列，避免每次查询窗口时重复对 reference 做字符到 bit 的转换。

典型流程如下：

1. `init_generate()`
2. 对每条染色体或参考串调用 `add_chromo(...)`
3. `finish_and_store("name")`
4. `init_load("name")`
5. 通过 `query(...)` 取出已经位对齐的 AVX 窗口

可以参考 `RefDBMain.cc` 中的最小示例。

## 仓库布局

- 顶层文件：项目自有算法、工具和基准代码
- `needleman_wunsch-0.3.5/`：vendored 的参考实现
- `bowtie2_augment/`、`bowtie2_skipED_augment/`：Bowtie2 集成实验副本
- `results/`：实验输出与脚本

## 维护说明

这是一个研究型代码库，不是工程化打磨完的通用库。使用时建议注意以下限制：

- 多处长度上限是写死的，例如 `128` 或 `256`
- 内存管理以手工分配 / 释放为主
- 驱动程序更偏向 benchmark harness，而不是面向最终用户的 CLI
- 一些命名和结构仍保留了原始实验代码的痕迹

本轮补充的注释和文档主要覆盖项目自有核心代码。第三方 vendored 目录为了避免无意义的大规模
diff，保持基本原样。

## TODO / Known Issues

- `vectorSHD_ED` 的 SIMD Levenshtein 路径与标量 `testLV` 还不一致；在 `testdata/smoke_input.txt` 上当前是 `3/6` 对 `5/6`
- `SIMD_ED` 的 Levenshtein correctness 还需要单独审计；`testNW 1 0 ...` 也会复现同类偏差
- `vectorSHD_ED` 的预编码 / SHD 快路径还没有做系统回归，当前不应作为 ground truth
- `needleman_wunsch-0.3.5/` 仍然是 vendored 老代码，现代编译器下还有不少 warning，没有做风格或结构清理
- 还没有完成正式的性能基线记录；目前只确认“能跑”，没有完成“与旧版本吞吐持平”的量化证明
- `-flto`、局部内存分配优化等低风险提速项还没开始做，需要在 correctness 固定后再评估
- `bowtie2_augment/` 和 `bowtie2_skipED_augment/` 这两套集成实验目录这一轮没有迁移
- 仓库当前会产生较多二进制和中间产物；清理规则和 `.gitignore` 还需要单独整理
