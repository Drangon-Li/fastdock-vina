# FastDock Vina

FastDock Vina 是基于 AutoDock Vina 的 CPU 优化版本，保留双精度计算和原有对接搜索设置，主要优化评分查表、对象传递和 BFGS 矩阵访问。

**A CPU-optimized AutoDock Vina derivative with double-precision calculations and unchanged docking search settings.**

本仓库发布选定的最终 D 版本源码、构建工具、最小运行示例和紧凑的基准证据。上游基线为 [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) 提交 `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645`，版本字符串为 `v1.2.7-27-g3c65c0b`。来源与修改记录见 [PROVENANCE](docs/PROVENANCE.md)。

## 加速效果

历史服务器实验在 19 个完整对接案例上，每个版本、每个案例重复运行 3 次。逐案例取墙钟时间中位数，再求和：

| 版本 | 19 个案例的中位数之和 | 相对上游耗时降低 |
| --- | ---: | ---: |
| 上游基线 | 3242.525 秒 | — |
| 第一阶段优化 v1 | 3074.107 秒 | 5.19% |
| 本仓库最终 D | 2989.231 秒 | **7.81%** |

最终 D 在 19/19 个完整对接案例上快于上游；相对 v1，总耗时再降低 **2.76%**，18/19 个案例更快，`basic-ad4-dock` 慢约 0.36%。这里的百分比来自上述中位数之和，不代表每个输入都能获得同样加速。

![逐案例墙钟加速比](docs/figures/wall-speedup.svg)

完整对接的峰值 RSS 相对上游降低 28.71%–73.27%，其中包含 v1 已有优化的贡献。D 相对 v1 的内存变化很小，且并非一致降低。实验环境、原始重复数据、正确性检查及适用边界见 [BENCHMARKS](docs/BENCHMARKS.md)。

## 如何优化

相对上游，核心修改集中在四个文件。D 在继承 v1 的基础上，只继续修改 BFGS 实现。

| 文件 | 修改及作用 |
| --- | --- |
| [`src/lib/precalculate.h`](src/lib/precalculate.h) | 相同的有序 XS 原子类型对共享查找表，减少重复构造和存储；保留按原子调用的重载与原来的操作数顺序。 |
| [`src/lib/cache.cpp`](src/lib/cache.cpp) | 缓存查找表指针，减少网格填充热循环中的重复查找。 |
| [`src/lib/vina.cpp`](src/lib/vina.cpp) | 用 `std::move` 转移已构造对象，避免深拷贝。 |
| [`src/lib/bfgs.h`](src/lib/bfgs.h) | BFGS 逆 Hessian 改用按列存储的完整对称矩阵，以双元素 CPU SIMD 同时处理相邻的独立行，并保持每行原有的 `j` 累加顺序。 |

BFGS 存储由 `n(n+1)/2` 个 double 增至 `n²` 个 double，并增加临时缓冲区，以换取更规则的访问和计算。SIMD 使用 GCC/Clang 的 `vector_size` 扩展。该实现使用 CPU，不依赖 GPU，也未通过降低精度、减少搜索量来提速。

## 构建

本次发布提供 Linux/macOS 的 C++ 命令行程序，未打包 Python binding。需要 Python 3.8+、支持 C++14 和 `vector_size` 的 GCC 或 Clang，以及 Boost 头文件与 thread、serialization、filesystem、program_options 库。对于 Boost 1.89 之前的版本，构建脚本还链接 system 库。

先获取代码，再从仓库根目录执行构建命令：

```bash
git clone https://github.com/Drangon-Li/fastdock-vina.git
cd fastdock-vina
```

Ubuntu / Debian：

```bash
sudo apt-get install build-essential python3 libboost-all-dev
python3 tools/build.py --jobs 4
```

macOS（已安装 Xcode Command Line Tools 和 Homebrew）：

```bash
brew install boost
python3 tools/build.py --boost-prefix "$(brew --prefix boost)" --jobs 4
```

自定义编译器、Boost 安装位置、输出目录，或同时构建 `vina_split`：

```bash
python3 tools/build.py --cxx g++ --boost-prefix /path/to/boost --jobs 4 --output-dir build --with-split
./build/vina --help
```

系统 Boost 可省略 `--boost-prefix`。构建脚本默认输出 `build/vina`，使用 `-O3 -DNDEBUG -std=c++14 -pthread -fvisibility=hidden` 和上述版本字符串，不额外启用浮点数学改写选项。历史基准所用的私有 Boost 安装目录不属于本仓库。

实测工具链为服务器 GCC 11.4 + Boost 1.82，以及发布检查时的 Apple Clang 21 + Boost 1.82。脚本读取当前 Boost 版本，并按照 [Boost 1.89 的变更说明](https://www.boost.org/releases/1.89.0/)在新版本中省略已移除的 System 占位库；本次没有对 Boost 1.89+ 进行实际构建和数值验证。

## 运行示例

仓库包含 1iep 受体、配体及搜索框。先执行只评分示例：

```bash
./build/vina \
  --receptor examples/basic/1iep_receptor.pdbqt \
  --ligand examples/basic/1iep_ligand.pdbqt \
  --config examples/basic/1iep_receptor.box.txt \
  --scoring vina --score_only --cpu 1 --seed 20261022
```

完整对接示例：

```bash
mkdir -p demo_out
./build/vina \
  --receptor examples/basic/1iep_receptor.pdbqt \
  --ligand examples/basic/1iep_ligand.pdbqt \
  --config examples/basic/1iep_receptor.box.txt \
  --scoring vina --exhaustiveness 8 --max_evals 0 --num_modes 9 \
  --cpu 1 --seed 20261022 \
  --out demo_out/1iep_out.pdbqt
```

这里的 `exhaustiveness=8` 用于运行演示；历史基准使用 `exhaustiveness=32`。示例运行不能代替基准复现。比较版本时，应固定输入、随机种子、CPU 数、搜索参数、编译环境及计时方法。

发布前已在 macOS arm64 上新编译上游与本版本，并运行上述两个示例：评分和对接日志一致，完整对接产生的 4 个构象 PDBQT 逐字节一致。检查参数与文件哈希见 [`release-validation.json`](benchmarks/release-validation.json)；这次启动检查不计入上面的服务器性能统计。

## 正确性与适用范围

历史检查覆盖 225 次 CLI 执行，其 PDBQT 输出逐字节一致、日志结果等价；另有 75 次 API 检查，其十六进制 double 能量值与可获取的坐标完全一致。这些检查支持已测输入和配置下的等价性，不能证明所有输入和平台上的普遍等价。

多配体、柔性受体及局部搜索的坐标 getter 覆盖不完整；Python binding 与 `--write_maps` 未测试。另有历史非有限值案例出现 NaN 符号差异，因此不宣称所有 IEEE 浮点特殊值都逐位一致。更多说明见 [BENCHMARKS](docs/BENCHMARKS.md)。

## 许可与致谢

本项目保留上游 [Apache License 2.0](LICENSE) 及源码署名。AutoDock Vina 的原始开发与持续维护工作来自 Oleg Trott、Arthur J. Olson、Stefano Forli 及 Scripps 团队等贡献者。本仓库的 CPU 优化基于这些工作；四个修改文件保留了修改说明。

使用本项目开展研究时，请引用 AutoDock Vina 原始论文及 1.2.0 论文，并注明本仓库所用提交：

- Trott O, Olson AJ. *AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading.* Journal of Computational Chemistry, 2010. [doi:10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334)
- Eberhardt J, Santos-Martins D, Tillack AF, Forli S. *AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings.* Journal of Chemical Information and Modeling, 2021. [doi:10.1021/acs.jcim.1c00203](https://doi.org/10.1021/acs.jcim.1c00203)
