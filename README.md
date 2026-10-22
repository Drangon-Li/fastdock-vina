# FastDock Vina

基于 [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina/tree/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645) 的 CPU 优化版本，保留双精度计算和原有搜索设置，不使用 GPU。

## 相比原版的效果

在 Intel Xeon Platinum 8362 上测试 19 个完整对接案例，每个案例重复 3 次，取各案例耗时中位数后求和：

| 版本 | 总耗时 | 耗时降低 |
| --- | ---: | ---: |
| 原始 AutoDock Vina | 3242.525 秒 | — |
| FastDock Vina | 2989.231 秒 | **7.81%** |

19 个案例均快于原版，峰值内存（RSS）降低 **28.71%–73.27%**。已测配置的对接 PDBQT 输出与原版逐字节一致；结果仅适用于已测环境和输入。详细数据与验证范围见 [实验记录](docs/BENCHMARKS.md)。

## 优化方法

- 共享相同原子类型对的评分查找表，减少重复计算与存储。
- 缓存热循环中的查表指针，通过移动对象避免深拷贝。
- 调整 BFGS 矩阵布局，使用 CPU SIMD 同时处理独立行，保留原有累加顺序。

## 编译

需要 Python 3.8+、支持 C++14 的 GCC/Clang 和 Boost（已验证 1.82）。

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

## 运行

使用仓库自带的 1iep 示例：

```bash
mkdir -p demo_out
./build/vina \
  --receptor examples/basic/1iep_receptor.pdbqt \
  --ligand examples/basic/1iep_ligand.pdbqt \
  --config examples/basic/1iep_receptor.box.txt \
  --scoring vina --cpu 1 --seed 20261022 \
  --exhaustiveness 8 --max_evals 0 --num_modes 9 \
  --out demo_out/1iep_out.pdbqt
```

结果保存至 `demo_out/1iep_out.pdbqt`。更多参数见 `./build/vina --help`。示例使用 `exhaustiveness=8`，上述性能测试使用 `32`。

## 许可与致谢

遵循 [Apache 2.0](LICENSE)，保留 AutoDock Vina 的版权与署名。研究使用请引用 [Vina 原始论文](https://doi.org/10.1002/jcc.21334)和 [Vina 1.2.0 论文](https://doi.org/10.1021/acs.jcim.1c00203)。源码来源见 [PROVENANCE](docs/PROVENANCE.md)。
