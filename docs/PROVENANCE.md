# 源码与数据来源

本仓库发布选定的 **D_column_simd2**，以 AutoDock Vina 提交
`3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645` 为原始基线。
这里仅保存一套最终计算源码；性能表中的 `baseline` 和 `v1` 是历史实验的对照标签，不是附带的另一套实现。

## 最终实现的身份

实验归档中的选定源码位于 `final-D/sources/v2`。归档外层的 `sources/v2` 是另一候选，未用于本次发布。
在冻结实验源码中，D 相对 v1 仅修改 `src/lib/bfgs.h`；相对上游原版共修改以下四个文件：

| 文件 | 最终版包含的修改 |
| --- | --- |
| [`precalculate.h`](../src/lib/precalculate.h) | 有序 XS 原子类型对共享预计算表，保留评分运算 |
| [`cache.cpp`](../src/lib/cache.cpp) | 网格计算复用类型对表指针 |
| [`vina.cpp`](../src/lib/vina.cpp) | 移动大对象，避免深拷贝 |
| [`bfgs.h`](../src/lib/bfgs.h) | 完整对称列存储和双 double CPU SIMD |

实验冻结的 `bfgs.h` SHA-256 为：

```text
3a60d72c0347a6edc3d642159ac984bb4491bebc65440c4fe394e25ef9f5c0c7
```

按上游 Apache 2.0 第 4(b) 条的要求，发布时在这四个文件中补充了明确的修改声明注释，保留上游版权与许可头。
除此之外，全部 63 个 C++ 源文件／头文件均来自选定冻结版本，计算代码没有改写。
[`source-manifest.json`](../benchmarks/source-manifest.json) 同时记录每个文件的**实验冻结哈希**和**发布哈希**，并标记四个新增声明的文件；不能把两类哈希混为一谈。
三个示例输入也逐字节取自同一上游快照。

新增的 `tools/build.py` 是独立构建入口，采用实验中的优化、精度和版本参数。
它不修改计算源码，也不打包 Boost 或二进制。CLI 的 `--version` 保留上游构建版本
`v1.2.7-27-g3c65c0b`；请通过本仓库 Git 提交和源码清单识别 FastDock 最终版。

## 测量数据

历史性能测量完成于 2026-10-15，采用独立 D 轮次中的原版、v1、D 三方结果。

| 文件 | 来源与范围 |
| --- | --- |
| [`benchmark.csv`](../benchmarks/benchmark.csv) | 原始 D 汇总 CSV 的逐字节副本；含 25 配置的中位数、范围与全部重复测量值 |
| [`runs.csv`](../benchmarks/runs.csv) | 从 225 条正式记录抽取计时、CPU、峰值 RSS、运行顺序、二进制与输出哈希；省略机器路径和进程环境 |
| [`summary.json`](../benchmarks/summary.json) | 从同一记录重算分组结果、回退项和 RSS 变化，附实验参数及原始汇总文件哈希 |
| [`wall-speedup.svg`](figures/wall-speedup.svg)、[`peak-rss.svg`](figures/peak-rss.svg) | 同一 D 轮次报告生成的图表副本 |

原始完整 `benchmark.json` 的 SHA-256：

```text
5ccf388d626c0018f58e047c351d934b89d13ed09c9ed6842cf054e6d8452e1b
```

导出时核对了 225 个 `(case, variant, repeat)` 唯一标识、成功状态、输出哈希一致性，以及每项中位数与逐次记录的对应关系。
这里发布的是紧凑证据集，未打包所有旧版本、服务器目录、完整 stdout/stderr、API 原始数值档案或失败候选。
哈希用于标识与核对文件；它本身不能替代重新运行全部科学验证。
测量方式、正确性覆盖和限制见 [BENCHMARKS.md](BENCHMARKS.md)。

## 发布时的启动检查

[`release-validation.json`](../benchmarks/release-validation.json) 记录本次发布包的独立启动检查。
它在同一台 Mac、同一编译器和 Boost 下，新建原版与发布版二进制，执行 README 的固定种子评分及完整对接示例，比较日志与完整 PDBQT。
此检查不计入历史服务器性能表，也不把不同机器的时间合并。

上游作者与引文见 [README](../README.md#许可与致谢)；完整许可见 [LICENSE](../LICENSE)。
