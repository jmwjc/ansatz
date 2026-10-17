# 快速开始

## 环境要求

- Rust 工具链（crate 使用 **edition 2024**，需要较新的 Rust 版本）；
- 一个 GMSH v4.1 ASCII 格式的 `.msh` 网格文件（2D/3D 示例所需的网格已随仓库提供在 `examples/msh/` 下）。

## 构建

```bash
cargo build
```

## 运行示例

仓库自带的端到端示例：1D 桁架（杆）在自重与右端集中拉力下的位移求解。

```bash
cargo run --example bar1d
```

输出大致如下：

```text
bar displacement: [2.0000000000000036e-7, 0.1950002000000004, 0.3800002000000008, ...]
Time elapsed: 2.616084ms
```

`bar1d` 不读网格文件：它直接按 11 个节点手工构造 10 个 `Seg2<1>` 单元。位移向量的最后一个分量约为 `1.5`，即杆右端在载荷作用下的位移，与解析解一致。

读取 GMSH 网格的 2D/3D 示例（需在**仓库根目录**运行，路径为相对路径）：

```bash
cargo run --example patchtest2d      # examples/msh/patchtest_tri3.msh、patchtest_quad4.msh
cargo run --example cantilever       # examples/msh/cantilever.msh
cargo run --release --example cook   # examples/msh/cook_{8,16,32,64}.msh
cargo run --example elasticity_cli -- examples/elasticity_cli.toml examples/msh/elasticity_cli.msh
```

## 运行测试

```bash
cargo test
```

## 性能剖析

`ansatz` 集成了 `dhat` 堆分配剖析，通过 feature 开关启用：

```bash
cargo run --release --example bar1d --features dhat-heap
```

运行后在仓库根目录生成 `dhat-heap.json`，可用 `dhat` 查看器打开。

## 静态检查

```bash
cargo clippy --all-targets
```

> 注意：当前 clippy 仍有若干历史警告，`PLAN.md` 的 M1 里程碑目标之一就是清零。
