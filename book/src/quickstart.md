# 快速开始

## 环境要求

- Rust 工具链（crate 使用 **edition 2024**，需要较新的 Rust 版本）；
- 一个 GMSH v4.1 ASCII 格式的 `.msh` 网格文件（示例所需的 `test/msh/*.msh` 已随仓库提供）。

## 构建

```bash
cargo build
```

## 运行示例

仓库自带一个端到端示例：1D 桁架（杆）在自重与右端集中拉力下的位移求解。

```bash
cargo run --example bar1d
```

输出大致如下：

```text
segments (Γᵍ): 40
triangles (Ω): 288
bar displacement: [2.0000000000000055e-7, 0.19500020000000057, ...]
Time elapsed: 5.102084ms
```

第一行表明示例读入了 `test/msh/patchtest.msh` 中的 40 个线单元（物理组 `Γᵍ`）与 288 个三角形单元（物理组 `Ω`）；位移向量的最后一个分量约为 `1.5`，这是杆右端在载荷作用下的位移，与解析解一致。

> **注意**：示例从工作目录读取 `test/msh/patchtest.msh`，因此须在**仓库根目录**运行。

## 运行测试

```bash
cargo test
```

## 性能剖析

`rustfem` 集成了 `dhat` 堆分配剖析，通过 feature 开关启用：

```bash
cargo run --release --example bar1d --features dhat-heap
```

运行后在仓库根目录生成 `dhat-heap.json`，可用 `dhat` 查看器打开。

## 静态检查

```bash
cargo clippy --all-targets
```

> 注意：当前 clippy 仍有若干历史警告，`PLAN.md` 的 M1 里程碑目标之一就是清零。
