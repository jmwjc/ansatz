# AGENTS.md

本文件面向 AI 编码代理，介绍 `rustfem` 项目的架构、构建方式与开发约定。项目暂无正式 README（仅有标题），`PLAN.md` 是主要的项目文档，本文件基于对实际代码的考察编写。

## 1. 项目概述

`rustfem` 是一个用 **Rust** 编写的通用有限元方法（FEM）库，目标是为结构力学/连续介质力学问题提供模块化的单元抽象、组装算子与稀疏线性求解。当前处于早期原型阶段：1D 桁架求解链路已可用（罚函数法施加本质边界条件 + faer 稀疏 Cholesky 求解），2D/3D 单元与积分方案框架已建立，物理问题扩展（平面应力/应变、线弹性等）仍在路线图上（见 `PLAN.md` 的 M1–M6 里程碑）。

### 技术栈与关键配置（`Cargo.toml`）

- 语言/工具链：Rust，**edition 2024**，crate 名 `rustfem`，版本 0.1.0。
- 核心依赖：
  - `faer = "0.24.4"` — 稠密/稀疏线性代数（`SparseColMat`、`Triplet`、`Col`、稀疏 Cholesky/LU）。版本升级需谨慎，PLAN.md 明确要求单独提交。
  - `dhat = "0.3.3"`（可选）— 堆分配性能剖析。
  - `rayon = "1.12.0"`（可选）— 并行。
- 开发依赖：`tempfile = "3"`（测试用）。
- `[profile.release] debug = 1` — release 构建保留调试符号，用于性能剖析。
- Features：
  - `default = ["sparse"]`；`sparse`、`std` 目前仅为占位开关；
  - `dhat-heap` / `dhat-ad-hoc` — 启用 dhat 剖析（`examples/bar1d.rs` 中已接入 `#[global_allocator]` 与 `Profiler`）；
  - `parallel` — 启用 rayon（暂未在代码中使用）。

## 2. 代码组织

库入口为 `src/lib.rs`。源码（`src/` 与 `examples/`）的注释和 rustdoc 统一使用英文。

```
src/
├── lib.rs            # 库入口，声明公共模块
├── node.rs           # Node<const D: usize> { id, coordinates: [f64; D] }
├── sparse_vector.rs  # Doublet（类比 faer 的 Triplet）与 try_new_from_doublet → Col<f64>
├── element/
│   ├── mod.rs        # 核心 trait 与积分规则参考数据（约 1100 行，含单元测试）
│   ├── seg2.rs       # Seg2<D>：2 节点线单元（1D/2D/3D）
│   ├── tri3.rs       # Tri3<D>：3 节点三角形
│   ├── quad4.rs      # Quad4<D>：4 节点四边形
│   └── poi1.rs       # Poi1<D>：1 节点点单元（集中力/约束作用点）
├── operation/
│   ├── mod.rs        # Variable<D, P, N> trait（对 &[T]/Vec<T>  blanket impl）
│   └── bar.rs        # 1D 桁架组装算子：stiffness / traction / body_force / displacement_penalty
└── io/
    ├── mod.rs        # （目前为空文件）
    └── gmsh.rs       # GMSH v4.1 网格读取：read_gmsh、GmshMesh、FromGmsh trait
test/msh/             # 测试用 GMSH 网格（patchtest.msh、patchtest1D.msh，含物理组 "Γᵍ"、"Ω"）
examples/
└── bar1d.rs          # 1D 桁架示例：GMSH 读取 + 组装 + 稀疏 Cholesky 求解（注释用英文）
```

### 核心抽象（理解代码前必读）

- **`Element<const D: usize, const P: usize, const N: usize>`**（`src/element/mod.rs:6`）：所有单元的核心 trait。`D` = 空间维数，`P` = 参数空间维数（0=点、1=线、2=面），`N` = 节点数。方法包括 `id()`、`shape()`、`derivative_shape()`、`vertices_coordinates()`、`coordinates()`、`jacobe()`（注意拼写为 `jacobe`，不是 `jacobian`）、`jacobe_mat()`。
- **`IntegrationScheme<const P: usize, const G: usize>`**：积分点与权重（关联类型风格，编译期常量）。`FullIntegration` / `ReducedIntegration` 是标记子 trait，用于区分完全积分与减缩积分。
- **积分规则是标记类型**（`GaussLegendre1D<N>`、`Dunavant<G>`、`TensorProductQuad<N>`、`Keast<G>`、`TensorProductHex<N>`、`Point`）：`IntegrationScheme` 的 impl 绑定到这些标记类型，以便同一 `(P, G)` 下区分不同拓扑。具体单元（如 `Seg2<D>`）再为自身实现 `IntegrationScheme` 并委托给对应标记类型。
- **`Variable<D, P, N>`**（`src/operation/mod.rs:5`）：组装算子的输入抽象，`type Item: Element<D, P, N>`，对 `&[T]` 和 `Vec<T>` 有 blanket impl。新增算子时照此模式接收一组单元。
- **`FromGmsh<D>`**（`src/io/gmsh.rs:19`）：单元类型与 GMSH element type 编号的桥接。新增单元接入 GMSH 读取时，在 io 层实现此 trait 即可（当前支持 Seg2=1、Tri3=2、Poi1=15）。
- **组装风格**：算子返回 `Vec<Triplet<usize, usize, f64>>`（刚度）或 `Vec<Doublet<usize, f64>>`（载荷），由调用方用 `faer::sparse::SparseColMat::try_new_from_triplets` 和 `sparse_vector::try_new_from_doublet` 合并。本质边界条件用**罚函数法**（`displacement_penalty`，罚数 α≈1e7），自然边界条件（体力/集中力）由 `body_force` / `traction` 组装。

`examples/bar1d.rs` 是目前唯一端到端示例：读取 `test/msh/patchtest.msh`，构造 10 个 `Seg2<1>` 组成的 1D 杆（左端固定、均布体力 f=1、右端拉力 P=1），罚函数加约束后 `sp_cholesky` 求解。

## 3. 构建与测试命令

```bash
cargo build              # 编译库
cargo run --example bar1d   # 运行 examples/bar1d.rs 中的 1D 桁架 + GMSH 读取示例
cargo test               # 运行全部单元测试（目前 13 个，全部内联在各模块的 #[cfg(test)] 中）
cargo clippy --all-targets  # 静态检查（注意：当前约 91 个警告，PLAN.md 的 M1 目标是清零，请勿新增）
cargo run --release --example bar1d --features dhat-heap   # 堆剖析；生成 dhat-heap.json（仓库根目录已有一份示例输出）
```

- 没有 `tests/` 目录、没有 benches；示例位于 `examples/` 目录（用 `cargo run --example <name>` 运行）；测试全部为模块内 `#[cfg(test)] mod tests`。
- 没有 CI（`.github/` 不存在）；PLAN.md M1 计划添加 GitHub Actions（fmt + clippy + test）。
- 没有部署/发布流程；M6 才计划 crates.io 发布。

## 4. 测试策略与约定

- **测试全部内联**：`src/element/mod.rs`（积分规则矩验证，如权重和、∫ξ²、Dunavant 33 点高精度）、各单元文件（形函数、雅可比、面积/体积）、`src/io/gmsh.rs`（用 `test/msh/` 下的真实网格文件做解析测试，这里用到 dev-dependency `tempfile`）。
- 浮点比较用近似相等（容差约 1e-12），项目里有现成的 `close(a, b)` 辅助函数可参考（`src/element/mod.rs:1041`）。
- 新增单元时应为形函数、导数、雅可比写单元测试；新增物理模型时应先写回归测试（PLAN.md 设计原则：测试驱动）。
- `test/msh/` 下的网格为 GMSH **v4.1** ASCII 格式；运行测试和示例依赖这些文件的相对路径，工作目录须为仓库根。

## 5. 代码风格与提交规范

- 源码（`src/` 与 `examples/`）的 rustdoc、注释统一使用**英文**；提交信息使用中文；示例 mesh 的物理组名使用 Unicode（如 `"Γᵍ"`、`"Ω"`）。
- 大量使用 `const` 泛型（`D/P/N/G`）和返回 `impl Iterator` 的 trait 方法（RPITIT）——新代码应保持这一风格，避免引入运行时动态分发。
- **提交信息遵循 [Conventional Commits](https://www.conventionalcommits.org/)**，详见 `PLAN.md` 第 4 节：`type(scope): subject` 格式，type 包括 `feat`/`fix`/`docs`/`style`/`refactor`/`perf`/`test`/`chore`/`build`/`ci`/`revert`，破坏性变更加 `!` 并在 footer 说明。
- 分支命名：`feat/<scope>-<description>`、`fix/<scope>-<description>`、`docs/<description>`。
- 设计原则（PLAN.md 第 5 节）：trait 优先、最小化修改、测试驱动、文档同步、性能可测。
- 注意：`src/io/mod.rs` 当前为空（0 字节），属正常状态；`jacobe` 等既有命名（含非常规拼写）是有意保留的，不要"修正"它们，以免破坏 API 一致性。

## 6. 已知缺口（与 PLAN.md 里程碑对照）

- M1：CI 未配置；`cargo clippy` 有约 91 个警告未清零。
- M2：`Problem` trait、`Truss1D`、求解器封装尚未建立；`examples/bar1d.rs` 是唯一示例。
- M3/M4：`Tri3`/`Quad4`/`Seg2` 已实现并测试，但尚未有 2D/3D 求解示例；`Material` trait、平面应力/应变问题未实现。
- M5：尚无 `Solver` trait；`parallel` feature 已声明但未使用；无基准测试。
- M6：README 未写、CHANGELOG 不存在。