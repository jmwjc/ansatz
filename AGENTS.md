# AGENTS.md

本文件面向 AI 编码代理，介绍 `ansatz` 项目的架构、构建方式与开发约定。项目暂无正式 README（仅有标题），`PLAN.md` 是主要的项目文档，本文件基于对实际代码的考察编写。

## 1. 项目概述

`ansatz` 是一个用 **Rust** 编写的**离散近似层**（不拥有求解器）：它把连续力学问题离散为稀疏代数系统——负责近似（形函数）、数值积分、组装，把 `Triplet`/`Doublet` 交给外部稀疏后端（默认 `faer`）。目标是同时容纳有限元与后续的无网格离散。

### 技术栈与关键配置（`Cargo.toml`）

- 语言/工具链：Rust，**edition 2024**，crate 名 `ansatz`，版本 0.3.0。
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
├── sparse_vector.rs  # Doublet（类比 faer 的 Triplet）与 try_new_from_doublet → Result<Col<f64>, DoubletError>
├── spatial/          # 独立 AABB / SpatialIndex（FEM 点定位与无网格邻域搜索共用）
├── element/
│   ├── mod.rs        # 核心 trait 与积分规则参考数据（约 1100 行，含单元测试）
│   ├── seg2.rs       # Seg2<D>：2 节点线单元（1D/2D/3D）
│   ├── tri3.rs       # Tri3<D>：3 节点三角形
│   ├── quad4.rs      # Quad4<D>：4 节点四边形
│   └── poi1.rs       # Poi1<D>：1 节点点单元（集中力/约束作用点）
├── operation/
│   ├── mod.rs        # Variable<D, C, P, N> trait（点查询/导数）与 Triplets/Doublets 类型别名
│   ├── bar.rs        # 1D 桁架：stiffness / point_load / body_force / point_displacement_penalty
│   ├── heat.rs       # 热传导：组装算子 + gradient / flux
│   └── elasticity.rs # 线弹性：组装算子 + strain / stress
└── io/
    ├── mod.rs        # （目前为空文件）
    └── gmsh.rs       # GMSH v4.1 读取：read_gmsh、GmshMesh、NodeTable、PhysicalGroup、FromGmsh
examples/
├── msh/              # 示例与测试用 GMSH 网格（bar1d、patchtest_tri3/quad4、cantilever、cook_*、elasticity_cli）
├── bar1d.rs          # 1D 桁架（手工构造节点）+ 罚函数 + 稀疏 Cholesky
├── patchtest2d.rs    # Tri3/Quad4 patch test
├── cantilever.rs     # 平面应力悬臂梁，与解析解对比
├── cook.rs           # Cook 膜平面应变，full/SRI 收敛对比
└── elasticity_cli.rs # 由 TOML 配置驱动的线弹性 CLI（+ elasticity_cli.toml）
book/                 # mdbook 用户手册（book.toml、src/*.md）
```

### 核心抽象（理解代码前必读）

- **`Element<const D: usize, const P: usize, const N: usize>`**（`src/element/mod.rs:6`）：所有单元的核心 trait。`D` = 空间维数，`P` = 参数空间维数（0=点、1=线、2=面），`N` = 节点数。方法包括 `id()`、`shape()`、`derivative_shape()`、`vertices_coordinates()`、`coordinates()`、`jacobe()`（注意拼写为 `jacobe`，不是 `jacobian`）、`jacobe_mat()`。
- **`IntegrationScheme<const P: usize, const G: usize>`**：积分点与权重（关联类型风格，编译期常量）。`FullIntegration` / `ReducedIntegration` 是标记子 trait，用于区分完全积分与减缩积分。
- **积分规则是标记类型**（`GaussLegendre1D<N>`、`Dunavant<G>`、`TensorProductQuad<N>`、`Keast<G>`、`TensorProductHex<N>`、`Point`）：`IntegrationScheme` 的 impl 绑定到这些标记类型，以便同一 `(P, G)` 下区分不同拓扑。具体单元（如 `Seg2<D>`）再为自身实现 `IntegrationScheme` 并委托给对应标记类型。
- **`Variable<D, C, P, N>`**（`src/operation/mod.rs`）：组装算子的输入抽象，`type Item: Element<D, P, N>`，对 `&[T]` 和 `Vec<T>` 有 blanket impl。`C` 是每节点自由度数/分量数；`Variable` 不编码张量阶，具体分量语义由问题模块解释。新增算子时照此模式接收一组单元。
- **`FromGmsh<D>`**（`src/io/gmsh.rs`）：单元类型与 GMSH element type 编号的桥接，声明 `GMSH_ELEMENT_TYPE`、`GMSH_NODE_COUNT`，并实现 `from_gmsh(ids, &NodeTable<D>) -> Result<Self, MeshError>`。`NodeTable` 是稠密节点坐标表（`id == index`），`GmshMesh::physical_groups()` 提供按 tag 升序的物理组迭代。新增单元接入 GMSH 读取时，在 io 层实现此 trait 即可（当前支持 Seg2=1、Tri3=2、Quad4=3、Poi1=15）。
- **组装风格**：算子返回 `Vec<Triplet<usize, usize, f64>>`（刚度）或 `Vec<Doublet<usize, f64>>`（载荷），由调用方用 `faer::sparse::SparseColMat::try_new_from_triplets` 和 `sparse_vector::try_new_from_doublet`（返回 `Result`）合并。`read_gmsh`、`GmshMesh::elements*`、`Variable::{value, derivative_value}` 也返回 `Result`，公开路径不使用 `unwrap`/`expect`/`panic!`。本质边界条件用**罚函数法**（`displacement_penalty`，罚数 α≈1e7），自然边界条件（体力/集中力）由 `body_force` / `traction` 组装。

`examples/` 下的示例均为端到端可运行：`bar1d`（手工构造 10 个 `Seg2<1>`，左端固定、均布体力 f=1、右端拉力 P=1，罚函数加约束后 `sp_cholesky` 求解）、`patchtest2d`（Tri3/Quad4 patch test）、`cantilever`（平面应力，与解析解比较 L²/能量范数误差）、`cook`（平面应变收敛对比）、`elasticity_cli`（TOML 配置驱动的 CLI）。

## 3. 构建与测试命令

```bash
cargo build              # 编译库
cargo run --example bar1d   # 1D 桁架示例（不需要网格文件）
cargo test               # 运行全部单元测试（目前 86 个，全部内联在各模块的 #[cfg(test)] 中）
cargo clippy --all-targets  # 静态检查（当前零警告，请勿新增）
cargo run --release --example bar1d --features dhat-heap   # 堆剖析；生成 dhat-heap.json（仓库根目录已有一份示例输出）
```

- 没有 `tests/` 目录、没有 benches；示例位于 `examples/` 目录（用 `cargo run --example <name>` 运行）；测试全部为模块内 `#[cfg(test)] mod tests`。
- 没有 CI（`.github/` 不存在）；PLAN.md M1 计划添加 GitHub Actions（fmt + clippy + test）。
- 没有部署/发布流程；M6 才计划 crates.io 发布。

## 4. 测试策略与约定

- **测试全部内联**：`src/element/mod.rs`（积分规则矩验证，如权重和、∫ξ²、Dunavant 33 点高精度）、各单元文件（形函数、雅可比、面积/体积）、`src/io/gmsh.rs`（用 `examples/msh/` 下的真实网格文件做解析测试，这里用到 dev-dependency `tempfile`）、`src/operation/{heat,elasticity}.rs`（patch test 与解析解回归）。
- 浮点比较用近似相等（容差约 1e-12），项目里有现成的 `close(a, b)` 辅助函数可参考（`src/element/mod.rs:1041`）。
- 新增单元时应为形函数、导数、雅可比写单元测试；新增物理模型时应先写回归测试（PLAN.md 设计原则：测试驱动）。
- `examples/msh/` 下的网格为 GMSH **v4.1** ASCII 格式；运行测试和示例依赖这些文件的相对路径，工作目录须为仓库根。

## 5. 代码风格与提交规范

- 源码（`src/` 与 `examples/`）的 rustdoc、注释统一使用**英文**；提交信息使用中文；示例 mesh 的物理组名使用 Unicode（如 `"Γᵍ"`、`"Ω"`）。
- 大量使用 `const` 泛型（`D/P/N/G`）和返回 `impl Iterator` 的 trait 方法（RPITIT）——新代码应保持这一风格，避免引入运行时动态分发。
- **提交信息遵循 [Conventional Commits](https://www.conventionalcommits.org/)**，详见 `PLAN.md` 第 4 节：`type(scope): subject` 格式，type 包括 `feat`/`fix`/`docs`/`style`/`refactor`/`perf`/`test`/`chore`/`build`/`ci`/`revert`，破坏性变更加 `!` 并在 footer 说明。
- 分支命名：`feat/<scope>-<description>`、`fix/<scope>-<description>`、`docs/<description>`。
- 设计原则（PLAN.md 第 5 节）：trait 优先、最小化修改、测试驱动、文档同步、性能可测。
- **几何量优先使用张量表示**：物理梯度、雅可比等几何量的推导与实现优先采用张量记号——协变切向量 `j_α = ∂x/∂ξ^α`、度量张量 `g_{αβ} = j_α · j_β`、逆度量 `g^{αβ}`，物理梯度即 `∇ₓN = j_α g^{αβ} ∂_β N`（对重复指标求和）。这样对任意空间维度 `D`（含高维嵌入，如空间三角形/四边形嵌入 3D）自然成立，无需引入运行时辅助函数或显式矩阵求逆。参考 `Tri3`/`Quad4` 的 `derivative_shape` 实现。
- **算子运算用张量形式（指标求和）**：组装算子内部的物理量运算（应变-位移、本构关系、刚度/载荷累加等）同样优先采用张量记号与指标求和，而不是显式拼 Voigt 矩阵。例如各向同性线弹性双线性形式写成 `ε(u) : C : ε(v) = λ (∂ₚNₐ)(∂qN_b) + μ [δ_pq (∂ₖNₐ)(∂ₖN_b) + (∂qNₐ)(∂ₚN_b)]`（对重复指标 `k` 求和），实现中直接对 `derivative_shape` 返回的物理梯度分量做显式求和，避免引入 B/D 矩阵等运行时辅助结构。参考 `src/operation/elasticity.rs` 的 `stiffness`。
- 注意：`src/io/mod.rs` 当前为空（0 字节），属正常状态；`jacobe` 等既有命名（含非常规拼写）是有意保留的，不要"修正"它们，以免破坏 API 一致性。

## 6. 已知缺口（与 PLAN.md 里程碑对照）

- M1：CI 未配置；`cargo clippy --all-targets` 当前零警告。
- M2：`Problem` trait、`Truss1D`、求解器封装尚未建立（按定位，求解器应保持为可选外部后端）。
- M3/M4：`Tri3`/`Quad4`/`Seg2` 已实现，`heat` 与 `elasticity`（平面应力/应变/3D）算子及 2D 算例已可用；`Material` trait 尚未抽象。
- M5：`parallel` feature 已声明但未使用；无基准测试。
- M6：`book/` 手册与 README 已建立，`CHANGELOG.md` 已建立（0.3.0）。