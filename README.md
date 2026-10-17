# ansatz

> 面向力学问题的离散近似与组装层：**你给出 ansatz，稀疏求解器负责剩下的。**

`ansatz` 把一个连续介质力学问题变成一个稀疏代数系统。*ansatz* 指"试探解的形式"：在有限元中它是网格上的分片多项式，在无网格方法中它是点云上的移动最小二乘 / 再生核近似。本库只对这一层负责——**近似（形函数）、数值积分、组装**——而不拥有线性求解器，解方程交给外部稀疏后端（默认 `faer`）。因此离散方式与线性代数相互独立。

```rust
use ansatz::element::quad4::Quad4;
use ansatz::operation::elasticity::{lame_constants, stiffness};
```

## 安装

> **当前状态**：尚未发布到 crates.io（计划见 `PLAN.md` 第 7 节），请以 git 或 path 依赖引入。

```toml
[dependencies]
# 方式一：git 依赖（推荐，跟随主干）
ansatz = { git = "https://github.com/jmwjc/ansatz" }

# 方式二：本地 path 依赖（在仓库内开发时）
# ansatz = { path = "../ansatz" }
```

发布到 crates.io 之后改为：

```bash
cargo add ansatz
```

### 环境要求

- Rust **edition 2024**（需要 Rust 1.85+，`cargo --version` 检查）；
- 线性求解后端 `faer` 由本库依赖自动引入，无需单独声明；
- 需要读取网格时，网格文件须为 **GMSH v4.1 ASCII** 格式。

### Feature 开关

| feature | 默认 | 作用 |
|---|---|---|
| `sparse` | ✅ | 稀疏组装路径（保留开关，当前仅为占位） |
| `std` | ❌ | 同上，占位开关 |
| `dhat-heap` | ❌ | 启用 `dhat` 堆剖析（示例中通过 `#[global_allocator]` 接入） |
| `dhat-ad-hoc` | ❌ | 启用 `dhat` ad-hoc 剖析 |
| `parallel` | ❌ | 启用 `rayon`（已声明，尚未接入组装循环） |

```bash
cargo run --release --example bar1d --features dhat-heap   # 堆剖析
```

## 设计

| 关注点 | 抽象 |
|---|---|
| 局部近似（形函数、逆映射 `local_coordinates`/`contains`） | `Element<D, P, N>` |
| 数值积分（Gauss–Legendre、Dunavant、Keast、张量积规则） | `IntegrationScheme<P, G>` + 标记类型 |
| 未知场（标量、向量或更高阶张量） | `Variable<D, O, P, N>` |
| `K` / `F` 的组装 | `operation::*` → `Triplet` / `Doublet` |
| 离散输入 | `io::gmsh` |
| 线性求解 | 外部：`faer` |

所有抽象都对空间维数 `D`、参数维数 `P`（`P < D` 覆盖嵌入高维空间的单元）、节点数 `N`、以及未知场的张量阶 `O`（0 = 标量，1 = 向量，…）泛化。每节点自由度数由 `C(D + O - 1, O)` 得出，因此标量场与向量场共用同一套组装骨架。

算子内部的物理量用**指标（求和）记号**书写，而不是显式拼出 Voigt `B`/`D` 矩阵——例如各向同性线弹性写作
`λ (∂ₚNₐ)(∂qN_b) + μ [δ_pq (∂ₖNₐ)(∂ₖN_b) + (∂qNₐ)(∂ₚN_b)]`，其中 `k` 对 `D` 个空间方向求和。因此同一个公式同时适用于 2D 平面应变、完整 3D，以及嵌入 3D 的壳类单元。

全部抽象在编译期完成：const 泛型 + 积分规则标记类型，无动态分发，也没有按维数复制的代码。

## 使用规则（初步）

以下约定是当前 API 的**初步契约**，在 1.0 之前可能随重构调整（破坏性变更会在 `CHANGELOG.md` 中标注）。

1. **场用 `Variable<D, O, P, N>` 描述**。算子的第一个参数是一"组"单元：`&[T]` 与 `Vec<T>` 均可（已有 blanket impl），不必构造额外的容器类型。
2. **节点编号从 0 开始且连续**。节点 `n` 拥有 `m = C(D + O - 1, O)` 个连续自由度，占据 `n*m .. (n+1)*m`；总自由度数为 `Variable::num_dofs()`，可直接用来创建 `n×n` 矩阵。
3. **刚度算子只产出下三角（含对角）**。因此组装后必须用 `Side::Lower` 做 `sp_cholesky`；不要假设拿到的是完整对称矩阵。
4. **本质边界条件用罚函数法**（`*_penalty` 算子），罚数取 `α ≈ 1e7`（相对材料刚度数量级选取）。罚函数同时贡献刚度与载荷，需分别 extend 到 `Triplets` 与 `Doublets`。
5. **自然边界条件独立组装**：体力（`body_force`）、面力/热流（`traction` / `heat_flux`）、对流（`convection`）、热源（`heat_source`）各自返回载荷条目，与刚度无关。体力与罚函数的分布量以**闭包**给出（如 `|x| [0.0, -rho * g]`），便于直接表达位置相关的载荷。
6. **求解器是外部组件**。库的公共 API 只返回 `Triplet`/`Doublet`（以及误差标量），从不返回矩阵类型；换 `faer` 的 LU / Cholesky 或换别的后端都不需要改动本库。
7. **单元须声明积分类别**：组装算子要求 `T::Item: FullIntegration<P, G>`（或 `ReducedIntegration<P, G>`），`G` 由调用方在类型参数上给出（如 `stiffness::<_, 2, 2, 4, 4>(&quads, λ, μ)`）。
8. **物理公式用指标求和，不引入 Voigt 矩阵**。新增算子时直接对 `derivative_shape` 返回的物理梯度分量求和，保持对任意 `D` 成立。
9. **维数约定 `P ≤ D`**：`P = D` 为体/面单元，`P < D` 为嵌入单元（空间杆、壳），`P = 0` 为点单元（集中力/约束点）。
10. **`jacobe` / `jacobe_mat` / `Doublet` 等命名是有意保留的 API**，不做"拼写修正"。
11. **单位制自洽**：库不做单位检查，材料常数、几何、载荷须使用同一套单位。

一条完整链路的样子：

```rust
// 1) 场：一组单元（从物理组取单元，或手工构造）
let units: Vec<Quad4<2>> = mesh.elements("Ω");
// 2) 组装：刚度（下三角）+ 载荷 + 罚函数约束（物理量以闭包给出）
let mut triplets = stiffness::<_, 2, 2, 4, 4>(&units, lambda, mu);
let mut doublets = body_force::<_, 2, 2, 4, 4>(&units, |_x| [0.0, -rho * g]);
let (pt, pd) = displacement_penalty(&fixed, 1e7, |_x| [0.0, 0.0]);
triplets.extend(pt);
doublets.extend(pd);
// 3) 求解：交给外部后端
let n = Variable::<2, 1, 2, 4>::num_dofs(&units);
let k = SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets).unwrap();
let f = try_new_from_doublet(n, &doublets);
let u = k.sp_cholesky(Side::Lower).unwrap().solve(&f);
```

## 当前状态

早期原型（0.1.0）。已可端到端运行：

- 单元：`Seg2<D>`、`Tri3<D>`、`Quad4<D>`、`Poi1<D>`；
- 积分：1D Gauss–Legendre，2D Dunavant / 张量积，3D Keast / 张量积；
- 算子：`bar`（桁架）、`heat`（导热 + 热容 + Robin/Dirichlet）、`elasticity`（平面应力/应变与 3D），本质边界条件用罚函数；
- GMSH v4.1 读取（`io::gmsh`），经 `faer` 稀疏 Cholesky/LU 求解；
- 示例：`bar1d`、`patchtest2d`、`cantilever`、`cook`、`elasticity_cli`。

## 构建与测试

```bash
cargo build
cargo test
cargo run --example bar1d      # 在仓库根目录运行（示例不读网格）
cargo run --example patchtest2d  # 读取 examples/msh/*.msh，需在仓库根目录运行
cargo run --release --example cook
```

## 文档

```bash
mdbook serve book      # 或：mdbook build book
```

- `book/` —— 用户手册、理论推导与算例讲解（中文）；
- `PLAN.md` —— 开发计划、里程碑与发布前安排；
- `AGENTS.md` —— 面向 AI 编码代理的架构与约定说明。

## 许可证

见 `LICENSE`（MIT）。
