#  ansatz

> 面向力学问题的离散近似与组装层：**你给出 ansatz，稀疏求解器负责剩下的。**

`ansatz` 把一个连续介质力学问题变成一个稀疏代数系统。*ansatz* 指"试探解的形式"：在有限元中它是网格上的分片多项式，在无网格方法中它是点云上的移动最小二乘 / 再生核近似。本库只对这一层负责——**近似（形函数）、数值积分、组装**——而不拥有线性求解器，解方程交给外部稀疏后端（默认 `faer`）。因此离散方式与线性代数相互独立。

```rust
use ansatz::element::quad4::Quad4;
use ansatz::operation::elasticity::{Lame, stiffness};
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
| 未知场（每节点 `C` 个自由度） | `Variable<D, C, P, N>` |
| `K` / `F` 的组装 | `operation::*` → `Triplet` / `Doublet` |
| 离散输入 | `io::gmsh` |
| 线性求解 | 外部：`faer` |

所有抽象都对空间维数 `D`、参数维数 `P`（`P < D` 覆盖嵌入高维空间的单元）、节点数 `N` 和每节点自由度数 `C` 泛化。`Variable` 不再编码张量阶，只描述“一组单元 + 每节点 `C` 个自由度”；具体分量的物理含义（标量、位移、应力等）由具体问题模块决定，因此同一个 `Variable` 骨架可以表示多变量自由度。

算子内部的物理量用**指标（求和）记号**书写，而不是显式拼出 Voigt `B`/`D` 矩阵——例如各向同性线弹性写作
`λ (∂ₚNₐ)(∂qN_b) + μ [δ_pq (∂ₖNₐ)(∂ₖN_b) + (∂qNₐ)(∂ₚN_b)]`，其中 `k` 对 `D` 个空间方向求和。因此同一个公式同时适用于 2D 平面应变、完整 3D，以及嵌入 3D 的壳类单元。

全部抽象在编译期完成：const 泛型 + 积分规则标记类型，无动态分发，也没有按维数复制的代码。

> **库的边界**：本库只提供近似、积分与组装原语。问题/求解器封装、约束施加策略由使用本库的主程序负责，不在本库范围内。

## 使用规则（初步）

以下约定是当前 API 的**初步契约**，在 1.0 之前可能随重构调整（破坏性变更会在 `CHANGELOG.md` 中标注）。

1. **场用 `Variable<D, C, P, N>` 描述**。`C` 是每节点自由度数/分量数：标量/温度 `C=1`、位移 `C=D`、对称二阶张量 `C=D(D+1)/2`、完整张量 `C=D^2` 等；`Variable` 不解释这些分量的张量含义。算子的第一个参数是一"组"单元：`&[T]` 与 `Vec<T>` 均可（已有 blanket impl），不必构造额外的容器类型。
2. **节点编号从 0 开始且连续，由前处理保证**。每个场使用独立网格，各自 renumber 到 `0..N-1`；节点 `n` 拥有 `m = C` 个连续自由度，占据 `n*m .. (n+1)*m`。矩阵维度为 `n_nodes * m`，其中 `n_nodes` 由前处理给出。本库不做节点重编号，也不提供多场 `DofMap`。
   > ⚠️ `Variable::num_dofs()` 只统计传入元素集合中的唯一节点 id，不能用于边界集合或子集网格的矩阵维度。非连续 id、孤立节点等问题应由前处理 renumber 和校验解决；详见下文「已知限制」第 1 条。
3. **刚度算子只产出下三角（含对角）**。因此组装后必须用 `Side::Lower` 做 `sp_cholesky`；不要假设拿到的是完整对称矩阵。
4. **本质边界条件用罚函数法**（`*_penalty` 算子），罚数取 `α ≈ 1e7`（相对材料刚度数量级选取）。罚函数同时贡献刚度与载荷，需分别 extend 到 `Triplets` 与 `Doublets`。
   > ⚠️ 罚函数是本库唯一提供的约束施加工具；精度与条件数受 `α` 影响。更高效的精确消元/约束施加策略由使用本库的主程序在组装后自行处理，本库不内置，也不计划提供全局施加流程。
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
let units: Vec<Quad4<2>> = mesh.elements("Ω").unwrap();
// 2) 组装：刚度（下三角）+ 载荷 + 罚函数约束（物理量以闭包给出）
let mut triplets = stiffness::<_, 2, 2, 4, 4>(&units, lambda, mu);
let mut doublets = body_force::<_, 2, 2, 4, 4>(&units, |_x| [0.0, -rho * g]);
let (pt, pd) = displacement_penalty(&fixed, 1e7, |_x| [0.0, 0.0]);
triplets.extend(pt);
doublets.extend(pd);
// 3) 求解：交给外部后端
let n = mesh.n_nodes() * 2; // 2D 位移场：n_nodes 由前处理 renumber 后给出
let k = SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets).unwrap();
let f = try_new_from_doublet(n, &doublets).unwrap();
let u = k.sp_cholesky(Side::Lower).unwrap().solve(&f);
```

## 已知限制与改进方向（使用评估）

> 本节记录对 0.2.0 使用体验的一次系统评估。数值内核（形函数、积分规则、指标求和组装）已较稳定：`cargo test` 86 项全过，`cargo clippy --all-targets` 零警告。以下问题主要影响**使用契约、健壮性与可扩展性**，按优先级排列，作为后续版本的输入。

### P0：结构性问题

#### 1. DOF 编号与节点 id 强耦合（不在本库解决，交由前处理）

所有算子直接用节点 id 生成全局自由度：

- 标量场：`dof = id`（`heat`、`bar`）；
- 向量场：`dof = id * D + p`（`elasticity`）；
- 高阶张量：`dof = id * m + k`，其中 `m = C`（每节点分量数）。

因此节点 id 必须稠密地落在 `0..N-1`。当前库不校验该条件，也不提供编号映射：

- 非连续 id 会让 `num_dofs()` 与 triplets 索引不一致；
- GMSH 读取只做 `tag - 1`，不保证 tag 连续；
- 示例使用 `max_id + 1`，与 `num_dofs()` 的口径不一致；孤立节点会引入零对角行；
- 如果不同场共用同一 id 空间会冲突（本库不处理，采用一场一网格规避）。

**决策**：本库不实现通用 `DofMap`，也不处理多物理场编号。后续采用“一场一网格、算子按场分别编写”的架构：

- 每个场（位移、转角、温度等）拥有独立网格，各自 renumber 到 `0..N_f-1`；
- 前处理负责节点编号、连通性重映射和校验；
- 库算子只按 `id*m+c` 组装，不做查表；
- 矩阵维度由前处理给出的 `n_nodes * m` 决定；
- 混合单元/非均质材料目前不纳入本库；如后期需要，由主程序按类型/区域新增不同的算子 `fn` 实现，不在库内引入 `enum`/`dyn` 运行时分发或 `Material` 抽象。

**前处理契约**（必须满足）：

1. 节点 id 为稠密的 `0..N-1`，`nodes[i].id == i`；
2. 该场的所有物理组（体、边界、点载荷、约束点）共用同一套映射；
3. 只保留被引用的节点，不要留下孤立节点；
4. 每个节点至少被一个“载体单元”（如体单元）覆盖，或显式约束，避免零对角行；
5. 保留原始 tag → dense id 的映射，用于输出与物理组查询；
6. `n_dofs = n_nodes * m`，并检查 `n_nodes * m` 不溢出；
7. 如果将来两个场要进入同一个全局矩阵，由驱动层加 `offset` 即可，不需要通用 `DofMap`。

**注意**：`Variable::num_dofs()` 只统计传入元素集合中的唯一 id，不能用于边界集合或子集网格的矩阵维度。矩阵维度统一使用前处理的 `n_nodes * m`。

### P1：健壮性与 API 细节

#### 2. GMSH 解析的存储与契约问题

当前状态：

- 节点 tag 减 1 后直接作为 id；非连续 tag 产生的 DOF 空洞由**前处理 gmsh renumbering** 解决，本库不再重复实现；
- tag 为 0 现在返回 `MeshError::InvalidNodeTag`，不再下溢 panic；
- `nodes` 仍是 `HashMap<usize, Node<D>>`，用于保留原始 tag 语义；
- `FromGmsh` 改为接收稠密 `NodeTable<D>`，元素构造不再在热路径查 `HashMap`；
- 物理组支持按 tag 升序迭代、按名称/标签查询，名称查找使用 O(1) 反查表。

##### 已实现：`NodeTable` + 零哈希 `FromGmsh`

- `NodeTable<D>` 只保存 `Vec<[f64; D]>`，满足 `id == index`，不保存冗余的 `Node.id`；
- 构造入口：
  - `NodeTable::new(Vec<[f64; D]>)`：直接接收前处理的稠密坐标；
  - `NodeTable::from_nodes(Vec<Node<D>>)`：前处理 renumber 后的稠密节点，校验 `node.id == index`；
  - `NodeTable::from_map(&HashMap<usize, Node<D>>)`：兼容 `GmshMesh` 的原始节点表，要求 id 为 `0..len`；
- `FromGmsh::from_gmsh(ids, &NodeTable<D>)`：用 O(1) 下标访问取坐标，不再依赖 `HashMap`；
- `GmshMesh::node_table()`：从当前节点表构造稠密表；
- `GmshMesh::elements_by_tag_with_table(tag, &table)`：复用同一张表构造多个物理组的单元，避免重复 O(N) 构建；
- `GmshMesh::elements` / `elements_by_tag` 保留为兼容 wrapper，内部构建一次 `NodeTable`；
- 新增错误：`NodeIndexOutOfRange`、`NodeIdMismatch`。

示例：

```rust
let table = mesh.node_table()?;
let boundary = mesh.physical_group_by_name("Γᵍ")?;
let segs: Vec<Seg2<2>> = mesh.elements_by_tag_with_table(boundary.tag, &table)?;
```

##### 已实现：物理组可迭代视图

```rust
pub struct PhysicalGroup<'a> {
    pub tag: i32,
    pub name: Option<&'a str>,
    pub element_type: i32,
    pub element_count: usize,
}
```

新增方法：

```rust
impl<const D: usize> GmshMesh<D> {
    pub fn physical_groups(&self) -> impl Iterator<Item = PhysicalGroup<'_>>;
    pub fn physical_group(&self, tag: i32) -> Result<PhysicalGroup<'_>, MeshError>;
    pub fn physical_group_by_name(&self, name: &str) -> Result<PhysicalGroup<'_>, MeshError>;
    pub fn physical_names(&self) -> impl Iterator<Item = (i32, &str)>;
    pub fn has_physical_group(&self, name: &str) -> bool;
}
```

- `physical_groups()` 按 tag 升序迭代，保证确定性；
- 未命名组（`name = None`）也包含在内；没有单元的命名组返回 `element_type = 0`、`element_count = 0`；
- 解析 `$PhysicalNames` 时建立 `name_to_tag` 反查表，`physical_group_by_name` 为 O(1)；
- 名称重复返回 `MeshError::DuplicatePhysicalName`；
- 未知名称返回 `UnknownPhysicalName`，未知 tag 返回 `UnknownPhysicalTag`。

示例：

```rust
for group in mesh.physical_groups() {
    println!(
        "tag={} name={:?} type={} n={}",
        group.tag, group.name, group.element_type, group.element_count
    );
}
```

##### 暂不处理

- **消除每单元坐标副本**：owned 元素仍保存 `[Node; N]`，`Quad4Ref` / 通用 `ElementRef` 借用视图留待后续；当前 `NodeTable` 只解决查表方式，不减少元素内部坐标副本；
- `GmshMesh.nodes` 仍为 `HashMap`，`elements_by_tag` 每次调用会构建一次 `NodeTable`（O(N)）；多次 materialize 时请使用 `node_table()` + `elements_by_tag_with_table` 复用。

> P1 中「类型派生与访问器」和「签名与命名一致性」已按方案 B 实现：`Node`/单元补齐 `Debug`/`Clone`/`PartialEq`，`GmshMesh` 提供摘要式 `Debug`；`bar::point_load`、`bar::point_displacement_penalty`、`heat::h1_semi_error`、`elasticity::Lame` 已替代旧名称。

### P2：性能与工程化

#### 3. 每个积分点重复计算 Jacobian 和形函数

算子目前普遍写成：

```rust
let dshape = elm.derivative_shape(ξg); // 内部算一次 J（Tri3/Quad4 还要求逆）
let jacobe = elm.jacobe(ξg);           // 又算一次 J
let fg = f(elm.coordinates(ξg));       // coordinates 内部又调一次 shape()
```

- `Tri3`/`Quad4` 的 `derivative_shape` 与 `jacobe` 各自调用 `jacobe_mat`；`Tri3` 的 J 是常量，却被每个积分点重复计算；
- `Seg2::derivative_shape` 内部调用 `jacobe` 求长度，算子又调用一次 `jacobe`；
- `coordinates` 重新调用 `shape`，于是体力/面力每点算两次形函数。

**改进方向**：给 `Element` 增加合并求值接口，一次返回 `shape`、`dshape`、`jacobe`、`coordinates`：

```rust
pub struct Evaluation<const D: usize, const P: usize, const N: usize> {
    pub shape: [f64; N],
    pub dshape: [[f64; D]; N],
    pub jacobe: f64,
    pub coordinates: [f64; D],
}

fn evaluate(&self, xi: [f64; P]) -> Evaluation<D, P, N>;
```

#### 4. `parallel` feature 未接入

`Cargo.toml` 声明了 `parallel` / `rayon`，但组装循环仍是单线程。单元级组装天然适合并行，建议按元素块 `par_iter` 收集，并保证结果顺序确定；同时用基准测试确认小网格上不会因并行反而变慢。

#### 5. Triplet 的内存与去重

当前返回 `Vec<Triplet>` 后由 faer 排序去重，内存峰值约为最终 nnz 的 1.5–2 倍。短期可在 `Assembler` 中提供 `with_capacity` 和按块 `reserve`；中期考虑直接组装 CSC 或“排序后合并”的可选路径。另需在 rustdoc 中明确写出“重复 triplet 会被求和”的隐式行为。

#### 6. GMSH I/O 的分配

`read_gmsh` 每行 `String::clear + read_line`，再 `split_whitespace().collect::<Vec<&str>>()`；节点表用 `HashMap`，连通性用 `Vec<Vec<usize>>`。大网格上分配次数和缓存不友好。建议改用 `read_until` + 字节解析，节点用 `Vec<Node<D>>`，连通性用扁平 `Vec<usize>` + offset。

#### 7. 文档、测试与工程化

- `cargo doc --no-deps` 目前有一个 broken intra-doc link（`CONTAINS_TOL` 为私有常量）；
- 没有 doctest，README 中的完整链路不会随 API 改动自动编译检查；
- 没有 `tests/` 集成测试；错误路径目前只在 `io::gmsh` / `sparse_vector` / `operation` 的单元测试中覆盖，尚未形成公开 API 的集成测试；
- `sparse`/`std` 是空 feature，`parallel` 未用；建议在 rustdoc 标注占位或直接移除；
- 没有 benchmark，PLAN M5 的性能目标无法回归验证；
- 建议 CI 增加 `cargo test --all-targets`（覆盖 examples）与 `cargo doc --no-deps`。

### 建议的落地顺序

| 阶段 | 内容 | 理由 |
|---|---|---|
| **0.3** | `Result` 错误类型 + GMSH 解析健壮性 + 前处理契约校验 | 先消除 panic 与静默少算 |
| **0.4** | 后处理（gradient/strain/stress + 节点磨光）+ `Variable` 批量求值/定位器 | 力学用户的核心工作流 |
| **0.5** | `evaluate()` 合并求值 + rayon 并行 + GMSH I/O 重写 + benchmark | 性能与大规模可用性 |
| **0.6** | 类型派生/访问器 + 文档/测试/feature 清理 | 收尾与发布准备 |

> P0.1 不在本库实现：节点 renumber 由前处理完成。若未来两个场需要进入同一全局矩阵，由驱动层加 `offset` 即可，不需要通用 `DofMap`。

如果只能先做一件事：**先做 `Result` 错误处理 + GMSH 解析健壮性**，避免 panic 和静默少算。问题/求解器封装、约束施加策略、混合单元/非均质材料均不在本库范畴：前两者由主程序负责，后者如后期需要则按类型/区域新增不同的算子 `fn` 实现。

## 后处理与空间索引

- 点查询：`Variable::{value, derivative_value}` 接受物理点 `x`；内部通过
  `Variable::locate(x)` 线性扫描单元，不存储空间索引。
- `value` 返回 `[f64; C]`，`derivative_value` 返回 `[[f64; C]; D]`
  （`derivative_value[i][k] = ∂u_k/∂x_i`）；均为固定大小数组，无堆分配。
- `gradient` / `divergence` / `curl` 不再定义在 `Variable` 上，而是由具体问题模块实现；
  当前 `heat::{gradient, flux}`、`elasticity::{strain, stress}` 提供局部量。
- `Variable::locate(x)` 返回 `Option<(&Item, [f64; P])>`，可用于在同一个点上继续调用
  问题相关的 `strain` / `stress` / `flux` 等函数。
- 本库**不输出节点导数或节点平均场**：有限元只有 C0 连续，节点处的导数在数学上不唯一；
  `derivative_value` / `strain` / `stress` / `flux` 仅作为误差估计与科研后处理的局部量。
- 空间索引：独立模块 `spatial`，提供 `Aabb<D>` 与后端无关的
  `SpatialIndex<T, D>`（当前为均匀网格，后续可替换 BVH）；
  供无网格邻域搜索或需要索引的调用方显式使用，FEM 点查询默认不做索引。
- `Aabb` 是计算几何/碰撞检测中的标准名称；公开 facade 使用
  `SpatialIndex`，避免把具体后端写进 API。AABB 只作为 broad phase，
  精确判定仍由 `local_coordinates` 或核函数完成。

参考：

- Christer Ericson, *Real-Time Collision Detection*, 2005；
- Ingo Wald, *On Fast Construction of SAH-based Bounding Volume Hierarchies*, 2007；
- Matthias Teschner et al., *Optimized Spatial Hashing for Collision Detection of Deformable Objects*, 2003；
- Mark de Berg et al., *Computational Geometry: Algorithms and Applications*；
- G. R. Liu, M. B. Liu, *Smoothed Particle Hydrodynamics: A Meshfree Particle Method*, 2003；
- W. K. Liu, S. Jun, Y. F. Zhang, *Reproducing kernel particle methods*, 1995。

## 当前状态

早期原型（0.3.0）。已可端到端运行：

- 单元：`Seg2<D>`、`Tri3<D>`、`Quad4<D>`、`Poi1<D>`；
- 积分：1D Gauss–Legendre，2D Dunavant / 张量积，3D Keast / 张量积；
- 算子：`bar`（桁架）、`heat`（导热 + 热容 + Robin/Dirichlet）、`elasticity`（平面应力/应变与 3D），本质边界条件用罚函数（限制见上文「已知限制」）；
- 后处理：`Variable::{value, derivative_value}`、`Variable::locate`、`spatial::{Aabb, SpatialIndex}`、`heat::{gradient, flux}`、`elasticity::{strain, stress}`；不提供节点导数/节点平均输出；
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
- `CHANGELOG.md` —— 版本变更记录（Keep a Changelog 格式）；
- `PLAN.md` —— 开发计划、里程碑与发布前安排；
- `AGENTS.md` —— 面向 AI 编码代理的架构与约定说明。

## 许可证

见 `LICENSE`（MIT）。
