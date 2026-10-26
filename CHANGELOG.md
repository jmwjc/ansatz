# Changelog

本文件记录 `ansatz` 的所有重要变更。

格式遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)，
版本号遵循 [语义化版本](https://semver.org/lang/zh-CN/)。
0.x 阶段的 MINOR 版本允许包含破坏性变更，破坏性变更会在对应版本的
「Breaking Changes」小节中单列并给出迁移方式。

## [Unreleased]

### Added

- `element::crtri3::CRTri3<D>`：Crouzeix–Raviart 非协调线性三角形，实现 `Element<D, 2, 3>`，
  自由度为三条边的平均值（基函数 `θᵢ = 1 − 2λᵢ`），2D 质量矩阵恰好对角 `(|T|/3) I`；
  `CRTri3::from_tri3(&tri, edge_id)` 由 `Tri3` 与"边 → 自由度"闭包构造（节点为边中点，
  `triangle_vertices()` 反射取回真实顶点），`FullIntegration<2, 3>` / `ReducedIntegration<2, 1>`；
- `element::tet4::Tet4<D>`：4 节点线性（协调）四面体，实现 `Element<D, 3, 4>`、`FromGmsh`
  （GMSH 类型 4），`volume()`、`FullIntegration<3, 4>` / `ReducedIntegration<3, 1>`；
- `element::crtet4::CRTet4<D>`：Crouzeix–Raviart 非协调线性四面体，实现 `Element<D, 3, 4>`，
  自由度为四个面的平均值（`θₐ = 1 − 3λₐ`），几何由面重心按 `ΣθₐCₐ = ΣλₐAₐ` 精确重现
  （`tetrahedron_vertices()` 反射取回顶点），质量为非对角（`2|T|/5` / `−|T|/20`，
  需 `FullIntegration<3, 4>`），`from_vertices([Node;4], face_id)` 由角点与"面 → 自由度"闭包构造；
- `element::crtri3::CRSeg3<D>` / `element::crtet4::CRTri4<D>`：CR 单元的 facet 迹单元（`CRSeg3` 为 2D 三角形的线段
  facet、3 个自由度；`CRTri4` 为 3D 四面体的三角 facet、4 个自由度）。自由度为父单元的全部 facet 自由度
  （`id() == 父单元 id()`），形函数为 `θₐ|facet`（2D：`shape[facet]=1`、`shape[lo]=ξ`、
  `shape[hi]=−ξ`），几何取 facet 顶点构成的 `Seg2`/`Tri3`（非等参，`jacobe = L/2` / `2×面积`）。
  于是一般边界算子（`heat_flux`、`convection`、`traction`、`temperature_penalty`、
  `displacement_penalty`、`*_lagrangian`）可直接用于 CR 的 Neumann/Robin/Dirichlet 处理；
  `CRSeg3::from_crtri3(&parent, facet)` / `CRTri4::from_crtet4(&parent, facet)` 按"对顶点"编号自动取迹；
- `GmshOptions::edges` 同时记录每条三角形边的父单元与局部 facet 编号；
  `GmshMesh::{crtri3_boundary, crtri3_boundary_by_tag}` 把边界线段组读成 `CRSeg3<D>`
  （自动确定父三角形与局部 facet，非网格边的线段返回 `MeshError::MissingEdge`）；
- `operation::heat::temperature_lagrangian(&temperature, &multiplier, T_D)`：用拉格朗日
  乘子精确施加温度本质边界条件，返回耦合块 $G_{aj} = \int_{\Gamma_D} \Lambda_a N_j\,d\Gamma$
  **全部** `M × N` 个三元组（`G` 位于鞍点系统 `[K Gᵀ; G 0]` 的非对角块，不能按下三角存储）
  与约束载荷 $g_a = \int_{\Gamma_D} \Lambda_a T_D\,d\Gamma$ 的 `Doublet` 列表；乘子即
  Dirichlet 边界热流 `λ = q_n`。温度与乘子可用**不同离散**（两组边界单元配对迭代：积分单元、
  积分点、雅可比取自温度单元，乘子单元只提供形函数与自由度 id），传同一个集合即得边界质量矩阵；
  鞍点系统对称不定，需 `sp_lu`；
- `operation::elasticity::displacement_lagrangian(&u, &lambda, ū)`：用（向量值）拉格朗日乘子
  精确施加位移本质边界条件，返回分量对角的耦合块
  $G_{(a,p),(b,q)} = \int_{\Gamma_u} \Lambda_a N_b\, \delta_{pq}\, d\Gamma$ **全部**元素与
  约束载荷 $g_{(a,p)} = \int_{\Gamma_u} \Lambda_a \bar u_p\, d\Gamma$；乘子即
  $-\boldsymbol\sigma\cdot\mathbf n$（负的边界反力密度）。同样支持位移与乘子的不同离散，
  需 `sp_lu` 解鞍点系统；
- `io::gmsh::EdgeTable`：Crouzeix–Raviart 的"边 → 自由度"编号表（`new(offset)`、`id(a, b)`
  按需分配且与端点顺序无关、`get`、`len`、`offset`）；
- `GmshOptions::faces`（默认 `false`）：读取时给所有四面体的面编号（`FaceTable`，`0..n_faces`），
  通过 `GmshMesh::face_numbering` 取回，供三维 CR 单元使用；
- `io::gmsh::FaceTable`：三维 CR 的"面 → 自由度"编号表（`new(offset)`、`id(a, b, c)` 与顶点顺序无关、
  `get`、`len`、`offset`）；
- `GmshMesh::{crtet4_elements, crtet4_elements_by_tag, crtet4_boundary, crtet4_boundary_by_tag}`：
  把 GMSH 四面体组读成 `CRTet4<D>`、把边界三角形组读成 `CRTri4<D>` 迹单元（未开 `faces` 时返回
  `MeshError::FaceNumberingMissing`，三角形不是任何四面体的面时返回 `MeshError::MissingFace`）；
- `GmshOptions::edges`（默认 `false`）：读取时给所有三角形建边编号（`0..n_edges`，
  与节点编号独立、允许重叠 —— 二者离散的是不同未知量），通过 `GmshMesh::edge_numbering`
  取回；与 `GmshOptions::renumber` 相互独立（建编号不需要重编号，但物化单元仍需稠密节点号段）；
- `GmshMesh::{crtri3_elements, crtri3_elements_by_tag}`：用读取时的边编号直接把 GMSH 三角形组读成
  `CRTri3<D>`（内部即 `Tri3` + `CRTri3::try_from_tri3`），无额外参数；相邻单元在公共边上共享同一
  自由度；未开 `edges` 时返回新增的 `MeshError::EdgeNumberingMissing`；
- `GmshMesh::{edge_table, crtri3_elements_with_table, crtri3_elements_by_tag_with_table}`：调用方
  自带 `EdgeTable` 的版本（`EdgeTable::new(offset)` 可把 CR 自由度整体平移到另一个场之后，
  适合混合系统；边缺失时返回新增的 `MeshError::MissingEdge`）；
- `operation::VariableError::ElementCountMismatch`：配对的两组单元数量不一致时由
  `temperature_lagrangian` / `displacement_lagrangian` 返回，避免静默截断；
- `io::gmsh::read_gmsh_with(path, GmshOptions)` 与 `io::gmsh::GmshOptions { renumber }`
  （默认 `false`，即与 `read_gmsh` 行为一致）：开启 `renumber` 时节点 id 按原始 id 升序压紧为
  `0..n_nodes` 并重映射单元连接关系，使 `node_table()` 必然无空洞、可在 `n_nodes..` 之后安全
  分配拉格朗日乘子自由度；连接关系引用缺失节点标签时在读取阶段返回 `MeshError::MissingNode`。

## [0.3.0] - 2026-10-24

### Breaking Changes

- `io::gmsh::read_gmsh` 返回 `Result<GmshMesh<D>, MeshError>`，不再在解析失败时
  `panic`/`expect`；
- `GmshMesh::{element_type, elements, elements_by_tag}` 与
  `FromGmsh::from_gmsh` 返回 `Result`；`FromGmsh` 新增 `GMSH_NODE_COUNT`，
  `from_gmsh` 的节点参数由 `&HashMap` 改为 `&NodeTable`，
  `elements` 会校验连通性长度，不再静默丢弃缺失节点的单元；
- `sparse_vector::try_new_from_doublet` 返回
  `Result<Col<f64>, DoubletError>`，越界索引不再触发 `panic`；
- `Variable` 去掉张量阶 `O`，保留 `C`（每节点自由度数），签名变为
  `Variable<D, C, P, N>`；`Variable` 不再表示张量，只表示一组单元上的多变量自由度；
- `Variable::gradient` / `divergence` / `curl` 从 trait 移除，改由具体问题模块实现；
  新增 `Variable::derivative_value`，返回 `Option<[[f64; C]; D]>`
  （`derivative_value[i][k] = ∂u_k/∂x_i`）；`Variable::value` 返回
  `Option<[f64; C]>`；
- `Variable::{value, derivative_value}` 为物理点查询；点定位通过
  `Variable::locate(x)` 线性扫描，不构建或存储空间索引；
- `bar::traction` 改名为 `bar::point_load`，并支持多维载荷 `[f64; D]`；
  `bar::displacement_penalty` 改名为 `bar::point_displacement_penalty`；
- `heat::h1_semierror` 改名为 `heat::h1_semi_error`；
- `elasticity::lame_constants` / `plane_stress_lambda` / `bulk_modulus`
  被具名结构体 `elasticity::Lame` 替代：
  `Lame::from_engineering` / `plane_stress` / `plane_strain` / `bulk_modulus`；
- `Node` 增加 `Copy`/`Debug`/`PartialEq`；单元类型增加 `Clone`/`Debug`/`PartialEq`
  和 `nodes()` / `node(index)` 访问器；`GmshMesh` 增加摘要式 `Debug`。

### 新增

- `io::gmsh::NodeTable<D>`：稠密节点坐标表（`id == index`），
  `GmshMesh::node_table()` 与 `elements_by_tag_with_table` 支持多组复用同一张表；
- `io::gmsh::PhysicalGroup<'a>` 与 `GmshMesh::physical_groups` /
  `physical_group` / `physical_group_by_name` / `physical_names` /
  `has_physical_group`；解析 `$PhysicalNames` 时建立 `name_to_tag` 反查表，
  物理组按 tag 升序迭代；
- `spatial` 模块：`Aabb<D>` 与通用 `SpatialIndex<T, D>`（AABB broad phase +
  均匀网格），供 FEM 点定位与无网格邻域搜索复用；
- 后处理：`Variable::{value, gradient, divergence, curl}`、`Variable::locate`、
  `heat::{gradient, flux}`、`elasticity::{strain, stress}`。
  不提供节点导数或节点平均输出（C0 连续下节点导数不唯一，属于调用方后处理）。

### 修复

- 新增 `MeshError` / `DoubletError` / `VariableError`，公开路径不再使用
  `unwrap`/`expect`/`panic!`，错误信息携带 path、section、line、物理组、
  单元与节点 tag 等上下文；
- GMSH 解析器改用带行号与 section 状态的 `MshReader`：修正
  `$PhysicalNames` 的字节偏移切分、节点坐标 `D.min` 静默补零、node tag 下溢，
  并校验节点/元素总数与结束标记；
- `elements()` 不再通过 `filter_map` 静默丢弃节点缺失的单元；
- 新增 `MeshError` 错误分支测试和 `try_new_from_doublet` / `Variable::value`
  越界测试。

## [0.2.0] - 2026-10-17

首个对外发布的版本。crate 由早期的 `rustfem` 更名为 `ansatz`，定位从"FEM 库"收窄为
**面向力学问题的离散近似与组装层**：负责近似（形函数）、数值积分、组装，线性求解交给
外部稀疏后端（默认 `faer`），为后续接入无网格离散留出空间。

### Breaking Changes

- crate 重命名为 `ansatz`，`use rustfem::…` 需全部改为 `use ansatz::…`（`402ef85`）；
- `Variable` 增加张量阶参数 `O`，签名由 `Variable<D, P, N>` 变为
  `Variable<D, O, P, N>`，每节点自由度数改为 `C(D + O - 1, O)`（`03e6891`）；
- `Element::derivative_shape` 改为返回**物理梯度**（张量表示，`[[f64; D]; N]`），
  不再返回参数空间导数（`3f0cec6`）。

### 新增

- 单元：`Seg2<D>`（1D/2D/3D 线单元）、`Tri3<D>`、`Quad4<D>`、`Poi1<D>`（点单元）
  （`1dabf8b`、`a4754ab`、`0e55c95`）；
- 积分规则：`GaussLegendre1D<N>`、`Dunavant<G>`、`TensorProductQuad<N>`、`Keast<G>`、
  `TensorProductHex<N>`、`Point`，并以 `FullIntegration` / `ReducedIntegration`
  区分完全积分与减缩积分（`f9c9b36`、`a4754ab`）；
- 热传导算子 `operation::heat`：`conduction_stiffness`、`capacity`、`heat_source`、
  `heat_flux`、`convection`、`temperature_penalty`、`l2_error`、`h1_semierror`（`3f0cec6`）；
- 线弹性算子 `operation::elasticity`：`stiffness`、`dilatational_stiffness`、
  `deviatoric_stiffness`、`body_force`、`traction`、`displacement_penalty`、
  `l2_error`、`energy_error`，以及 `lame_constants`、`plane_stress_lambda`、
  `bulk_modulus`（`511bf1a`、`5599db0`）；
- `Variable::value` 支持在任意物理点插值求值，`Element::local_coordinates` /
  `Element::contains` 提供逆映射与点包含判定（`03e6891`）；
- `GmshMesh::element_type` 支持在运行时按物理组选择单元类型（`53ba80d`）；
- 端到端示例：`bar1d`、`patchtest2d`、`cantilever`、`cook`、`elasticity_cli`
  （由 TOML 配置驱动，`d593699`、`702ede0`、`f200b5c`）；
- mdBook 用户手册（`book/`），含用法、理论与算例三部分（`4d3487f`、`509c5f9`、`8346bbf`）；
- `dhat-heap` / `dhat-ad-hoc` feature 与全局分配器接入，用于堆分配剖析。

### 变更

- 组装算子改为返回稀疏条目：`Triplets` / `Doublets` / `TripletsAndDoublets`
  类型别名，由调用方合并为 `faer` 稀疏矩阵与向量（`6d282c5`）；
- 线弹性刚度改为独立实现，体积项 `dilatational_stiffness` 使用减缩积分（`5599db0`）；
- 算子内部物理量统一改写为指标（求和）记号，移除显式 Voigt `B`/`D` 矩阵，
  使同一公式适用于任意空间维数 `D`（含嵌入单元）（`3f0cec6`、`5599db0`）；
- 示例统一迁移到 `examples/`，网格迁移到 `examples/msh/`（`4d3487f`）；
- 项目定位改写为离散近似层，README 补充安装方式与 11 条初步使用规则（`2f34a7a`）；
- `Cargo.toml` 补齐发布元数据（`description`、`keywords`、`categories`、`license`、
  `readme`、`repository`、`documentation`、`exclude`）并声明
  `rust-version = "1.85.0"`，`cargo package` 内容由 69 个文件收敛到 40 个（`9dda107`）。

### 修复

- 清零 `cargo clippy --all-targets` 警告（`6d282c5`）；
- 修正文档中失效的网格路径（`test/msh/` → `examples/msh/`）与过期的状态描述
  （示例清单、测试数量、`src/main.rs` 等）（`1c5d7ec`）。

### 文档

- 新增 `PLAN.md` 第 7 节「首次正式发布前的安排」：发布门槛、时间线、元数据检查、
  版本策略、发布 Runbook 与发布后收尾（`1c5d7ec`）；
- `AGENTS.md` 全面刷新（模块树、`Variable<D, O, P, N>`、算子清单、已知缺口）（`1c5d7ec`）；
- 补充线弹性理论、悬臂梁与 Cook 膜算例文档（`509c5f9`、`8346bbf`）。

### 测试

- 58 个内联单元测试：积分规则矩验证、形函数/雅可比、GMSH 解析、
  热传导与线弹性的解析解回归、`Tri3`/`Quad4` patch test；
- 已用 Rust 1.85.0 实测 `cargo check --all-targets` 与 `cargo test` 通过，锚定 MSRV。

## [0.1.0] - 未发布

早期内部原型（crate 名 `rustfem`），未发布到 crates.io，保留此条目仅用于说明版本沿革。

[Unreleased]: https://github.com/jmwjc/ansatz/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/jmwjc/ansatz/releases/tag/v0.3.0
[0.2.0]: https://github.com/jmwjc/ansatz/releases/tag/v0.2.0
