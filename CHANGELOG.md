# Changelog

本文件记录 `ansatz` 的所有重要变更。

格式遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)，
版本号遵循 [语义化版本](https://semver.org/lang/zh-CN/)。
0.x 阶段的 MINOR 版本允许包含破坏性变更，破坏性变更会在对应版本的
「Breaking Changes」小节中单列并给出迁移方式。

## [Unreleased]

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
