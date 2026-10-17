# Changelog

本文件记录 `ansatz` 的所有重要变更。

格式遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)，
版本号遵循 [语义化版本](https://semver.org/lang/zh-CN/)。
0.x 阶段的 MINOR 版本允许包含破坏性变更，破坏性变更会在对应版本的
「Breaking Changes」小节中单列并给出迁移方式。

## [Unreleased]

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

[Unreleased]: https://github.com/jmwjc/ansatz/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/jmwjc/ansatz/releases/tag/v0.2.0
