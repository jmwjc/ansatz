# ansatz 开发计划（PLAN.md）

> 说明：README 目前仅包含项目标题，缺少正式文档，本计划是主要项目文档。本文基于对现有代码的考察（2026-10 更新），反映**当前实际架构**：`src/element`（单元与积分规则）、`src/operation`（组装算子）、`src/io`（GMSH 读取）三大模块，1D 桁架端到端求解链路已可用。历史计划中提到的 `src/problem`、`src/approximation.rs`、`src/operations.rs`、`src/solver.rs`、`GaussSeg2` 等均已不存在或被取代，以本文为准。

## 1. 项目目标与定位

`ansatz` 是一个用 Rust 编写的**离散近似层**：它把一个连续介质力学问题变成一个稀疏代数系统。

*ansatz*（试探解的形式）在有限元中是在网格上的分片多项式，在无网格方法中是在点云上的移动最小二乘 / 再生核近似。本库只对这一层负责——**近似（形函数）、数值积分、组装**——而不拥有线性求解器：解方程交给外部稀疏后端（默认为 `faer`）。因此离散方式（FEM，未来还有无网格）与线性代数相互独立。目标是为结构力学/连续介质力学问题提供：

- 模块化的近似抽象（`Element<D, P, N>`，const 泛型，零成本抽象）与积分规则（`IntegrationScheme<P, G>`）；
- 可扩展的组装算子（`Variable<D, O, P, N>` 输入抽象，返回 `Triplet`/`Doublet` 由调用方合并为稀疏矩阵/向量）；
- 线性求解交给外部稀疏后端（默认 `faer`：`SparseColMat`、稀疏 Cholesky/LU），本库不绑定求解器；
- 1D/2D/3D 桁架、弹性体等常见问题的求解示例；
- 清晰的 API、完整文档与回归测试。

## 2. 当前状态（已完成）

截至 2026-10，以下工作已完成并有测试覆盖（58 个单元测试全部通过）：

- **库结构**：`src/lib.rs` 库入口（无 `src/main.rs`）；模块划分为 `element/`、`operation/`、`io/`，示例均在 `examples/`。
- **单元与积分规则**（`src/element/`）：
  - `Seg2<D>`（2 节点线单元，支持 1D/2D/3D）、`Tri3<D>`、`Quad4<D>`、`Poi1<D>`（点单元）；
  - 积分规则标记类型：`GaussLegendre1D<N>`、`Dunavant<G>`（2D 三角）、`TensorProductQuad<N>`、`Keast<G>`（3D 四面体）、`TensorProductHex<N>`、`Point`；`FullIntegration`/`ReducedIntegration` 标记子 trait 区分完全/减缩积分；
  - 形函数、导数、雅可比（`jacobe`/`jacobe_mat`，该拼写为有意保留的既有命名）均有单元测试。
- **组装算子**（`src/operation/`）：`bar`（1D 桁架）、`heat`（热传导：刚度/热容/热源/热流/对流/温度罚函数/误差范数）、`elasticity`（线弹性：刚度/体积-偏量分解/体力/面力/位移罚函数/L²/能量误差）——均以指标求和实现，不拼 Voigt 矩阵。
- **稀疏工具**（`src/sparse_vector.rs`）：`Doublet` 与 `try_new_from_doublet` → `faer::Col<f64>`。
- **GMSH v4.1 读取**（`src/io/gmsh.rs`）：`read_gmsh`、`GmshMesh`、`FromGmsh` trait（已支持 Seg2=1、Tri3=2、Poi1=15），用 `examples/msh/` 真实网格做解析测试。
- **端到端示例**（`examples/`）：`bar1d`（10 个 `Seg2<1>` 杆单元，左端固定、均布体力 f=1、右端拉力 P=1，罚函数加约束后 `sp_cholesky` 求解）、`patchtest2d`（Tri3/Quad4 patch test）、`cantilever`（平面应力解析解对比）、`cook`（平面应变收敛）、`elasticity_cli`（TOML 驱动的 CLI）。
- **性能剖析设施**：`dhat-heap`/`dhat-ad-hoc` feature 与 `#[global_allocator]` 已接入。

## 3. 里程碑（Milestones）

| 里程碑 | 目标 | 主要交付物 | 验收标准 |
|--------|------|------------|----------|
| **M1：代码质量与 CI** | clippy 清零，建立持续集成 | clippy 警告修复、`.github/workflows/ci.yml`（fmt + clippy + test）、`main.rs` 精简为纯示例 | `cargo clippy --all-targets` 零警告；GitHub Actions 全绿 |
| **M2：问题抽象与求解器封装** | 把 `main.rs` 中的求解链路提炼为库 API | `Problem` trait、`Truss1D` 问题类型、`Solver` trait（封装 faer Cholesky/LU）、`examples/` 目录、`Element` 类型 GMSH 注册便捷接口 | 1D 杆端受拉示例通过库 API 调用，与当前 `main.rs` 结果一致（回归测试） |
| **M3：2D 求解链路** | 打通三角形/四边形网格上的首个 2D 物理问题 | `Tri3`/`Quad4` 的 B 矩阵与单元刚度组装算子、2D 泊松或平面应力问题、patch test 验证 | `examples/msh/` 下 2D 网格 patch test 数值结果与理论值在 1e-8 内一致 |
| **M4：材料与物理问题扩展** | 从单物理模型扩展到材料库 | `Material` trait（`Truss`、`LinearElastic`）、平面应力/平面应变、3D 线弹性骨架、`Quad4` 减缩积分验证 | 至少新增 2 个可运行示例（平面应力悬臂梁等） |
| **M5：性能与并行** | 中大网格下的性能优化 | `dhat` 剖析与热点消除、`rayon` 接入单元组装并行化、基准测试 | 中等规模网格（>1e4 DOF）release 下组装+求解 5 秒内完成 |
| **M6：文档与发布** | 完成用户文档并准备 crates.io 发布 | 完整 README、公共 API rustdoc 全覆盖、`CHANGELOG.md`、`docs/` 设计文档、`cargo publish --dry-run` 通过（执行清单见第 7 节） | 示例可一键运行；docs.rs 可编译 |

## 4. 任务分解与时间线（建议）

每个里程碑预计 1–2 周，视投入时间而定。M1–M3 原有目标中单元框架部分已提前完成（见第 2 节），以下按剩余工作排序：

### M1：代码质量与 CI（第 1 周）

- [ ] 清零 `cargo clippy --all-targets` 警告（当前约 91 个；修复时遵循"最小化修改"，不引入行为变化）；
- [ ] 配置 `.github/workflows/ci.yml`（rustfmt --check + clippy --all-targets -- -D warnings + cargo test）；
- [x] 精简 `main.rs`：求解逻辑全部下沉到库，`main.rs` 只保留示例调用；
- [x] 为端到端 1D 桁架示例补回归测试（锚定当前数值结果）。

### M2：问题抽象与求解器封装（第 2 周）

- [ ] 定义 `Problem` trait（几何、材料、边界、载荷、组装、求解、后处理）；
- [ ] 定义 `Solver` trait（`factorize`/`solve`），用 faer 稀疏 Cholesky/LU 实现默认后端；
- [ ] 实现 `Truss1D` 问题类型，封装 `bar` 算子的调用与罚函数边界条件；
- [ ] 抽象 `BoundaryCondition`（本质/自然边界）与 `Load`（体力/面力/集中力）；
- [ ] 新建 `examples/truss1d.rs`，迁移 `main.rs` 示例；
- [ ] 完善 `FromGmsh` 的单元注册机制，使新单元接入 GMSH 读取只需在 io 层实现 trait。

### M3：2D 求解链路（第 3–4 周）

- [ ] 为 `Tri3<D>`/`Quad4<D>` 实现 B 矩阵（应变-位移矩阵）与单元刚度组装算子；
- [ ] 实现首个 2D 物理问题（建议先做标量泊松问题，再做平面应力）；
- [ ] 用 `examples/msh/patchtest_tri3.msh` 做 patch test：常应变场下数值解应与理论解一致；
- [ ] 验证 `Quad4` 的完全积分与减缩积分（`FullIntegration`/`ReducedIntegration`）行为差异；
- [ ] 为 2D 问题编写单元测试与回归测试。

### M4：材料与物理问题扩展（第 5–6 周）

- [ ] 引入 `Material` trait（`Truss`、`LinearElastic`（平面应力/平面应变/3D））；
- [ ] 实现平面应力/平面应变问题类型；
- [ ] 实现 3D 线弹性问题骨架（复用 `Keast`/`TensorProductHex` 积分规则）；
- [ ] 添加 `examples/poisson2d.rs`、`examples/elasticity2d.rs`（平面应力悬臂梁与解析解对比）；
- [ ] 评估梁单元（Euler-Bernoulli/Timoshenko）需求，决定是否纳入 M4 或单列里程碑。

### M5：性能与并行（第 7–8 周）

- [ ] 用 `dhat-heap` 剖析组装热点（Triplet 分配、单元循环），消除冗余分配；
- [ ] 用 `rayon` 并行化单元级组装（element loop → par_iter 收集 Triplet），验证 `parallel` feature；
- [ ] 生成中等规模网格（>1e4 DOF），建立基准测试（可用 `cargo bench` 或简单计时脚本）；
- [ ] 验证 release 模式（`debug = 1` 保留符号便于剖析）下组装+求解性能目标。

### M6：文档与发布（第 9–10 周）

- [ ] 重写 `README.md`（简介、安装、快速开始、API 概览、示例链接）；
- [ ] 补齐所有公共 API 的 rustdoc（中文）；
- [ ] 编写 `docs/element.md`、`docs/assembly.md`、`docs/solver.md` 设计文档；
- [x] 维护 `CHANGELOG.md`（基于 Conventional Commits 历史，已建立 0.2.0 条目）；
- [x] 升级到 `0.2.0` 并执行 `cargo publish --dry-run` 验证（已于 0.2.0 打 tag 时完成）。

## 5. Git Commit 规范化

本项目采用 [Conventional Commits](https://www.conventionalcommits.org/) 作为提交信息规范，便于自动生成 CHANGELOG 和版本号决策。

### 5.1 提交信息格式

```text
<type>(<scope>): <subject>

<body>

<footer>
```

- **type**（必填）：提交类别，见下表；
- **scope**（可选）：模块或功能范围，如 `element`、`operation`、`io`、`solver`、`truss1d`、`ci`；
- **subject**（必填）：简短描述，使用祈使句，首字母小写，结尾不加句号；
- **body**（可选）：详细说明动机、实现细节、设计取舍；
- **footer**（可选）：关联 issue、破坏性变更说明等。

### 5.2 类型说明

| 类型 | 含义 | 示例 |
|------|------|------|
| `feat` | 新功能 | `feat(element): add Tri3 shape functions` |
| `fix` | 修复 bug | `fix(operation): correct boundary condition indexing` |
| `docs` | 文档变更（README、rustdoc、设计文档） | `docs: add quickstart guide` |
| `style` | 代码格式调整（不影响逻辑） | `style: apply rustfmt` |
| `refactor` | 重构（不新增功能也不修复 bug） | `refactor(io): unify gmsh parsing` |
| `perf` | 性能优化 | `perf(sparse): reduce triplet allocations` |
| `test` | 测试相关 | `test(truss1d): add regression test` |
| `chore` | 构建/工具/杂项 | `chore(ci): add GitHub Actions workflow` |
| `build` | 构建系统、依赖变更 | `build: bump faer to 0.25` |
| `ci` | CI/CD 配置 | `ci: run clippy on PR` |
| `revert` | 回滚提交 | `revert: feat(element): add Tri3 shape functions` |

### 5.3 破坏性变更

若提交引入不兼容变更，在 `type` 后加 `!`，并在 `footer` 中说明迁移方式：

```text
feat(element)!: replace Element trait with Approximation

BREAKING CHANGE: The `Element<D, P, N>` trait is removed. Use `Approximation<DIM, NDOF>` instead.
```

### 5.4 示例提交

```text
feat(solver): wrap faer sparse Cholesky in Solver trait

- Add `Solver` trait with `solve` and `factorize` methods.
- Implement `FaerCholeskySolver` for symmetric positive definite matrices.
- Use `Side::Lower` as default fill-reducing ordering.

Refs: #12
```

### 5.5 推荐分支与合并策略

- 功能开发：`feat/<scope>-<description>`（如 `feat/element-tri3`）；
- Bug 修复：`fix/<scope>-<description>`；
- 文档：`docs/<description>`；
- 每个 PR/MR 尽量只对应一个 `type`，标题按规范书写；
- 合并前要求 CI 通过、至少一名 reviewer 通过。

## 6. 设计原则

- **最小化修改**：每次提交只解决一个明确问题；
- **trait 优先**：用 Rust trait 抽象单元、材料、问题、求解器，方便扩展；保持 const 泛型 + `impl Iterator`（RPITIT）的零成本抽象风格，避免运行时动态分发；
- **测试驱动**：新增物理模型或单元时，先写单元测试与回归测试（patch test 优先）；
- **文档同步**：每次 API 变更同步更新 rustdoc 与 README（中文）；
- **性能可测**：关键路径使用 `dhat` 或自定义 benchmark 验证；
- **命名稳定**：既有命名（如 `jacobe`、`Doublet`）是有意保留的 API 一部分，不随个人喜好"修正"。

## 7. 首次正式发布（v0.2.0）前的安排

> 本节是 M6「文档与发布」的细化与执行清单。0.2.0 的目标是：**把当前已验证的 API 面（`element` / `operation` / `io` + `bar` / `heat` / `elasticity`）冻结成一个可用、可安装、可引用的版本**，而不在发布前引入新特性。

### 7.1 发布门槛（全部满足才允许 `cargo publish`）

- [ ] `cargo fmt --check` 通过；
- [ ] `cargo clippy --all-targets -- -D warnings` 零警告（M1 的清零结果）；
- [ ] `cargo test` 全绿（当前 58 个测试）；
- [ ] `cargo doc --no-deps` 无警告，公共 API rustdoc 全覆盖；
- [ ] CI（GitHub Actions：fmt + clippy + test + doc）在 `main` 上全绿，且开启分支保护；
- [ ] 5 个示例（`bar1d`/`patchtest2d`/`cantilever`/`cook`/`elasticity_cli`）全部可运行，README 中的命令逐条实跑验证；
- [x] `CHANGELOG.md` 覆盖本次发布的全部变更；
- [ ] 未稳定的 API 已标注（`#[doc(hidden)]` 或移入 `unstable` 模块），公开面经过一次逐项审查；
- [ ] `cargo package --list` 与 `cargo publish --dry-run` 通过，且打包内容仅含预期文件。

### 7.2 时间线（T = 发布日）

| 时间 | 工作重点 |
|------|----------|
| **T−4 周** | 功能冻结：只接受 bugfix 与文档；开始整理 `CHANGELOG.md`；确定 0.2.0 的公开 API 清单 |
| **T−3 周** | API 审查：命名一致性（`jacobe` 等既定拼写保持不变）、泛型参数顺序（`<T, D, P, N, G>`）、错误处理（`unwrap` 是否应在库内出现）、用 `cargo public-api` 导出基线快照 |
| **T−2 周** | 文档周：README/`book/` 补齐并校对所有代码片段；`rustdoc` 示例与文档测试能编译；检查 `book/` 中所有路径真实存在 |
| **T−1 周** | 发布演练：在干净目录做 `cargo package` + 试装（`cargo new` 后以 path 依赖引入）；`--dry-run` 修元数据；确认 MSRV |
| **T 日** | 打 tag → `cargo publish` → 验证 docs.rs 构建 → 发 announcement |
| **T+1 周** | 收集反馈，补 issue 模板/讨论区；如需修缺陷则发 0.2.1（不撤回已发布版本） |

### 7.3 元数据检查（`Cargo.toml`）

- [x] `name = "ansatz"`（已在 crates.io 确认未被占用）、`version = "0.2.0"`、`edition = "2024"`；
- [x] `description`：一句英文、不含“in Rust”、点明“离散近似/组装层”；
- [x] `keywords`：最多 5 个且为小写（当前 `fem`/`meshfree`/`discretization`/`tensor`/`assembly`）；
- [x] `categories`：`science`、`mathematics`、`simulation`（需与 crates.io 合法分类一致）；
- [x] `license = "MIT"` 且仓库内存在 `LICENSE`；
- [x] `readme = "README.md"`、`repository`、`documentation`（docs.rs 链接）、`homepage`（可选）；
- [x] `rust-version = "1.85.0"`（MSRV）：edition 2024 的下限；已用 `cargo +1.85.0 check --all-targets` 与 `cargo +1.85.0 test` 验证，仍需在 CI 中加一个固定 1.85 job 守住；
- [x] `exclude`：已排除 `book/`、`dhat-heap.json`、`PLAN.md`、`AGENTS.md`（`cargo package --list` 现为 40 个文件）；`target/` 由 cargo 自动排除；
- [ ] `[features]` 中 `sparse`/`std` 若仍为占位开关，在 rustdoc 中注明，或直接移除；
- [ ] `[profile.release] debug = 1` 保留（用于 `dhat` 剖析），并在 README 中说明。

### 7.4 版本与兼容性策略

- **0.x 语义**：`0.MINOR.PATCH`。MINOR 允许破坏性变更，PATCH 仅修复；
- **破坏性变更**必须：提交信息用 `type!:` + footer `BREAKING CHANGE:`，`CHANGELOG.md` 单列「Breaking」小节并给出迁移示例；
- **依赖升级**（尤其 `faer`）：单独 `build:` 提交，升级后留一周观察期，并重跑全部 patch test 与算例数值对比；
- **1.0.0 的判定条件**：公共 API 冻结 ≥1 个 MINOR 周期且无破坏性变更；无网格离散（`ansatz` 抽象的真正检验）落地并稳定；MSRV 承诺明确；已有至少一个外部使用者在真实项目中可用。

### 7.5 发布操作步骤（Runbook）

```bash
# 0) 工作区必须干净
git status --short

# 1) 全量校验
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test
cargo doc --no-deps

# 2) 打包内容审查（确认无 target/、dhat-heap.json）
cargo package --list

# 3) 演练发布
cargo publish --dry-run

# 4) 打 tag 并发布
git tag -a v0.2.0 -m "ansatz v0.2.0"
git push origin v0.2.0
cargo publish

# 5) 发布后验证：docs.rs 构建成功；在新目录试装
cargo new /tmp/ansatz-smoke && cd /tmp/ansatz-smoke && cargo add ansatz
```

- **撤回**：仅在严重缺陷时 `cargo yank --version 0.2.0`；yank 不删除代码，必须同时发布修复版本；
- **发布后 24 小时内**：检查 docs.rs 构建结果、README 渲染、`cargo add` 是否可拉取、示例能否编译。

### 7.6 发布后收尾

- [ ] GitHub Release 基于 `v0.2.0` tag，正文摘录 `CHANGELOG.md`；
- [ ] 把 README 中的「尚未发布到 crates.io」提示升级为「安装」；
- [ ] 同步更新 `book/src/quickstart.md` 的安装方式（git 依赖 → `cargo add ansatz`）；
- [ ] 添加 issue 模板（bug / feature）与 `CONTRIBUTING.md`（可选）；
- [ ] 在 `docs.rs` 配置中启用全部可选 feature 的文档构建（如 `all-features = true`）。

## 8. 风险与依赖

| 风险 | 应对 |
|------|------|
| `faer` 版本升级可能带来 API 不兼容 | 锁定版本并在升级时单独提 `build` 类型提交 |
| `Quad4` 减缩积分可能引入沙漏模式 | M3 中通过 patch test 与 `FullIntegration`/`ReducedIntegration` 对比验证 |
| 2D/3D 问题组装逻辑与 1D 差异大，算子接口可能需演进 | `Variable` trait 保持最小接口，必要时以 `!` 破坏性变更单独提交 |
| 缺少文档导致协作困难 | 每个里程碑结束时更新 README 与设计文档 |
| 性能瓶颈未及时发现 | M5 引入基准测试与 `dhat` 剖析 |

## 9. 结论

`ansatz` 已完成 1D 桁架端到端求解链路与单元/积分规则框架，处于"框架就绪、待抽象与扩展"阶段。本计划把剩余工作组织为 M1（代码质量与 CI）→ M2（问题与求解器抽象）→ M3（2D 求解链路）→ M4（材料与物理扩展）→ M5（性能与并行）→ M6（文档与发布）六个里程碑，每个里程碑有明确的交付物和验收标准，并严格遵守 Conventional Commits 规范。
