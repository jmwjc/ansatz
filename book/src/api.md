# API 参考

`ansatz` 的公共 API 按模块组织。以下为速查清单；各条目的完整说明见库内 rustdoc。

## `node`

- `Node<const D: usize>`：节点，字段 `id: usize`、`coordinates: [f64; D]`。
  - `Node::new(id, coordinates) -> Self`。

## `element`

核心 trait 与积分规则（`src/element/mod.rs`）：

- `Element<const D: usize, const P: usize, const N: usize>`：
  - `id() -> [usize; N]`
  - `shape([f64; P]) -> [f64; N]`
  - `derivative_shape([f64; P]) -> [f64; N]`
  - `vertices_coordinates() -> [[f64; D]; N]`
  - `coordinates([f64; P]) -> [f64; D]`
  - `jacobe([f64; P]) -> f64`
  - `jacobe_mat([f64; P]) -> [[f64; P]; D]`
- `IntegrationScheme<const P: usize, const G: usize>`：`points() -> [[f64; P]; G]`、`weights() -> [f64; G]`。
- `FullIntegration<P, G>` / `ReducedIntegration<P, G>`：标记子 trait。
- 积分规则标记类型：`GaussLegendre1D<N>`、`Dunavant<G>`、`TensorProductQuad<N>`、`Keast<G>`、`TensorProductHex<N>`、`Point`。

子模块：

- `element::seg2::Seg2<D>`：`new(n1, n2)`，2 节点线单元；
- `element::tri3::Tri3<D>`：`new(n1, n2, n3)`，3 节点三角形；
- `element::quad4::Quad4<D>`：`new(n1, n2, n3, n4)`，4 节点四边形；
- `element::poi1::Poi1<D>`：`new(node)`，1 节点点单元。

## `operation`

- `operation::Variable<const D: usize, const C: usize, const P: usize, const N: usize>`：`C` 为每节点自由度数/分量数；提供 `iter`、`num_dofs`、`locate` / `value` / `derivative_value`；对 `&[T]` 与 `Vec<T>` 有 blanket impl。
- `Variable::value(x, coefficients) -> Result<Option<[f64; C]>, VariableError>`：在物理点上求值，点外返回 `Ok(None)`；内部用 `locate(x)` 线性扫描单元。
- `Variable::derivative_value(x, coefficients) -> Result<Option<[[f64; C]; D]>, VariableError>`：返回 `derivative_value[i][k] = ∂u_k / ∂x_i`。
- `gradient` / `divergence` / `curl` 由具体问题模块实现，不在 `Variable` 上。
- `spatial`：`Aabb<D>`、`SpatialIndex<T, D>`（AABB broad phase + 均匀网格，可复用于无网格邻域搜索）。
- `heat::gradient` / `heat::flux`；`elasticity::strain` / `elasticity::stress`。
- `elasticity::Lame::{from_engineering, three_dimensional, plane_stress, plane_strain, bulk_modulus}`：具名 Lamé 常量。
- `operation::bar`（1D 桁架组装算子）：
  - `stiffness(v, ea) -> Vec<Triplet<usize, usize, f64>>`
  - `point_load(v, [f64; D]) -> Vec<Doublet<usize, f64>>`
  - `body_force(v, f) -> Vec<Doublet<usize, f64>>`
  - `point_displacement_penalty(v, penalty, u) -> (Vec<Triplet>, Vec<Doublet>)`

## `io::gmsh`

- `read_gmsh<const D: usize>(path) -> Result<GmshMesh<D>, MeshError>`
- `GmshMesh<D>`：字段 `nodes`、`physical_names`；方法 `element_type(name)`、`elements<T>(name)`、`elements_by_tag<T>(tag)`、`node_table()`、`elements_by_tag_with_table<T>(tag, &NodeTable)`、`physical_groups()`、`physical_group(tag)`、`physical_group_by_name(name)`、`physical_names()`、`has_physical_group(name)`。
- `NodeTable<D>`：稠密节点坐标表（`id == index`），`from_nodes`、`from_map`、`get`、`coords`。
- `PhysicalGroup<'a>`：`tag`、`name`、`element_type`、`element_count` 的只读视图。
- `FromGmsh<const D: usize>`：`GMSH_ELEMENT_TYPE`、`GMSH_NODE_COUNT` 常量与 `from_gmsh(ids, &NodeTable) -> Result<Self, MeshError>`。

## `sparse_vector`

- `Doublet<I, T>`：`(index, value)` 对，`Doublet::new(index, value)`。
- `try_new_from_doublet(n, &[Doublet]) -> Result<faer::Col<f64>, DoubletError>`：把载荷条目累加为稠密向量。

## 常用外部依赖（faer）

组装产物经 faer 合并求解：

- `faer::sparse::SparseColMat::try_new_from_triplets(nrows, ncols, &triplets)`;
- `k.sp_cholesky(faer::Side::Lower)` → `llt.solve(&f)`。

## 最小示例

完整可运行示例见 `examples/bar1d.rs`（`cargo run --example bar1d`）；「组装与求解」一章也给出了完整链路片段。
