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

- `operation::Variable<const D, const P, const N>`：`type Item`、`iter()`；对 `&[T]` 与 `Vec<T>` 有 blanket impl。
- `operation::bar`（1D 桁架组装算子）：
  - `stiffness(v, ea) -> Vec<Triplet<usize, usize, f64>>`
  - `traction(v, p) -> Vec<Doublet<usize, f64>>`
  - `body_force(v, f) -> Vec<Doublet<usize, f64>>`
  - `displacement_penalty(v, penalty, u) -> (Vec<Triplet>, Vec<Doublet>)`

## `io::gmsh`

- `read_gmsh<const D: usize>(path) -> GmshMesh<D>`
- `GmshMesh<D>`：字段 `nodes`、`physical_names`；方法 `elements<T>(name)`、`elements_by_tag<T>(tag)`。
- `FromGmsh<const D: usize>`：`GMSH_ELEMENT_TYPE` 常量与 `from_gmsh(...)`。

## `sparse_vector`

- `Doublet<I, T>`：`(index, value)` 对，`Doublet::new(index, value)`。
- `try_new_from_doublet(n, &[Doublet]) -> faer::Col<f64>`：把载荷条目累加为稠密向量。

## 常用外部依赖（faer）

组装产物经 faer 合并求解：

- `faer::sparse::SparseColMat::try_new_from_triplets(nrows, ncols, &triplets)`;
- `k.sp_cholesky(faer::Side::Lower)` → `llt.solve(&f)`。

## 最小示例

完整可运行示例见 `examples/bar1d.rs`（`cargo run --example bar1d`）；「组装与求解」一章也给出了完整链路片段。
