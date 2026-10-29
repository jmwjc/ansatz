# 核心抽象

`ansatz` 的公共 API 由少数几个 trait 组织而成。理解这些 trait 是使用与扩展库的关键。

## 单元 trait：`ElementBase` / `Element` / `HdivElement`

单元抽象按「几何」与「形函数」拆成三个 trait，三个 const 泛型参数都是编译期常量：

- `D`：**空间维数**（物理空间，1/2/3）；
- `P`：**参数空间维数**（0 = 点、1 = 线、2 = 面）；
- `N`：**自由度/节点数**。

### `ElementBase<D, P, N>`：所有单元共享的几何与自由度

```rust
pub trait ElementBase<const D: usize, const P: usize, const N: usize> {
    fn id(&self) -> [usize; N];
    fn vertices_coordinates(&self) -> [[f64; D]; N];
    fn coordinates(&self, parametric_coordinates: [f64; P]) -> [f64; D];
    fn jacobe(&self, parametric_coordinates: [f64; P]) -> f64;
    fn jacobe_mat(&self, parametric_coordinates: [f64; P]) -> [[f64; D]; P];
    fn local_coordinates(&self, x: [f64; D]) -> Option<[f64; P]>;
    fn contains(&self, x: [f64; D]) -> bool { /* local_coordinates(x).is_some() */ }
}
```

它不涉及"每个自由度是标量还是向量"，因此标量单元与向量值单元共享同一套定位、雅可比、度量与边界盒逻辑。

### `Element<D, P, N>`：标量单元

```rust
pub trait Element<const D: usize, const P: usize, const N: usize>: ElementBase<D, P, N> {
    fn shape(&self, parametric_coordinates: [f64; P]) -> [f64; N];
    fn derivative_shape(&self, parametric_coordinates: [f64; P]) -> [[f64; D]; N];
}
```

每个自由度对应一个**标量**形函数，这是 `Seg2`、`Tri3`、`Quad4`、`Tet4`、`Poi1` 与 CR 单元族实现的 trait。`ScalarElement` 是它的可读别名，用于强调与 `HdivElement` 的区别。

### `HdivElement<D, P, N>`：$H(\mathrm{div})$ 向量值单元

```rust
pub trait HdivElement<const D: usize, const P: usize, const N: usize>: ElementBase<D, P, N> {
    fn shape(&self, parametric_coordinates: [f64; P]) -> [[f64; D]; N];
    fn divergence_shape(&self, parametric_coordinates: [f64; P]) -> [f64; N];
}
```

每个自由度对应一个 **$D$ 维向量**形函数（例如 RT0 的边通量基），因此不能由 `Element::shape` 表达。最低阶 Raviart–Thomas 单元 `RTTri3` 即实现 `ElementBase + HdivElement`。

> `HdivElement::shape` 与 `Element::shape` **同名**，但返回 `[[f64; D]; N]`。本库的单元不会同时实现两个 trait，因此普通方法调用不会歧义；若某个类型两者都实现，需用 UFCS（`Element::shape` / `HdivElement::shape`）消歧。

### 方法语义

- `id()`：单元的 `N` 个全局自由度编号，用于组装时定位稀疏矩阵的行列；
- `vertices_coordinates()`：`N` 个节点的坐标（CR 是边中点/面重心，迹单元是父单元顶点）；
- `Element::shape(ξ)`：参数坐标 $ξ$ 处的 `N` 个标量形函数值；
- `HdivElement::shape(ξ)`：参数坐标 $ξ$ 处的 `N` 个向量形函数（每个一个 $D$ 维向量）；
- `derivative_shape(ξ)`：标量形函数的**物理**梯度 `[[f64; D]; N]`（`[node][spatial]`）；
- `divergence_shape(ξ)`：向量形函数的散度（仅 `HdivElement`）；
- `coordinates(ξ)`：等参映射，把参数坐标映射回物理坐标；
- `jacobe(ξ)`：面积/体积比（`|det J|` 的推广），数值积分时作为被积函数的度量因子；
- `jacobe_mat(ξ)`：雅可比矩阵 `J[a][i] = ∂xᵢ/∂ξₐ`（`[parametric][spatial]`）；
- `local_coordinates(x)` / `contains(x)`：逆映射与点定位。

> 通过 `T: Element` 或 `T: HdivElement` 的泛型约束可以直接调用 `ElementBase` 的方法；但对**具体类型**（如 `tri.jacobe(...)`）调用时，需要把 `ElementBase` 引入作用域。

> 注意拼写：方法名是 `jacobe` 而非 `jacobian`，这是有意保留的既有命名，属于公开 API 的一部分，请勿"修正"。

## `IntegrationScheme<P, G>`：积分规则

```rust
pub trait IntegrationScheme<const P: usize, const G: usize> {
    fn points() -> [[f64; P]; G];
    fn weights() -> [f64; G];
}
```

`G` 是积分点个数。`points()` 返回参考单元上的 `G` 个积分点坐标，`weights()` 返回对应权重。

### 完全积分与减缩积分

```rust
pub trait FullIntegration<const P: usize, const G: usize>: IntegrationScheme<P, G> {}
pub trait ReducedIntegration<const P: usize, const G: usize>: IntegrationScheme<P, G> {}
```

`FullIntegration` / `ReducedIntegration` 是**标记子 trait**：它们不增加任何方法，只为"同一套 `(P, G)` 下的某条积分规则属于完全还是减缩"提供编译期证据。组装算子（如 `stiffness`）以 `FullIntegration<P, G>` 作为约束，从而在编译期保证不会误用减缩积分。

### 标记类型

不同拓扑可能在相同的 `(P, G)` 下共享积分点个数，因此积分规则通过**标记类型**区分：

| 标记类型 | 拓扑 | 空间 |
|---|---|---|
| `GaussLegendre1D<N>` | 参考区间 `[-1, 1]` | 1D |
| `Dunavant<G>` | 参考三角形 $\{(0,0),(1,0),(0,1)\}$ | 2D |
| `TensorProductQuad<N>` | 参考四边形 `[-1, 1]²` | 2D |
| `Keast<G>` | 参考四面体 | 3D |
| `TensorProductHex<N>` | 参考立方体 `[-1, 1]³` | 3D |
| `Point` | 0D 单点 | 0D |

具体单元再为自身实现 `IntegrationScheme`，并委托给对应标记类型。例如 `Seg2` 的完全积分委托给 `GaussLegendre1D<2>`，`Tri3` 的完全积分委托给 `Dunavant<3>`。

## `Variable<D, C, P, N>`：算子输入与后处理抽象

```rust
pub trait Variable<
    const D: usize,
    const C: usize,
    const P: usize,
    const N: usize,
> {
    type Item: Element<D, P, N>;
    fn iter(&self) -> impl Iterator<Item = &Self::Item>;
}
```

`C` 是每节点自由度数/分量数；`Variable` 不再编码张量阶，分量的物理含义由具体问题模块解释。`Variable` 提供 `num_dofs`、点查询 `locate(x)`、`value`、`derivative_value`；点查询内部线性扫描单元，不存储空间索引。

组装算子（`stiffness`、`body_force` 等）通过 `Variable` 接收"一组单元"，而无需关心调用方持有的是切片还是 `Vec`。库为 `&[T]` 与 `Vec<T>` 提供了 blanket impl，因此两种容器都能直接传入：

```rust
use ansatz::element::seg2::Seg2;
use ansatz::operation::bar::stiffness;

let bars: Vec<Seg2<1>> = /* ... */;
let k1 = stiffness(&bars, ea);   // Vec<T> → 借用后按切片处理

let slice: &[Seg2<1>] = &bars;
let k2 = stiffness(slice, ea);   // 同样合法
```

## `FromGmsh<D>`：GMSH 桥接

```rust
pub trait FromGmsh<const D: usize>: Sized {
    const GMSH_ELEMENT_TYPE: i32;
    const GMSH_NODE_COUNT: usize;
    fn from_gmsh(node_ids: &[usize], nodes: &NodeTable<D>) -> Result<Self, MeshError>;
}
```

单元类型实现 `FromGmsh<D>` 后，即可通过 `GmshMesh::elements` 按物理组名构造类型化单元向量。当前已支持 `Seg2`（类型号 1）、`Tri3`（类型号 2）、`Quad4`（类型号 3）、`Poi1`（类型号 15）。

## 组装风格：`Triplet` 与 `Doublet`

组装算子不直接写矩阵/向量，而是返回**稀疏条目**，由调用方统一合并：

- 刚度项 → `Vec<Triplet<usize, usize, f64>>`（faer 的三元组），用 `SparseColMat::try_new_from_triplets` 合并；
- 载荷项 → `Vec<Doublet<usize, f64>>`（本库定义的 `(index, value)` 对），用 `try_new_from_doublet` 合并。

这种"收集条目、最后合并"的风格把"每单元贡献"与"全局数据结构"解耦，既便于并行化（M5 里程碑方向），也让稀疏矩阵的填充方式对调用方透明。
