# 单元类型

本库目前提供十二种单元，均以 `Node<D>` 作为节点存储；对协调单元（`Seg2`/`Tri3`/`Quad4`/`Poi1`/`Tet4`）节点就是顶点；非协调的 `CRTri3` 节点是边中点、`CRTet4` 节点是面重心；分片常数 `P0Tri3` 的"节点"是单元重心，只用来携带单元自由度；边界迹单元 `CRSeg3`/`CRTri4`/`RTSeg3` 的"节点"是其父单元的顶点；`RTTri3` 是唯一**向量值**单元（$H(\mathrm{div})$ 协调），其"节点"同样是边中点，存放的是边通量自由度。

## `Seg2<D>`：2 节点线单元

```rust
pub struct Seg2<const D: usize> { /* ... */ }

impl<const D: usize> Seg2<D> {
    pub fn new(n1: Node<D>, n2: Node<D>) -> Self;
}
```

参考单元为 1D 区间 `[-1, 1]`。形函数为线性拉格朗日插值：

$$
N_1(\xi) = \frac{1-\xi}{2}, \qquad N_2(\xi) = \frac{1+\xi}{2}
$$

`Seg2<D>` 可嵌入 1D/2D/3D 空间（`D` 任意），例如：

```rust
use ansatz::element::seg2::Seg2;
use ansatz::node::Node;

// 1D 杆单元，端点 x = 0 与 x = 0.1
let e = Seg2::new(Node::new(0, [0.0]), Node::new(1, [0.1]));

// 2D 空间中的杆单元（如桁架）
let e2 = Seg2::new(Node::new(0, [0.0, 0.0]), Node::new(1, [1.0, 1.0]));
```

积分方案：

- 完全积分 `FullIntegration<1, 2>`：`GaussLegendre1D<2>`（2 点，代数精度 3）；
- 减缩积分 `ReducedIntegration<1, 1>`：`GaussLegendre1D<1>`（1 点，代数精度 1）。

## `Tri3<D>`：3 节点三角形

```rust
pub struct Tri3<const D: usize> { /* ... */ }

impl<const D: usize> Tri3<D> {
    pub fn new(n1: Node<D>, n2: Node<D>, n3: Node<D>) -> Self;
}
```

参考三角形顶点为 `(0, 0)`、`(1, 0)`、`(0, 1)`。形函数是**面积坐标**：

$$
N_1(\xi,\eta) = 1 - \xi - \eta, \qquad N_2(\xi,\eta) = \xi, \qquad N_3(\xi,\eta) = \eta
$$

积分方案：

- 完全积分 `FullIntegration<2, 3>`：`Dunavant<3>`（3 点，代数精度 2）；
- 减缩积分 `ReducedIntegration<2, 1>`：`Dunavant<1>`（1 点，形心）。

## `Tet4<D>`：4 节点线性四面体

```rust
pub struct Tet4<const D: usize> { /* ... */ }

impl<const D: usize> Tet4<D> {
    pub fn new(n1: Node<D>, n2: Node<D>, n3: Node<D>, n4: Node<D>) -> Self;
    pub fn volume(&self) -> f64;
}
```

参考四面体为 `(0,0,0)`、`(1,0,0)`、`(0,1,0)`、`(0,0,1)`，形函数就是体积（重心）坐标：

$$
N_1 = 1 - \xi - \eta - \zeta,\qquad N_2 = \xi,\qquad N_3 = \eta,\qquad N_4 = \zeta
$$

物理梯度与逆映射都通过度量张量 $g_{\alpha\beta} = J_\alpha\cdot J_\beta$ 计算（与 `Tri3` 同一套公式），因此 `jacobe` 为 $6\times$ 体积、对嵌入高维（`D > 3`）也成立。实现了 `FromGmsh`（GMSH 类型 4）。

积分方案：

- 完全积分 `FullIntegration<3, 4>`：`Keast<4>`（4 点，代数精度 2，协调质量矩阵需要）；
- 减缩积分 `ReducedIntegration<3, 1>`：`Keast<1>`（1 点，线性四面体刚度被积为常数）。

## `CRTri3<D>`：Crouzeix–Raviart 非协调线性三角形

```rust
pub struct CRTri3<const D: usize> { /* ... */ }

impl<const D: usize> CRTri3<D> {
    pub fn new(m1: Node<D>, m2: Node<D>, m3: Node<D>) -> Self;
    pub fn from_tri3(tri: &Tri3<D>, edge_id: impl FnMut(usize, usize) -> usize) -> Self;
    pub fn triangle_vertices(&self) -> [[f64; D]; 3];
}
```

三个自由度是**边的平均值**（对线性函数即边中点值）：`nodes[i]` 存"顶点 `i` 的对边"的中点，`Node::id` 是该边的全局自由度编号。以面积坐标 $\lambda_i$ 表示，基函数为

$$
\theta_i = 1 - 2\lambda_i \quad (i = 1,2,3),
$$

满足 $\frac{1}{|e_j|}\int_{e_j}\theta_i\,ds = \delta_{ij}$。它**不是插值型**单元（在顶点取 $\pm1$），全局空间只保证跨边**平均**连续，所以 $V_h \not\subset H^1$；相应地 2D 质量矩阵天然对角（$|T|/3$）。

自由度挂在边上，网格顶点 id 不能直接当自由度编号，构造单元时需给出**边 → 自由度**映射：

```rust
use std::collections::HashMap;
use ansatz::element::crtri3::CRTri3;

let n_vertices = 3;
let mut edges: HashMap<(usize, usize), usize> = HashMap::new();
let cr = CRTri3::from_tri3(&tri, |a, b| {
    let next = edges.len();
    *edges.entry((a, b)).or_insert(next + n_vertices) // 边自由度接在顶点之后
});
```

`CRTri3` 没有 `FromGmsh` 实现：GMSH 的三角形只给出顶点 id，边到自由度的编号必须另外维护。`io::gmsh` 为此提供了 `GmshMesh::edge_table()` 与 `GmshMesh::crtri3_elements(name, &mut edges)`（见「GMSH 网格读取」章的「读取 Crouzeix–Raviart（CR）单元」）。

积分方案：

- 完全积分 `FullIntegration<2, 3>`：`Dunavant<3>`（质量矩阵被积函数为二次，必须 3 点）；
- 减缩积分 `ReducedIntegration<2, 1>`：`Dunavant<1>`（刚度被积函数为常数，1 点即精确）。

## `CRTet4<D>`：Crouzeix–Raviart 非协调线性四面体

```rust
pub struct CRTet4<const D: usize> { /* ... */ }

impl<const D: usize> CRTet4<D> {
    pub fn new(f1: Node<D>, f2: Node<D>, f3: Node<D>, f4: Node<D>) -> Self;
    pub fn from_tet4(tet: &Tet4<D>, face_id: impl FnMut([usize; 3]) -> usize) -> Self;
    pub fn from_vertices(vertices: [Node<D>; 4], face_id: impl FnMut([usize; 3]) -> usize) -> Self;
    pub fn tetrahedron_vertices(&self) -> [[f64; D]; 4];
    pub fn volume(&self) -> f64;
}
```

`CRTri3` 的三维版本，实现 `Element<D, 3, 4>`：四个自由度是**面的平均值**（线性函数上等于面重心值），`nodes[a]` 存"顶点 `a` 的对面"的重心，`Node::id` 是该面的全局自由度编号。基函数

$$
\theta_a = 1 - 3\lambda_a \quad (a = 1,\dots,4),
$$

满足 $\frac{1}{|f_b|}\int_{f_b}\theta_a\,dA = \delta_{ab}$；`facet` 编号同样是"对顶点"。几何由四个面重心按 $\sum_a\theta_a C_a = \sum_a\lambda_a A_a$ 精确重现，真实顶点用 `CRTet4::tetrahedron_vertices()`（$A_a = \sum_b C_b - 3C_a$）取回。

与二维不同，三维 CR 的**质量矩阵不再对角**（$M_{aa} = 2|T|/5$、$M_{ab} = -|T|/20$），所以要用 4 点 Keast 规则（`FullIntegration<3, 4>`）；刚度被积函数为常数，1 点（`ReducedIntegration<3, 1>`）即精确。面的全局编号由调用方（或用 `GmshOptions { faces: true }` 时的 gmsh 读取器）维护：`from_tet4`/`from_vertices` 的闭包接收一个面的三个角点 id（升序）并返回该面的自由度编号；无法失败的查表用 `from_*`，需要传播错误的用 `try_from_tet4`/`try_from_vertices`。

## `CRSeg3<D>` / `CRTri4<D>`：CR 的 facet 迹单元（边界条件专用）

```rust
pub struct CRSeg3<const D: usize> { /* ... */ }  // 2D：线段 facet，3 个自由度
pub struct CRTri4<const D: usize> { /* ... */ }  // 3D：三角 facet，4 个自由度

impl<const D: usize> CRSeg3<D> {
    pub fn new(facets: [usize; 3], vertices: [Node<D>; 3], facet: usize) -> Option<Self>;
    pub fn from_crtri3(parent: &CRTri3<D>, facet: usize) -> Option<Self>;  // CRSeg3
    // CRTri4::from_crtet4(&CRTet4<D>, facet)
    pub fn local_facet(&self) -> usize;
}
```

CR 在一条边上的迹依赖相邻三角形的**全部三个**边自由度，所以边界积分（Neumann 通量、Robin、Dirichlet 罚函数/乘子）需要这个"父单元的迹"单元，而不是 `Seg2`：

- 命名沿袭 `Seg2`/`Tri3` 的习惯："**自身形态 + 自由度个数**"：`CRSeg3` 是 3 自由度的 CR 线段（`CRTri3` 的迹），`CRTri4` 是 4 自由度的 CR **三角形**（`CRTet4` 的面迹）。名字里的 `Tri` 说的是单元自身的形态，四面体的面本来就是三角形；`Tet` 留给 3D 内部单元；
- `facet = k` 表示"父顶点 $k$ 不在其上的 facet"，与 `CRTri3` 的自由度顺序一致，因此 `id()` 就是父单元的 `id()`；
- 迹基（2D，`lo < hi` 为 facet 两顶点下标）：`shape[facet] = 1`、`shape[lo] = ξ`、`shape[hi] = −ξ`；
- 几何取 facet 顶点构成的 `Seg2`/`Tri3`（`jacobe = L/2` 或 `2×面积`），是**非等参**的；
- 积分：`CRSeg3` 用 `FullIntegration<1, 2>` / `ReducedIntegration<1, 1>`，`CRTri4` 用 `FullIntegration<2, 3>` / `ReducedIntegration<2, 1>`。

## `RTTri3<D>`：Raviart–Thomas（RT0）三角形

```rust
pub struct RTTri3<const D: usize> { /* ... */ }

impl<const D: usize> HdivElement<D, 2, 3> for RTTri3<D> {
    fn shape(&self, xi: [f64; 2]) -> [[f64; D]; 3];
    fn divergence_shape(&self, xi: [f64; 2]) -> [f64; 3];
}

impl<const D: usize> RTTri3<D> {
    pub fn new(m1: Node<D>, m2: Node<D>, m3: Node<D>) -> Self;
    pub fn from_tri3(tri: &Tri3<D>, edge_id: impl FnMut(usize, usize) -> usize) -> Self;
    pub fn triangle_vertices(&self) -> [[f64; D]; 3];
    pub fn area(&self) -> f64;
    pub fn reference_vector_shape(&self, xi: [f64; 2]) -> [[f64; 2]; 3];
    pub fn mass_matrix(&self) -> [[f64; 3]; 3];
    pub fn flux_dof(&self, edge: usize, q: impl Fn([f64; D]) -> [f64; D]) -> Option<f64>;
}
```

最低阶 Raviart–Thomas 单元（RT0），局部空间是仿射向量场 $\boldsymbol q(\boldsymbol x)=\boldsymbol a+c\,\boldsymbol x$，三个自由度是三条边上的**法向通量**

$$
\Phi_e(\boldsymbol q)=\int_{\Gamma_e}\boldsymbol q\cdot\boldsymbol n_e\,ds ,
$$

对应基函数 $\boldsymbol q_I=(\boldsymbol x-\boldsymbol x_I)/(2|T|)$（$\boldsymbol x_I$ 为边 $\Gamma_I$ 的对顶点，$\boldsymbol n_I$ 为外法向），满足 $\Phi_J(\boldsymbol q_I)=\delta_{IJ}$、$\nabla\cdot\boldsymbol q_I=1/|T|$。全局空间 $H(\mathrm{div})$ 协调：切向可以跳跃，法向通量跨边单值。

`RTTri3` 与 `CRTri3` 的数据布局完全相同（`nodes[i]` 是边 $i$ 的中点，`Node::id` 是该边的全局自由度），因此可以共用 `EdgeTable`；但它的基底是**向量**，无法用 `Element::shape`（每个自由度一个标量）表达，所以它不实现标量 `Element`，而实现 `ElementBase<D, 2, 3>`（几何与自由度）与 `HdivElement<D, 2, 3>`（`HdivElement::shape` / `divergence_shape`）；上面代码列出的是剩下不进入通用抽象的固有工具方法。质量矩阵是**稠密**的（与 CR 的对角质量矩阵不同），用 `FullIntegration<2, 3>` / `Dunavant<3>`；散度为常数，`ReducedIntegration<2, 1>` 即精确。理论、与 CR 的等价（杂交化 RT0 凝聚后即 CR 刚度）见「混合格式与 $H(\mathrm{div})$：Raviart–Thomas（RT0）」章。 形函数与分量的可视化见 `cargo run --example rt0`（输出 `target/rt0_shape_functions.png` 与 `target/rt0_shape_components.png`）。

## `P0Tri3<D>`：分片常数（P0）三角形

```rust
pub struct P0Tri3<const D: usize> { /* ... */ }

impl<const D: usize> P0Tri3<D> {
    pub fn new(vertices: [Node<D>; 3], dof_id: usize) -> Self;
    pub fn from_tri3(tri: &Tri3<D>, dof_id: usize) -> Self;
    pub fn vertices(&self) -> &[Node<D>; 3];
    pub fn centroid(&self) -> &Node<D>;
}
```

局部空间为常数，`P = 2`、`N = 1`：形函数恒为 1、梯度为 0，`jacobe = 2×面积`，几何（逆映射、`contains`）委托给三个角点构成的 `Tri3`。唯一的"节点"是重心，只承载全局**单元自由度** id（因此不能按顶点编号）。实现 `Element<D, 2, 1>`，积分 `FullIntegration<2, 3>` / `ReducedIntegration<2, 1>`。

典型用途是混合（Hellinger–Reissner）热传导中的温度空间（与 `RTTri3` 组成 RT0/P0 稳定对），见「热传导方程」§7 与 `heat::mixed_divergence`。

## `RTSeg3<D>`：RT0 的边界迹

```rust
pub struct RTSeg3<const D: usize> { /* ... */ }

impl<const D: usize> RTSeg3<D> {
    pub fn new(facets: [usize; 3], vertices: [Node<D>; 3], facet: usize) -> Option<Self>;
    pub fn from_rttri3(parent: &RTTri3<D>, facet: usize) -> Option<Self>;
    pub fn local_facet(&self) -> usize;
    pub fn normal_flux_shape(&self, xi: [f64; 1]) -> [f64; 3];
}
```

`RTTri3` 的一条边上的法向迹。三个自由度是**父三角形的全部三个边通量自由度**（`id() == 父单元 id()`），`facet = k` 表示"父顶点 $k$ 的对边"；几何取该边两端点构成的 `Seg2`（`jacobe = L/2`），是**非等参**的。`normal_flux_shape` 返回 $\boldsymbol q_i\cdot\boldsymbol n$（$i$ 为父单元的三个通量 dof），RT0 下为 $\delta_{i,\mathrm{facet}}/\ell_{\mathrm{facet}}$。

它实现 `ElementBase<D, 1, 3>` 与 `IntegrationScheme<1, 1>` / `<1, 2>`，但不实现标量 `Element`。用途是混合格式中 Dirichlet 数据的自然边界项 $\ell_i=-\int_{\Gamma_D}T_D\,\boldsymbol q_i\cdot\boldsymbol n\,ds$（见 `heat::mixed_dirichlet`）。

## `Quad4<D>`：4 节点四边形

```rust
pub struct Quad4<const D: usize> { /* ... */ }

impl<const D: usize> Quad4<D> {
    pub fn new(n1: Node<D>, n2: Node<D>, n3: Node<D>, n4: Node<D>) -> Self;
}
```

参考四边形为 `[-1, 1]²`，节点 1–4 依次对应 `(-1,-1)`、`(1,-1)`、`(1,1)`、`(-1,1)`。形函数为双线性：

$$
N_1 = \frac{(1-\xi)(1-\eta)}{4}, \quad
N_2 = \frac{(1+\xi)(1-\eta)}{4}, \quad
N_3 = \frac{(1+\xi)(1+\eta)}{4}, \quad
N_4 = \frac{(1-\xi)(1+\eta)}{4}
$$

积分方案：

- 完全积分 `FullIntegration<2, 4>`：`TensorProductQuad<2>`（2×2 = 4 点，代数精度 3）；
- 减缩积分 `ReducedIntegration<2, 1>`：`TensorProductQuad<1>`（1 点，形心）。

## `Poi1<D>`：1 节点点单元

```rust
pub struct Poi1<const D: usize> { /* ... */ }

impl<const D: usize> Poi1<D> {
    pub fn new(node: Node<D>) -> Self;
}
```

点单元的参数空间维数为 0（$P = 0$），形函数恒等于 1，雅可比恒等于 1。它不产生刚度，只用于**定位**——代表一个约束点（`displacement_penalty`）或一个集中力作用点（`traction`）：

```rust
use ansatz::element::poi1::Poi1;
use ansatz::node::Node;

// 表示"节点 0 处施加约束/载荷"
let left = Poi1::new(Node::new(0, [0.0]));
```

## 关于 `derivative_shape` 的一个设计细节

`Element` trait 的 `derivative_shape` 返回 `[f64; N]`，对 2D 单元（`Tri3`/`Quad4`）而言，这个长度不足以容纳完整梯度（对 $ξ$、$η$ 各一个偏导）。因此：

- 对 1D 单元（`Seg2`），`derivative_shape` 直接返回**物理导数** `dN/dx`（已除以雅可比）；
- 对 2D 单元，`derivative_shape` 只返回对第一个参数坐标 $ξ$ 的偏导 $∂N/∂ξ$；对 $η$ 的偏导以及物理梯度需结合 `jacobe_mat` 另行计算。

这是当前 1D 求解链路先行、2D 组装尚未接入的直接体现；2D 的 B 矩阵与单元刚度组装属于 `PLAN.md` 的 M3 里程碑。
