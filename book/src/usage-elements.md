# 单元类型

本库目前提供四种单元，均以 `Node<D>` 作为顶点存储。

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
use rustfem::element::seg2::Seg2;
use rustfem::node::Node;

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
use rustfem::element::poi1::Poi1;
use rustfem::node::Node;

// 表示"节点 0 处施加约束/载荷"
let left = Poi1::new(Node::new(0, [0.0]));
```

## 关于 `derivative_shape` 的一个设计细节

`Element` trait 的 `derivative_shape` 返回 `[f64; N]`，对 2D 单元（`Tri3`/`Quad4`）而言，这个长度不足以容纳完整梯度（对 $ξ$、$η$ 各一个偏导）。因此：

- 对 1D 单元（`Seg2`），`derivative_shape` 直接返回**物理导数** `dN/dx`（已除以雅可比）；
- 对 2D 单元，`derivative_shape` 只返回对第一个参数坐标 $ξ$ 的偏导 $∂N/∂ξ$；对 $η$ 的偏导以及物理梯度需结合 `jacobe_mat` 另行计算。

这是当前 1D 求解链路先行、2D 组装尚未接入的直接体现；2D 的 B 矩阵与单元刚度组装属于 `PLAN.md` 的 M3 里程碑。
