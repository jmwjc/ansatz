# 组装与求解

`ansatz` 把"算贡献"与"建矩阵"分开：组装算子返回稀疏条目，调用方统一合并为稀疏矩阵/向量后求解。

## 组装算子

1D 桁架的全部组装算子位于 `ansatz::operation::bar`：

### `stiffness` — 单元刚度

```rust
pub fn stiffness<T, const P: usize, const N: usize, const G: usize>(
    v: &T,
    ea: f64,
) -> Vec<Triplet<usize, usize, f64>>
where
    T: Variable<1, 1, P, N>,
    T::Item: Element<1, P, N> + FullIntegration<P, G>,
```

- `ea`：轴向刚度（弹性模量 × 截面积）；
- 返回**下三角**（含对角线）的稀疏三元组。

### `body_force` — 体力（分布载荷）

```rust
pub fn body_force<T, const P: usize, const N: usize, const G: usize>(
    v: &T,
    f: impl Fn([f64; 1]) -> f64,
) -> Vec<Doublet<usize, f64>>
```

- `f`：体力密度（单位长度上的轴向力），按单元积分 $Fᵢ = ∫ Nᵢ f dx$。

### `traction` — 集中力

```rust
pub fn traction<T, const D: usize>(v: &T, p: f64) -> Vec<Doublet<usize, f64>>
```

- `v`：一组点单元（`Poi1`），每个点单元对应一个加载点；
- `p`：力的大小，施加 $Fᵢ += Nᵢ P$。

### `displacement_penalty` — 位移约束（罚函数法）

```rust
pub fn displacement_penalty<T, const D: usize>(
    v: &T,
    penalty: f64,
    u: f64,
) -> (Vec<Triplet<usize, usize, f64>>, Vec<Doublet<usize, f64>>)
```

- `v`：一组点单元（`Poi1`），每个点单元对应一个被约束点；
- `penalty`：罚数 $α$（示例取 `1e7`）；
- `u`：指定位移 $ū$；
- 返回一对 `(刚度三元组, 载荷双元组)`：$Kᵢⱼ += α Nᵢ Nⱼ$，$Fᵢ += α Nᵢ ū$。

## 合并为稀疏矩阵/向量

### 载荷向量

```rust
use ansatz::sparse_vector::try_new_from_doublet;

let f = try_new_from_doublet(np, &doublets).unwrap();   // → faer::Col<f64>
```

`try_new_from_doublet` 把 `Doublet` 累加到长度为 `np` 的稠密向量中（重复下标自然累加）。

### 刚度矩阵

```rust
use faer::sparse::SparseColMat;

let k = SparseColMat::<usize, f64>::try_new_from_triplets(np, np, &triplets).unwrap();
```

`np` 是总自由度（节点）数。

## 求解

刚度矩阵对称正定，使用稀疏 Cholesky 分解：

```rust
use faer::prelude::*;
use faer::Side;

let llt = k.sp_cholesky(Side::Lower).unwrap();
let d = llt.solve(&f);
```

## 完整链路

```rust
use faer::Side;
use faer::prelude::*;
use faer::sparse::SparseColMat;
use ansatz::element::poi1::Poi1;
use ansatz::element::seg2::Seg2;
use ansatz::node::Node;
use ansatz::operation::bar::{
    body_force, point_displacement_penalty, point_load, stiffness,
};
use ansatz::sparse_vector::try_new_from_doublet;

// 1. 构造单元
let bars: Vec<Seg2<1>> = /* ... */;
let left  = Poi1::new(Node::new(0, [0.0]));       // 左端约束点
let right = Poi1::new(Node::new(10, [1.0]));      // 右端加载点

// 2. 刚度：单元 + 罚函数约束
let mut triplets = stiffness(&bars, 1.0);
let (pen_k, pen_f) = point_displacement_penalty(&vec![left], 1e7, 0.0);
triplets.extend(pen_k);

// 3. 载荷：体力 + 集中力 + 罚函数载荷
let mut doublets = body_force(&bars, |_| 1.0);
doublets.extend(point_load(&vec![right], [1.0]));
doublets.extend(pen_f);

// 4. 合并并求解
let k = SparseColMat::<usize, f64>::try_new_from_triplets(11, 11, &triplets).unwrap();
let f = try_new_from_doublet(11, &doublets).unwrap();
let llt = k.sp_cholesky(Side::Lower).unwrap();
let d = llt.solve(&f);
```

> 说明：$np = 11$ 是节点总数；`stiffness` 与 `body_force` 中 `FullIntegration` 约束由 `Seg2<1>` 的完全积分（2 点 Gauss）满足。
