# 单元刚度与全局组装

弱形式的全局积分可**按单元分解**，逐单元算出**单元刚度矩阵**，再累加为全局稀疏矩阵。本章讲通用的组装机制；具体单元的刚度闭式结果见算例部分。

## 单元分解

弱形式中的积分在单元剖分 $\Omega = \bigcup_e \Omega^e$ 上分解：

$$
a(u,v) = \sum_e a^e(u,v), \qquad
l(v) = \sum_e l^e(v)
$$

每个单元上的双线性形式 $a^e$ 只涉及该单元的节点。在单元局部，场变量由该单元的 $N$ 个形函数插值：

$$
u^e = \sum_{j=1}^{N} u_j^e N_j
$$

于是**单元刚度矩阵**与**单元载荷向量**为

$$
K^e_{ij} = a^e(N_j, N_i), \qquad
F^e_i = l^e(N_i), \qquad i,j = 1,\dots,N
$$

## 数值积分计算单元刚度

$K^e_{ij}$ 是参考单元上的积分，通常用高斯求积计算（见「数值积分」章）：

$$
K^e_{ij} \approx \sum_{g} \left[ \text{被积函数在积分点 } \xi_g \text{ 的取值} \right] \cdot |J| \cdot w_g
$$

`ansatz` 的组装算子正是按此结构写的：对每个单元、每个积分点，用形函数导数、雅可比 $|J|$ 与权重 $w_g$ 累加出 `Triplet`/`Doublet` 条目。

## 全局组装

单元刚度定义在**局部自由度**上，需映射到**全局自由度**。单元的 $id()$ 返回其 $N$ 个全局节点号：

```rust
let id = elm.id();   // 例如 [3, 7]
```

单元贡献 $K^e_{ij}$ 累加到全局矩阵的 $(id_i, id_j)$ 位置：

$$
K_{id_i, id_j} \mathrel{+}= K^e_{ij}
$$

调用方用 faer 把三元组合并为稀疏矩阵：

```rust
let k = SparseColMat::<usize, f64>::try_new_from_triplets(np, np, &triplets).unwrap();
```

重复下标自动累加，不同单元共享节点处的贡献自然叠加。

## 对称性与下三角存储

当双线性形式对称（$K_{ij} = K_{ji}$）时，只需生成**下三角**（含对角线）三元组：

```rust
for i in 0..N {
    for j in 0..i + 1 {   // 仅 j ≤ i
        // ...
    }
}
```

这使条目减半，也契合 `SparseColMat` 按列稀疏存储。

## 稀疏性

每个节点只与相邻节点通过单元相连，因此全局矩阵非零元极少，是**稀疏**的。这是有限元能求解大规模问题的关键：稀疏矩阵只存储非零元，稀疏求解器（如 Cholesky/LU）利用其结构把复杂度从 $O(n^3)$ 降到接近 $O(n^{1.5})$ 量级。

## 载荷组装

载荷项同样逐单元累加为 `Doublet`（$(index, value)$ 对）：

- 体力/面力：按单元数值积分 $\int N_i f\, d\Omega \approx \sum_g N_i(\xi_g)\, f\, |J|\, w_g$；
- 集中力：直接加在对应节点自由度上；
- 罚函数载荷：见「边界条件」章。

```rust
let mut doublets = body_force(&bars, |_| 1.0);   // 体力
doublets.extend(traction(&vec![right], 1.0));     // 集中力
let f = try_new_from_doublet(np, &doublets);      // 合并为 Col<f64>
```

最终 $K$ 与 $F$ 交给稀疏求解器（对称正定时用 Cholesky）求解。
