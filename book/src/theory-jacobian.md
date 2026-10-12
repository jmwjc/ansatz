# 雅可比与坐标映射

形函数定义在参考单元（参数坐标 $ξ$）上，但积分与求导发生在物理空间（坐标 `x`）。连接两者的桥梁是**雅可比矩阵**。

## 雅可比矩阵

等参映射 $x = x(ξ)$ 的雅可比矩阵为

$$
J = \frac{\partial x}{\partial \xi}
$$

- 1D：`J` 是标量 $dx/dξ$；
- 2D：`J` 是 `2 × 2` 矩阵 $J_{ij} = \partial x_i / \partial \xi_j$。

`rustfem` 的 `Element::jacobe_mat` 返回的正是它：$J[i][d] = ∂x_i/∂ξ_d$。

## 链式法则

形函数的物理导数通过雅可比求得：

$$
\frac{dN_a}{dx} = \frac{dN_a}{d\xi} \left( \frac{dx}{d\xi} \right)^{-1} = \frac{dN_a/d\xi}{J}
$$

`Seg2` 正是这么做的：先算 `jacobe`（即 $J = dx/dξ$），再返回 $dN/dx = (dN/dξ) / J$。对线性形函数 $N_1 = (1-ξ)/2$、$N_2 = (1+ξ)/2$，有

$$
\frac{dN_1}{d\xi} = -\frac12, \quad \frac{dN_2}{d\xi} = +\frac12
$$

而等参映射在 1D 下的 $J = dx/dξ = (x_2 - x_1)/2 = L/2$，故

$$
\frac{dN_1}{dx} = -\frac{1}{L}, \qquad \frac{dN_2}{dx} = +\frac{1}{L}
$$

与 `Seg2::derivative_shape` 返回的 `[-1/L, 1/L]` 一致。

## 积分中的度量因子：`jacobe`

做变量替换 $x \to ξ$ 时，微元变换为

$$
dx = |J|\, d\xi
$$

因此面积/体积比 `|J|` 是积分被积函数必须乘上的度量因子。`Element::jacobe` 返回的正是这个量：

- 1D（`Seg2`）：$jacobe = L/2$；
- 2D 平面单元（`Tri3`/`Quad4` 在 2D 空间）：$jacobe = |det J|$。

### 高维嵌入时的推广

当单元嵌入更高维空间（如 2D 空间的三角形嵌入 3D 空间），`J` 不再是方阵，`det J` 无定义。此时用 $J^T J$ 的行列式的平方根：

$$
|J| = \sqrt{\det(J^T J)}
$$

它在 `J` 为方阵时退化为 `|det J|`。`Tri3` 与 `Quad4` 的 `jacobe` 正是按此实现（对 $D = 2$ 即 `|det J|`，对 $D = 3$ 即两切向量叉积的模）。

## 具体值

- `Seg2` 单元（端点 $x_1, x_2$）：$J = (x_2 - x_1)/2$，$jacobe = L/2$；
- `Tri3` 三角形（2D）：$jacobe = |det J| = 2 \times \text{面积}$（测试中用 $2 \times \text{面积} = 12$ 验证）；
- `Quad4` 单位正方形 $[0,1]^2$：在形心 `(0,0)` 处 $jacobe = 1/4$（参考正方形面积 4 与物理面积 1 之比）；
- `Poi1`：0 维，$jacobe = 1$。

## 为什么用 `jacobe` 而非 `det`

方法名 `jacobe` 是有意保留的既有命名（API 稳定性的一部分）。它的语义是"雅可比度量因子"，即积分变量替换中出现的 `|J|`，而不是雅可比矩阵本身（矩阵由 `jacobe_mat` 返回）。
