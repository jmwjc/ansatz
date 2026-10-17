# 雅可比与坐标映射

形函数定义在参考单元（参数坐标 $ξ$）上，但积分与求导发生在物理空间（坐标 `x`）。连接两者的桥梁是**雅可比矩阵**。

## 雅可比矩阵

等参映射 $x = x(ξ)$ 的雅可比矩阵为

$$
J = \frac{\partial x}{\partial \xi}
$$

- 1D 单元（`Seg2` 嵌入任意维度）：`J` 是一个 $D$ 维切向量 $\partial x/\partial \xi$；
- 2D 单元（`Tri3`/`Quad4`）：`J` 由两个 $D$ 维切向量组成（$\partial x/\partial \xi$ 与 $\partial x/\partial \eta$）。

`ansatz` 的 `Element::jacobe_mat` 返回 `[[f64; D]; P]`，索引为 `[parametric][spatial]`：$J[d][i] = \partial x_i/\partial \xi_d$，即每个参数方向对应一个 $D$ 维切向量。

## 链式法则

形函数的物理导数通过雅可比求得：

$$
\frac{dN_a}{dx} = \frac{dN_a}{d\xi} \left( \frac{dx}{d\xi} \right)^{-1} = \frac{dN_a/d\xi}{J}
$$

`Seg2` 正是这么做的：先算 `jacobe`（即 $L/2$），再返回物理梯度 $\nabla_x N = \pm \hat t / L$，其中 $\hat t = (x_2 - x_1)/L$ 是沿单元切向的单位向量。对线性形函数 $N_1 = (1-ξ)/2$、$N_2 = (1+ξ)/2$，有

$$
\frac{dN_1}{d\xi} = -\frac12, \quad \frac{dN_2}{d\xi} = +\frac12
$$

而等参映射在 1D 下的 $dx/dξ = (x_2 - x_1)/2$，$L = |x_2 - x_1|$，故物理梯度（$D$ 维向量）

$$
\nabla_x N_1 = -\frac{x_2 - x_1}{L^2}, \qquad \nabla_x N_2 = +\frac{x_2 - x_1}{L^2}
$$

在 1D（$D = 1$）下退化为 $\mp 1/L$，与 `Seg2::derivative_shape` 一致。

## 三角形单元的物理梯度（张量表示）

2D/3D 单元的物理梯度用**张量记号**推导最清晰。以三角形 `Tri3` 为例：三个顶点 $x_1, x_2, x_3 \in \mathbb{R}^D$（$D = 2$ 或 $3$），等参映射

$$
x(\xi, \eta) = N_1 x_1 + N_2 x_2 + N_3 x_3,
\qquad
N_1 = 1 - \xi - \eta,\; N_2 = \xi,\; N_3 = \eta
$$

下面希腊字母 $\alpha, \beta \in \{\xi, \eta\}$ 遍历参数方向，拉丁字母 $i = 1, \dots, D$ 遍历空间方向。

### 协变切向量与度量张量

两个协变切向量（各为 $D$ 维）

$$
j_\xi = \frac{\partial x}{\partial \xi} = x_2 - x_1,
\qquad
j_\eta = \frac{\partial x}{\partial \eta} = x_3 - x_1
$$

张成三角形的切平面。**度量张量**（第一基本形式）$g_{\alpha\beta} = j_\alpha \cdot j_\beta$ 是 $2 \times 2$ 矩阵

$$
g =
\begin{bmatrix}
j_\xi \cdot j_\xi & j_\xi \cdot j_\eta \\
j_\xi \cdot j_\eta & j_\eta \cdot j_\eta
\end{bmatrix}
$$

其行列式 $\det g = (j_\xi \cdot j_\xi)(j_\eta \cdot j_\eta) - (j_\xi \cdot j_\eta)^2$，**逆度量** $g^{\alpha\beta}$ 为

$$
g^{\xi\xi} = \frac{j_\eta \cdot j_\eta}{\det g},
\qquad
g^{\xi\eta} = g^{\eta\xi} = -\frac{j_\xi \cdot j_\eta}{\det g},
\qquad
g^{\eta\eta} = \frac{j_\xi \cdot j_\xi}{\det g}
$$

### 物理梯度通式

链式法则 $\partial_\alpha N = (\nabla_x N) \cdot j_\alpha$ 的解是

$$
\nabla_x N_a = j_\alpha\, g^{\alpha\beta}\, \partial_\beta N_a
$$

（对重复指标 $\alpha, \beta$ 求和）。展开为

$$
\nabla_x N_a =
j_\xi \left(g^{\xi\xi} \partial_\xi N_a + g^{\xi\eta} \partial_\eta N_a\right)
+
j_\eta \left(g^{\eta\xi} \partial_\xi N_a + g^{\eta\eta} \partial_\eta N_a\right)
$$

三角形的参数导数是常数，$(\partial_\xi N_a, \partial_\eta N_a)$ 依次为 $(-1,-1)$、$(1,0)$、$(0,1)$。

### D = 2 与 D = 3

- $D = 2$：$j_\xi, j_\eta$ 是平面内向量，$g$ 满秩，上式即 $\nabla_x N = J^{-\mathsf T} \nabla_\xi N$；
- $D = 3$：$j_\xi, j_\eta$ 是空间向量，$J$（$3 \times 2$）非方阵，$g = J^\mathsf T J$ 是切平面上的度量张量，上式给出**切平面内的梯度**（最小二乘解），对曲面三角形依然成立。

`Tri3::derivative_shape` 正是按此通式直接实现，不依赖任何辅助函数，对 $D = 2$ 与 $D = 3$ 均成立。

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
