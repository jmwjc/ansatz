# 悬臂梁的力学模型（平面应力）

悬臂梁受端部横向荷载是线弹性有限元最经典的验证算例：它的**精确解**在平面应力假设下存在闭式表达（Saint-Venant 半逆解法），因此可以直接用数值解与解析解对比，量化离散误差。本章给出完整的强形式 → 精确解 → 应力/应变 → 边界条件自洽性推导；对应的可运行算例与误差验证见「悬臂梁算例（cantilever）」。

> 记法说明：梁长 $L$、梁高 $H$（沿 $y$ 方向）、宽度 $b$（沿 $z$ 方向）。本章的 $H$ 即原始笔记里的 $D$，这里改用 $H$ 以免与空间维数 $D$ 混淆。

## 1. 问题设定

考虑矩形域

$$
\Omega = [0, L] \times \left[-\frac{H}{2}, \frac{H}{2}\right] \subset \mathbb{R}^2
$$

左端 $x = 0$ 固定，右端 $x = L$ 承受沿 $y$ 方向的横向集中力 $P$（按 Saint-Venant 原理等效为端部剪力分布）。材料为各向同性线弹性体，杨氏模量 $E$、泊松比 $\nu$，按**平面应力**假设处理（薄板：$\sigma_{zz} = \sigma_{xz} = \sigma_{yz} = 0$，见「线弹性力学」章的平面应力/平面应变一节）。

截面惯性矩（宽 $b = 1$）为

$$
I = \frac{b H^3}{12}.
$$

问题无体力，$\mathbf f = \mathbf 0$。

## 2. 精确解

### 应力函数

取 Airy 应力函数

$$
\phi(x, y) = \frac{P}{2I}\left[\frac{(x - L)\,y^3}{3} - \frac{H^2 x y}{4}\right],
$$

由 $\sigma_{xx} = \partial^2\phi/\partial y^2$、$\sigma_{yy} = \partial^2\phi/\partial x^2$、$\sigma_{xy} = -\partial^2\phi/\partial x\partial y$ 得应力场

$$
\begin{Bmatrix}
\sigma_{xx} \\ \sigma_{yy} \\ \sigma_{xy}
\end{Bmatrix}
=
\begin{Bmatrix}
-\dfrac{P(L-x)\,y}{I} \\[2mm]
0 \\[2mm]
\dfrac{P}{2I}\left(\dfrac{H^2}{4} - y^2\right)
\end{Bmatrix}.
$$

该应力场自动满足无体力的平衡方程（见第 4 节），并在上下表面 $y = \pm H/2$ 上给出零面力。

### 位移场

在平面应力本构下积分几何方程，并用 $x = 0$ 处 $u = 0$、$y = 0$ 处 $v = 0$ 等条件固定刚体位移，得精确位移场

$$
\left\{
\begin{array}{l}
u(x, y) = -\dfrac{P y}{6 E I}\left[(6L - 3x)x + (2 + \nu)\left(y^2 - \dfrac{H^2}{4}\right)\right] \\[3mm]
v(x, y) = \dfrac{P}{6 E I}\left[3\nu y^2 (L - x) + (4 + 5\nu)\dfrac{H^2 x}{4} + (3L - x)x^2\right]
\end{array}
\right.
$$

### 位移梯度

对位移直接求偏导：

$$
\left\{
\begin{array}{l}
u_{,x} = -\dfrac{P y}{E I}(L - x) \\[3mm]
u_{,y} = \dfrac{P}{6 E I}\left[-(6L - 3x)x - (2 + \nu)\left(3y^2 - \dfrac{H^2}{4}\right)\right] \\[3mm]
v_{,x} = \dfrac{P}{6 E I}\left[(6L - 3x)x - 3\nu y^2 + (4 + 5\nu)\dfrac{H^2}{4}\right] \\[3mm]
v_{,y} = \nu\,\dfrac{P y}{E I}(L - x)
\end{array}
\right.
$$

## 3. 应变与应力

把位移梯度代入几何方程，得应变场（采用工程剪应变 $2\varepsilon_{xy}$）

$$
\begin{Bmatrix}
\varepsilon_{xx} \\ \varepsilon_{yy} \\ 2\varepsilon_{xy}
\end{Bmatrix}
=
\begin{Bmatrix}
u_{,x} \\ v_{,y} \\ u_{,y} + v_{,x}
\end{Bmatrix}
=
\frac{P}{E I}
\begin{Bmatrix}
-(L - x)\,y \\[1mm]
\nu (L - x)\,y \\[1mm]
(1 + \nu)\left(\dfrac{H^2}{4} - y^2\right)
\end{Bmatrix}.
$$

代入平面应力本构关系（见「线弹性力学」章）

$$
\mathbf D_{\text{plane stress}} = \frac{E}{1-\nu^2}
\begin{bmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & \dfrac{1-\nu}{2}
\end{bmatrix}
$$

即得第 2 节的应力场：

$$
\begin{Bmatrix}
\sigma_{xx} \\ \sigma_{yy} \\ \sigma_{xy}
\end{Bmatrix}
=
\frac{E}{1-\nu^2}
\begin{bmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & \dfrac{1-\nu}{2}
\end{bmatrix}
\begin{Bmatrix}
\varepsilon_{xx} \\ \varepsilon_{yy} \\ 2\varepsilon_{xy}
\end{Bmatrix}
=
\begin{Bmatrix}
-\dfrac{P(L-x)y}{I} \\[2mm]
0 \\[2mm]
\dfrac{P}{2I}\left(\dfrac{H^2}{4} - y^2\right)
\end{Bmatrix}.
$$

## 4. 平衡与边界条件自洽性

### 无体力

将应力代入平衡方程（$\mathbf b$ 为体力）：

$$
\mathbf b = -
\begin{Bmatrix}
\sigma_{xx,x} + \sigma_{xy,y} \\
\sigma_{yx,x} + \sigma_{yy,y}
\end{Bmatrix}
=
-
\begin{Bmatrix}
\dfrac{P y}{I} + \left(-\dfrac{P y}{I}\right) \\[2mm]
0
\end{Bmatrix}
= \mathbf 0,
$$

故问题确为无体力，与设定一致。

### 自然边界条件（上、下表面与右端）

- **上、下表面** $y = \pm H/2$：外法向 $\mathbf n = (0, \pm 1)$，面力

  $$
  \mathbf t = \boldsymbol\sigma \cdot \mathbf n
  = \begin{Bmatrix} \sigma_{xy} \\ \sigma_{yy} \end{Bmatrix}
  = \begin{Bmatrix} \dfrac{P}{2I}\left(\dfrac{H^2}{4} - y^2\right) \\ 0 \end{Bmatrix}\Bigg|_{y=\pm H/2}
  = \mathbf 0,
  $$

  即上下表面为**自由面**（齐次 Neumann，等效为自然边界条件，组装时无需处理）。

- **右端** $x = L$：外法向 $\mathbf n = (1, 0)$，面力

  $$
  \mathbf t = \boldsymbol\sigma \cdot \mathbf n
  = \begin{Bmatrix} \sigma_{xx} \\ \sigma_{xy} \end{Bmatrix}\Bigg|_{x=L}
  = \begin{Bmatrix} 0 \\ \dfrac{P}{2I}\left(\dfrac{H^2}{4} - y^2\right) \end{Bmatrix}.
  $$

  其合力恰为剪力 $P$：

  $$
  \int_{-H/2}^{H/2} \sigma_{xy}\, dy
  = \frac{P}{2I}\left(\frac{H^2}{4}H - \frac{H^3}{12}\right)
  = \frac{P H^3}{12 I}
  = P \quad (\text{因 } I = H^3/12).
  $$

  因此在右端施加**抛物线形剪力分布**，等价于端部集中力 $P$。

### 本质边界条件（左端）

左端 $x = 0$ 不再取"位移为零"的固支，而是**施加精确位移**作为 Dirichlet 条件：

$$
u(0, y) = -\frac{P y}{6 E I}(2 + \nu)\left(y^2 - \frac{H^2}{4}\right), \qquad
v(0, y) = \frac{P \nu L}{2 E I} y^2.
$$

这样精确解同时满足域内方程与全部边界条件，成为该定解问题的真正解——有限元解应向它收敛。这与"完全固支梁"（$u = v = 0$ 于 $x=0$）不同：后者的应力场在端部附近与 Saint-Venant 解有差异。

## 5. 与有限元离散的对应

把精确解作为验证基准时，离散格式与边界条件一一落位如下：

| 连续问题 | 有限元实现 |
|---|---|
| 平面应力本构 $\mathbf D_{\text{plane stress}}$ | `stiffness(v, λ*, μ)`，其中 $\lambda^* = E\nu/(1-\nu^2)$、$\mu = E/(2(1+\nu))$ |
| 左端 Dirichlet 位移 | `displacement_penalty`（罚函数法） |
| 右端抛物线剪力 $\mathbf t$ | `traction`（自然边界条件） |
| 上下表面 $\mathbf t = 0$ | 不施加，自动满足 |
| 无体力 | 不施加 `body_force` |

由于精确位移场含三次项，双线性四边形（`Quad4`）无法精确表示，数值解存在离散误差。误差用「线弹性力学」章的两个算子度量：

$$
\|u_h - u\|_{L^2} = \left(\int_\Omega |u_h - u|^2\,d\Omega\right)^{1/2},
\qquad
\|u_h - u\|_E = \left(\int_\Omega (\varepsilon_h - \varepsilon) : \mathbf C : (\varepsilon_h - \varepsilon)\,d\Omega\right)^{1/2}
$$

分别由 `l2_error` 与 `energy_error` 计算。随网格加密，两种误差均应单调下降——这构成对弹性算子的端到端验证。

## 6. 数值验证结果

取参数 $E = 1000$、$\nu = 0.3$、$P = 1$、$L = 8$、$H = 1$（宽 $b = 1$，故 $I = 1/12$）。左端施加本节的精确位移（罚函数法），右端施加本节的抛物线剪力，网格由 `examples/msh/cantilever.geo` 生成。

### 与解析解的对比（$80\times 8$ 网格，640 个 `Quad4`，1458 个自由度）

| 量 | 解析解 | 有限元解 | 相对误差 |
|---|---|---|---|
| 端部挠度 $v(L, 0)$ | $2.070000$ | $2.058744$ | $0.54\%$ |
| 端部横向位移 $u(L, 0)$ | $0$ | $5.07\times10^{-13}$ | — |
| $L^2$ 位移误差 $\|u_h-u\|_{L^2}$ | — | $1.553\times10^{-2}$ | — |
| 能量范数误差 $\|u_h-u\|_E$ | — | $1.064\times10^{-1}$ | — |

端部横向位移 $u(L,0)$ 在机器精度内为零（理论值恰为 $0$），挠度相对误差约 $0.5\%$，量级与网格分辨率相符。

### 网格收敛性

精确位移场含三次项，双线性 `Quad4` 无法精确表示，故误差非机器零，但随网格加密稳定下降（同一几何文件调整 `nL`、`nH`）：

| 网格 | 单元数 | $\|u_h-u\|_{L^2}$ | $\|u_h-u\|_E$ | 端部挠度 $v(L,0)$ |
|---|---|---|---|---|
| $40 \times 4$ | 160 | $6.14\times10^{-2}$ | $2.11\times10^{-1}$ | 2.0255 |
| $80 \times 8$ | 640 | $1.55\times10^{-2}$ | $1.06\times10^{-1}$ | 2.0587 |
| $160 \times 16$ | 2560 | $3.81\times10^{-3}$ | $5.33\times10^{-2}$ | 2.0672 |
| $240 \times 24$ | 5760 | $1.63\times10^{-3}$ | $3.56\times10^{-2}$ | 2.0688 |
| 精确解 | — | 0 | 0 | 2.0700 |

$L^2$ 位移误差约按 $O(h^2)$ 下降，端部挠度单调收敛到解析值 $2.07$。这张表同时验证了 `stiffness`、`displacement_penalty`、`traction` 三个组装算子与 `l2_error`、`energy_error` 两个误差算子的正确性。

### 应力分布

由第 3 节的解析应力场，$\sigma_{xx} = -P(L-x)y/I$ 在固定端上、下表面取极值 $\pm P L H/(2I)$，$\sigma_{xy} = \tfrac{P}{2I}(H^2/4 - y^2)$ 仅在端部随 $y$ 抛物变化、与 $x$ 无关。有限元解的 von Mises 应力云图（真实变形，`scale = 1`）如下：

![悬臂梁变形 von Mises 应力云图](./cantilever_stress.png)

固定端应力最高（红），自由端最低（蓝）；中性轴附近 $\sigma_{xx}\to 0$，von Mises 应力由剪力 $\sigma_{xy}$ 主导而偏低，与解析应力分布一致。

> 上述数值均可由「悬臂梁算例（cantilever）」中的 `cargo run --example cantilever` 复现；该示例还会输出端部挠度、两种误差，并渲染上图的 `target/cantilever_stress.png`。
