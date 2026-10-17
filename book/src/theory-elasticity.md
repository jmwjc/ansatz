# 线弹性力学（Linear Elasticity）

线弹性力学是结构/连续介质力学的核心物理问题，也是 `rustfem` 在 M3/M4 里程碑要扩展的方向（`Material` trait、平面应力/平面应变、3D 线弹性骨架）。与拉压杆、热传导不同，弹性体中的主变量是**向量场**（位移 $\mathbf u$），应力与应变是**二阶张量**，因此需要先厘清场量之间的基本算子，再走"强形式 → 弱形式 → 伽辽金离散 → 单元组装"的同一条有限元骨架。

> 说明：弹性组装算子已实现于 `src/operation/elasticity.rs`（`stiffness`、`body_force`、`traction`、`displacement_penalty`、`lame_constants`、`plane_stress_lambda`）；`Material` trait、平面应力/应变与 3D 线弹性问题类型仍在路线图上（`PLAN.md` M3/M4）。

## 1. 场量与基础微分算子

### 位移、应变、应力

在域 $\Omega \subset \mathbb{R}^D$（$D = 2$ 或 $3$）上定义三个基本场：

- **位移场** $\mathbf u(\mathbf x)$：向量场，每点有 $D$ 个分量 $u_i(\mathbf x)$；
- **应变张量** $\boldsymbol\varepsilon(\mathbf x)$：二阶对称张量，$\varepsilon_{ij} = \varepsilon_{ji}$，描述局部变形；
- **应力张量** $\boldsymbol\sigma(\mathbf x)$：二阶对称张量，$\sigma_{ij} = \sigma_{ji}$，描述内部单位面积上的内力。

下标 $i, j = 1, \dots, D$ 遍历空间方向，重复下标按爱因斯坦求和约定求和。

### 梯度与对称梯度（应变–位移算子）

位移的**梯度** $\nabla \mathbf u$ 是二阶张量，分量 $\partial u_i/\partial x_j$。它可分解为对称部分与反对称部分：

$$
\nabla \mathbf u = \underbrace{\tfrac12\big(\nabla \mathbf u + (\nabla \mathbf u)^\mathsf T\big)}_{\boldsymbol\varepsilon} + \underbrace{\tfrac12\big(\nabla \mathbf u - (\nabla \mathbf u)^\mathsf T\big)}_{\boldsymbol\omega}
$$

- **对称部分**即小变形下的**应变张量**（几何方程），也是本章的**应变–位移算子**：

$$
\varepsilon_{ij}(\mathbf u) = \frac12 \left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i} \right)
$$

- **反对称部分** $\boldsymbol\omega$ 是**旋转张量**（刚体转动），在线弹性（应力只依赖应变）中不产生应力，故本构关系中不出现。

工程上常把该算子记为 $\varepsilon = \mathcal{S}\mathbf u$（$\mathcal{S}$ 即对称梯度算子 $\operatorname{sym}\nabla$），离散化后就是"应变–位移矩阵" $B$（见第 4 节）。

### 散度与平衡算子

应力张量的**散度**是向量，分量 $\partial \sigma_{ij}/\partial x_j$（对第二个下标求和）。平衡方程要求散度（内力净流出）+ 体力 = 0：

$$
\nabla \cdot \boldsymbol\sigma + \mathbf f = 0
\qquad\Longleftrightarrow\qquad
\frac{\partial \sigma_{ij}}{\partial x_j} + f_i = 0
$$

$(\nabla \cdot \boldsymbol\sigma)_i = \partial_j \sigma_{ij}$ 是**平衡算子**，它与对称梯度算子互为（分部积分意义下的）伴随，这正是弱形式能"把导数从 $\mathbf u$ 转移到权函数 $\mathbf v$"的算子根源。

### 迹、球量与偏量分解

对称张量可分解为**球量（体积变化）**与**偏量（形状变化）**两部分：

$$
\boldsymbol\varepsilon = \underbrace{\tfrac{1}{D}(\operatorname{tr}\boldsymbol\varepsilon)\, \mathbf I}_{\text{球量}} + \underbrace{\left(\boldsymbol\varepsilon - \tfrac{1}{D}(\operatorname{tr}\boldsymbol\varepsilon)\,\mathbf I\right)}_{\text{偏量 }\boldsymbol\varepsilon'}
$$

其中 $\operatorname{tr}\boldsymbol\varepsilon = \varepsilon_{kk}$ 是**体积应变**（体膨胀率），$\mathbf I$ 是单位张量（$\delta_{ij}$）。应力的球量部分 $p = -\sigma_{kk}/D$ 是**静水压**（负号约定压力为正），偏量部分 $\boldsymbol\sigma'$ 是**剪应力**。体积变化与形状变化的解耦在线弹性本构中体现为两个独立的材料常数（见下）。

## 2. 三大类方程（弹性力学方程组）

弹性力学方程组由三类方程构成：**平衡（动量）方程**、**几何（运动学）方程**、**本构方程**，外加边界条件。它们分别来自牛顿定律、变形定义、材料性质。

### 平衡方程（静力）

忽略惯性项，域内每点满足力平衡：

$$
\frac{\partial \sigma_{ij}}{\partial x_j} + f_i = 0
\qquad\Longleftrightarrow\qquad
\nabla \cdot \boldsymbol\sigma + \mathbf f = 0
$$

其中 $\mathbf f$ 是单位体积体力（如重力 $\mathbf f = \rho \mathbf g$）。这是 $D$ 个方程。

### 几何（运动学）方程

小变形假设下，应变是位移梯度的对称部分：

$$
\varepsilon_{ij} = \frac12 \left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i} \right)
$$

它保证应变场由连续位移场导出（自动满足相容性条件）。这是 $D(D+1)/2$ 个独立分量。

### 本构方程（广义胡克定律）

线弹性假设下应力与应变成线性关系（**广义胡克定律**）：

$$
\sigma_{ij} = C_{ijkl}\, \varepsilon_{kl}
$$

$C_{ijkl}$ 是四阶**弹性张量**，含 $D^4$ 个分量，但由对称性（$\sigma$、$\varepsilon$ 对称，且存在应变能密度）约束后独立常数大幅减少。对**各向同性**材料只剩两个独立常数：

$$
\sigma_{ij} = \lambda\, \varepsilon_{kk}\, \delta_{ij} + 2\mu\, \varepsilon_{ij}
$$

其中 $\lambda$、$\mu$ 是**拉梅常数**。$\mu$ 即**剪切模量** $G$，对应偏量（形状变化）；组合 $K = \lambda + \tfrac{2}{3}\mu$ 是**体积模量**，对应球量（体积变化）。

### 拉梅常数与工程常数

工程上更常用杨氏模量 $E$ 与泊松比 $\nu$，与拉梅常数的换算为：

$$
\lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}, \qquad
\mu = G = \frac{E}{2(1+\nu)}
$$

物理含义：$E$ 是单轴拉伸刚度，$\nu$ 是横向收缩比（$0 \le \nu < 0.5$，不可压缩材料 $\nu \to 0.5$）；$\mu$ 是剪切刚度，$K$ 是体积压缩刚度：

$$
K = \lambda + \frac{2}{3}\mu = \frac{E}{3(1-2\nu)}
$$

可见 $\nu \to 0.5$ 时 $K \to \infty$（不可压缩），这会带来有限元中的**体积锁定**问题，是减缩积分/混合格式的动机之一。

### 位移法：Navier–Cauchy 方程

把几何方程代入本构方程再代入平衡方程，消去应力与应变，得到只含位移的强形式（**Navier–Cauchy 方程**）：

$$
\mu \nabla^2 \mathbf u + (\lambda + \mu)\, \nabla(\nabla \cdot \mathbf u) + \mathbf f = 0
$$

分量形式：

$$
\mu \frac{\partial^2 u_i}{\partial x_j \partial x_j} + (\lambda + \mu) \frac{\partial^2 u_j}{\partial x_i \partial x_j} + f_i = 0
$$

这是**位移法**求解的出发点，也是有限元弱形式背后要离散的偏微分算子。

## 3. 边界条件

边界 $\partial\Omega = \Gamma_u \cup \Gamma_t$ 分为互不相交的两部分：

- **位移（Dirichlet）边界** $\Gamma_u$：给定位移 $\mathbf u = \bar{\mathbf u}$（本质边界条件）；
- **面力（Neumann）边界** $\Gamma_t$：给定面力 $\boldsymbol\sigma \cdot \mathbf n = \bar{\mathbf t}$（自然边界条件），其中 $\mathbf n$ 是外法向，分量 $\sigma_{ij} n_j = \bar t_i$。

与拉压杆、热传导完全对应：位移边界需显式强加（罚函数法），面力边界在分部积分时自动进入载荷项。

## 4. Voigt 记号与矩阵算子

因 $\boldsymbol\varepsilon$、$\boldsymbol\sigma$ 对称，独立分量只有 $D(D+1)/2$ 个，工程上把它们排成列向量（**Voigt 记号**），把四阶张量 $C$ 展成矩阵 $D$。

### 向量化

2D（平面）取 $\boldsymbol\varepsilon = [\varepsilon_{xx},\ \varepsilon_{yy},\ \gamma_{xy}]^\mathsf T$，其中 $\gamma_{xy} = 2\varepsilon_{xy}$ 是**工程剪应变**；3D 取

$$
\boldsymbol\varepsilon = [\varepsilon_{xx},\ \varepsilon_{yy},\ \varepsilon_{zz},\ \gamma_{xy},\ \gamma_{yz},\ \gamma_{zx}]^\mathsf T, \qquad
\boldsymbol\sigma = [\sigma_{xx},\ \sigma_{yy},\ \sigma_{zz},\ \sigma_{xy},\ \sigma_{yz},\ \sigma_{zx}]^\mathsf T
$$

采用工程剪应变后，本构关系 $\boldsymbol\sigma = \mathbf D \boldsymbol\varepsilon$ 中的 $\mathbf D$ 保持对称，且切应力–剪应变系数恰为 $\mu$（无需 2 倍修正）。

### B 矩阵（应变–位移矩阵）

位移插值 $\mathbf u = \sum_a N_a \mathbf u_a$ 代入几何方程，应变与节点位移线性相关：

$$
\boldsymbol\varepsilon = \sum_a B_a\, \mathbf u_a
$$

2D（平面）下，节点 $a$ 的 $B_a$（$3 \times 2$）为

$$
B_a =
\begin{bmatrix}
\partial N_a/\partial x & 0 \\
0 & \partial N_a/\partial y \\
\partial N_a/\partial y & \partial N_a/\partial x
\end{bmatrix}
$$

3D 下，节点 $a$ 的 $B_a$（$6 \times 3$）为

$$
B_a =
\begin{bmatrix}
\partial N_a/\partial x & 0 & 0 \\
0 & \partial N_a/\partial y & 0 \\
0 & 0 & \partial N_a/\partial z \\
\partial N_a/\partial y & \partial N_a/\partial x & 0 \\
0 & \partial N_a/\partial z & \partial N_a/\partial y \\
\partial N_a/\partial z & 0 & \partial N_a/\partial x
\end{bmatrix}
$$

$B$ 矩阵即离散化后的**应变–位移算子** $\mathcal{S}$，其元素是形函数的物理梯度分量（见「雅可比与坐标映射」章）。

### D 矩阵（弹性矩阵）

各向同性 3D 弹性矩阵（按上述 Voigt 排列）为

$$
\mathbf D =
\begin{bmatrix}
\lambda + 2\mu & \lambda & \lambda & 0 & 0 & 0 \\
\lambda & \lambda + 2\mu & \lambda & 0 & 0 & 0 \\
\lambda & \lambda & \lambda + 2\mu & 0 & 0 & 0 \\
0 & 0 & 0 & \mu & 0 & 0 \\
0 & 0 & 0 & 0 & \mu & 0 \\
0 & 0 & 0 & 0 & 0 & \mu
\end{bmatrix}
$$

## 5. 平面应力与平面应变

实际工程中，许多问题可降维为平面问题，分两种情况，区别在于哪一个面外分量为零。

### 平面应力

适用于**薄板**面内受力：$\sigma_{zz} = \sigma_{xz} = \sigma_{yz} = 0$（面外无应力），而 $\varepsilon_{zz} \ne 0$（板自由变薄）。由 $\sigma_{zz} = 0$ 解出 $\varepsilon_{zz}$ 代入后得到

$$
\mathbf D_{\text{plane stress}} = \frac{E}{1-\nu^2}
\begin{bmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & (1-\nu)/2
\end{bmatrix}
$$

### 平面应变

适用于**长柱体**（轴向无限长、荷载沿轴向不变）：$\varepsilon_{zz} = \varepsilon_{xz} = \varepsilon_{yz} = 0$（面外无应变），而 $\sigma_{zz} \ne 0$（约束产生的面外应力）。此时

$$
\mathbf D_{\text{plane strain}} = \frac{E}{(1+\nu)(1-2\nu)}
\begin{bmatrix}
1-\nu & \nu & 0 \\
\nu & 1-\nu & 0 \\
0 & 0 & (1-2\nu)/2
\end{bmatrix}
= \begin{bmatrix}
\lambda + 2\mu & \lambda & 0 \\
\lambda & \lambda + 2\mu & 0 \\
0 & 0 & \mu
\end{bmatrix}
$$

### 对比

| | 平面应力 | 平面应变 |
|---|---|---|
| 适用对象 | 薄板 | 长柱体 / 坝体 |
| 面外假设 | $\sigma_{zz} = 0$ | $\varepsilon_{zz} = 0$ |
| 面外应变 $\varepsilon_{zz}$ | 非零（自由收缩） | 零 |
| 面外应力 $\sigma_{zz}$ | 零 | $\nu(\sigma_{xx}+\sigma_{yy})$ |
| $\mathbf D$ 中 $(1,1)$ 元 | $E/(1-\nu^2)$ | $E(1-\nu)/[(1+\nu)(1-2\nu)]$ |

两者的 $\mathbf D$ 都是 $3\times3$ 对称正定矩阵，因此平面应力/平面应变只需一套 B 矩阵组装代码、一个 $3\times3$ 的材料矩阵，差异仅体现在 $\mathbf D$ 的取值——这正是 `Material` trait 要抽象的点。

## 6. 最小势能原理与弱形式

### 总势能

线弹性问题的总势能是应变能减去外力功：

$$
\Pi(\mathbf u) = \frac12 \int_\Omega \boldsymbol\varepsilon(\mathbf u) : \boldsymbol\sigma\, d\Omega
- \int_\Omega \mathbf f \cdot \mathbf u\, d\Omega
- \int_{\Gamma_t} \bar{\mathbf t} \cdot \mathbf u\, d\Gamma
$$

其中 $\boldsymbol\varepsilon : \boldsymbol\sigma = \varepsilon_{ij}\sigma_{ij}$ 是双点积（分量全收缩）。把 $\boldsymbol\sigma = \mathbf C : \boldsymbol\varepsilon$ 代入，应变能密度为 $\tfrac12 \boldsymbol\varepsilon : \mathbf C : \boldsymbol\varepsilon$。

### 变分与弱形式

平衡解使 $\Pi$ 取驻值（对线弹性为最小值）。对 $\Pi$ 取变分，令 $\delta\Pi = 0$，取 $\delta\mathbf u = \mathbf v$（满足 $\mathbf v = 0$ 于 $\Gamma_u$），得**弱形式**：

$$
\int_\Omega \boldsymbol\varepsilon(\mathbf v) : \mathbf C : \boldsymbol\varepsilon(\mathbf u)\, d\Omega
= \int_\Omega \mathbf f \cdot \mathbf v\, d\Omega
+ \int_{\Gamma_t} \bar{\mathbf t} \cdot \mathbf v\, d\Gamma
$$

即

$$
a(\mathbf u, \mathbf v) = l(\mathbf v)
$$

其中

$$
a(\mathbf u, \mathbf v) = \int_\Omega \boldsymbol\varepsilon(\mathbf u) : \mathbf C : \boldsymbol\varepsilon(\mathbf v)\, d\Omega,
\qquad
l(\mathbf v) = \int_\Omega \mathbf f \cdot \mathbf v\, d\Omega + \int_{\Gamma_t} \bar{\mathbf t} \cdot \mathbf v\, d\Gamma
$$

该弱形式即**虚功原理**：内虚功 = 外虚功。注意 $a$ 关于 $\mathbf u,\mathbf v$ 对称（$\mathbf C$ 的对称性保证），且强制（应变能非负），故离散后的刚度矩阵对称正定。

## 7. 有限元离散

### 位移插值

每个单元内用形函数插值位移的每个分量：

$$
\mathbf u^e(\mathbf x) = \sum_{a=1}^{N} N_a(\mathbf x)\, \mathbf u_a^e, \qquad
\mathbf u_a^e = [u_{a,1}, \dots, u_{a,D}]^\mathsf T
$$

每个节点贡献 $D$ 个自由度，单元自由度总数为 $D \cdot N$。应变由 $B_a$ 矩阵组装：

$$
\boldsymbol\varepsilon^e = \sum_a B_a\, \mathbf u_a^e
$$

### 单元刚度

代入双线性形式，单元刚度按节点分块：

$$
K^e_{ab} = \int_{\Omega^e} B_a^\mathsf T\, \mathbf D\, B_b\, d\Omega
\qquad (a,b = 1,\dots,N)
$$

每个块是 $D\times D$ 矩阵。数值积分计算：

$$
K^e_{ab} \approx \sum_{g} B_a(\xi_g)^\mathsf T\, \mathbf D\, B_b(\xi_g)\, |J|\, w_g
$$

这是标量问题（拉压杆 $K^e_{ij} = \int EA\, dN_i/dx\, dN_j/dx$）的直接推广：标量梯度点积换成 $B^\mathsf T \mathbf D B$。

### 张量（指标）形式的单元刚度：`stiffness` 的精确表达式

`operation::elasticity::stiffness` 不拼 Voigt 矩阵，而是直接把各向同性双线性形式展开成分量求和（见第 6 节）：

$$
a(\mathbf u, \mathbf v) = \int_\Omega \left[ \lambda\,(\nabla\cdot\mathbf u)(\nabla\cdot\mathbf v) + \mu\left(\nabla\mathbf u : \nabla\mathbf v + \nabla\mathbf u : (\nabla\mathbf v)^\mathsf T\right) \right] d\Omega
$$

其中 $\nabla\mathbf u : \nabla\mathbf v = \partial_i u_j\,\partial_i v_j$、$\nabla\mathbf u : (\nabla\mathbf v)^\mathsf T = \partial_i u_j\,\partial_j v_i$（对重复指标求和）。该恒等式由 $\boldsymbol\varepsilon(\mathbf u):\mathbf C:\boldsymbol\varepsilon(\mathbf v) = \lambda\,\operatorname{tr}\boldsymbol\varepsilon(\mathbf u)\,\operatorname{tr}\boldsymbol\varepsilon(\mathbf v) + 2\mu\,\boldsymbol\varepsilon(\mathbf u):\boldsymbol\varepsilon(\mathbf v)$ 与 $\operatorname{tr}\boldsymbol\varepsilon = \nabla\cdot\mathbf u$ 推出。

取试探函数 $\mathbf u = N_a\mathbf e_p$、权函数 $\mathbf v = N_b\mathbf e_q$（即节点 $a$ 的第 $p$ 个分量与节点 $b$ 的第 $q$ 个分量），逐项代入得单元刚度条目

$$
K_{(a,p),(b,q)} = \int_\Omega \left[ \lambda\,\partial_p N_a\,\partial_q N_b + \mu\left( \delta_{pq}\,\partial_k N_a\,\partial_k N_b + \partial_q N_a\,\partial_p N_b \right) \right] d\Omega
$$

（记 $\partial_p N_a \equiv \partial N_a/\partial x_p$；$\partial_k N_a\,\partial_k N_b$ 对 $k$ 求和，即 $\nabla N_a\cdot\nabla N_b$。）三个括号项的来源：

| 项 | 连续形式 | $u = N_a\mathbf e_p,\ v = N_b\mathbf e_q$ 时 |
|---|---|---|
| 体积项 | $\lambda(\nabla\cdot u)(\nabla\cdot v)$ | $\lambda\,\partial_p N_a\,\partial_q N_b$ |
| 剪切项 ① | $\mu\,\nabla u:\nabla v$ | $\mu\,\delta_{pq}\,\partial_k N_a\,\partial_k N_b$ |
| 剪切项 ② | $\mu\,\nabla u:(\nabla v)^\mathsf T$ | $\mu\,\partial_q N_a\,\partial_p N_b$ |

数值积分后（`derivative_shape` 返回的物理梯度记为 `dshape`，$|J|$、$w_g$ 见「数值积分」章）：

$$
K_{(a,p),(b,q)} \approx \sum_g \Big[ \lambda\,(\text{dshape}[a][p])(\text{dshape}[b][q]) + \mu\big( \delta_{pq}\,(\text{dshape}[a][k]\,\text{dshape}[b][k]) + (\text{dshape}[a][q])(\text{dshape}[b][p]) \big) \Big] |J|\, w_g
$$

与代码的逐项对应（代码中节点下标为 `i, j`、分量下标为 `a, b`，即数学记号里的 $a,b\to i,j$、$p,q\to a,b$）：

| 数学记号 | `stiffness` 代码 |
|---|---|
| $\partial_k N_a\,\partial_k N_b$ | `grad_dot = (0..D).map(|k| dshape[i][k] * dshape[j][k]).sum()` |
| $\delta_{pq}$ | `delta = if a == b { 1.0 } else { 0.0 }` |
| $\lambda\,\partial_p N_a\,\partial_q N_b + \mu(\delta_{pq}\,\partial_k N_a\,\partial_k N_b + \partial_q N_a\,\partial_p N_b)$ | `value = lambda * dshape[i][a] * dshape[j][b] + mu * (delta * grad_dot + dshape[i][b] * dshape[j][a])` |

自由度编号采用节点优先（node-major）：全局自由度 $= \mathrm{id}[a]\cdot D + p$。因 $K$ 对称，只生成下三角（含对角），条目写入 `(max, min)` 位置后由 faer 累加合并。

该分量形式与 Voigt 形式 $B_a^\mathsf T \mathbf D B_b$ 是**同一个双线性形式**的两种表示，二者都等于 $a(N_a\mathbf e_p, N_b\mathbf e_q)$。平面应力时把 $\lambda$ 换成 $\lambda^* = E\nu/(1-\nu^2)$（见第 5 节），$\mu$ 不变。

### 单元载荷

- **体力**：$F_a = \int_{\Omega^e} N_a\, \mathbf f\, d\Omega \approx \sum_g N_a(\xi_g)\, \mathbf f\, |J|\, w_g$（$D$ 维向量）；
- **面力**：$F_a = \int_{\Gamma_t^e} N_a\, \bar{\mathbf t}\, d\Gamma$（沿边界单元积分，`Poi1` 点单元退化为集中力）。

### 全局组装与求解

把各单元的 $K^e_{ab}$ 按节点自由度映射到全局自由度上累加，得到全局方程

$$
\mathbf K \mathbf U = \mathbf F
$$

每个节点的 $D$ 个自由度编号连续（或交错），`id()` 返回的全局节点号需展开为 $D$ 个自由度。$\mathbf K$ 对称正定，用稀疏 Cholesky 求解；位移边界条件用罚函数法强加（与拉压杆的 `displacement_penalty` 同构，只是约束向量化）。

## 8. 与 rustfem 抽象的对应

### 物理梯度与 B 矩阵

`Element::derivative_shape` 已返回形函数的**物理梯度** $\nabla_x N_a$（$D$ 维向量，见「形函数与等参单元」「雅可比与坐标映射」章）。B 矩阵无需额外计算雅可比，直接用这些分量填表：

- 2D：$B_a$ 的第 1 行取 $(\partial N_a/\partial x, 0)$，第 2 行取 $(0, \partial N_a/\partial y)$，第 3 行取 $(\partial N_a/\partial y, \partial N_a/\partial x)$；
- 3D：按 4.2 节的 $6\times3$ 模板填。

因此弹性算子的实现只需：拿到每个积分点的 $\nabla_x N_a$、$\mathbf D$、$|J|$、$w_g$，做 $B_a^\mathsf T \mathbf D B_b$ 的矩阵乘法后累加。

### 弹性算子（`src/operation/elasticity.rs`，对照 bar.rs / heat）

| 现有（bar.rs） | 弹性算子（`operation::elasticity`） | 物理含义 |
|---|---|---|
| `stiffness(v, ea)` | `stiffness(v, λ, μ)` | $K_{(a,p),(b,q)} = \sum_g [\lambda\,\partial_p N_a\,\partial_q N_b + \mu(\delta_{pq}\,\partial_k N_a\,\partial_k N_b + \partial_q N_a\,\partial_p N_b)]\, |J|\, w_g$（见上文「张量（指标）形式的单元刚度」小节） |
| `body_force(v, f)` | `body_force(v, f)` | $F_{a,p} = \sum_g N_a\, f_p\, |J|\, w_g$（向量载荷） |
| `traction(v, p)` | `traction(v, t)` | $F_{a,p} = \int_{\Gamma_t} N_a\, \bar t_p\, d\Gamma$（面力，向量） |
| `displacement_penalty(v, α, ū)` | `displacement_penalty(v, α, ū)` | 罚函数强加向量位移（$\delta_{pq}$，每个分量各罚一次） |
| —（新增） | `lame_constants(E, ν)` / `plane_stress_lambda(E, ν)` | 由 $E,\nu$ 生成 $\lambda,\mu$ 或平面应力 $\lambda^*$ |

要点：

- `elasticity::stiffness` 与 `bar::stiffness` 同骨架，仅把标量 $EA$ 换成 $\lambda,\mu$、把标量导数点积换成张量分量求和（见上文「张量（指标）形式的单元刚度」小节）：二者都在同一套 `Variable<D, O, P, N>` + `IntegrationScheme` 骨架上累加 `Triplet`；
- 位移是向量，每个节点的 $D$ 个自由度对应 $D\times D$ 的刚度分块与 $D$ 维载荷分量，全局自由度编号须把节点号映射为 $D$ 个连续编号；
- **张量表示**贯穿始终：$\boldsymbol\varepsilon = \operatorname{sym}\nabla\mathbf u$、$\boldsymbol\sigma = \mathbf C:\boldsymbol\varepsilon$ 用指标记号写出，而 $\nabla_x N_a$ 由 `derivative_shape` 按 $\nabla_x N = j_\alpha g^{\alpha\beta}\partial_\beta N$ 直接给出——这对任意空间维度 $D$（含曲面单元嵌入 3D）自然成立，无需运行时矩阵求逆（见「雅可比与坐标映射」章）。

### 常应变单元与 patch test

`Tri3` 的形函数梯度是常数，故 $B_a$ 也是常数矩阵，单元内应变为常数（**常应变三角形 CST**），恰好能精确表示常应变场——这是 patch test 的理论依据：对常应力/常应变问题，数值解应与解析解在机器精度内一致（`PLAN.md` M3 验收标准）。`Quad4` 则因双线性项存在非零的"沙漏"变形，需结合完全积分/减缩积分对比验证。
