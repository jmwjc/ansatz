# 混合格式与 $H(\mathrm{div})$：Raviart–Thomas（RT0）单元

前面各章的单元（`Seg2`、`Tri3`、`Quad4`）都是 $H^1$ 协调的标量单元；「非协调单元：Crouzeix–Raviart（CR）」章则把自由度从顶点搬到边平均，得到 $H^1$ 非协调的标量单元。本章介绍第三个方向：**Raviart–Thomas 单元**（Raviart & Thomas, 1977）——它是**向量值**单元，自由度不再是任何函数的取值，而是每条边上的**法向通量**，全局空间是 $H(\mathrm{div})$ 而不是 $H^1$。

RT 与 CR 是一对"对偶"的离散：两者的自由度**都挂在边上、每条边一个**（RT0 的通量矩正好是 CR 相容性余项里出现的通量矩），并且**把 RT0 杂交化（hybridization）后消去通量，得到的正是 CR 刚度矩阵**（§6）。因此两章可以共用同一套边编号与组装数据。

> 实现：`ansatz::element::rttri3::RTTri3<D>`（最低阶 RT0 三角形，3 个边通量自由度）；见末节「与 ansatz 抽象的对应」。

## 0. 记号

本章沿用 CR 章的记号，只增加通量场 $\boldsymbol q$：

| 记号 | 含义 | 出现位置 |
|---|---|---|
| $\boldsymbol x_I$、$A_T$、逆时针 | 三角形顶点、面积、定向 | CR 章 §2 |
| $\Gamma_I$、$\ell_I$、$\boldsymbol n_I$、$\boldsymbol t_I$ | 顶点 $I$ 的**对边**、边长、**外法向**、切向 | 本章 |
| $\boldsymbol D_I = \ell_I\boldsymbol n_I$ | 边法向量 | CR 章 §3 |
| $\boldsymbol M_I$、$\lambda_I$ | 边中点、面积（重心）坐标 | CR 章 §2 |
| $m_{e,n} = \frac1{\ell_e}\int_{\Gamma_e} u\,L_n\,ds$ | **试探**（势）矩，$n=0$ 即边平均 | CR 章 §2 |
| $\mu_{e,n} = \frac1{\ell_e}\int_{\Gamma_e} \partial_{\boldsymbol n}u\,L_n\,ds$ | **通量**矩 | 本章 |
| $\Phi_e = \int_{\Gamma_e} \boldsymbol q\cdot\boldsymbol n_e\,ds$ | 边通量（RT 的自由度） | 本章 |

其中 $L_n$ 是边上的 Legendre 多项式（$L_0=1$、$L_1(\xi)=\xi$），$\xi = 2s/\ell-1$ 是边坐标。RT0 的自由度就是 $n=0$ 的通量矩 $\Phi_e=\ell_e\mu_{e,0}$。

## 1. 为什么需要另一个空间：混合格式

把模型问题写成一阶系统，通量记为 $\boldsymbol q := -\nabla u$：

$$
\boldsymbol q + \nabla u = \boldsymbol 0,\qquad
\nabla\cdot\boldsymbol q = f \ \text{in }\Omega,\qquad
u = 0\ \text{on }\Gamma_D,\qquad
\boldsymbol q\cdot\boldsymbol n = \bar t\ \text{on }\Gamma_N .
$$

由此立刻得到两条结论：

1. 第一个方程的弱形式要求 $\boldsymbol q$ 具有**散度**，即
   $$
   \boldsymbol q\in H(\mathrm{div};\Omega)
   :=\bigl\{\boldsymbol v\in L^2(\Omega)^2:\ \nabla\cdot\boldsymbol v\in L^2(\Omega)\bigr\}.
   $$
   这个空间里的向量场**不连续**，但它的**法向迹是单值的**：$H(\mathrm{div})$ 正是"$\boldsymbol v\cdot\boldsymbol n$ 在每条内部边上良定义"的空间。切向分量可以任意跳跃，法向分量必须连续。
2. 因此边上的自然未知量不是任何函数的取值，而是**法向通量**。这正是 §0 中 $\{\mu_{e,n}\}$ 那一族量。

RT 单元的定义就是"**自由度恰好是通量矩**"的单元。不要指望它有节点值。

## 2. 局部空间

$$
\mathcal{RT}_k(T) := \mathcal P_k(T)^2 \oplus \boldsymbol x\,\widetilde{\mathcal P}_k(T),
\qquad
\dim\mathcal{RT}_k(T) = (k+1)(k+3),
$$

其中 $\widetilde{\mathcal P}_k$ 是 $k$ 次**齐次**多项式空间。本章只用到最低阶 $k=0$：

| $k$ | 局部空间 | 每个三角形 $\dim$ | 自由度 |
|---|---|---|---|
| $0$ | $\boldsymbol q(\boldsymbol x) = \boldsymbol a + c\,\boldsymbol x$ | $3$ | 三条边各一个 $\Phi_e$ |
| $1$ | 再加 $\boldsymbol a_1x+\boldsymbol a_2y+c\lvert\boldsymbol x\rvert^2$ | $8$ | 三条边各 $\Phi_e,\Psi_e$（6 个）+ 2 个内部 |

通量自由度有两种等价归一化：

$$
\boxed{\ \Phi_e(\boldsymbol q) := \int_{\Gamma_e}\boldsymbol q\cdot\boldsymbol n_e\,ds
= \ell_e\,\mu_{e,0}(\boldsymbol q)\ }
\qquad
\boxed{\ \Psi_e(\boldsymbol q) := \int_{\Gamma_e}\boldsymbol q\cdot\boldsymbol n_e\,L_1(\xi)\,ds
= \ell_e\,\mu_{e,1}(\boldsymbol q)\ }
$$

> **归一化提醒**：CR 章使用的是**平均** $\mu_{e,n}$（量纲为通量），而经典 RT 定义使用**积分** $\Phi_e,\Psi_e$（量纲为通量 × 长度）。二者相差因子 $\ell_e$，只有积分归一化才使基函数满足 §3.1 的 Kronecker 性质。不要混用：$\Phi_e = \ell_e\mu_{e,0}$。

端点/中点上的取值对 $\mathcal{RT}_0$ 没有意义：边上的法向通量在 RT0 中是**常数**（RT1 才允许线性变化），这正是"一个边一个自由度"的原因。

## 3. 最低阶单元（RT0）的显式构造

### 3.1 基函数与 Kronecker 性质

约定 $\boldsymbol x_I$ 是边 $\Gamma_I$ 的**对顶点**（与 CR 章"对边/对顶点"的编号一致），$\boldsymbol n_I$ 为该边的**外**法向，则边 $I$ 的 RT0 基函数为

$$
\boxed{\ \boldsymbol q_I(\boldsymbol x) = \frac{1}{2A_T}\bigl(\boldsymbol x - \boldsymbol x_I\bigr)\ }
\qquad\Longrightarrow\qquad
\Phi_J(\boldsymbol q_I) = \delta_{IJ}.
$$

**证明**：在 $\Gamma_I$ 上，$\boldsymbol x-\boldsymbol x_I$ 沿法向，长度为高 $h_I$，且 $\tfrac12\ell_I h_I = A_T$，故

$$
\boldsymbol q_I\cdot\boldsymbol n_I
= \frac{1}{2A_T}\cdot\frac{2A_T}{\ell_I} = \frac{1}{\ell_I}
$$

在 $\Gamma_I$ 上是**常数**，于是 $\Phi_I(\boldsymbol q_I) = \ell_I\cdot\frac{1}{\ell_I} = 1$。对 $J\ne I$，边 $\Gamma_J$ 连接 $\boldsymbol x_I$ 与第三个顶点，故 $\boldsymbol M_J-\boldsymbol x_I$ 与 $\Gamma_J$ 相切，$\int_{\Gamma_J}\boldsymbol q_I\cdot\boldsymbol n_J\,ds = \ell_J(\boldsymbol M_J-\boldsymbol x_I)\cdot\boldsymbol n_J = 0$。$\blacksquare$

和 CR 一样，这里起作用的是"边矩"而不是"点值"：$\Phi_J(\boldsymbol q_I)=\delta_{IJ}$ 是**边通量矩**上的 Kronecker 性质。

### 3.2 散度

$$
\boxed{\ \nabla\cdot\boldsymbol q_I = \frac{1}{A_T}\ }
\qquad\Longrightarrow\qquad
\boxed{\ \nabla\cdot\boldsymbol q_h
= \frac{1}{A_T}\sum_I \Phi_I
= \frac{1}{A_T}\sum_I \ell_I\,\mu_{I,0}\ }
$$

第二式与 CR 的离散梯度 $\nabla u^h = \frac{1}{A_T}\sum_I\boldsymbol D_I m_I$ 结构完全平行：**离散散度是零阶通量矩的 $\ell$ 加权和，离散梯度是试探矩的 $\boldsymbol D$ 加权和**。

### 3.3 离散散度定理（精确）

对任意 $\boldsymbol q_h\in\mathcal{RT}_0$ 和任意在 $T$ 上为常数的 $w$：

$$
\int_T \nabla\cdot\boldsymbol q_h\,w\,d\Omega
= \sum_I \Phi_I\,w
= \oint_{\partial T}(\boldsymbol q_h\cdot\boldsymbol n)\,w\,ds .
$$

即**散度定理在单元上精确成立，边自由度就是边通量**。这是 CR 恒等式 $\sum_I\boldsymbol D_I m_I = A_T\nabla u^h$ 在 RT 一侧的对应物。

### 3.4 单元矩阵

| 矩阵 | RT0 | CR（对照） |
|---|---|---|
| 质量 $\int_T\boldsymbol q_I\cdot\boldsymbol q_J$ | **稠密** | $\boldsymbol M = \frac{A_T}{3}\boldsymbol I$（对角） |
| "$\boldsymbol B$" $\int_T\nabla\cdot\boldsymbol q_J\,w$ | $\boldsymbol B=(1,1,1)$，**与几何无关** | $\boldsymbol K = \frac{1}{A_T}\boldsymbol D^{\mathsf T}\boldsymbol D$ |
| 刚度 | 不存在（混合格式没有单元刚度） | $\ell_i\ell_j\boldsymbol n_i\cdot\boldsymbol n_j/A_T$ |

RT0 质量矩阵有闭式

$$
\int_T\boldsymbol q_I\cdot\boldsymbol q_J\,d\Omega
= \frac{1}{48A_T}\Bigl[
\sum_{K=1}^{3}(\boldsymbol x_K-\boldsymbol x_I)\cdot(\boldsymbol x_K-\boldsymbol x_J)
+ 9\,(\boldsymbol c-\boldsymbol x_I)\cdot(\boldsymbol c-\boldsymbol x_J)\Bigr],
$$

其中 $\boldsymbol c$ 是形心（由 $\int_T x_ix_j = \frac{A_T}{12}\bigl(\sum_K x_{K,i}x_{K,j}+9c_ic_j\bigr)$ 得到；与 3 点积分逐项相同）。质量矩阵稠密是 $H(\mathrm{div})$ 协调性的代价——CR 的质量矩阵对角，恰恰因为 CR **不**是 $H(\mathrm{div})$ 协调的。

## 4. 参数坐标与 Piola 映射

### 4.1 三种等价写法

取参考三角形上的面积坐标 $\lambda=(\lambda_1,\lambda_2,\lambda_3)$，$\lambda_1=1-\xi-\eta$，$\lambda_2=\xi$，$\lambda_3=\eta$，$\boldsymbol x=\sum_I\lambda_I\boldsymbol x_I$。RT0 基函数有三种等价表达：

| 形式 | 表达式 | 用途 |
|---|---|---|
| 物理（§3.1） | $\boldsymbol q_I=\dfrac{1}{2A_T}\bigl(\boldsymbol x-\boldsymbol x_I\bigr)$ | 逐点恒等式、证明 |
| **重心坐标** | $\boldsymbol q_I=\dfrac{1}{2A_T}\Bigl[\lambda_J\bigl(\boldsymbol x_J-\boldsymbol x_I\bigr)+\lambda_K\bigl(\boldsymbol x_K-\boldsymbol x_I\bigr)\Bigr]$ | 不依赖内部点、纯边/几何组装 |
| **参考 + Piola** | $\hat{\boldsymbol q}_1=\binom{\xi}{\eta},\ \hat{\boldsymbol q}_2=\binom{\xi-1}{\eta},\ \hat{\boldsymbol q}_3=\binom{\xi}{\eta-1}$，$\boldsymbol q_I=\dfrac{1}{\lvert\det\boldsymbol J\rvert}\,\boldsymbol J\,\hat{\boldsymbol q}_I(\xi,\eta)$ | 数值积分、曲边/一般单元 |

后两种都来自同一个恒等式 $\boldsymbol x-\boldsymbol x_I=\sum_L\lambda_L(\boldsymbol x_L-\boldsymbol x_I)$（$L=I$ 项为零）以及 $\lvert\det\boldsymbol J\rvert = 2A_T$、$\boldsymbol J(\xi,\eta)=\xi(\boldsymbol x_2-\boldsymbol x_1)+\eta(\boldsymbol x_3-\boldsymbol x_1)$。

关键点是这是**逆变（contravariant）Piola 映射**：它保持法向通量，因而边通量矩 $\Phi_e$ 是天然自由度（§2）；同时也保证 $D=3$ 中嵌入的曲/平三角形同样成立（此时 $\lvert\det\boldsymbol J\rvert$ 取 $\sqrt{\det(\boldsymbol J^{\mathsf T}\boldsymbol J)}$，即 `jacobe`）。

### 4.2 重心坐标形式的散度

利用 CR 章的法向记号 $\boldsymbol D_I=\ell_I\boldsymbol n_I$（$\boldsymbol n_I$ 外法向）有

$$
\nabla\lambda_I = -\frac{1}{2A_T}\boldsymbol D_I ,
$$

以及两个恒等式 $\boldsymbol D_J\cdot(\boldsymbol x_J-\boldsymbol x_I)=\boldsymbol D_K\cdot(\boldsymbol x_K-\boldsymbol x_I)=-2A_T$（对顶点到对边直线的距离为 $h_L=2A_T/\ell_L$，取内侧），于是从参数形式直接得到 §3.2：

$$
\nabla\cdot\boldsymbol q_I
= \frac{1}{2A_T}\sum_{L\ne I}\nabla\lambda_L\cdot\bigl(\boldsymbol x_L-\boldsymbol x_I\bigr)
= \frac{1}{2A_T}\cdot\frac{-1}{2A_T}\bigl(-4A_T\bigr)
= \frac{1}{A_T}.
$$

### 4.3 边上的迹

在 $\Gamma_I$ 上 $\lambda_I=0$；若 $s\in[0,\ell_I]$ 从 $\boldsymbol x_J$ 起量，则 $\lambda_J=1-s/\ell_I$、$\lambda_K=s/\ell_I$，于是

$$
\boldsymbol q_I\big|_{\Gamma_I}
= \frac{1}{2A_T}\Bigl[\bigl(1-\tfrac{s}{\ell_I}\bigr)(\boldsymbol x_J-\boldsymbol x_I)
+\tfrac{s}{\ell_I}(\boldsymbol x_K-\boldsymbol x_I)\Bigr],
\qquad
\boldsymbol q_I\cdot\boldsymbol n_I = \frac{1}{\ell_I}\ \text{（常数）},
$$

而 $\boldsymbol q_J\cdot\boldsymbol n_I$、$\boldsymbol q_K\cdot\boldsymbol n_I$ 沿边**仿射**且在**中点**（$\lambda_J=\lambda_K=\tfrac12$）处为零——这正是 §3.1 中 $\Phi_I(\boldsymbol q_I)=\delta_{IJ}$ 的参数解释。

### 4.4 与 Legendre 矩的字典

边上的 Legendre 坐标是 $\xi_{\text{edge}}=2s/\ell_I-1=2\lambda_K-1$，所以

$$
\mu_{I,0}=\frac{1}{\ell_I}\int_{\Gamma_I}\boldsymbol q\cdot\boldsymbol n_I\,ds,
\qquad
\mu_{I,1}=\frac{1}{\ell_I}\int_{\Gamma_I}\boldsymbol q\cdot\boldsymbol n_I\,(2\lambda_K-1)\,ds,
$$

即**边上的 Legendre 基就是"对顶点面积坐标"的多项式基**，RT1 只需把 $\{1,\,2\lambda_K-1\}$ 加进边自由度。

> 实用提示：代码里通常用"参考 + Piola"形式求值，因为积分点与逆映射都在 $(\xi,\eta)$ 上；而"重心坐标"形式适合只由 $\{\boldsymbol x_I\}$、$\{\ell_I\}$、$\{\boldsymbol n_I\}$ 描述几何的边/通量组装（不需要单元内部点）。

## 5. 全局空间与泊松混合格式

$$
\mathcal{RT}_0(\mathcal T_h):=\bigl\{\boldsymbol v_h:\ \boldsymbol v_h|_T\in\mathcal{RT}_0(T),\
\boldsymbol v_h\cdot\boldsymbol n \text{ 在每条内部边上单值}\bigr\},
$$

**每条边一个自由度**——与 CR 空间的数据布局完全相同（CR 用边表编号自由度，RT0 可以原样复用同一张表）。泊松方程的混合格式是

$$
\text{求 }(\boldsymbol q_h,u_h)\in\mathcal{RT}_0\times\mathcal P_0^d:\quad
\begin{cases}
\displaystyle\int_\Omega\boldsymbol q_h\cdot\boldsymbol r\,d\Omega
=\int_\Omega u_h\,\nabla\cdot\boldsymbol r\,d\Omega,
&\forall\boldsymbol r\in\mathcal{RT}_0,\\[2mm]
\displaystyle\int_\Omega\nabla\cdot\boldsymbol q_h\,w\,d\Omega
=\int_\Omega f\,w\,d\Omega,
&\forall w\in\mathcal P_0^d,
\end{cases}
$$

即鞍点系统

$$
\begin{bmatrix}\boldsymbol A & \boldsymbol B^{\mathsf T}\\ \boldsymbol B & \boldsymbol 0\end{bmatrix}
\begin{bmatrix}\boldsymbol q\\ \boldsymbol u\end{bmatrix}
=\begin{bmatrix}\boldsymbol 0\\ \boldsymbol F\end{bmatrix},
\qquad
\boldsymbol A=\text{（RT0 质量矩阵）},\quad
\boldsymbol B=\text{（散度关联矩阵）}.
$$

性质：

1. **局部守恒精确**：取 $w=\mathbb 1_T$，第二式给出 $\int_T\nabla\cdot\boldsymbol q_h=\int_T f$，即 $\nabla\cdot\boldsymbol q_h|_T = (\Pi_{P_0}f)|_T$ **逐单元成立**。散度误差只剩 $f$ 的 $P_0$ 插值误差，没有相容性缺陷。
2. **协调性**：$\Phi_e$ 单值意味着法向通量跨边连续，误差分析中**没有跳跃项**（对比 CR：其全部误差来自势的跳跃）。
3. 标准估计（Brezzi–Babuška 理论，$\boldsymbol q\in H^1$、$f\in H^1$）：
   $$
   \lVert\boldsymbol q-\boldsymbol q_h\rVert_0\le Ch\lvert\boldsymbol q\rvert_1,
   \qquad
   \lVert u-u_h\rVert_0\le Ch\lvert u\rvert_1 .
   $$
4. 系统对称**不定**，求解需 `sp_lu`（与 CR 的对称正定系统不同）。

## 6. 与 CR 的等价：杂交化的 RT0 就是 CR 单元

在 $\partial T$ 上引入拉格朗日乘子 $\lambda$（即**势的迹**）来强制法向通量连续，逐单元写出杂交问题：

$$
\int_T\boldsymbol q_h\cdot\boldsymbol r\,d\Omega-\int_Tu_h\,\nabla\cdot\boldsymbol r\,d\Omega
=-\oint_{\partial T}\lambda_h\,(\boldsymbol r\cdot\boldsymbol n)\,ds,
\qquad
\int_T\nabla\cdot\boldsymbol q_h\,w\,d\Omega=\int_T f w\,d\Omega .
$$

一个三角形的边界上的乘子空间是 $\mathcal P_0(\partial T)$——**逐边常数**。而逐边常数的迹正是 CR 的边平均自由度 $m_e=\frac{1}{\ell_e}\int_{\Gamma_e}u\,ds$。所以**乘子空间就是 CR 空间**。

### 6.1 定理（静态凝聚）

局部消去 $(\boldsymbol q_h,u_h)$（static condensation）。记 $\lambda_e$ 为边常数、$\Phi_e$ 为流出 $T$ 的通量，则凝聚后的算子**正是 CR 单元刚度矩阵差一个负号**：

$$
\boxed{\ \Phi_e = -\sum_j K^{CR}_{ej}\lambda_j,
\qquad
K^{CR}_{ej}=\frac{\ell_e\ell_j}{A_T}\,\boldsymbol n_e\cdot\boldsymbol n_j
=\Bigl(\frac{1}{A_T}\boldsymbol D^{\mathsf T}\boldsymbol D\Bigr)_{ej}\ }
$$

（即 CR 章 §2 的刚度矩阵）。因此对内部边施加通量连续，$\sum_{T\ni e}\Phi_e^{(T)}=0$（再加上 Dirichlet 数据），就得到

$$
\sum_j\Bigl(\sum_{T\ni e}K^{CR,(T)}_{ej}\Bigr)\lambda_j=\text{（边界数据）},
$$

这正是 CR 的整体方程 $\boldsymbol K^{CR}\boldsymbol m=\boldsymbol f$，其中 $\boldsymbol m=\lambda$。于是

$$
\boxed{\ \text{CR（非协调，}H^1\text{）}\;=\;\text{RT0（混合，}H(\mathrm{div})\text{），消去通量后。}\ }
$$

这是二阶椭圆问题"混合格式与非协调格式等价"的三角形情形（Arnold–Brezzi, *Mixed and nonconforming finite element methods*, M²AN 19 (1985)），上面的推导是它在单元层面的证明。ansatz 的 `rttri3_cr_gradient_flux_dofs_are_cr_stiffness_rows` 测试用线性场的后处理恒等式数值验证了同一结论。

## 7. 后处理：把 CR 解变成 RT 通量

若已有 CR 解 $u_h$，最自然的 RT0 通量就是它的逐单元梯度 $\nabla_h u_h$。其通量自由度正是 CR 最小二乘/后处理小节里出现的那些数：

$$
\Phi_e(\nabla_h u_h)
=\ell_e\,\mu_{e,0}(\nabla_h u_h)
=\ell_e\,\boldsymbol n_e\cdot\nabla u_h
=\frac{\ell_e}{A_T}\sum_j\ell_j(\boldsymbol n_e\cdot\boldsymbol n_j)m_j
=\sum_j K^{CR}_{ej}\,m_j ,
$$

即 **CR 梯度的 RT 通量自由度恰好是 CR 刚度矩阵的行**——这是 §6 恒等式反过来读。这解释了为什么两个族可以共用同一套"每条边一个自由度"的组装数据。

## 8. 与 ansatz 抽象的对应

### 8.1 为什么需要拆分 `Element`

`Element<D, P, N>` 的 `shape` 返回 `[f64; N]`——每个自由度一个**标量**。RT0 的基函数是**向量**，一个自由度对应一个 $D$ 维向量。为了让 RT0 也进入单元抽象，`ansatz` 把 `Element` 拆成三层：

- `ElementBase<D, P, N>`：自由度 id、节点坐标与等参几何（`coordinates` / `jacobe` / `jacobe_mat` / `local_coordinates` / `contains`）；
- `Element<D, P, N>: ElementBase`：标量形函数 `shape` / `derivative_shape`（现有全部标量单元，可读别名 `ScalarElement`）；
- `HdivElement<D, P, N>: ElementBase`：向量形函数 `shape` / `divergence_shape`。

`RTTri3` 实现 `ElementBase + HdivElement`，不再是游离在 `element` 外的类型；但它仍不进入标量 `Variable`/刚度算子链路（RT 的混合组装需要通量质量与散度关联，属于后续工作）。

```rust
pub struct RTTri3<const D: usize> { /* nodes: 三个边中点 */ }

impl<const D: usize> ElementBase<D, 2, 3> for RTTri3<D> {
    fn id(&self) -> [usize; 3];                       // 三条边的通量自由度
    fn vertices_coordinates(&self) -> [[f64; D]; 3];  // 三个边中点
    fn coordinates(&self, xi: [f64; 2]) -> [f64; D];
    fn jacobe(&self, xi: [f64; 2]) -> f64;
    fn jacobe_mat(&self, xi: [f64; 2]) -> [[f64; D]; 2];
    fn local_coordinates(&self, x: [f64; D]) -> Option<[f64; 2]>;
}

impl<const D: usize> HdivElement<D, 2, 3> for RTTri3<D> {
    fn shape(&self, xi: [f64; 2]) -> [[f64; D]; 3];   // qᵢ = J q̂ᵢ / jacobe
    fn divergence_shape(&self, xi: [f64; 2]) -> [f64; 3];    // 1/|T|
}

impl<const D: usize> RTTri3<D> {
    // 构造：与 CRTri3 完全相同的"边中点 + 边自由度 id"布局
    pub fn new(m1: Node<D>, m2: Node<D>, m3: Node<D>) -> Self;
    pub fn from_tri3(tri: &Tri3<D>, edge_id: impl FnMut(usize, usize) -> usize) -> Self;
    pub fn try_from_tri3<E>(tri: &Tri3<D>, edge_id: impl FnMut(usize, usize) -> Result<usize, E>)
        -> Result<Self, E>;

    // RT 专属工具
    pub fn triangle_vertices(&self) -> [[f64; D]; 3];
    pub fn area(&self) -> f64;
    pub fn reference_vector_shape(&self, xi: [f64; 2]) -> [[f64; 2]; 3];
    pub fn divergence(&self, coefficients: [f64; 3]) -> f64;
    pub fn mass_matrix(&self) -> [[f64; 3]; 3];
    pub fn edge_normal(&self, edge: usize) -> Option<[f64; D]>;
    pub fn flux_dof(&self, edge: usize, q: impl Fn([f64; D]) -> [f64; D]) -> Option<f64>;
}
```

### 8.2 积分规则

- 质量矩阵被积函数为二次，用 3 点 `Dunavant<3>`（`FullIntegration<2, 3>`）；
- 散度/通量被积函数为常数或仿射，1 点 `Dunavant<1>`（`ReducedIntegration<2, 1>`）即精确；
- `flux_dof` 在边上用 2 点 Gauss（`GaussLegendre1D<2>`），对 $\le 3$ 次多项式精确。

### 8.3 测试

`src/element/rttri3.rs` 的单元测试覆盖了本章的全部结论：

| 测试 | 验证内容 |
|---|---|
| `rttri3_basis_forms_agree` | 物理、重心、参考 + Piola 三种基函数形式逐点一致（$\lesssim10^{-13}$） |
| `rttri3_flux_dofs_are_kronecker_deltas` | $\Phi_J(\boldsymbol q_I)=\delta_{IJ}$ |
| `rttri3_divergence_is_unit_and_the_divergence_theorem_holds` | $\nabla\cdot\boldsymbol q_I=1/A_T$，$\oint\Phi_I=\int_T\nabla\cdot\boldsymbol q_h$ |
| `rttri3_mass_matrix_is_dense_and_matches_the_closed_form` | §3.4 的闭式；矩阵对称正定且**非对角** |
| `rttri3_cr_gradient_flux_dofs_are_cr_stiffness_rows` | §7 的恒等式，即 §6 定理的反向表述 |
| `rttri3_from_tri3_assigns_edge_ids_and_midpoints` | 边编号、边中点、反射顶点、面积 |
| `rttri3_works_for_a_triangle_embedded_in_3d` | $D=3$ 嵌入三角形：Piola、外法向、Kronecker 性质 |
| `elementbase_is_shared_with_scalar_elements` | `ElementBase` 的泛型几何接口同时接受 `Tri3`（标量）与 `RTTri3`（向量值） |

可视化示例 `examples/rt0.rs`（`cargo run --example rt0`）在参考三角形上画出三个基函数：

- `target/rt0_shape_functions.png`：三个基函数的**向量场图**（箭头，颜色表示 $|\boldsymbol q_i|$）；
- `target/rt0_shape_components.png`：每个基函数的 $x$、$y$ **分量图**（热图 + 仿射分量的等值线 + 色条）。

## 9. 小结与开放项

- RT0 是 $H(\mathrm{div})$ 协调的**向量值**单元，局部空间为仿射向量场，**每条边一个自由度**——法向通量 $\Phi_e$；
- 基函数 $\boldsymbol q_I=(\boldsymbol x-\boldsymbol x_I)/(2A_T)$，满足 $\Phi_J(\boldsymbol q_I)=\delta_{IJ}$、$\nabla\cdot\boldsymbol q_I=1/A_T$，散度定理逐单元精确；
- 质量矩阵**稠密**（与 CR 的对角质量矩阵恰成对照），"散度矩阵" $\boldsymbol B=(1,1,1)$ 与几何无关；混合格式局部守恒精确，但系统对称不定；
- **杂交化消去通量后得到 CR 刚度矩阵**：两个族是同一个对象的两面；CR 解的后处理通量自由度就是 CR 刚度矩阵的行；
- 实现上 RT0 与 CR 共享"边中点节点 + 边 id"的布局；RT0 实现 `ElementBase + HdivElement`（向量值），暂不进入现有标量组装算子。

开放项：

1. **整体（非逐单元）形式的 §6 定理**，含 Dirichlet 数据与乘子符号约定的完整推导；
2. **RT1**：每条边加 $\Psi_e=\ell_e\mu_{e,1}$ 加两个内部矩；内部自由度使 RT1 不能再被凝聚成 CR 型势方法；
3. **组装层**：`heat::mixed_flux_mass` / `mixed_divergence` / `mixed_dirichlet` 已提供 RT0/P0 的混合弱形式；通用的向量值 `Variable` 与基于网格的 RT0 组装尚未建立；
4. **边界迹单元**：`RTSeg3` 已提供（父单元三个边通量、`normal_flux_shape`），`heat::mixed_dirichlet` 用它组装 Dirichlet 自然项；尚未接入 GMSH 读取与全局边方向符号。
