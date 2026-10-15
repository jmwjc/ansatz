# 热传导方程（Heat Equation）

热传导（扩散）方程是继拉压杆之后的第二个完整物理问题。它与拉压杆共享同一套有限元骨架（强形式 → 弱形式 → 伽辽金离散 → 单元组装 → 边界条件），但新增了两个关键维度：**时间导数项**（瞬态问题）与**更丰富的边界条件**（Dirichlet / Neumann / Robin / Cauchy / 混合等）。本章给出从强形式到离散求解的完整理论，并在最后与 `rustfem` 的抽象逐项对应。

> 说明：当前库尚未实现热算子（`PLAN.md` M3 的"2D 泊松问题"即稳态热传导的等价形式）。本章给出理论框架与抽象映射，为后续实现铺路。

## 1. 强形式

### 守恒形式

温度场 $T(x, t)$ 满足能量守恒：

$$
\rho c_p \frac{\partial T}{\partial t} = \nabla \cdot (k \nabla T) + Q
$$

- $\rho$：密度；
- $c_p$：比热容；
- $k$：热导率（可随位置变化，各向异性时为张量 $k_{ij}$）；
- $Q$：单位体积的发热率（体积热源）。

物理含义：单位体积内能的变化率 = 净流入的热传导通量散度 + 内部产热。式中 $q = -k\nabla T$ 是**傅里叶定律**——热从高温区流向低温区。

### 常系数简化

当 $\rho, c_p, k$ 为常数，两边除以 $\rho c_p$ 得

$$
\frac{\partial T}{\partial t} = \alpha \nabla^2 T + \frac{Q}{\rho c_p}
$$

其中

$$
\alpha = \frac{k}{\rho c_p}
$$

是**热扩散率**（单位 $\mathrm{m^2/s}$），衡量温度场扩散的快慢。

### 稳态与泊松方程

当 $\partial T/\partial t = 0$，退化为稳态热传导（泊松方程）：

$$
-\nabla \cdot (k \nabla T) = Q
$$

它与拉压杆平衡方程 $-\frac{d}{dx}\left(EA\frac{du}{dx}\right) = f$ **结构完全相同**：温度 $T$ 对应位移 $u$，热导率 $k$ 对应轴向刚度 $EA$，热源 $Q$ 对应体力 $f$。这一对应关系意味着热算子可以复用桁架算子的全部骨架。

## 2. 边界条件

热方程对时间是**一阶**（初值问题）、对空间是**二阶**（边值问题），因此求解需要**一个初始条件**，以及边界每点**一个边界条件**。按类型，边界条件分三类经典条件，外加若干组合与特殊情形。

### 2.1 初始条件（时间）

给定初始温度分布：

$$
T(x, 0) = T_0(x), \qquad x \in \Omega
$$

时间一阶只需一个初始条件；空间二阶要求边界每点各有一个边界条件，二者共同构成定解问题。

### 2.2 第一类：Dirichlet（温度给定）

在边界 $\Gamma_D$ 上给定**温度**：

$$
T(x, t) = T_D(x, t), \qquad x \in \Gamma_D
$$

这是**本质边界条件**：直接约束试探空间，离散后必须显式强加（罚函数法或消元法），否则刚度矩阵奇异。物理上对应"边界与恒温热源接触"。

### 2.3 第二类：Neumann（热通量给定）

在边界 $\Gamma_N$ 上给定**热通量**（沿外法向流出为正）：

$$
-k \frac{\partial T}{\partial n} = q_n(x, t), \qquad x \in \Gamma_N
$$

其中 $\partial T/\partial n = \nabla T \cdot n$ 是外法向方向导数。这是**自然边界条件**：在分部积分时自动进入载荷项，无需显式处理。

- $q_n > 0$：热量从域内向外流出（散热）；
- $q_n = 0$：**绝热边界**（齐次 Neumann），常见于对称面或理想隔热面；
- 由 $q = -k\nabla T$，通量向外为正对应 $\partial T/\partial n < 0$。

### 2.4 第三类：Robin（对流换热 / 牛顿冷却）

在边界 $\Gamma_R$ 上给定**温度与通量的线性组合**：

$$
-k \frac{\partial T}{\partial n} = h\big(T - T_\infty\big), \qquad x \in \Gamma_R
$$

- $h$：对流换热系数（$\mathrm{W/(m^2\cdot K)}$）；
- $T_\infty$：环境（流体）温度。

这是**牛顿冷却定律**：边界与环境的换热率正比于温差。它既非纯本质也非纯自然——离散后 $h$ 项**同时进入刚度矩阵与载荷向量**（见弱形式节）。Robin 条件是 Dirichlet（$h \to \infty$）与绝热（$h \to 0$）之间的中间情形。

### 2.5 Cauchy：Dirichlet 与 Neumann 同点给定

Cauchy 边界条件指在同一边界点上**同时给定函数值与其法向导数**：

$$
T = T_D, \qquad -k\frac{\partial T}{\partial n} = q_n \quad \text{（同一 } x \in \Gamma_C \text{）}
$$

它主要出现在常微分方程初值问题的语境；对椭圆/抛物型偏微分方程，信息通常过约束、**不适定**，仅在少数反问题或耦合模型中采用。此处列出以保证"所有边界条件"的完整性，实际热传导问题很少使用。

### 2.6 混合边界条件

边界 $\partial\Omega$ 按物理段拆分为互不相交的并：

$$
\partial\Omega = \Gamma_D \cup \Gamma_N \cup \Gamma_R
$$

各段施加不同类型条件。实际工程问题几乎都是混合的：例如一块散热器，底面定温（Dirichlet）、侧面绝热（Neumann）、顶面对流（Robin）。

### 2.7 其他情形（简述）

- **周期性边界条件**：$T|_{\Gamma_1} = T|_{\Gamma_2}$ 且通量连续，用于代表单元法（RVE）建模；
- **对称边界条件**：对称面上 $\partial T/\partial n = 0$，等价于绝热（齐次 Neumann），用于取半域建模；
- **辐射边界（非线性）**：$q = \varepsilon\sigma(T^4 - T_\infty^4)$（Stefan–Boltzmann 定律），是 $T$ 的四次非线性项，需 Newton 迭代求解，本章不展开。

## 3. 弱形式（空间半离散）

### 加权残差与分部积分

用权函数 $v$ 乘守恒方程并在 $\Omega$ 上积分，对扩散项做分部积分（格林第一恒等式）：

$$
\int_\Omega v\,\nabla\cdot(k\nabla T)\,d\Omega
= \int_{\partial\Omega} v\,k\frac{\partial T}{\partial n}\,d\Gamma
- \int_\Omega k\,\nabla v \cdot \nabla T\,d\Omega
$$

代入并整理：

$$
\int_\Omega \rho c_p\, v \frac{\partial T}{\partial t}\,d\Omega
+ \int_\Omega k\,\nabla v\cdot\nabla T\,d\Omega
= \int_\Omega Q\, v\,d\Omega
+ \int_{\partial\Omega} v\, k\frac{\partial T}{\partial n}\,d\Gamma
$$

### 各边界条件的落点

用边界条件替换边界项 $k\frac{\partial T}{\partial n}$：

- 在 $\Gamma_D$ 上：温度给定，令 $v = 0$（本质约束），边界项贡献为零；
- 在 $\Gamma_N$ 上：$k\frac{\partial T}{\partial n} = -q_n$，贡献 $- \int_{\Gamma_N} v\, q_n\, d\Gamma$（自然，进载荷）；
- 在 $\Gamma_R$ 上：$k\frac{\partial T}{\partial n} = -h(T - T_\infty)$，贡献 $- \int_{\Gamma_R} h\,v\,T\,d\Gamma + \int_{\Gamma_R} h\,v\,T_\infty\,d\Gamma$（**一半进刚度、一半进载荷**）。

代入得完整弱形式：

$$
\int_\Omega \rho c_p\,v\,\dot T\,d\Omega
+ \int_\Omega k\,\nabla v\cdot\nabla T\,d\Omega
+ \int_{\Gamma_R} h\,v\,T\,d\Gamma
= \int_\Omega Q\,v\,d\Omega
- \int_{\Gamma_N} v\,q_n\,d\Gamma
+ \int_{\Gamma_R} h\,v\,T_\infty\,d\Gamma
$$

其中 $\dot T = \partial T/\partial t$。

### 半离散形式

空间离散 $T_h = \sum_j T_j(t) N_j$，取 $v = N_i$，得到**半离散**（时间仍连续）的常微分方程组：

$$
M \dot{\mathbf T} + (K + K_R)\,\mathbf T = \mathbf F
$$

各项定义：

- **热容（质量）矩阵** $M_{ij} = \int_\Omega \rho c_p\, N_i N_j\, d\Omega$；
- **热传导（刚度）矩阵** $K_{ij} = \int_\Omega k\,\nabla N_i \cdot \nabla N_j\, d\Omega$；
- **对流（Robin）附加刚度** $K_{R,ij} = \int_{\Gamma_R} h\, N_i N_j\, d\Gamma$；
- **载荷向量** $F_i = \int_\Omega Q\,N_i\,d\Omega - \int_{\Gamma_N} q_n\,N_i\,d\Gamma + \int_{\Gamma_R} h\,T_\infty\,N_i\,d\Gamma$。

## 4. 时间离散：θ 方法

把时间区间离散为步长 $\Delta t$，记 $\bar K = K + K_R$。对 $M\dot{\mathbf T} = \mathbf F - \bar K\mathbf T$ 应用单步 **θ 方法**：

$$
M\frac{\mathbf T^{n+1} - \mathbf T^n}{\Delta t}
= \theta\big(\mathbf F^{n+1} - \bar K\mathbf T^{n+1}\big)
+ (1-\theta)\big(\mathbf F^n - \bar K\mathbf T^n\big)
$$

整理为每步的线性方程组：

$$
\left(\frac{M}{\Delta t} + \theta \bar K\right)\mathbf T^{n+1}
= \left(\frac{M}{\Delta t} - (1-\theta)\bar K\right)\mathbf T^n
+ \theta\mathbf F^{n+1} + (1-\theta)\mathbf F^n
$$

三种特例：

| θ | 名称 | 精度 | 稳定性 |
|---|---|---|---|
| 0 | 前向欧拉（显式） | 一阶 | 条件稳定，$\Delta t \le \Delta t_{\text{cr}}$ |
| 1/2 | Crank–Nicolson | 二阶 | 无条件稳定，但大 $\Delta t$ 会振荡 |
| 1 | 后向欧拉（隐式） | 一阶 | 无条件稳定 |

- **隐式**（$\theta = 1/2$ 或 $1$）每步需求解线性方程组；左端矩阵 $M/\Delta t + \theta\bar K$ 对称正定，可用稀疏 Cholesky——正是 `rustfem` 当前的求解路径；
- **显式**（$\theta = 0$）无需解方程组，但稳定步长受 $\Delta t_{\text{cr}} = O(h^2/\alpha)$ 限制（网格越密，允许步长越小）。

### 质量矩阵与集中化

1D `Seg2` 上的一致（协调）质量矩阵与**集中质量**矩阵：

$$
M^e_{\text{consistent}} = \frac{\rho c_p L^e}{6}
\begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}, \qquad
M^e_{\text{lumped}} = \frac{\rho c_p L^e}{2}
\begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}
$$

集中化把各行元素加到对角上，使 $M$ 成为对角阵，是显式格式的常见配套；一致质量与隐式格式配合精度更好。

## 5. 稳态热传导：与拉压杆的逐项对应

稳态（$\dot{\mathbf T} = 0$）退化为纯空间问题 $\bar K\mathbf T = \mathbf F$，与拉压杆完全同构：

| 热传导 | 拉压杆 | 角色 |
|---|---|---|
| 温度 $T$ | 轴向位移 $u$ | 主变量 |
| 热导率 $k$ | $EA$ | 材料刚度 |
| 热源 $Q$ | 体力 $f$ | 体载荷 |
| 边界热通量 $q_n$ | 集中力/面力 $P$ | 自然边界 |
| 给定温度 $T_D$ | 给定位移 $\bar u$ | 本质边界 |
| $K^e = \dfrac{k}{L^e}\begin{bmatrix}1&-1\\-1&1\end{bmatrix}$ | $K^e = \dfrac{EA}{L^e}\begin{bmatrix}1&-1\\-1&1\end{bmatrix}$ | 单元刚度 |

### 一维稳态解析解（对照基准）

取 $k$、$Q$ 为常数，对 $k\,T'' = -Q$ 积分两次：

$$
T(x) = -\frac{Q}{2k}x^2 + C_1 x + C_2
$$

常数 $C_1, C_2$ 由边界条件确定。例如左端定温 $T(0) = T_L$、右端对流 $-k\,T'(L) = h\big(T(L)-T_\infty\big)$，可解出闭式解作为算例验证基准——其结构与「拉压杆的力学模型」中 $u(x) = -\frac{f}{2EA}x^2 + \cdots$ 完全一致。

## 6. 与 rustfem 抽象的对应

热传导完全复用现有抽象（`Element`、`Variable`、`IntegrationScheme`），只需新增一组热算子。参照 `src/operation/bar.rs` 的四个函数：

| 现有（bar.rs） | 热算子（规划） | 物理含义 |
|---|---|---|
| `stiffness(v, ea)` | `conduction_stiffness(v, k)` | $K^e_{ij} = \sum_g k\,\nabla N_i\cdot\nabla N_j\,|J|\,w_g$ |
| `body_force(v, f)` | `heat_source(v, Q)` | $F_i = \sum_g Q\,N_i\,|J|\,w_g$ |
| `traction(v, p)` | `heat_flux(v, q_n)` | $F_i = -\int_{\Gamma_N} N_i\, q_n\, d\Gamma$（注意符号，沿边界积分） |
| `displacement_penalty(v, α, ū)` | `temperature_penalty(v, α, T_D)` | 罚函数强加定温（沿边界积分） |
| —（新增） | `convection(v, h, T_∞)` | Robin 项：刚度 $h\int N_i N_j\,d\Gamma$ + 载荷 $h T_\infty\int N_i\,d\Gamma$（沿边界积分） |
| —（新增） | `capacity(v, ρc_p)` | 质量矩阵 $M_{ij} = \sum_g \rho c_p N_i N_j\,|J|\,w_g$ |

要点：

- `conduction_stiffness` 与 `stiffness` 仅差一个系数（$k$ 代替 $ea$）：二者都计算物理梯度 $\nabla_x N$ 的点积 $\nabla N_i \cdot \nabla N_j$（`derivative_shape` 现返回关于物理坐标的完整梯度），线单元上 $\nabla_x N = \pm \hat t/L$ 沿切向，退化为 $dN_i/dx \cdot dN_j/dx$；
- `heat_flux` 与 `traction` 的差别是**符号**：Neumann 通量 $q_n$ 进入载荷时带负号（$-\int N_i q_n\,d\Gamma$），源自约定 $k\frac{\partial T}{\partial n} = -q_n$；二者均作用于任意维度的边界单元（`Poi1` 点单元退化为集中通量）；
- 温度 Dirichlet 条件用**罚函数法**，与 `displacement_penalty` 完全同构：对任意维度边界单元生成 $K_{ij} \mathrel{+}= \alpha\int N_i N_j\,d\Gamma$、$F_i \mathrel{+}= \alpha\int N_i T_D\,d\Gamma$（`Poi1` 点单元退化为 $K_{ii}\mathrel{+}=\alpha$、$F_i\mathrel{+}=\alpha T_D$；见「边界条件与罚函数法」章）。

### 实现要点

- **稳态**热传导只需 `conduction_stiffness` + `heat_source` + `temperature_penalty`（+ 可选的 `heat_flux` / `convection`），组装后走稀疏 Cholesky，与 `bar1d` 同一条求解链路；
- **瞬态**需额外 `capacity` 与时间循环，每步解 $\left(\frac{M}{\Delta t} + \theta\bar K\right)\mathbf T^{n+1} = \cdots$，隐式格式仍保持对称正定；
- 2D/3D 上 `Tri3`/`Quad4` 的 $\nabla N_i$ 由 `derivative_shape` 直接返回（已结合雅可比矩阵求物理梯度，见「雅可比与坐标映射」章），单元刚度即 $B^\mathsf T k B$ 形式；`conduction_stiffness`、`capacity`、`heat_source` 均已泛化到任意空间维度 `D`，对任意维度单元物理正确。
