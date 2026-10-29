# 热传导方程（Heat Equation）

热传导（扩散）方程是继拉压杆之后的第二个完整物理问题。它与拉压杆共享同一套有限元骨架（强形式 → 弱形式 → 伽辽金离散 → 单元组装 → 边界条件），但新增了两个关键维度：**时间导数项**（瞬态问题）与**更丰富的边界条件**（Dirichlet / Neumann / Robin / Cauchy / 混合等）。本章给出从强形式到离散求解的完整理论，并在最后与 `ansatz` 的抽象逐项对应。

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

这是**本质边界条件**：直接约束试探空间，离散后必须显式强加（罚函数法或消元法），否则刚度矩阵奇异。物理上对应"边界与恒温热源接触"。两种强加方式的完整理论见「边界条件与罚函数法」章与本章 §6（拉格朗日乘子法）。

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

- **隐式**（$\theta = 1/2$ 或 $1$）每步需求解线性方程组；左端矩阵 $M/\Delta t + \theta\bar K$ 对称正定，可用稀疏 Cholesky——正是 `ansatz` 当前的求解路径；
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

## 6. 本质边界条件的强加：拉格朗日乘子法

### 为什么另起一套做法

§2.2 的 Dirichlet 条件 $T = T_D$ 在离散后必须**显式强加**。库当前采用罚函数法（「边界条件与罚函数法」章），以"近似满足 + 需调参 $\alpha$"为代价换取矩阵对称正定。热传导问题还有一个更"数学正确"的选择：**拉格朗日乘子法**——它把本质条件当作**等式约束**用乘子精确施加，并且**顺带给出 Dirichlet 边界上的热流**。热工程中 $q_n$ 往往正是最关心的输出（散热功率、热流密度分布），所以这一优点在热问题上比在弹性力学中更突出。

### 约束问题与 Lagrangian

记

$$
a(T,v) = \int_\Omega k\,\nabla T\cdot\nabla v\,d\Omega + \int_{\Gamma_R} h\,T\,v\,d\Gamma,\qquad
l(v) = \int_\Omega Q\,v\,d\Omega - \int_{\Gamma_N} q_n\,v\,d\Gamma + \int_{\Gamma_R} h\,T_\infty\,v\,d\Gamma
$$

（稳态情形；瞬态时 $a$ 还需加质量项，见"瞬态：半离散 DAE"）。试探空间与检验空间为

$$
V_D = \bigl\{T \in H^1(\Omega) : T = T_D \ \text{on}\ \Gamma_D\bigr\},\qquad
V_0 = \bigl\{v \in H^1(\Omega) : v = 0 \ \text{on}\ \Gamma_D\bigr\},
$$

弱形式（§3）即

$$
\text{求 } T \in V_D:\quad a(T,v) = l(v),\qquad \forall v \in V_0 .
$$

这等价于能量泛函 $\Pi(T) = \tfrac12 a(T,T) - l(T)$ 在**仿射流形** $V_D$ 上取极小。把"$T = T_D$ on $\Gamma_D$"看成等式约束、引入乘子 $\lambda \in \Lambda = L^2(\Gamma_D)$，构造 **Lagrangian**

$$
\mathcal{L}(T,\lambda) = \Pi(T) + \int_{\Gamma_D}\lambda\,\bigl(T - T_D\bigr)\,d\Gamma .
$$

鞍点条件（分别对 $T$ 与 $\lambda$ 取变分）给出

$$
\begin{cases}
 a(T,v) + \displaystyle\int_{\Gamma_D}\lambda\,v\,d\Gamma = l(v), & \forall v \in V,\\[2mm]
 \displaystyle\int_{\Gamma_D}\mu\,\bigl(T - T_D\bigr)\,d\Gamma = 0, & \forall \mu \in \Lambda .
\end{cases}
$$

第一个方程的检验空间从 $V_0$ **放宽**到 $V$：约束被"搬进"方程，$v$ 在 $\Gamma_D$ 上不再要求为零，代价是多出一个未知场 $\lambda$。

### 乘子的物理意义

对照原始弱形式中的边界项 $\int_{\partial\Omega} v\,k\frac{\partial T}{\partial n}d\Gamma$ 与上式多出的 $\int_{\Gamma_D}\lambda\,v\,d\Gamma$，可见

$$
\lambda = -k\frac{\partial T}{\partial n}\bigg|_{\Gamma_D} = q_n\big|_{\Gamma_D} ,
$$

即**乘子就是 Dirichlet 边界上未知的向外热流**，符号与 §2.3 的约定（向外为正）一致。这一点很实用：罚函数法只能用 $\alpha(T_h - T_D)$ 事后近似这个量，而乘子法是它的直接解。

### 离散：鞍点系统

取 $T_h = \sum_j T_j N_j$，乘子用 $\Gamma_D$ 上的**迹空间**离散：

$$
\lambda_h = \sum_{a} \lambda_a \, N_a\big|_{\Gamma_D},
$$

自由度取在 $\Gamma_D$ 的边界节点上（与 $T_h$ 的节点自由度一一对应）。记

$$
G_{aj} = \int_{\Gamma_D} N_a N_j\,d\Gamma,\qquad
g_a = \int_{\Gamma_D} N_a T_D\,d\Gamma ,
$$

$G$ 只在 $\Gamma_D$ 的边界节点间非零，是**边界质量矩阵**：一侧耦合边界节点与相邻内部节点（$G_{aj}$，$j \notin \Gamma_D$），另一侧给出边界节点自身的边界质量矩阵块（$G_{ab}$，$a,b \in \Gamma_D$，对称正定）。$F_j = l(N_j)$。离散鞍点系统为

$$
\begin{bmatrix} K & G^{\mathsf T} \\ G & 0 \end{bmatrix}
\begin{bmatrix} \mathbf T \\ \boldsymbol\lambda \end{bmatrix}
=
\begin{bmatrix} \mathbf F \\ \mathbf g \end{bmatrix},
$$

其中 $K$ 泛指空间算子矩阵（含 Robin 刚度）。按自由节点 $F$ 与约束节点 $D$（$\Gamma_D$ 上的节点）分块，并把已知的 $\mathbf T_D$ 移到右端：

$$
\begin{bmatrix} K_{FF} & G_F^{\mathsf T} \\ G_F & 0 \end{bmatrix}
\begin{bmatrix} \mathbf T_F \\ \boldsymbol\lambda \end{bmatrix}
=
\begin{bmatrix} \mathbf F_F - K_{FD}\mathbf T_D \\ \mathbf g - G_D\mathbf T_D \end{bmatrix}.
$$

要点：

- **无需罚数**：约束不是"近似满足"，因此不存在 $\alpha$ 取小则约束不满足、取大则病态的两难；
- **约束行是积分意义的**：$\int_{\Gamma_D} N_a (T_h - T_D)\,d\Gamma = 0$ 即 $G(\mathbf T - \mathbf T_D) = \mathbf 0$。因 $G_D$ 是边界质量矩阵（可逆），它蕴含 $T_h|_{\Gamma_D} = I_h T_D$ **逐节点精确成立**（$I_h$ 为节点插值）；
- **矩阵对称但不定**：对角块为零、鞍点结构，**不能**用 Cholesky，需带主元的稀疏 $LDL^{\mathsf T}$ / LU（`faer` 的稀疏 LU），这是与罚函数法（保持对称正定、可走 `sp_cholesky`）的根本差别；
- 额外未知数只在 $\Gamma_D$ 上，规模远小于总自由度，通常不是瓶颈。

### 乘子空间的选取与混合离散

以上取乘子为 $T_h$ 的**迹空间**（$\lambda_h$ 与 $T_h$ 共用边界节点），此时 $G$ 是方形的边界质量矩阵、约束逐点成立。乘子空间其实可以不同：只要 $\lambda$ 与 $T$ 使用**同一套积分单元与积分点**（同一几何、同一参数域 $\xi$），形函数算法可以各自独立。常见搭配：

| $\lambda$ 空间 | $G$ 的形状 | 约束的含义 | 备注 |
|---|---|---|---|
| $T_h$ 的迹空间（等阶） | $N\times N$ 方阵，即边界质量矩阵 | 逐节点精确 $T_h = I_h T_D$ | 最常用，$G_D$ 可逆 |
| 每个边界单元一个常值（$P_0$） | 长方形 $N_e \times N$ | 积分平均 $\int_e (T_h - T_D)\,d\Gamma = 0$ | 自由度更少，只"弱"满足；适合复杂几何或要与其它耦合场匹配时 |
| 高于 $T_h$ 的迹空间 | 行数更多（超定） | —— | 通常破坏 inf-sup，不建议 |

要点：混合离散下 $G$ 不再是方阵，$G\mathbf T = \mathbf g$ 的行数（乘子自由度）可以多于或少于约束节点数；可解性由离散 inf-sup 条件保证——等阶与"每个边界单元一个常数"两种常见搭配都满足（$G$ 行满秩），而比 $T_h$ 迹空间更粗或更细的乘子都需要单独检查。实现上只需分别给出 $T$ 与 $\lambda$ 的形函数在共享积分点处的取值，因此 `ansatz` 的 `temperature_lagrangian` 接收**两组边界单元**（配对迭代，几何、积分规则取温度单元，乘子单元只提供形函数与自由度编号），并返回**完整**的耦合块 $G$（$M\times N$ 全部元素）与 $\mathbf g$，因为 $G$ 位于鞍点系统的非对角块，不能按对称三角存储。

### Schur 补与静力凝聚

由第一个块行解出 $\mathbf T_F$ 再代入约束行，得关于乘子的**对偶（Schur 补）系统**：

$$
S\,\boldsymbol\lambda = G_F K_{FF}^{-1}\bigl(\mathbf F_F - K_{FD}\mathbf T_D\bigr) - \bigl(\mathbf g - G_D\mathbf T_D\bigr),
\qquad S = G_F K_{FF}^{-1} G_F^{\mathsf T} .
$$

$S$ 对称半正定，$G_F$ 列满秩时正定，故乘子方程本身可用共轭梯度；但每次矩阵-向量乘都要解一次 $K_{FF}$，因此

- 只需解一次（稳态）：直接解鞍点系统（稀疏 LU）最省事；
- 需反复求解（瞬态、非线性）：预分解 $K_{FF}$，再在乘子空间上迭代（Uzawa / 对偶迭代），可把 $\Gamma_D$ 上的小规模求解与内部大规模求解解耦。

### 瞬态：半离散 DAE

瞬态时约束在**每一时刻**都成立，而质量矩阵在约束方向上为零（温度值被强制给定，不能演化），故半离散系统是**指标 2 的微分代数方程（DAE）**：

$$
\begin{bmatrix} M_{FF} & 0 \\ 0 & 0 \end{bmatrix}
\begin{bmatrix} \dot{\mathbf T}_F \\ \dot{\boldsymbol\lambda} \end{bmatrix}
+
\begin{bmatrix} K_{FF} & G_F^{\mathsf T} \\ G_F & 0 \end{bmatrix}
\begin{bmatrix} \mathbf T_F \\ \boldsymbol\lambda \end{bmatrix}
=
\begin{bmatrix} \mathbf F_F - K_{FD}\mathbf T_D \\ \mathbf g - G_D\mathbf T_D \end{bmatrix}.
$$

第二块行没有时间导数：$\boldsymbol\lambda$ 不参与时间推进，只是每个时刻的"瞬时反力"（边界热流）。θ 方法只作用在第一个块行上。以 $\theta = 1$（后向欧拉）为例，每步求解

$$
\begin{bmatrix} \dfrac{M_{FF}}{\Delta t} + K_{FF} & G_F^{\mathsf T} \\[1mm] G_F & 0 \end{bmatrix}
\begin{bmatrix} \mathbf T_F^{n+1} \\ \boldsymbol\lambda^{n+1} \end{bmatrix}
=
\begin{bmatrix} \dfrac{M_{FF}}{\Delta t}\mathbf T_F^{n} + \mathbf F_F^{n+1} - K_{FD}\mathbf T_D^{n+1} \\[1mm] \mathbf g^{n+1} - G_D\mathbf T_D^{n+1} \end{bmatrix}.
$$

一般 $\theta$ 的情形把惯性项写成残差形式更清楚：对第一个块行要求

$$
\frac{M_{FF}}{\Delta t}\bigl(\mathbf T_F^{n+1} - \mathbf T_F^{n}\bigr)
+ \theta\,\mathbf r^{n+1} + (1-\theta)\,\mathbf r^{n} = \mathbf 0,\qquad
\mathbf r^{n} \coloneqq K_{FF}\mathbf T_F^{n} + G_F^{\mathsf T}\boldsymbol\lambda^{n} + K_{FD}\mathbf T_D^{n} - \mathbf F_F^{n},
$$

即每步求解

$$
\begin{bmatrix} \dfrac{M_{FF}}{\Delta t} + \theta K_{FF} & \theta G_F^{\mathsf T} \\[1mm] G_F & 0 \end{bmatrix}
\begin{bmatrix} \mathbf T_F^{n+1} \\ \boldsymbol\lambda^{n+1} \end{bmatrix}
=
\begin{bmatrix} \dfrac{M_{FF}}{\Delta t}\mathbf T_F^{n} - (1-\theta)\bigl(K_{FF}\mathbf T_F^{n} + G_F^{\mathsf T}\boldsymbol\lambda^{n} + K_{FD}\mathbf T_D^{n} - \mathbf F_F^{n}\bigr) + \theta\mathbf F_F^{n+1} - \theta K_{FD}\mathbf T_D^{n+1} \\[1mm] \mathbf g^{n+1} - G_D\mathbf T_D^{n+1} \end{bmatrix}.
$$

注意两点：

1. 右端出现 $\boldsymbol\lambda^{n}$ **不是**乘子被"积分"——乘子始终是代数量（约束瞬时成立），只是 θ 加权需要上一时刻的完整残差；若取 $\theta = 1$ 该项消失，无需保存历史乘子；
2. $\theta \ne 1$ 时耦合块为 $\theta G_F^{\mathsf T}$，矩阵不再严格对称；把乘子重新标度 $\tilde{\boldsymbol\lambda} = \theta\boldsymbol\lambda$ 可恢复到 $[M/\Delta t + \theta K_{FF},\ G_F^{\mathsf T};\ G_F,\ 0]$ 的对称形式（右端相应变为 $\tilde{\boldsymbol\lambda}^{n}/\theta$），否则直接用稀疏 LU。工程上常见做法是瞬态用拉格朗日乘子（$\theta = 1$、一阶但无条件稳定），而 Crank–Nicolson 配合罚函数法或消元法。

### 适定性与 inf-sup 条件

乘子空间不能随意取。离散鞍点问题的可解性要求**离散 inf-sup（LBB）条件**：

$$
\inf_{\mu_h \in \Lambda_h}\ \sup_{v_h \in V_h}
\frac{\int_{\Gamma_D}\mu_h\,v_h\,d\Gamma}{\lVert v_h \rVert_{1} \lVert \mu_h \rVert_{0,\Gamma_D}} \ge \beta > 0 .
$$

本节"乘子取 $T_h$ 在 $\Gamma_D$ 上的迹空间"的做法满足该条件：约束块 $G_D$ 就是边界质量矩阵，在乘子空间上强制正定。相比之下，若在同一节点上放多个乘子、或乘子自由度多于约束所能分辨的模态，$G_D$ 奇异、鞍点系统不可解。

**退化情形**：当 $\Gamma_D$ 退化为孤立节点（`Poi1` 点单元，例如 2D 中只固定一个角点的温度），迹空间维数就是节点数，$G$ 退化为单位阵，$\lambda_a$ 变为该节点的**总热流**（单位 W，而非 W/m²）。

### 与罚函数法的关系

罚函数法可以看成把乘子**消去**后的正则化版本。在 Lagrangian 上加罚项得**增广拉格朗日**

$$
\mathcal{L}_\alpha(T,\lambda) = \mathcal{L}(T,\lambda) + \frac{\alpha}{2}\int_{\Gamma_D}\bigl(T - T_D\bigr)^2 d\Gamma ,
$$

若在乘子上做 Uzawa 迭代 $\lambda \leftarrow \lambda + \alpha(T - T_D)$，就是把精确约束逐步"罚"出来的乘子实现。经典罚函数法相当于只做一次、把乘子初值取零并**丢弃乘子**：它的罚刚度项 $\alpha\int_{\Gamma_D} N_i N_j\,d\Gamma$ 与罚载荷项 $\alpha\int_{\Gamma_D} N_i T_D\,d\Gamma$ 正是本节 $G$ 与 $\mathbf g$ 的 $\alpha$ 倍，即

$$
K \;\longrightarrow\; K + \alpha G,\qquad \mathbf F \;\longrightarrow\; \mathbf F + \alpha\mathbf g
$$

（只改 $D$ 相关行列）。因此 $\alpha \to \infty$ 时罚函数解收敛到鞍点解的 $\mathbf T$ 部分，而 $\alpha(T_h - T_D) \to \lambda$。这解释了两点工程经验：

- 罚函数法**没有**额外的约束精度控制手段，只能靠加大 $\alpha$，且以条件数为代价；
- 想直接得到边界热流，用拉格朗日乘子法更干净；若已用罚函数法，反力 $\alpha(T_h - T_D)$ 精度较差，更可靠的做法是用能量（热平衡）校验来反推总热流。

### 能量（热平衡）校验

无论用哪种方法，稳态解都应满足全局热平衡。稳态时取 $v \equiv 1$：

$$
\sum_{a} \lambda_a = \int_\Omega Q\,d\Omega - \int_{\Gamma_N} q_n\,d\Gamma + \int_{\Gamma_R} h\bigl(T_\infty - T\bigr)\,d\Gamma ,
$$

即"净产热扣除 Neumann 与 Robin 的流出后，剩下的全部经 Dirichlet 边界流出"。这条恒等式既是拉格朗日乘子法算例的现成单元测试断言，也是罚函数法结果的自检手段（用 $K\mathbf T - \mathbf F$ 在约束节点上的残差近似 $\lambda$）。

## 7. 混合（Hellinger–Reissner）形式：通量与温度独立

前两节把温度作为唯一未知量；混合形式改为把**通量** $\boldsymbol q$ 与**温度** $T$ 都作为独立未知量，一阶系统为（约定 $\boldsymbol q=-k\nabla T$，$\boldsymbol n$ 外法向）

$$
\boldsymbol q+k\nabla T=\boldsymbol 0,\qquad \nabla\cdot\boldsymbol q+f=0 .
$$

对第一个方程乘测试函数 $\boldsymbol r$ 并**逐单元分部积分**，第二个方程乘 $w$，得到 §4.2 的形式 B：

$$
\text{(B1)}\quad \int_\Omega\frac1k\boldsymbol q\cdot\boldsymbol r\,d\Omega
-\int_\Omega T\,\nabla\cdot\boldsymbol r\,d\Omega
=-\int_{\Gamma_D}T_D\,\boldsymbol r\cdot\boldsymbol n\,ds ,
$$

$$
\text{(B2)}\quad -\int_\Omega\nabla\cdot\boldsymbol q\,w\,d\Omega
=\int_\Omega f\,w\,d\Omega .
$$

与原始形式相比有两条关键变化：温度成为 $L^2$ 乘子（允许跨单元间断），Dirichlet 数据 $T_D$ 变成**自然边界条件**（进入 (B1) 右端），而 Neumann 数据 $\boldsymbol q\cdot\boldsymbol n=q_n$ 变成通量空间上的**本质约束**。

离散取 $\boldsymbol q_h\in\mathcal{RT}_0$、$T_h\in\mathcal{P}_0$，记 $\boldsymbol A_{ij}=\int\frac1k\boldsymbol q_i\cdot\boldsymbol q_j$、$C_{ce}=-\int N_c\nabla\cdot\boldsymbol q_e$、$\ell_i=-\int_{\Gamma_D}T_D\,\boldsymbol q_i\cdot\boldsymbol n$、$m_c=\int fN_c$，得到鞍点系统

$$
\begin{bmatrix}\boldsymbol A & \boldsymbol C^{\mathsf T}\\ \boldsymbol C & \boldsymbol 0\end{bmatrix}
\begin{bmatrix}\boldsymbol q\\ \boldsymbol T\end{bmatrix}
=\begin{bmatrix}\boldsymbol\ell\\ \boldsymbol m\end{bmatrix},
\qquad \boldsymbol A\ \text{对称正定},\ \boldsymbol C\ \text{矩形}.
$$

它对称但**不定**，必须用 `sp_lu`（或 $LDL^{\mathsf T}$），不能走 Cholesky。若消去 $\boldsymbol q$，得到 $\boldsymbol T$ 上的 Schur 补 $-\boldsymbol C\boldsymbol A^{-1}\boldsymbol C^{\mathsf T}$；局部守恒是精确的：取 $w=1$ 的 (B2) 给出 $\nabla\cdot\boldsymbol q_h=(\Pi_0f)|_T$。

### 与 ansatz 的对应

| 弱形式项 | 算子 | 说明 |
|---|---|---|
| $\boldsymbol A$ | `heat::mixed_flux_mass(&flux, k)` | 返回**完整** $\boldsymbol A$（鞍点矩阵走 LU） |
| $\boldsymbol C$ | `heat::mixed_divergence(&flux, &temperature)` | 返回完整 $M\times N$ 块，行温度、列通量；$\boldsymbol C^{\mathsf T}$ 由行列互换 |
| $\boldsymbol\ell$ | `heat::mixed_dirichlet(&traces, T_D)` | 边界迹单元为 `RTSeg3`；含**全部**三个父单元通量自由度 |
| $\boldsymbol m$ | `heat::heat_source(&temperature, f)` | 对 $\mathcal P_0$ 即 $\int_T f$ |

对应的单元：通量用 `element::rttri3::RTTri3<D>`（实现 `HdivElement<D,2,3>`），温度用 `element::p0tri3::P0Tri3<D>`（单元一个自由度，`Element<D,2,1>`），Dirichlet 边界积分用 `element::rttri3::RTSeg3<D>`（`RTSeg3::from_rttri3(&rt, facet)`，`normal_flux_shape` 给出 $\boldsymbol q_i\cdot\boldsymbol n$）。

**边界条件**：Dirichlet 用 `mixed_dirichlet` 的自然项；Neumann 是本质条件，直接把边界边的通量自由度固定为 $\Phi_e=\ell_e q_n$，**不要**再加 `heat_flux` 载荷。Robin 条件在 (B1) 的边界项中代入 $\boldsymbol r\cdot\boldsymbol n=h(T-T_\infty)$，会同时产生 $T$ 的质量型边界块与常数载荷；目前尚未提供专用的混合 Robin 算子。

**方向符号**：`RTTri3` 的基函数按单元**外法向**定向。相邻两个单元在公共边上的外法向相反，因此整体组装时必须对每个（单元，边）引入 $\sigma_{e,T}=\pm1$（以全局边法向为准），并把它乘到相应的 $\boldsymbol A$ 行/列、$\boldsymbol C$ 列与 $\boldsymbol\ell$ 分量上。单元素算例（测试 `mixed_heat_saddle_point_reproduces_a_linear_temperature`）中 $\sigma\equiv1$，无需处理。

求解一个单元上的 RT0/P0 混合系统，$k=1$、$f=0$、$T=2+x+2y$ 在三条边上给定，可精确重现 $\boldsymbol q=(-1,-2)$ 与单元平均温度 $17/3$，该测试同时验证了 $\boldsymbol A$、$\boldsymbol C$、$\boldsymbol\ell$ 的符号。

## 8. 与 ansatz 抽象的对应

热传导完全复用现有抽象（`Element`、`Variable`、`IntegrationScheme`），只需新增一组热算子。参照 `src/operation/bar.rs` 的四个函数：

| 现有（bar.rs） | 热算子（规划） | 物理含义 |
|---|---|---|
| `stiffness(v, ea)` | `conduction_stiffness(v, k)` | $K^e_{ij} = \sum_g k\,\nabla N_i\cdot\nabla N_j\,|J|\,w_g$ |
| `body_force(v, f)` | `heat_source(v, Q)` | $F_i = \sum_g Q\,N_i\,|J|\,w_g$ |
| `traction(v, p)` | `heat_flux(v, q_n)` | $F_i = -\int_{\Gamma_N} N_i\, q_n\, d\Gamma$（注意符号，沿边界积分） |
| `displacement_penalty(v, α, ū)` | `temperature_penalty(v, α, T_D)` | 罚函数强加定温（沿边界积分） |
| —（已实现） | `temperature_lagrangian(T, λ, T_D)` | 拉格朗日乘子：返回完整耦合块 $G_{aj} = \int_{\Gamma_D} \Lambda_a N_j\,d\Gamma$（$M\times N$ 全部元素）与 $g_a = \int_{\Gamma_D} \Lambda_a T_D\,d\Gamma$；$T$ 与 $\lambda$ 可用不同形函数，需解鞍点系统（见 §6） |
| —（新增） | `convection(v, h, T_∞)` | Robin 项：刚度 $h\int N_i N_j\,d\Gamma$ + 载荷 $h T_\infty\int N_i\,d\Gamma$（沿边界积分） |
| —（新增） | `capacity(v, ρc_p)` | 质量矩阵 $M_{ij} = \sum_g \rho c_p N_i N_j\,|J|\,w_g$ |
| —（新增） | `mixed_flux_mass(flux, k)` | 混合形式通量块 $A_{ij}=\int (1/k)\boldsymbol q_i\cdot\boldsymbol q_j$（RT0，完整矩阵） |
| —（新增） | `mixed_divergence(flux, T)` | 混合形式耦合块 $C_{ce}=-\int N_c\nabla\cdot\boldsymbol q_e$（完整矩形块） |
| —（新增） | `mixed_dirichlet(traces, T_D)` | 混合形式 Dirichlet 自然载荷 $\ell_i=-\int_{\Gamma_D}T_D\,\boldsymbol q_i\cdot\boldsymbol n$（`RTSeg3` 边界迹） |

要点：

- `conduction_stiffness` 与 `stiffness` 仅差一个系数（$k$ 代替 $ea$）：二者都计算物理梯度 $\nabla_x N$ 的点积 $\nabla N_i \cdot \nabla N_j$（`derivative_shape` 现返回关于物理坐标的完整梯度），线单元上 $\nabla_x N = \pm \hat t/L$ 沿切向，退化为 $dN_i/dx \cdot dN_j/dx$；
- `heat_flux` 与 `traction` 的差别是**符号**：Neumann 通量 $q_n$ 进入载荷时带负号（$-\int N_i q_n\,d\Gamma$），源自约定 $k\frac{\partial T}{\partial n} = -q_n$；二者均作用于任意维度的边界单元（`Poi1` 点单元退化为集中通量）；
- 温度 Dirichlet 条件用**罚函数法**，与 `displacement_penalty` 完全同构：对任意维度边界单元生成 $K_{ij} \mathrel{+}= \alpha\int N_i N_j\,d\Gamma$、$F_i \mathrel{+}= \alpha\int N_i T_D\,d\Gamma$（`Poi1` 点单元退化为 $K_{ii}\mathrel{+}=\alpha$、$F_i\mathrel{+}=\alpha T_D$；见「边界条件与罚函数法」章）；
- 若要**精确**满足定温约束并直接输出 Dirichlet 边界热流 $\lambda = q_n|_{\Gamma_D}$，则改用 §6 的拉格朗日乘子法：`temperature_lagrangian(&T, &λ, T_D)` 组装耦合块 $G_{aj} = \int_{\Gamma_D} \Lambda_a N_j\,d\Gamma$ 与 $g_a = \int_{\Gamma_D} \Lambda_a T_D\,d\Gamma$（两组边界单元配对迭代：积分单元/积分点/雅可比取温度单元，乘子单元只提供形函数与自由度 id；传同一个集合即得边界质量矩阵）。它返回**完整**的 $G$ 而非下三角，转置块 $G^{\mathsf T}$ 由行列互换得到；求解器须支持对称不定系统（`sp_lu` 稀疏 LU / $LDL^{\mathsf T}$），不能走当前的 `sp_cholesky` 链路；且 LU 读取整矩阵，故 `conduction_stiffness` 等返回的下三角三元组也须先对称展开。两组单元数量不一致时返回 `VariableError::ElementCountMismatch`。

### 实现要点

- **稳态**热传导只需 `conduction_stiffness` + `heat_source` + `temperature_penalty`（+ 可选的 `heat_flux` / `convection`），组装后走稀疏 Cholesky，与 `bar1d` 同一条求解链路；
- **瞬态**需额外 `capacity` 与时间循环，每步解 $\left(\frac{M}{\Delta t} + \theta\bar K\right)\mathbf T^{n+1} = \cdots$，隐式格式仍保持对称正定；
- 2D/3D 上 `Tri3`/`Quad4` 的 $\nabla N_i$ 由 `derivative_shape` 直接返回（已结合雅可比矩阵求物理梯度，见「雅可比与坐标映射」章），单元刚度即 $B^\mathsf T k B$ 形式；`conduction_stiffness`、`capacity`、`heat_source` 均已泛化到任意空间维度 `D`，对任意维度单元物理正确。
