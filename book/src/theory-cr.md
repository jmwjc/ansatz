# 非协调单元：Crouzeix–Raviart（CR）单元

前面各章的单元（`Seg2`、`Tri3`、`Quad4`）都是**协调**（conforming）单元：相邻单元在公共边界上逐点连续，因而全局空间是 $H^1$ 的子空间，可以直接套用伽辽金理论。CR 单元（Crouzeix & Raviart, 1973）主动**放弃逐点连续**，只保留公共边界上的**积分平均连续**。它用极小的代价换来两个重要特性：一个好的不可压缩（Stokes/近不可压弹性）速度空间，以及**对角的质量矩阵**。

> 实现：`ansatz::element::crtri3::CRTri3<D>`（其中一个自由度 = 一条边的平均值），边界条件用 `CRSeg3<D>` 迹单元；见末节「与 ansatz 抽象的对应」。

## 1. 动机：为什么放弃协调性

协调 $P_1$ 单元的自由度取在**顶点**，等价于要求函数在顶点、进而沿整个界面逐点连续。这在两个场合成为障碍：

1. **不可压缩问题**：Stokes 方程要求速度和压力满足离散 inf-sup（LBB）条件。$P_1/P_1$（同阶）不满足，会出现棋盘状压力；需要一个"骑在边上"的速度空间。CR 速度 + 分片常数压力 $P_1^{\text{nc}}/P_0$ 是经典稳定对，且比 Taylor–Hood（$P_2/P_1$）的自由度更少、质量矩阵更简单。
2. **弱连续性已经足够**：很多问题的弱形式只通过界面**通量/平均值**耦合，逐点连续是多余的强要求。用边平均作自由度，可以少一层网格层的约束（自由度挂在边上，而不是顶点）。

于是我们把自由度换成：**函数在每条边上的平均值**。相邻两个单元共享同一条边、同一个边自由度，它们的界面值仍可能有跳跃，但跳跃的**平均值**为零。

## 2. 2D 三角形 CR 单元

### 自由度与局部空间

在三角形 $T$（顶点 $A_1,A_2,A_3$，边 $e_i$ 是顶点 $A_i$ 的**对边**）上取

- 局部空间：$\mathcal P_1(T)$（线性函数），$\dim = 3$；
- 自由度：三条边的平均值
  $$
  \ell_i(v) = \frac{1}{|e_i|}\int_{e_i} v\, ds, \qquad i = 1,2,3,
  $$
  对线性函数而言，边上的平均值恰好等于**边中点**处的函数值。

自由度个数与空间维数一致（都是 3），且下面的基函数显式构造证明了**唯一可解性**（unisolvence），因此它们定义一个合法的有限元。

### 基函数

记 $\lambda_i$ 为面积（重心）坐标，$M_i$ 为边 $e_i$ 的中点。取

$$
\theta_i = 1 - 2\lambda_i, \qquad i = 1,2,3 .
$$

验证：在 $e_i$ 上 $\lambda_i \equiv 0$，故 $\theta_i \equiv 1$，其边平均为 1；在 $e_j$（$j\ne i$）上 $\lambda_i$ 由 1 线性变到 0，平均为 $1/2$，故 $\theta_i$ 的边平均为 $1-2\cdot\tfrac12 = 0$。即

$$
\frac{1}{|e_j|}\int_{e_j}\theta_i\,ds = \delta_{ij},
$$

这正是 Kronecker δ 性质——只是发生在**边平均**上，而不是在节点值上。

| 量 | 表达式 | 说明 |
|---|---|---|
| 基函数 | $\theta_i = 1 - 2\lambda_i$ | 分片线性 |
| 顶点处取值 | $\theta_i(A_i) = -1$，$\theta_i(A_j) = +1\ (j\ne i)$ | **不是** Kronecker δ |
| 边中点取值 | $\theta_i(M_i) = 1$，$\theta_i(M_j) = 0\ (j\ne i)$ | 中点值即该边平均 |
| 单位分解 | $\sum_i \theta_i = 1$ | 常数场可精确表示 |
| 几何再现 | $\sum_i \theta_i M_i = \mathbf x$ | 见下文等参映射 |
| 梯度 | $\nabla\theta_i = -2\nabla\lambda_i$ | 单元内为**常数** |

（沿边 $e_j$ $(j\ne i)$，$\lambda_i$ 从顶点 $A_i$ 处的 1 线性降到顶点 $A_k$ 处的 0，故 $\theta_i$ 从 $-1$ 变到 $+1$，其中点值、亦即边平均正好为 0。）

有两点值得强调：

- **不是插值型单元**：$\theta_i$ 在顶点上取 $\pm1$ 而非 $\delta_{ij}$，所以"节点值"没有意义，只有边平均有意义。后处理、边界条件、误差评估都必须按"边平均"来理解。
- **几何再现**：$\sum_i\theta_i(\mathbf x)M_i$ 是线性映射，且在三个中点处与恒等映射一致，故恒等于 $\mathbf x$。这意味着用中点作"节点"、$\theta_i$ 作形函数的等参映射**恰好重现仿射映射**——直线三角形上几何不会被扭曲。

### 等参映射与雅可比

三角形几何与 `Tri3` 完全相同，只是基函数换成了 $\theta_i$：

$$
\mathbf x(\xi) = \sum_{i=1}^{3}\theta_i(\xi)\, M_i ,
\qquad
J = \frac{\partial\mathbf x}{\partial\xi}
$$

$J$ 是常数矩阵，$|J| = 2|T|$（参考三角形面积为 $1/2$）。物理梯度同样按张量公式给出：

$$
\nabla_x\theta_i = j_\alpha\, g^{\alpha\beta}\, \partial_\beta\theta_i
= -2\, j_\alpha\, g^{\alpha\beta}\, \partial_\beta\lambda_i ,
$$

其中 $j_\alpha = \partial\mathbf x/\partial\xi^\alpha$、$g_{\alpha\beta} = j_\alpha\cdot j_\beta$，与「雅可比与坐标映射」章完全一致。

### 单元刚度矩阵

泊松型双线性形式 $a(v,w)=\int_T \nabla v\cdot\nabla w\,d\Omega$ 的被积函数在单元内是**常数**，于是

$$
K_{ij} = \int_T \nabla\theta_i\cdot\nabla\theta_j\,d\Omega
= 4\,|T|\;\nabla\lambda_i\cdot\nabla\lambda_j .
$$

因为它只是 $P_1$ 刚度矩阵乘 4，CR 的刚度矩阵与 `Tri3` 的仅差一个系数。以参考三角形（顶点 $(0,0),(1,0),(0,1)$，$|T|=1/2$）为例，$\nabla\lambda_1=(-1,-1)$、$\nabla\lambda_2=(1,0)$、$\nabla\lambda_3=(0,1)$，得

$$
K = \begin{bmatrix}
 4 & -2 & -2\\
-2 &  2 &  0\\
-2 &  0 &  2
\end{bmatrix},
$$

每一行之和为零：常数在核中，正是"刚体位移"（零能量模式）的离散对应。

### 单元质量矩阵：对角！

对线性三角形有解析公式 $\int_T \lambda_1^{a}\lambda_2^{b}\lambda_3^{c}\,d\Omega = 2|T|\,a!\,b!\,c!/(a+b+c+2)!$。代入 $\theta_i = 1-2\lambda_i$：

$$
\int_T\theta_i\theta_j\,d\Omega
= |T|\left[1 - \frac{4}{3} + \frac{1+\delta_{ij}}{3}\right]
= \frac{|T|}{3}\,\delta_{ij},
$$

即

$$
M_{ij} = \int_T \theta_i\theta_j\,d\Omega = \frac{|T|}{3}\,\delta_{ij},
\qquad
M = \frac{|T|}{3}\,\mathbf I .
$$

**质量矩阵天然对角**（且 $\frac{|T|}{3}=\frac14\sum_j M_{ij}$ 与集中质量一致）。这是 CR 单元最著名的性质：无需"集中化"操作就得到集中质量矩阵，因而在显式动力学、质量集中、以及需要 $M^{-1}$ 的算法（如压力 Poisson 投影）中极具吸引力。

> 注意：这一性质是**2D 三角形专属**的。3D 四面体 CR 的质量矩阵不再对角（见 §5）。

### 积分规则

- 刚度：被积函数为常数 → **1 点**（重心）积分即精确，对应 `ReducedIntegration<2, 1>` / `Dunavant<1>`；
- 质量：被积函数为二次 → 需 **3 点**（`Dunavant<3>`）才精确；但由于解析值已知（$\frac{|T|}{3}\delta_{ij}$），也可以直接使用该式而不做积分。

## 3. 全局空间、非协调性与收敛性

### 全局空间

$$
V_h = \Bigl\{ v_h \in L^2(\Omega) : v_h|_T \in \mathcal P_1(T)\ \forall T,\
\text{且 } \int_e \llbracket v_h\rrbracket\, ds = 0 \ \forall\, \text{内部边 } e \Bigr\},
$$

其中 $\llbracket\cdot\rrbracket$ 是跨边跳跃。每个内部边只有一个自由度（被相邻两个单元共享），边界边各有一个，所以 $V_h$ 的维数等于**边的总数**，而不是顶点数。关键点是

$$
V_h \not\subset H^1(\Omega),
$$

因为 $v_h$ 在界面处允许有（均值为零的）跳跃。这就是"**非协调**"（nonconforming）的含义。

跳跃 $\llbracket v_h\rrbracket$ 沿边是线性函数且均值为零，于是它一定在边的两端反号：跳跃是"翻转"型的，而不是任意函数。这一点是误差分析的全部来源。

### 误差分析：Strang 第二引理

对 $V_h \not\subset H^1$ 不能直接用 Céa 引理，需要分解为**逼近误差 + 相容性误差**：

$$
\|u - u_h\|_{1,h}
\;\le\; C\left(
\inf_{v_h\in V_h}\|u - v_h\|_{1,h}
\;+\;
\sup_{w_h\in V_h}\frac{|a_h(u,w_h) - l(w_h)|}{\|w_h\|_{1,h}}
\right),
$$

其中 $a_h(\cdot,\cdot)=\sum_T a_T(\cdot,\cdot)$、$\|\cdot\|_{1,h}$ 是**分片（broken）** $H^1$ 范数 $\bigl(\sum_T\|\cdot\|_{H^1(T)}^2\bigr)^{1/2}$。

- **逼近误差**：CR 插值算子 $\Pi_h$ 由边平均定义，$\int_{e}\Pi_h u\,ds = \int_e u\,ds$。由于单元上仍是分片线性，$\|u-\Pi_h u\|_{1,h}\le C h|u|_{H^2}$。
- **相容性误差**：$u$ 在单元内是光滑的，逐单元分部积分得
  $$
  a_h(u,w_h)-l(w_h) = \sum_{e}\int_e \llbracket w_h\rrbracket\,\bigl(\nabla u\cdot \mathbf n\bigr)\,ds .
  $$
  因为 $\int_e\llbracket w_h\rrbracket\,ds = 0$，可以把 $\nabla u\cdot\mathbf n$ 换成它在该边上的**均值**：只有通量的**波动部分**进入误差，
  $$
  \Bigl|\sum_e\int_e \llbracket w_h\rrbracket(\nabla u\cdot\mathbf n - \overline{\nabla u\cdot\mathbf n})\,ds\Bigr|
  \le C\sum_e h_e^{1/2}|u|_{H^2(T)}\,\|\llbracket w_h\rrbracket\|_{L^2(e)}
  \le C h\,|u|_{H^2}\,\|w_h\|_{1,h},
  $$
  即相容性误差是 $O(h)$ 的。**"跳跃均值为零"这一条弱连续性恰好把 $O(1)$ 的相容性误差降到了 $O(h)$**——这正是 CR 以放弃逐点连续为代价、仍能拿到最优收敛阶的机制。

合并两项，对 $u\in H^2(\Omega)$：

$$
\|u - u_h\|_{1,h} \le C h\,|u|_{H^2(\Omega)},
\qquad
\|u - u_h\|_{L^2(\Omega)} \le C h^2\,|u|_{H^2(\Omega)} .
$$

与协调 $P_1$ 的阶数**完全相同**，只是 $H^1$ 误差要在分片范数下度量。

### 广义 patch test

CR 能精确重现**线性场**：对仿射 $u$，$\Pi_h u = u$ 逐单元成立，且 $a_h(u,w_h)=0=l(w_h)$（因跳跃均值为零），所以 $u_h = u$。这就是非协调单元的"广义 patch test"（在边平均意义下精确，而不是在顶点值意义下）。它也是上面的相容性分析能给出 $O(h)$ 的根本原因：常应变/常通量场不会产生任何相容性误差。

### 本质边界条件

因为**没有顶点自由度**，Dirichlet 条件 $u = \bar u$ 应施加在边的平均值上：

$$
\frac{1}{|e|}\int_e u_h\,ds = \frac{1}{|e|}\int_e \bar u\,ds \qquad (\text{边界边 } e).
$$

用罚函数法时，直接在边界边上积分 $\alpha\int_e (u_h-\bar u)^2\,ds$（即 `displacement_penalty` / `temperature_penalty` 的形式）与 CR 的自由度天然匹配：它惩罚的正是边上的 $L^2$ 偏差，而不会引入顶点层面的虚假约束；$\alpha\to\infty$ 时该罚项把**迹**在 $L^2$ 意义下压向 $\bar u$，对连续数据等价于两端点插值，比"只约束边平均"略强，仍属标准罚函数法。用乘子法时约束行 $\int_e \Lambda(u_h-\bar u) = 0$ 亦然。实现这些边界积分需要一个"迹单元"，见 §6 的 `CRSeg3`。

## 4. 典型用途

- **Stokes 方程**：速度用 $(\text{CR})^D$、压力用 $P_0$。该对满足离散 inf-sup 条件（原始文献的出发点），速度分片线性、压力分片常数，自由度数目与 $P_1/P_1$ 相当，但压力不会出现棋盘失稳；
- **近不可压弹性 / 大体积模量**：$(\text{CR})^D/P_0$ 是位移–压力（混合）格式，用于规避体积锁定（与「线弹性力学」的体积–偏量分解、`Cook 膜`算例中的 SRI 是同一目的的两条路线）；
- **广义/弱连续性界面问题**：接触、粘结、多材料界面，只要求界面量的平均匹配；
- **显式动力学与需要 $M^{-1}$ 的算法**：2D 质量矩阵对角，$M^{-1}$ 即逐自由度常数；
- **后处理**：CR 与最低阶 Raviart–Thomas（RT0）互为旋转对偶（2D 中把 RT0 的向量场旋转 90° 得到 CR），可用于通量恢复/超收敛后处理；RT0 的边通量与 CR 的刚度矩阵逐行对应，详见「混合格式与 $H(\mathrm{div})$：Raviart–Thomas（RT0）」章。

与常见单元的对比：

| 单元 | 空间 | 自由度位置 | 跨单元连续 | 质量矩阵 | 典型用途 |
|---|---|---|---|---|---|
| 协调 $P_1$（`Tri3`） | $\mathcal P_1$ | 顶点值 | 逐点 | 非对角 | 泊松、弹性（协调） |
| **CR**（$\mathcal P_1^{\text{nc}}$） | $\mathcal P_1$ | 边平均 | 平均 | **对角**（2D） | Stokes、近不可压、显式 |
| RT0 | $\mathcal P_0^2 \oplus \mathbf x\,\mathcal P_0$ | 边通量 | 法向逐点 | 非对角 | 混合泊松、通量 |

## 5. 3D：四面体 CR 单元

四面体有 4 个面，取自由度为该面上的平均值（重心处的值），即

$$
\theta_i = 1 - 3\lambda_i, \qquad i = 1,\dots,4,
$$

$D$ 维单纯形的一般形式为 $\theta_i = 1 - D\lambda_i$，自由度数 $= D+1$（面数），每个面一个自由度。刚度矩阵仍为

$$
K_{ij} = D^2\,|T|\,\nabla\lambda_i\cdot\nabla\lambda_j ,
$$

但**质量矩阵不再对角**。用一般单纯形积分公式（$\int_T\lambda_i = |T|/(D+1)$、$\int_T\lambda_i\lambda_j = |T|/((D+1)(D+2))$，$i\ne j$）可得

$$
\int_T\theta_i^2 = \frac{2}{5}|T|,
\qquad
\int_T\theta_i\theta_j = -\frac{1}{20}|T| \quad (i\ne j),
$$

行和 $=\frac{2}{5}|T| - 3\cdot\frac{1}{20}|T| = \frac14|T|$，与 $\int_T\theta_i = |T|/4$ 一致。非对角元为负：3D CR 的质量矩阵虽仍对称正定，但**不再对角**，2D 的对角性质无法推广。

## 6. 与 ansatz 抽象的对应

CR 单元可以直接塞进现有的 `Element` / `Variable` / `IntegrationScheme` 骨架，唯一的新问题是**全局编号**。

### 单元类型

`ansatz::element::crtri3::CRTri3<D>` 实现 `Element<D, 2, 3>`：参数维 $P=2$、基函数个数 $N=3$，与 `Tri3` 相同，差别只在形函数、其梯度以及"节点"的含义（三个边中点）：

```rust
// 参考三角形 (0,0), (1,0), (0,1)：面积坐标 λ = [1 - ξ - η, ξ, η]
fn shape(&self, [ξ, η]: [f64; 2]) -> [f64; 3] {
    let lambda = [1.0 - ξ - η, ξ, η];
    std::array::from_fn(|i| 1.0 - 2.0 * lambda[i])
}

// 物理梯度为常数：∇θ_i = -2 ∇λ_i
fn derivative_shape(&self, xi: [f64; 2]) -> [[f64; D]; 3] { /* -2 × ∇λ_i */ }
```

几何部分（`jacobe`、`jacobe_mat`、`local_coordinates`）由三个中点反射出的三角形顶点
$A_i = M_j + M_k - M_i$ 委托给 `Tri3` 完成，因此与 `Tri3` 逐位一致；`coordinates` 仍按
$\sum_i\theta_i M_i$ 计算，结果与 $\sum_i\lambda_i A_i$ 相同。

集成规则：刚度用 `Dunavant<1>`（`ReducedIntegration<2, 1>`，被积函数为常数），质量用
`Dunavant<3>`（`FullIntegration<2, 3>`）。

### 自由度挂在"边"上

`ansatz` 用 `Element::id()` 返回**节点 id** 作为全局自由度编号，而 CR 的自由度属于**边**。`CRTri3<D>` 采用**中点节点法**：它的三个"节点"就是三条边的中点，`Node::id` 存放该边自由度的全局编号、`Node::coordinates` 存放中点坐标。由 $\sum_i\theta_i M_i=\mathbf x$ 可知等参映射仍精确重现三角形几何（见 §2），故这种做法自洽；真实的三角形顶点可用 `CRTri3::triangle_vertices()`（反射公式 $A_i = M_j+M_k-M_i$）取回。

因此构造 CR 单元的核心是**给每条边分配唯一 id**。`CRTri3::from_tri3(&tri, edge_id)` 为此而设：`edge_id(a, b)` 接收一条边的两个端点顶点 id（升序），返回该边的全局自由度编号。GMSH 网格的读取链路已封装好这一步：

```rust
use ansatz::element::crtri3::CRTri3;
use ansatz::io::gmsh::{GmshOptions, read_gmsh_with};

let mesh = read_gmsh_with::<2>(
    path,
    GmshOptions {
        renumber: true, // 顶点压紧为 0..n_vertices
        edges: true,    // 读取时给边编号，接在节点之后
    },
)?;
let crs: Vec<Cr3<2>> = mesh.crtri3_elements("Ω")?;
```

`EdgeTable` 内部就是"顶点对 → 边编号"的哈希表，与端点顺序无关，因此相邻单元在公共边上拿到同一个自由度。若手工构造网格，等价的写法是

```rust
let mut edges = ansatz::io::gmsh::EdgeTable::new(n_vertices);
let crs: Vec<CRTri3<2>> = tris.iter().map(|t| CRTri3::from_tri3(t, |a, b| edges.id(a, b))).collect();
```

编号布局：

```text
顶点 id:  0 .. n_vertices   (几何用；CR 解不使用这些自由度)
边   id:  0 .. n_edges      (CR 的未知数)
```

两套编号**相互独立、允许重叠**：它们索引的是不同的未知量，所以不需要像"节点 + 乘子"那样错开。只有当 CR 与另一个场要放进同一个全局方程组时，才需要用 `EdgeTable::new(offset)` 把边自由度整体平移。

「GMSH 网格读取」章为此提供了两个预处理选项：`GmshOptions { edges: true }` 在读取时为所有三角形的边建好编号（`0..n_edges`，跨所有三角形组统一，故子域界面共享同一个自由度），`GmshOptions { renumber: true }` 只负责把顶点压紧为 `0..n_vertices`，好让 `node_table()` 能构造单元。

### 边界条件的离散：facet 迹单元 `CRSeg3` / `CRTri4`

CR 没有顶点自由度，而**一条边上的迹**依赖相邻三角形的**全部三个**边自由度：

$$
u_h\big|_e(\xi) = \theta_1\big|_e(\xi)\,u_1 + \theta_2\big|_e(\xi)\,u_2 + \theta_3\big|_e(\xi)\,u_3 ,
$$

其中只有 $u_{e}$ 属于这条边。因此边界积分不能用 `Seg2`（它的自由度是**顶点值**），必须用"父单元的迹"作为边界单元：

```rust
// 二维：父三角形 3 个 facet（边），facet 自身是线段（P = 1，3 个自由度）
pub struct CRSeg3<const D: usize> { /* ... */ }   // impl Element<D, 1, 3>
// 三维：父四面体 4 个 facet（面），facet 自身是三角形（P = 2，4 个自由度）
pub struct CRTri4<const D: usize> { /* ... */ }   // impl Element<D, 2, 4>
```

二维的迹用 `CRSeg3::from_crtri3(&parent, facet)` 构造，三维用 `CRTri4::from_crtet4(&parent, facet)`，父单元分别是 `CRTri3` 与 `CRTet4`；`CRTet4` 自身由 `CRTet4::from_tet4(&tet4, face_id)` 从协调四面体 `Tet4` 生成（面自由度编号由 `FaceTable` 或调用方维护）。

命名沿袭 `Seg2`/`Tri3` 的"拓扑 + 自由度个数"习惯：`CRSeg3` = CR 的、3 个自由度的线段，`CRTri4` = CR 的、4 个自由度的三角形。

**局部 facet 编号**：`facet = k` 表示"父顶点 $k$ 不在其上的那个 facet"，与父单元的自由度顺序（`nodes[i]` = 顶点 $i$ 的对面/对边）完全一致。因此

$$
\texttt{CRSeg3::id()} \equiv \texttt{CRTri3::id()} ,
$$

迹单元与父单元共用同一套全局自由度编号，边界条件和内部组装天然对齐。facet 自身的参数方向由"**facet 顶点按父顶点下标升序排列**"固定，所以 `shape` 与网格读取器给出的连线方向无关。

二维迹基函数（`lo < hi` 为 facet 的两个顶点下标）：

| `facet` | facet 顶点 | `shape(ξ)` |
|---|---|---|
| 0 | (1, 2) | `[1, ξ, −ξ]` |
| 1 | (0, 2) | `[ξ, 1, −ξ]` |
| 2 | (0, 1) | `[ξ, −ξ, 1]` |

即 `shape[facet] = 1`、`shape[lo] = ξ`、`shape[hi] = −ξ`（三者之和恒为 1）；三维同构：`shapeₐ = 1 − 3λₐ`，其中 $\lambda$ 由 facet 的重心坐标散播进父单元的四分量数组（对顶点分量为 0），此时 `shape[facet] ≡ 1`。

几何是 $P_1$（或三维的 $P_2$）的，直接委托给 facet 顶点构成的 `Seg2`/`Tri3`：`jacobe = L/2`（二维）或 `2×面积`（三维）。注意它是**非等参**的（$\sum_a\texttt{shape}_a x_a \ne x$），边界算子只用到 `shape`、`jacobe` 与 `coordinates`，因此无碍。

有了迹单元，所有边界算子都能**原样复用**：

```rust
let neumann: Vec<CRSeg3<2>> = mesh.crtri3_boundary("Γⁿ")?;
let dirichlet: Vec<CRSeg3<2>> = mesh.crtri3_boundary("Γᵈ")?;

let mut k = conduction_stiffness(&tris, k_cond);
let mut f = heat_flux(&neumann, q_n);                          // ∫ θₐ q dΓ（严格非协调载荷）
let (pk, pf) = temperature_penalty(&dirichlet, 1e7, ū);        // α ∫ (θₐuₐ − ū)² dΓ
k.extend(pk);
f.extend(pf);
```

gmsh 侧只需打开读取期预处理，父单元与局部 facet 编号由读取器自动确定（用顶点对匹配三角形局部的三条边）：

```rust
// 2D
let mesh = read_gmsh_with::<2>(path, GmshOptions { renumber: true, edges: true, faces: false })?;
let facets: Vec<CRSeg3<2>> = mesh.crtri3_boundary("Γᵈ")?;

// 3D
let mesh = read_gmsh_with::<3>(path, GmshOptions { renumber: true, edges: false, faces: true })?;
let faces: Vec<CRTri4<3>> = mesh.crtet4_boundary("Γᵈ")?;
```

线性精确解的 patch test（`crtri3_patchtest_with_dirichlet_penalty`）给出机器精度量级的吻合（误差 < $10^{-6}$），验证了"读取 → 迹基 → 罚函数 → 求解"整条链。

### 组装与后处理

- 刚度组装与 `Tri3` 完全相同（只是 `shape`/`derivative_shape` 不同），`conduction_stiffness`、`stiffness` 等算子**无需修改**即可处理 CR 单元（`CRTri3` 的传导刚度为同网格 `Tri3` 的 4 倍，测试 `crtri3_assembles_with_the_heat_conduction_operator` 即验证这一点）；
- 质量矩阵在 2D 对角，`capacity` 逐单元组装后自然得到对角矩阵（若用 1 点积分只会在形心处取 $\theta_i=\frac13$，得到 $\frac{|T|}{9}\mathbf 1\mathbf 1^\mathsf T$ 而**不是**对角——务必用 3 点积分或解析值）；
- 由于 $V_h\not\subset H^1$，误差必须在 broken 范数下度量：`l2_error`（逐单元 $L^2$）直接适用；`h1_semi_error` 逐单元累加 $\int_T|\nabla u_h-\nabla u|^2$，得到的是**分片** $H^1$ 半范——正是 §3 中理论估计所用的范数，无需修改；
- 结果场的云图应逐单元绘制（跨边不连续），不要假设相邻单元共享轨迹。

## 小结

- CR 用**边平均**替代顶点值：局部仍是 $\mathcal P_1$，基函数 $\theta_i = 1 - D\lambda_i$，靠"跳跃均值为零"的弱连续性工作；
- 刚度与 $P_1$ 只差系数 $D^2$，**2D 质量矩阵恰好对角**（$|T|/3$），3D 不再对角；
- 非协调性把相容性误差降到 $O(h)$（Strang 第二引理），收敛阶与协调 $P_1$ 相同（$H^1$ 分片一阶、$L^2$ 二阶），并精确重现线性场（广义 patch test）；
- 经典舞台是 Stokes/近不可压问题（CR–$P_0$ 稳定）与需要集中质量的显式算法；
- 引入成本主要在**边的全局编号**，其余与 `Tri3` 完全同构。
