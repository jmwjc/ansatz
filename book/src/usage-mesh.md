# GMSH 网格读取

`ansatz` 的 `io` 模块支持读取 **GMSH v4.1 ASCII** 格式的 `.msh` 文件。

## 入口函数

```rust
use ansatz::io::gmsh::read_gmsh;

let mesh = read_gmsh::<2>("examples/msh/patchtest_tri3.msh").unwrap();
```

`read_gmsh::<D>` 解析文件并返回中间产物 `GmshMesh<D>`，其中 `D` 是期望的空间维数。它**只读文件、不构造具体单元类型**；类型化单元向量按需由 `GmshMesh::elements` 构造。

### 预处理选项：节点重编号

`read_gmsh` 等价于"不启用任何预处理"。需要重编号时用 `read_gmsh_with`：

```rust
use ansatz::io::gmsh::{read_gmsh_with, GmshOptions};

let mesh = read_gmsh_with::<2>(
    "examples/msh/patchtest_tri3.msh",
    GmshOptions { renumber: true },
)
.unwrap();
```

`GmshOptions { renumber: true, edges: false, faces: false }`：三个开关分别控制节点重编号、三角形边编号（2D CR）、四面体面编号（3D CR），`Default` 全为 `false`（与 `read_gmsh` 完全相同，节点 id 保持 GMSH 原始编号 `tag - 1`）。`renumber` 设为 `true` 时：

- 节点 id 按**原始 id 升序**压紧为 `0..n_nodes`，单元连接关系同步重映射；
- 因此 `node_table()` 必然成功（无空洞），`NodeTable::len()` 就是节点总数；
- 连接关系引用了 `$Nodes` 中不存在的节点标签时，**读取阶段**即返回 `MeshError::MissingNode`（否则要等到物化单元时才报错）。

**为什么拉格朗日乘子需要它**：乘子自由度要接在主变量自由度之后，而乘子的编号是按节点 id 计算的（`id·D + p`）。只有网格节点占据连续的 `0..n` 段，才能安全地把乘子节点放在 `n..n+m` 而不会撞号：

```rust
use ansatz::element::Element;
use ansatz::element::poi1::Poi1;
use ansatz::io::gmsh::{GmshOptions, read_gmsh_with};
use ansatz::node::Node;

let mesh =
    read_gmsh_with::<1>("examples/msh/bar1d.msh", GmshOptions { renumber: true }).unwrap();
let table = mesh.node_table().unwrap();
let n = table.len(); // 网格节点占据 0..n

// 定温/固定的边界点（物理标签 0 的点单元）→ 每个点一个乘子节点，编号 n..n+m。
let ends: Vec<Poi1<1>> = mesh.elements_by_tag(0).unwrap();
let multipliers: Vec<Poi1<1>> = ends
    .iter()
    .map(|end| {
        let id = end.id()[0];
        Poi1::new(Node::new(n + id, table.get(id).unwrap()))
    })
    .collect();
```

这样 `temperature_lagrangian` / `displacement_lagrangian` 直接产出的乘子行列就是全局自由度编号，无需再手动平移。

## `GmshMesh<D>`

```rust
pub struct GmshMesh<const D: usize> {
    pub nodes: HashMap<usize, Node<D>>,
    pub physical_names: HashMap<i32, String>,
    // blocks: HashMap<i32, ElementBlock>  (私有)
}
```

关键方法：

```rust
// 按物理组名构造类型化单元向量
pub fn elements<T: FromGmsh<D>>(&self, physical_name: &str) -> Result<Vec<T>, MeshError>;

// 按物理标签构造
pub fn elements_by_tag<T: FromGmsh<D>>(&self, tag: i32) -> Result<Vec<T>, MeshError>;

// 复用同一张稠密节点表构造多个组
pub fn node_table(&self) -> Result<NodeTable<D>, MeshError>;
pub fn elements_by_tag_with_table<T: FromGmsh<D>>(
    &self,
    tag: i32,
    nodes: &NodeTable<D>,
) -> Result<Vec<T>, MeshError>;

// 物理组只读视图（按 tag 升序）
pub fn physical_groups(&self) -> impl Iterator<Item = PhysicalGroup<'_>>;
pub fn physical_group(&self, tag: i32) -> Result<PhysicalGroup<'_>, MeshError>;
pub fn physical_group_by_name(&self, name: &str) -> Result<PhysicalGroup<'_>, MeshError>;
pub fn physical_names(&self) -> impl Iterator<Item = (i32, &str)>;
pub fn has_physical_group(&self, name: &str) -> bool;
```

一个物理组内的所有单元共享同一种 GMSH 元素类型，因此 `Vec<T>` 中的 `T` 是单一类型。

`NodeTable<D>` 是稠密节点坐标表（`id == index`）。`GmshMesh::elements` 内部会构建一次；需要同时 materialize 多个物理组时，先调用 `node_table()` 再用 `elements_by_tag_with_table` 复用。

## 用法示例

```rust
use ansatz::element::poi1::Poi1;
use ansatz::element::seg2::Seg2;
use ansatz::element::tri3::Tri3;
use ansatz::io::gmsh::read_gmsh;

let mesh = read_gmsh::<2>("examples/msh/patchtest_tri3.msh").unwrap();

let table = mesh.node_table().unwrap();
let boundary = mesh.physical_group_by_name("Γᵍ").unwrap();
let region = mesh.physical_group_by_name("Ω").unwrap();

let segs: Vec<Seg2<2>> = mesh
    .elements_by_tag_with_table(boundary.tag, &table)
    .unwrap();
let tris: Vec<Tri3<2>> = mesh
    .elements_by_tag_with_table(region.tag, &table)
    .unwrap();

for group in mesh.physical_groups() {
    println!(
        "tag={} name={:?} type={} n={}",
        group.tag, group.name, group.element_type, group.element_count
    );
}
println!("segments: {}", segs.len());
println!("triangles: {}", tris.len());
```

输出：

```text
segments: 40
triangles: 288
```

## 读取 Crouzeix–Raviart（CR）单元

CR 单元的自由度挂在**边**上，GMSH 只给出三角形的顶点，所以 `CRTri3` 不能像 `Tri3` 那样用一个
`FromGmsh` 实现直接读出。`io::gmsh` 为此提供两样东西：`EdgeTable`（边 → 自由度编号）与
`GmshOptions::edges`（**读取时**就把边编号建好）。打开该选项后即可直接取到 CR 单元：

```rust
use ansatz::element::crtri3::CRTri3;
use ansatz::io::gmsh::{GmshOptions, read_gmsh_with};

let mesh = read_gmsh_with::<2>(
    "examples/msh/patchtest_tri3.msh",
    GmshOptions {
        renumber: true, // 节点压紧为 0..n_nodes
        edges: true,    // 读取时给所有三角形的边编号
    },
)?;
let crs: Vec<CRTri3<2>> = mesh.crtri3_elements("Ω")?;

let edges = mesh.edge_numbering().unwrap();
let nodes = mesh.node_table()?;
println!(
    "nodes = {}, cr3 = {}, edges = {} (dofs {}..{})",
    nodes.len(),
    crs.len(),
    edges.len(),
    edges.offset(),
    edges.offset() + edges.len()
);
```

输出：

```text
nodes = 165, cr3 = 288, edges = 452 (dofs 0..452)
```

### 边编号与节点编号是两套独立的编号

CR 的自由度离散的是与顶点场**不同的未知量**，所以边编号与节点编号**不需要错开**：
`GmshOptions { edges: true }` 得到的边 id 就是 `0..n_edges`，与节点 id `0..n_nodes` 可以重叠，
不会冲突。只有当你把 CR 与另一个场（例如协调的 `Tri3` 位移、或拉格朗日乘子）放进**同一个全局
方程组**时，才需要用 `EdgeTable::new(offset)` 把 CR 的自由度整体平移（见下面的「两套接口」）。

需要 `renumber` 的理由与边编号无关：物化任何类型化单元（`elements`、`crtri3_elements`）都要用稠密
的 `node_table()` 取顶点坐标，而它要求节点 id 恰为 `0..n_nodes`。所以对标签不连续的网格仍需
`renumber: true`，否则返回 `MeshError::NodeIndexOutOfRange`（测试
`edge_numbering_works_without_renumbering` 覆盖了这一情形：只开 `edges` 时边编号照常建好，
但 `crtri3_elements` 会因节点表不稠密而失败）。

### 两套接口

| 接口 | 边编号来源 | 用途 |
|---|---|---|
| `mesh.edge_table()` + `crtri3_elements_with_table(name, &nodes, &mut edges)` | 调用方 | 想自己控制编号、或与其它网格/乘子自由度共享一张表 |
| `GmshOptions { edges: true }` + `crtri3_elements(name)` | 读取时（跨所有三角形组统一编号） | 常规用法，`crtri3_elements` 无额外参数 |

其余要点：

- `EdgeTable` 对同一条边只分配一个编号，且与端点顺序无关，所以相邻两个三角形在公共边上拿到
  **同一个** CR 自由度；读取时的编号跨所有三角形组，因此子域之间的界面也是共享的；
- 边界组与体组**共用同一个 `EdgeTable`** 才能保持一致：用 `elements_by_tag_with_table` 式的
  复用路径（`crtri3_elements_by_tag_with_table(tag, &nodes, &mut edges)`）读多个组；
- 组内不是三角形（GMSH 类型 2）时返回 `MeshError::ElementTypeMismatch`，与 `elements::<Tri3>` 一致；
- 未开 `edges` 却调用 `crtri3_elements` 返回 `MeshError::EdgeNumberingMissing`；手工表若缺少某条边
  的编号返回 `MeshError::MissingEdge`；
- 底层实现就是 `Tri3` + `CRTri3::try_from_tri3`（`from_tri3` 的失败可传播版本）。

### 读取 CR 的边界（facet）单元

CR 的边界条件不能用 `Seg2`：`Seg2` 的自由度是顶点值，而 CR 在边上的迹由相邻三角形的**三个**边自由度决定。读取器把边界线段与其所属三角形配对，返回 `CRSeg3<D>`（迹单元）：

```rust
use ansatz::element::crtri3::CRSeg3;
use ansatz::io::gmsh::{GmshOptions, read_gmsh_with};

let mesh = read_gmsh_with::<2>(
    path,
    GmshOptions {
        renumber: true,
        edges: true, // 读取时记录每条边的父三角形
    },
)?;
let boundary: Vec<CRSeg3<2>> = mesh.crtri3_boundary("Γᵈ")?;

// 直接交给现有边界算子：Neumann 用 heat_flux，Dirichlet 用 temperature_penalty
```

要点：

- `CRSeg3::id()` 返回**父三角形的三个边自由度**（与 `crtri3_elements` 的编号一致），因此边界项与内部刚度共享同一套自由度；
- `local_facet()` 是局部 facet 编号：**父顶点 k 不在其上**的那条边（与 `CRTri3` 的自由度顺序相同），参数方向由 facet 顶点按下标升序固定，与 GMSH 段的连线方向无关；
- `jacobe` 是 facet 的长度缩放（`L/2`），所以 `heat_flux`/`temperature_penalty` 的积分面积为真实边界长度；
- 错误：未开 `edges` → `EdgeNumberingMissing`；组内不是线段 → `ElementTypeMismatch`；
  线段不是任何三角形的边（例如 2D 网格里的孤立 1D 组）→ `MissingEdge`；
- 端到端验证见 `crtri3_patchtest_with_dirichlet_penalty`（线性精确解的 patch test）。

### 读取三维 CR（四面体）单元

三维的配套是 `Tet4` → `CRTet4`（内部）与 `CRTri4`（边界迹），对应预处理开关是 `faces`：

```rust
use ansatz::element::crtet4::{CRTet4, CRTri4};
use ansatz::io::gmsh::{GmshOptions, read_gmsh_with};

let mesh = read_gmsh_with::<3>(
    path,
    GmshOptions {
        renumber: true,
        edges: false,
        faces: true, // 读取时给所有四面体的面编号
    },
)?;
let tets: Vec<CRTet4<3>> = mesh.crtet4_elements("Ω")?;      // 4 个面自由度
let boundary: Vec<CRTri4<3>> = mesh.crtet4_boundary("Γ")?;  // 面迹单元，id 与父四面体相同

let faces = mesh.face_numbering().unwrap(); // 0..n_faces
```

要点：

- `FaceTable` 按"三个顶点 id（升序）"登记面，因此共享面在两侧拿到**同一个**自由度；编号跨所有四面体组，与节点/边编号互不相干；
- `crtet4_elements` 需要 `GmshOptions::faces`（否则 `FaceNumberingMissing`），并且只接受 GMSH 类型 4 的组；
- `crtet4_boundary` 把边界三角形（类型 2）与所属四面体配对，返回 `CRTri4`（局部 facet 编号 = 对顶点），不是网格边的线段；
- 2D 与 3D 互不干扰：`crtri3_*` 用 `edges`，`crtet4_*` 用 `faces`；
- `Tet4` 本身（协调单元）可以直接 `mesh.elements::<Tet4<3>>("Ω")` 读取，不需要任何预处理。

## 支持的单元类型## 支持的单元类型

通过 `FromGmsh<D>` trait 桥接 GMSH 元素类型号与本库单元类型：

| GMSH 类型号 | 本库单元 | 拓扑 |
|---|---|---|
| 1 | `Seg2<D>` | 2 节点线 |
| 2 | `Tri3<D>` | 3 节点三角形 |
| 4 | `Tet4<D>` | 4 节点四面体 |
| 15 | `Poi1<D>` | 1 节点点 |

## 新增单元接入

要让一个新单元类型支持 GMSH 读取，只需在 `io` 层为该类型实现 `FromGmsh<D>`：声明 `GMSH_ELEMENT_TYPE` 与 `GMSH_NODE_COUNT` 常量，并实现 `from_gmsh`（根据节点编号列表与节点表构造单元，返回 `Result<Self, MeshError>`）。此后无需改动 `read_gmsh` 本身。

## 文件格式要求

- 必须为 GMSH **v4**（版本号在 `[4.0, 5.0)` 区间），二进制格式不受支持；
- 必须包含 `$Entities` 与 `$Nodes` 段；
- 物理组名可选（用于 `elements("名称")`），缺失时只能通过物理标签访问。

示例与测试用网格位于 `examples/msh/`：`bar1d.msh`（1D）、`patchtest_tri3.msh` / `patchtest_quad4.msh`（2D，含物理组 `Γᵍ`、`Ω`）、`cantilever.msh`、`cook_{8,16,32,64}.msh`、`elasticity_cli.msh`。
