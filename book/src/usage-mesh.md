# GMSH 网格读取

`ansatz` 的 `io` 模块支持读取 **GMSH v4.1 ASCII** 格式的 `.msh` 文件。

## 入口函数

```rust
use ansatz::io::gmsh::read_gmsh;

let mesh = read_gmsh::<2>("examples/msh/patchtest_tri3.msh").unwrap();
```

`read_gmsh::<D>` 解析文件并返回中间产物 `GmshMesh<D>`，其中 `D` 是期望的空间维数。它**只读文件、不构造具体单元类型**；类型化单元向量按需由 `GmshMesh::elements` 构造。

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

## 支持的单元类型

通过 `FromGmsh<D>` trait 桥接 GMSH 元素类型号与本库单元类型：

| GMSH 类型号 | 本库单元 | 拓扑 |
|---|---|---|
| 1 | `Seg2<D>` | 2 节点线 |
| 2 | `Tri3<D>` | 3 节点三角形 |
| 15 | `Poi1<D>` | 1 节点点 |

## 新增单元接入

要让一个新单元类型支持 GMSH 读取，只需在 `io` 层为该类型实现 `FromGmsh<D>`：声明 `GMSH_ELEMENT_TYPE` 与 `GMSH_NODE_COUNT` 常量，并实现 `from_gmsh`（根据节点编号列表与节点表构造单元，返回 `Result<Self, MeshError>`）。此后无需改动 `read_gmsh` 本身。

## 文件格式要求

- 必须为 GMSH **v4**（版本号在 `[4.0, 5.0)` 区间），二进制格式不受支持；
- 必须包含 `$Entities` 与 `$Nodes` 段；
- 物理组名可选（用于 `elements("名称")`），缺失时只能通过物理标签访问。

示例与测试用网格位于 `examples/msh/`：`bar1d.msh`（1D）、`patchtest_tri3.msh` / `patchtest_quad4.msh`（2D，含物理组 `Γᵍ`、`Ω`）、`cantilever.msh`、`cook_{8,16,32,64}.msh`、`elasticity_cli.msh`。
