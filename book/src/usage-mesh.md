# GMSH 网格读取

`rustfem` 的 `io` 模块支持读取 **GMSH v4.1 ASCII** 格式的 `.msh` 文件。

## 入口函数

```rust
use rustfem::io::gmsh::read_gmsh;

let mesh = read_gmsh::<2>("test/msh/patchtest.msh");
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
pub fn elements<T: FromGmsh<D>>(&self, physical_name: &str) -> Vec<T>;

// 按物理标签构造
pub fn elements_by_tag<T: FromGmsh<D>>(&self, tag: i32) -> Vec<T>;
```

一个物理组内的所有单元共享同一种 GMSH 元素类型，因此 `Vec<T>` 中的 `T` 是单一类型。

## 用法示例

```rust
use rustfem::element::poi1::Poi1;
use rustfem::element::seg2::Seg2;
use rustfem::element::tri3::Tri3;
use rustfem::io::gmsh::read_gmsh;

let mesh = read_gmsh::<2>("test/msh/patchtest.msh");

let segs: Vec<Seg2<2>> = mesh.elements("Γᵍ");   // 线单元（物理组名 Γᵍ）
let tris: Vec<Tri3<2>> = mesh.elements("Ω");    // 三角形单元（物理组名 Ω）

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

要让一个新单元类型支持 GMSH 读取，只需在 `io` 层为该类型实现 `FromGmsh<D>`：声明 `GMSH_ELEMENT_TYPE` 常量，并实现 `from_gmsh`（根据节点编号列表与节点表构造单元，找不到节点时返回 `None`）。此后无需改动 `read_gmsh` 本身。

## 文件格式要求

- 必须为 GMSH **v4**（版本号在 `[4.0, 5.0)` 区间），二进制格式不受支持；
- 必须包含 `$Entities` 与 `$Nodes` 段；
- 物理组名可选（用于 `elements("名称")`），缺失时只能通过物理标签访问。

测试网格位于 `test/msh/`：`patchtest.msh`（2D，含物理组 `Γᵍ`、`Ω`）与 `patchtest1D.msh`（1D）。
