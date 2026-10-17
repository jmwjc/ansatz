# 1D 桁架杆（bar1d）

`examples/bar1d.rs` 是库的最小端到端示例：手工构造一个 11 节点 / 10 单元的 1D 杆，施加约束与载荷，罚函数法加边界条件后稀疏 Cholesky 求解。下面逐段讲解。（需要读取 GMSH 网格的示例见 `patchtest2d`、`cantilever`、`cook`、`elasticity_cli`。）

## 1. 性能剖析的全局分配器

```rust
#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;
```

仅在 `dhat-heap` feature 启用时生效，替换全局分配器以记录堆分配。

## 2. 构造 1D 杆

`bar1d` 不依赖网格文件，直接按 `xᵢ = i/10` 生成节点与 `Seg2<1>` 单元（GMSH 读取的用法见「GMSH 网格读取」一章与 `patchtest2d` 示例）。

```rust
let np = 11;          // 11 个节点
let ne = np - 1;      // 10 个单元
let bars: Vec<Seg2<1>> = (0..ne)
    .map(|i| {
        Seg2::new(
            Node::new(i, [i as f64 / ne as f64]),
            Node::new(i + 1, [(i + 1) as f64 / ne as f64]),
        )
    })
    .collect();
```

把区间 `[0, 1]` 均分为 10 个 `Seg2<1>` 单元，节点坐标 $x_i = i/10$。

```rust
let left  = Poi1::new(Node::new(0, [0.0]));     // 左端 x = 0
let right = Poi1::new(Node::new(np - 1, [1.0])); // 右端 x = 1
```

两个点单元分别标记左端（约束点）与右端（加载点）。

## 3. 刚度组装 + 罚函数约束

```rust
let mut triplets = stiffness(&bars, 1.0);       // EA = 1

let alpha: f64 = 1e7;
let (penalty_triplets, penalty_doublets) =
    displacement_penalty(&vec![left], alpha, 0.0);
triplets.extend(penalty_triplets);

let k = SparseColMat::<usize, f64>::try_new_from_triplets(np, np, &triplets).unwrap();
```

- `stiffness(&bars, 1.0)`：组装单元刚度（$EA = 1$）；
- `displacement_penalty`：左端施加固定约束 $u(0) = 0$，罚数 $α = 1e7$；
- 罚刚度并入总刚度后，`try_new_from_triplets` 合并为 `11 × 11` 稀疏矩阵。

## 4. 载荷组装

```rust
let mut doublets = body_force(&bars, |_| 1.0);   // 自重：均布体力 f = 1
doublets.extend(traction(&vec![right], 1.0));     // 右端集中拉力 P = 1
doublets.extend(penalty_doublets);                // 罚载荷
let f = try_new_from_doublet(np, &doublets);
```

三种载荷依次叠加：

1. 均布体力（自重），密度 $f = 1$；
2. 右端集中拉力 $P = 1$；
3. 罚函数法的等效载荷项（本例 $ū = 0$，故该项为零，但代码仍显式叠加）。

## 5. 求解

```rust
let llt = k.sp_cholesky(Side::Lower).unwrap();
let d = llt.solve(&f);

println!("bar displacement: {:?}", d);
```

刚度矩阵对称正定，做稀疏 Cholesky 分解后回代，得到位移向量。

## 结果解读

输出位移向量最后一个分量约 `1.5`。由理论部分可知，该杆（$EA = 1$、$L = 1$、均布体力 $f = 1$、右端集中力 $P = 1$）的解析解为

$$
u(x) = -\frac{x^2}{2} + 2x
$$

故 $u(1) = 1.5$，与数值结果一致。
