// 1D bar example
//
// Run with: cargo run --example bar1d

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use faer::Side;
use faer::prelude::*;
use faer::sparse::SparseColMat;

use rustfem::element::poi1::Poi1;
use rustfem::element::seg2::Seg2;
use rustfem::node::Node;
use rustfem::operation::bar::{body_force, displacement_penalty, stiffness, traction};
use rustfem::sparse_vector::try_new_from_doublet;
use std::time::Instant;

fn main() {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let start = Instant::now();

    // ---------- 1D truss assembly + solve ----------
    // Build a 1D bar made of 10 Seg2<1> elements, fixed at the left end,
    // subject to its own weight (uniform body force f = 1) and a concentrated
    // tensile load P = 1 at the right end.
    let np = 11;
    let ne = np - 1;
    let bars: Vec<Seg2<1>> = (0..ne)
        .map(|i| {
            Seg2::new(
                Node::new(i, [i as f64 / ne as f64]),
                Node::new(i + 1, [(i + 1) as f64 / ne as f64]),
            )
        })
        .collect();

    // Constraint point at the left end and load point at the right end (point elements).
    let left = Poi1::new(Node::new(0, [0.0]));
    let right = Poi1::new(Node::new(np - 1, [1.0]));

    let mut triplets = stiffness(&bars, 1.0);

    // Left end fixed: impose the displacement constraint u(0) = 0 via the penalty method.
    let alpha: f64 = 1e7;
    let (penalty_triplets, penalty_doublets) = displacement_penalty(&vec![left], alpha, 0.0);
    triplets.extend(penalty_triplets);

    let k = SparseColMat::<usize, f64>::try_new_from_triplets(np, np, &triplets).unwrap();

    // Uniform body force (self-weight): natural boundary condition (body force).
    let mut doublets = body_force(&bars, |_| 1.0);
    // Unit tensile load at the right end: natural boundary condition (concentrated force).
    doublets.extend(traction(&vec![right], 1.0));
    // Penalty load.
    doublets.extend(penalty_doublets);
    let f = try_new_from_doublet(np, &doublets);

    let llt = k.sp_cholesky(Side::Lower).unwrap();
    let d = llt.solve(&f);

    println!("bar displacement: {d:?}");

    println!("Time elapsed: {:?}", start.elapsed());
}
