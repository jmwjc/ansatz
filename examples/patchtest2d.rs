// 2D steady-state heat-conduction patch test
//
// Run with: cargo run --example patchtest2d
//
// The patch test verifies that the finite-element assembly exactly reproduces
// a linear temperature field T(x, y) = x + y on the unit square [0, 1]², for
// both a triangular mesh and a quadrilateral mesh read directly from GMSH.
// In each mesh, the physical group `Ω` is the interior and `Γᵍ` the boundary,
// where the Dirichlet temperature is prescribed. Accuracy is measured in the
// L² and H¹ (semi-norm) error.

use faer::prelude::*;
use faer::sparse::SparseColMat;
use faer::Side;

use ansatz::element::quad4::Quad4;
use ansatz::element::seg2::Seg2;
use ansatz::element::tri3::Tri3;
use ansatz::element::{Element, FullIntegration};
use ansatz::io::gmsh::{read_gmsh, FromGmsh};
use ansatz::operation::heat::{conduction_stiffness, h1_semierror, l2_error, temperature_penalty};
use ansatz::sparse_vector::try_new_from_doublet;

/// Runs a steady-state patch test on a mesh read from GMSH: assembles the
/// conduction stiffness, imposes the prescribed boundary temperature (linear
/// field `T = x + y`) by the penalty method, solves, and measures the L² and
/// H¹ (semi-norm) error against the exact field.
///
/// Degrees of freedom are numbered by the 0-based node ids assigned by
/// [`read_gmsh`] (GMSH node ids decremented by one).
fn run_patch_test<T, const N: usize, const G: usize>(path: &str, label: &str)
where
    T: FromGmsh<2> + Element<2, 2, N> + FullIntegration<2, G>,
{
    let mesh = read_gmsh::<2>(path);
    let units: Vec<T> = mesh.elements("Ω");
    let segs: Vec<Seg2<2>> = mesh.elements("Γᵍ");

    let max_id = mesh.nodes.keys().copied().max().unwrap();
    let n_dofs = max_id + 1;

    let exact = |x: [f64; 2]| x[0] + x[1];
    let exact_grad = |_: [f64; 2]| [1.0, 1.0];
    let alpha = 1e7;

    let mut triplets = conduction_stiffness(&units, 1.0);
    let mut doublets = Vec::new();

    // Prescribed boundary temperature, integrated over the boundary segments.
    let (ptt, ptd) = temperature_penalty(&segs, alpha, exact);
    triplets.extend(ptt);
    doublets.extend(ptd);

    let kk = SparseColMat::<usize, f64>::try_new_from_triplets(n_dofs, n_dofs, &triplets).unwrap();
    let f = try_new_from_doublet(n_dofs, &doublets);
    let llt = kk.sp_cholesky(Side::Lower).unwrap();
    let d = llt.solve(&f);

    let l2 = l2_error(&units, |id| d[id], &exact);
    let h1 = h1_semierror(&units, |id| d[id], &exact_grad);

    println!(
        "  {}: {} elements, L2 error: {:.3e}, H1 semi-norm error: {:.3e}",
        label,
        units.len(),
        l2,
        h1
    );
    assert!(
        l2 < 1e-6 && h1 < 1e-6,
        "{label} patch test failed: L2 = {l2:.3e}, H1 = {h1:.3e}"
    );
}

fn main() {
    println!("== patch test: triangular mesh (GMSH) ==");
    run_patch_test::<Tri3<2>, 3, 3>("examples/msh/patchtest_tri3.msh", "triangles");

    println!("== patch test: quadrilateral mesh (GMSH) ==");
    run_patch_test::<Quad4<2>, 4, 4>("examples/msh/patchtest_quad4.msh", "quads");

    println!("patch test passed for both triangular and quadrilateral meshes.");
}
