// Cook's membrane (plane strain): full integration vs selective reduced
// integration (SRI), for a compressible and a nearly incompressible material.
//
// Run with: cargo run --example cook
//
// The trapezoidal membrane is clamped on the left edge and loaded by a uniform
// shear traction on the right edge (resultant 1 N). The reference quantity is
// the vertical displacement at the top-right corner (48, 60).
//
// `full` assembles the Q4 stiffness with full 2x2 Gauss integration. `SRI`
// (selective reduced integration) assembles the dilatational part with reduced
// integration and the deviatoric part with full integration, which mitigates
// volumetric locking. Two material cases are studied:
//   * nu = 1/3      (compressible)         benchmark reference 23.9
//   * nu = 0.49999  (nearly incompressible) converged reference ~19.4
// The script prints convergence tables and renders a plot for each case.
// See `book/src/example-cook.md`.

use faer::Side;
use faer::prelude::*;
use faer::sparse::SparseColMat;
use plotters::prelude::*;

use ansatz::element::quad4::Quad4;
use ansatz::element::seg2::Seg2;
use ansatz::io::gmsh::read_gmsh;
use ansatz::operation::elasticity::{
    bulk_modulus, deviatoric_stiffness, dilatational_stiffness, displacement_penalty,
    lame_constants, stiffness, traction,
};
use ansatz::sparse_vector::try_new_from_doublet;

const E: f64 = 1.0;
const TOTAL_LOAD: f64 = 1.0;
const EDGE_LEN: f64 = 16.0; // loaded right edge length
const REF_POINT: [f64; 2] = [48.0, 60.0];
const ALPHA: f64 = 1e7;

/// Solves the Cook membrane on one mesh for a given Poisson ratio and returns
/// `(reference vertical displacement, elements, DOFs)`.
///
/// `sri = false` uses the fully integrated Q4 stiffness; `sri = true` uses
/// selective reduced integration (dilatational reduced + deviatoric full).
fn solve(path: &str, sri: bool, nu: f64) -> (f64, usize, usize) {
    let mesh = read_gmsh::<2>(path);
    let units: Vec<Quad4<2>> = mesh.elements("Ω");
    let left: Vec<Seg2<2>> = mesh.elements("Γu");
    let right: Vec<Seg2<2>> = mesh.elements("Γt");

    let max_id = mesh.nodes.keys().copied().max().unwrap();
    let n_dofs = (max_id + 1) * 2;

    // Plane strain: λ, μ from the 3D Lamé constants.
    let (lambda, mu) = lame_constants(E, nu);

    let mut triplets = if sri {
        let bulk = bulk_modulus::<2>(lambda, mu);
        let mut t = dilatational_stiffness(&units, bulk);
        t.extend(deviatoric_stiffness(&units, mu));
        t
    } else {
        stiffness(&units, lambda, mu)
    };

    // Left edge clamped: u = v = 0 (penalty method).
    let (pt, mut doublets) = displacement_penalty(&left, ALPHA, |_| [0.0, 0.0]);
    triplets.extend(pt);

    // Right edge: uniform shear traction in +y, resultant 1 N.
    doublets.extend(traction(&right, |_| [0.0, TOTAL_LOAD / EDGE_LEN]));

    let k = SparseColMat::<usize, f64>::try_new_from_triplets(n_dofs, n_dofs, &triplets).unwrap();
    let f = try_new_from_doublet(n_dofs, &doublets);
    let llt = k.sp_cholesky(Side::Lower).unwrap();
    let d = llt.solve(&f);

    let ref_id = mesh
        .nodes
        .iter()
        .find(|(_, node)| {
            (node.coordinates[0] - REF_POINT[0]).abs() < 1e-9
                && (node.coordinates[1] - REF_POINT[1]).abs() < 1e-9
        })
        .map(|(id, _)| *id)
        .expect("reference node (48, 60) not found");

    (d[ref_id * 2 + 1], units.len(), n_dofs)
}

/// Renders a full-vs-SRI convergence plot to a PNG file.
fn plot_convergence(
    full: &[(f64, f64)],
    sri: &[(f64, f64)],
    title: &str,
    y_bounds: (f64, f64),
    reference: Option<(f64, String)>,
    out: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let (x_lo, x_hi) = (0.0_f64, 64.0_f64); // element count per edge
    let (y_lo, y_hi) = y_bounds;
    let root = BitMapBackend::new(out, (1000, 640)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 24))
        .margin(20)
        .x_label_area_size(55)
        .y_label_area_size(70)
        .build_cartesian_2d(x_lo..x_hi, y_lo..y_hi)?;

    chart
        .configure_mesh()
        .x_desc("elements per edge n")
        .y_desc("vertical displacement")
        .x_labels(9)
        .x_label_formatter(&|x| format!("{x:.0}"))
        .y_label_formatter(&|y| format!("{y:.1}"))
        .draw()?;

    // Optional reference line.
    if let Some((value, label)) = reference {
        chart
            .draw_series(std::iter::once(PathElement::new(
                vec![(x_lo, value), (x_hi, value)],
                RGBColor(120, 120, 120).stroke_width(2),
            )))?
            .label(label)
            .legend(|(x, y)| {
                PathElement::new(
                    vec![(x, y), (x + 24, y)],
                    RGBColor(120, 120, 120).stroke_width(2),
                )
            });
    }

    let line = |points: &[(f64, f64)], color: RGBColor| {
        LineSeries::new(points.iter().copied(), color.stroke_width(2))
    };

    chart
        .draw_series(line(full, BLUE))?
        .label("full integration")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 24, y)], BLUE.stroke_width(2)));
    chart.draw_series(full.iter().map(|p| Circle::new(*p, 6, BLUE.filled())))?;

    chart
        .draw_series(line(sri, RED))?
        .label("SRI (dilatational reduced)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 24, y)], RED.stroke_width(2)));
    chart.draw_series(sri.iter().map(|p| Circle::new(*p, 6, RED.filled())))?;

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.85))
        .border_style(BLACK)
        .draw()?;

    root.present()?;
    Ok(())
}

fn main() {
    std::fs::create_dir_all("target").ok();

    let meshes = [
        ("examples/msh/cook_8.msh", 8u32),
        ("examples/msh/cook_16.msh", 16u32),
        ("examples/msh/cook_32.msh", 32u32),
        ("examples/msh/cook_64.msh", 64u32),
    ];

    // (label, ν, y-range, reference line, output file)
    let cases = [
        (
            "nu = 1/3 (compressible)",
            1.0 / 3.0,
            (19.5, 24.5),
            Some((23.9, "benchmark 23.9".to_string())),
            "target/cook_convergence.png",
        ),
        (
            "nu = 0.49999 (nearly incompressible)",
            0.49999,
            (4.0, 21.0),
            Some((19.4, "converged ~ 19.4".to_string())),
            "target/cook_incompressible.png",
        ),
    ];

    for (label, nu, y_bounds, reference, out) in cases {
        println!("== Cook's membrane (plane strain, E = {E}, {label}) ==");
        println!(
            "{:>4} {:>9} {:>7} {:>10} {:>10}",
            "n", "elements", "DOFs", "v_full", "v_sri"
        );

        let mut full_points = Vec::new();
        let mut sri_points = Vec::new();
        for (path, n) in meshes {
            let (v_full, ne, nd) = solve(path, false, nu);
            let (v_sri, _, _) = solve(path, true, nu);
            println!("{n:>4} {ne:>9} {nd:>7} {v_full:>10.4} {v_sri:>10.4}");
            full_points.push((n as f64, v_full));
            sri_points.push((n as f64, v_sri));
        }

        let title = format!("Cook's membrane ({label}): vertical displacement at (48, 60)");
        plot_convergence(&full_points, &sri_points, &title, y_bounds, reference, out)
            .expect("failed to render the convergence plot");
        println!("  convergence plot: {out}\n");
    }
}
