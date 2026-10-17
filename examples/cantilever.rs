// Plane-stress cantilever beam under a transverse tip load
//
// Run with: cargo run --example cantilever
//
// The domain is the rectangle [0, L] x [-H/2, H/2] (width b = 1). The left end
// is clamped by the analytical displacement field (essential boundary
// condition, imposed with the penalty method); the right end carries the
// analytically equilibrating parabolic shear traction (natural boundary
// condition); the top and bottom edges are traction free. See
// `book/src/example-cantilever-model.md` for the derivation of the exact
// solution and `book/src/usage-example-cantilever.md` for a walkthrough.
//
// The finite-element solution is compared against the exact solution through
// the L² and energy-norm error operators in `ansatz::operation::elasticity`,
// and a deformed-shape von Mises stress cloud is rendered with `plotters` to
// `target/cantilever_stress.png` (a second CLI argument overrides the path).

use faer::Side;
use faer::prelude::*;
use faer::sparse::SparseColMat;
use plotters::prelude::*;
use plotters::style::text_anchor::{HPos, Pos, VPos};

use ansatz::element::Element;
use ansatz::element::quad4::Quad4;
use ansatz::element::seg2::Seg2;
use ansatz::io::gmsh::read_gmsh;
use ansatz::operation::elasticity::{
    displacement_penalty, energy_error, lame_constants, l2_error, plane_stress_lambda, stiffness,
    traction,
};
use ansatz::sparse_vector::try_new_from_doublet;

// Geometry (width b = 1).
const L: f64 = 8.0; // beam length
const H: f64 = 1.0; // beam height

// Material (plane stress).
const E: f64 = 1000.0;
const NU: f64 = 0.3;

// Tip shear load.
const P: f64 = 1.0;

/// Second moment of area of the rectangular cross-section,
/// `I = b H³ / 12` with width `b = 1`.
fn second_moment() -> f64 {
    H.powi(3) / 12.0
}

/// Analytical displacement field `(u, v)` of the cantilever beam.
fn exact_displacement(x: [f64; 2]) -> [f64; 2] {
    let (xx, y) = (x[0], x[1]);
    let i = second_moment();
    let u = -(P * y / (6.0 * E * i)) * ((6.0 * L - 3.0 * xx) * xx + (2.0 + NU) * (y * y - H * H / 4.0));
    let v = (P / (6.0 * E * i))
        * (3.0 * NU * y * y * (L - xx) + (4.0 + 5.0 * NU) * H * H * xx / 4.0 + (3.0 * L - xx) * xx * xx);
    [u, v]
}

/// Analytical (plane-stress) strain tensor field at `x`.
fn exact_strain(x: [f64; 2]) -> [[f64; 2]; 2] {
    let (xx, y) = (x[0], x[1]);
    let i = second_moment();
    let exx = -(P / (E * i)) * (L - xx) * y;
    let eyy = NU * (P / (E * i)) * (L - xx) * y;
    // 2 ε_xy = u,_y + v,_x = (P / EI)(1 + ν)(H²/4 − y²).
    let exy = 0.5 * (P / (E * i)) * (1.0 + NU) * (H * H / 4.0 - y * y);
    [[exx, exy], [exy, eyy]]
}

/// Analytical traction on the loaded right end `x = L`, whose outward normal
/// is `(1, 0)`: `t = σ·n = (σ_xx, σ_xy)|_{x=L} = (0, P/(2I)(H²/4 − y²))`.
fn exact_traction(x: [f64; 2]) -> [f64; 2] {
    let y = x[1];
    let sxy = P / (2.0 * second_moment()) * (H * H / 4.0 - y * y);
    [0.0, sxy]
}

// ---------------------------------------------------------------------------
// Visualization with plotters: deformed-shape von Mises stress cloud.
// ---------------------------------------------------------------------------

/// von Mises stress for a plane-stress state.
fn von_mises(sxx: f64, syy: f64, sxy: f64) -> f64 {
    (sxx * sxx - sxx * syy + syy * syy + 3.0 * sxy * sxy).sqrt()
}

/// A simple blue → cyan → green → yellow → red colour map for `t ∈ [0, 1]`.
fn jet(t: f64) -> RGBColor {
    let t = t.clamp(0.0, 1.0);
    let (r, g, b) = if t < 0.25 {
        (0.0, 4.0 * t, 1.0)
    } else if t < 0.5 {
        (0.0, 1.0, 1.0 - 4.0 * (t - 0.25))
    } else if t < 0.75 {
        (4.0 * (t - 0.5), 1.0, 0.0)
    } else {
        (1.0, 1.0 - 4.0 * (t - 0.75), 0.0)
    };
    RGBColor((r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8)
}

/// Plane-stress tensor `(σ_xx, σ_yy, σ_xy)` at a parametric point of an
/// element, from the nodal displacement vector `d`.
fn stress_at(
    elm: &Quad4<2>,
    xi: [f64; 2],
    d: &Col<f64>,
    lambda: f64,
    mu: f64,
) -> (f64, f64, f64) {
    let id = elm.id();
    let g = elm.derivative_shape(xi);
    let (mut exx, mut eyy, mut exy) = (0.0, 0.0, 0.0);
    for a in 0..4 {
        let (ux, uy) = (d[id[a] * 2], d[id[a] * 2 + 1]);
        exx += g[a][0] * ux;
        eyy += g[a][1] * uy;
        exy += 0.5 * (g[a][1] * ux + g[a][0] * uy);
    }
    let sxx = (lambda + 2.0 * mu) * exx + lambda * eyy;
    let syy = lambda * exx + (lambda + 2.0 * mu) * eyy;
    let sxy = 2.0 * mu * exy;
    (sxx, syy, sxy)
}

/// Deformed coordinates of a parametric point: `x + scale · u_h`.
fn deformed_point(elm: &Quad4<2>, xi: [f64; 2], d: &Col<f64>, scale: f64) -> (f64, f64) {
    let id = elm.id();
    let verts = elm.vertices_coordinates();
    let shape = elm.shape(xi);
    let mut x = 0.0;
    let mut y = 0.0;
    for a in 0..4 {
        x += shape[a] * (verts[a][0] + scale * d[id[a] * 2]);
        y += shape[a] * (verts[a][1] + scale * d[id[a] * 2 + 1]);
    }
    (x, y)
}

/// Renders a deformed-shape von Mises stress cloud to a PNG file with plotters.
fn plot_deformed_stress(
    units: &[Quad4<2>],
    d: &Col<f64>,
    lambda: f64,
    mu: f64,
    out: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let scale = 1.0; // deformation exaggeration factor (true deformation)
    let k = 4; // subdivisions per element edge, for a smooth cloud

    // Sub-cell polygons (on the deformed shape) coloured by the von Mises
    // stress at their centre.
    struct Cell {
        pts: [(f64, f64); 4],
        value: f64,
    }
    let mut cells: Vec<Cell> = Vec::with_capacity(units.len() * k * k);
    let (mut xmin, mut xmax) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut ymin, mut ymax) = (f64::INFINITY, f64::NEG_INFINITY);

    for elm in units {
        for i in 0..k {
            for j in 0..k {
                let a0 = -1.0 + 2.0 * i as f64 / k as f64;
                let a1 = -1.0 + 2.0 * (i + 1) as f64 / k as f64;
                let b0 = -1.0 + 2.0 * j as f64 / k as f64;
                let b1 = -1.0 + 2.0 * (j + 1) as f64 / k as f64;
                let corners = [[a0, b0], [a1, b0], [a1, b1], [a0, b1]];
                let mut pts = [(0.0, 0.0); 4];
                for (c, xi) in corners.iter().enumerate() {
                    pts[c] = deformed_point(elm, *xi, d, scale);
                    xmin = xmin.min(pts[c].0);
                    xmax = xmax.max(pts[c].0);
                    ymin = ymin.min(pts[c].1);
                    ymax = ymax.max(pts[c].1);
                }
                let (sxx, syy, sxy) =
                    stress_at(elm, [(a0 + a1) / 2.0, (b0 + b1) / 2.0], d, lambda, mu);
                cells.push(Cell {
                    pts,
                    value: von_mises(sxx, syy, sxy),
                });
            }
        }
    }

    let vmin = cells.iter().map(|c| c.value).fold(f64::INFINITY, f64::min);
    let vmax = cells.iter().map(|c| c.value).fold(f64::NEG_INFINITY, f64::max);
    let vspan = (vmax - vmin).max(f64::MIN_POSITIVE);

    // Data ranges: geometry plus padding, with the right margin reserved for
    // the colour bar.
    let xspan = xmax - xmin;
    let yspan = ymax - ymin;
    let pad_x = 0.03 * xspan;
    let pad_y = 0.06 * yspan;
    let x_lo = xmin - pad_x;
    let x_hi = xmax + 0.18 * xspan;
    let y_lo = ymin - pad_y;
    let y_hi = ymax + 0.14 * yspan; // extra headroom for the title / colour-bar label
    let cb0 = xmax + 0.07 * xspan;
    let cb1 = xmax + 0.11 * xspan;

    // Equal-aspect canvas: the plotting area is proportional to the data spans.
    let plot_w = 1100u32;
    let plot_h = (plot_w as f64 * (y_hi - y_lo) / (x_hi - x_lo)).max(150.0) as u32;

    let root = BitMapBackend::new(out, (plot_w, plot_h)).into_drawing_area();
    root.fill(&WHITE)?;

    {
        let mut chart = ChartBuilder::on(&root).build_cartesian_2d(x_lo..x_hi, y_lo..y_hi)?;

        // Undeformed beam outline, for reference.
        let undef = [
            (0.0, -H / 2.0),
            (L, -H / 2.0),
            (L, H / 2.0),
            (0.0, H / 2.0),
            (0.0, -H / 2.0),
        ];
        chart.draw_series(std::iter::once(PathElement::new(
            undef.to_vec(),
            RGBColor(175, 175, 175).stroke_width(2),
        )))?;

        // Coloured stress cloud (filled sub-cell polygons).
        chart.draw_series(cells.iter().map(|c| {
            let t = (c.value - vmin) / vspan;
            Polygon::new(c.pts.to_vec(), jet(t))
        }))?;

        // Deformed mesh edges.
        for elm in units {
            let verts = elm.vertices_coordinates();
            let id = elm.id();
            let mut edges: Vec<(f64, f64)> = verts
                .iter()
                .enumerate()
                .map(|(a, v)| (v[0] + scale * d[id[a] * 2], v[1] + scale * d[id[a] * 2 + 1]))
                .collect();
            edges.push(edges[0]);
            chart.draw_series(std::iter::once(PathElement::new(
                edges,
                RGBColor(90, 90, 90).stroke_width(1),
            )))?;
        }

        // Frame around the deformed geometry.
        let frame = [
            (xmin - pad_x, ymin - pad_y),
            (xmax + pad_x, ymin - pad_y),
            (xmax + pad_x, ymax + pad_y),
            (xmin - pad_x, ymax + pad_y),
            (xmin - pad_x, ymin - pad_y),
        ];
        chart.draw_series(std::iter::once(PathElement::new(
            frame.to_vec(),
            BLACK.stroke_width(2),
        )))?;

        // Colour bar.
        let n = 120;
        chart.draw_series((0..n).map(|i| {
            let t0 = i as f64 / n as f64;
            let t1 = (i + 1) as f64 / n as f64;
            Rectangle::new([(cb0, ymin + t0 * yspan), (cb1, ymin + t1 * yspan)], jet(t0))
        }))?;
        chart.draw_series(std::iter::once(PathElement::new(
            vec![
                (cb0, ymin),
                (cb1, ymin),
                (cb1, ymax),
                (cb0, ymax),
                (cb0, ymin),
            ],
            BLACK.stroke_width(1),
        )))?;

        // Colour-bar labels: min / mid / max.
        let bar_style = |anchor: Pos| {
            ("sans-serif", 22)
                .into_text_style(&root)
                .pos(anchor)
                .color(&BLACK)
        };
        for (value, y) in [
            (vmin, ymin),
            (0.5 * (vmin + vmax), 0.5 * (ymin + ymax)),
            (vmax, ymax),
        ] {
            chart.draw_series(std::iter::once(Text::new(
                format!("{value:.3}"),
                (cb1 + 0.008 * xspan, y),
                bar_style(Pos::new(HPos::Left, VPos::Center)),
            )))?;
        }
        chart.draw_series(std::iter::once(Text::new(
            "σ_vm",
            (0.5 * (cb0 + cb1), ymax + 0.02 * yspan),
            bar_style(Pos::new(HPos::Center, VPos::Bottom)),
        )))?;

        // Title.
        chart.draw_series(std::iter::once(Text::new(
            format!(
                "Plane-stress cantilever — von Mises stress (deformation ×{scale:.0}, {} elements)",
                units.len()
            ),
            (0.5 * (xmin + xmax), y_hi - 0.005 * yspan),
            bar_style(Pos::new(HPos::Center, VPos::Top)),
        )))?;
    }

    root.present()?;
    Ok(())
}

fn main() {
    // Optional first argument overrides the mesh path (useful for refinement
    // studies); defaults to the committed mesh.
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "examples/msh/cantilever.msh".to_string());
    let mesh = read_gmsh::<2>(&path);
    let units: Vec<Quad4<2>> = mesh.elements("Ω");
    let left: Vec<Seg2<2>> = mesh.elements("Γu");
    let right: Vec<Seg2<2>> = mesh.elements("Γt");

    // 0-based node ids from the GMSH reader; each node owns D = 2 DOFs.
    let max_id = mesh.nodes.keys().copied().max().unwrap();
    let n_dofs = (max_id + 1) * 2;

    // Plane-stress material: use λ* = Eν/(1 − ν²) for the in-plane volumetric
    // coupling, while μ is unchanged.
    let lambda = plane_stress_lambda(E, NU);
    let (_, mu) = lame_constants(E, NU);

    // ---- assembly ----
    let mut triplets = stiffness(&units, lambda, mu);
    let mut doublets = Vec::new();

    // Left end: exact displacement (essential boundary condition, penalty).
    let alpha = 1e7;
    let (pt, pd) = displacement_penalty(&left, alpha, exact_displacement);
    triplets.extend(pt);
    doublets.extend(pd);

    // Right end: exact shear traction (natural boundary condition).
    doublets.extend(traction(&right, exact_traction));

    // ---- solve ----
    let k = SparseColMat::<usize, f64>::try_new_from_triplets(n_dofs, n_dofs, &triplets).unwrap();
    let f = try_new_from_doublet(n_dofs, &doublets);
    let llt = k.sp_cholesky(Side::Lower).unwrap();
    let d = llt.solve(&f);

    // ---- error against the analytical solution ----
    let l2 = l2_error(&units, |id, p| d[id * 2 + p], exact_displacement);
    let energy = energy_error(
        &units,
        |id, p| d[id * 2 + p],
        exact_strain,
        lambda,
        mu,
    );

    // Tip deflection at (L, 0) for reference.
    let tip_id = mesh
        .nodes
        .iter()
        .find(|(_, node)| (node.coordinates[0] - L).abs() < 1e-9 && node.coordinates[1].abs() < 1e-9)
        .map(|(id, _)| *id)
        .expect("tip node at (L, 0) not found");
    let tip = [d[tip_id * 2], d[tip_id * 2 + 1]];
    let tip_exact = exact_displacement([L, 0.0]);

    println!("== plane-stress cantilever beam ==");
    println!("  mesh: {} quadrilateral elements, {} DOFs", units.len(), n_dofs);
    println!(
        "  tip deflection  FE = ({:.6e}, {:.6e}), exact = ({:.6e}, {:.6e})",
        tip[0], tip[1], tip_exact[0], tip_exact[1]
    );
    println!("  L2 error        : {l2:.6e}");
    println!("  energy-norm error: {energy:.6e}");

    // The exact field is cubic, so it is not representable by bilinear
    // elements and the errors are not machine-zero; they decrease under mesh
    // refinement (see `book/src/usage-example-cantilever.md`).
    assert!(l2 < 5e-2, "L2 error too large: {l2:e}");
    assert!(energy < 2e-1, "energy-norm error too large: {energy:e}");

    // ---- deformed-shape von Mises stress cloud (plotters) ----
    let out = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "target/cantilever_stress.png".to_string());
    if let Some(parent) = std::path::Path::new(&out).parent() {
        std::fs::create_dir_all(parent).ok();
    }
    plot_deformed_stress(&units, &d, lambda, mu, &out).expect("failed to render the stress cloud");
    println!("  stress cloud    : {out}");
}
