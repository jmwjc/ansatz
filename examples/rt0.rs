// Raviart–Thomas (RT0) shape functions on the reference triangle.
//
// Run with: cargo run --example rt0
//
// The RT0 basis on the unit triangle (0,0), (1,0), (0,1) is
//
//   q1 = (x, y),   q2 = (x - 1, y),   q3 = (x, y - 1),
//
// i.e. q_i = (x - A_i) / (2 |T|) with A_i the vertex opposite edge i. The dof of
// edge j is the normal flux Φ_j(q) = ∫_{Γ_j} q·n_j ds, so Φ_j(q_i) = δ_ij, and
// every basis field has the constant divergence ∇·q_i = 1/|T| = 2.
//
// Two figures are written:
//   * target/rt0_shape_functions.png — one panel per basis field: the vector
//     field drawn as arrows (colour = |q_i|);
//   * target/rt0_shape_components.png — one panel per basis field and component
//     (x and y), with a heat map, iso-lines of the affine component, and a
//     per-panel colour bar.
//
// See `book/src/theory-rt.md` and `src/element/rttri3.rs`.

use plotters::prelude::*;
use plotters::style::text_anchor::{HPos, Pos, VPos};

use ansatz::element::HdivElement;
use ansatz::element::rttri3::RTTri3;
use ansatz::element::tri3::Tri3;
use ansatz::node::Node;

/// Vertices of the reference triangle, indexed like the RT0 edge dofs (edge `i`
/// is opposite vertex `i`).
const TRI: [(f64, f64); 3] = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];

fn main() {
    std::fs::create_dir_all("target").ok();

    let rt = reference_rt();
    let area = rt.area();
    println!(
        "RT0 on the reference triangle (area = {area}, 1/|T| = {})",
        1.0 / area
    );
    println!("  q1 = (x, y)      q2 = (x - 1, y)      q3 = (x, y - 1)");
    println!("  Φ_j(q_i) = δ_ij,   ∇·q_i = 1/|T| = {}", 1.0 / area);

    // Numerical check of the defining property of the dofs.
    let mut dof_error = 0.0_f64;
    for (i, opposite) in rt.triangle_vertices().iter().enumerate() {
        let q_i = |x: [f64; 2]| {
            [
                (x[0] - opposite[0]) / (2.0 * area),
                (x[1] - opposite[1]) / (2.0 * area),
            ]
        };
        for j in 0..3 {
            let phi = rt.flux_dof(j, q_i).expect("valid edge index");
            let expected = if i == j { 1.0 } else { 0.0 };
            dof_error = dof_error.max((phi - expected).abs());
        }
    }
    println!("  max |Φ_j(q_i) - δ_ij| = {dof_error:.1e}\n");

    plot_shape_functions(&rt, "target/rt0_shape_functions.png")
        .expect("failed to render the RT0 shape functions");
    println!("  vector field figure:   target/rt0_shape_functions.png");

    plot_components(&rt, "target/rt0_shape_components.png")
        .expect("failed to render the RT0 shape function components");
    println!("  component figure:      target/rt0_shape_components.png");
}

/// RT0 element of the reference triangle, with edge dofs 0, 1, 2.
fn reference_rt() -> RTTri3<2> {
    let tri = Tri3::new(
        Node::new(0, [TRI[0].0, TRI[0].1]),
        Node::new(1, [TRI[1].0, TRI[1].1]),
        Node::new(2, [TRI[2].0, TRI[2].1]),
    );
    let mut next = 0;
    RTTri3::from_tri3(&tri, |_, _| {
        let id = next;
        next += 1;
        id
    })
}

/// Formula of basis field `index`, for the panel titles.
fn basis_label(index: usize) -> &'static str {
    ["q1 = (x, y)", "q2 = (x - 1, y)", "q3 = (x, y - 1)"][index]
}

/// Formula of one component of a basis field, for the panel titles.
fn component_label(index: usize, component: usize) -> &'static str {
    [
        ["q1_x = x", "q1_y = y"],
        ["q2_x = x - 1", "q2_y = y"],
        ["q3_x = x", "q3_y = y - 1"],
    ][index][component]
}

/// Value of component `component` of basis field `index` at `point`.
fn component(rt: &RTTri3<2>, index: usize, component: usize, point: (f64, f64)) -> f64 {
    rt.shape([point.0, point.1])[index][component]
}

/// Closed outline of the reference triangle.
fn triangle_outline() -> Vec<(f64, f64)> {
    vec![TRI[0], TRI[1], TRI[2], TRI[0]]
}

/// Sequential colour map: dark blue → teal → sand → dark red for `t ∈ [0, 1]`.
fn sequential(t: f64) -> RGBColor {
    let t = t.clamp(0.0, 1.0);
    let stops: [(f64, (u8, u8, u8)); 4] = [
        (0.0, (31, 51, 153)),
        (0.35, (60, 170, 170)),
        (0.70, (245, 214, 96)),
        (1.0, (192, 45, 38)),
    ];
    for pair in stops.windows(2) {
        let ((t0, c0), (t1, c1)) = (pair[0], pair[1]);
        if t <= t1 {
            let s = (t - t0) / (t1 - t0);
            let mix = |a: u8, b: u8| (a as f64 + s * (b as f64 - a as f64)).round() as u8;
            return RGBColor(mix(c0.0, c1.0), mix(c0.1, c1.1), mix(c0.2, c1.2));
        }
    }
    let (_, c) = stops[3];
    RGBColor(c.0, c.1, c.2)
}

/// Intersection of the level set `component = level` with the triangle.
///
/// The component is affine on the triangle, so its level sets are straight
/// lines; intersecting with the three edges yields the segment to draw.
fn contour_segment(values: [f64; 3], level: f64) -> Option<((f64, f64), (f64, f64))> {
    let mut points: Vec<(f64, f64)> = Vec::new();
    for (a, b) in [(0, 1), (1, 2), (2, 0)] {
        let (va, vb) = (values[a], values[b]);
        if va != vb && (va - level) * (vb - level) <= 0.0 {
            let t = (level - va) / (vb - va);
            let point = (
                TRI[a].0 + t * (TRI[b].0 - TRI[a].0),
                TRI[a].1 + t * (TRI[b].1 - TRI[a].1),
            );
            if !points
                .iter()
                .any(|q| (q.0 - point.0).abs() < 1e-9 && (q.1 - point.1).abs() < 1e-9)
            {
                points.push(point);
            }
        }
    }
    (points.len() >= 2).then(|| (points[0], points[1]))
}

/// One panel per basis field: arrows of `q_i`, coloured by `|q_i|`.
fn plot_shape_functions(rt: &RTTri3<2>, out: &str) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(out, (1620, 540)).into_drawing_area();
    root.fill(&WHITE)?;
    let panels = root.split_evenly((1, 3));

    let max_magnitude = TRI
        .iter()
        .flat_map(|p| rt.shape([p.0, p.1]))
        .map(|q| (q[0] * q[0] + q[1] * q[1]).sqrt())
        .fold(0.0_f64, f64::max);
    let arrow_scale = 0.20 / max_magnitude;

    for (index, cell) in panels.iter().enumerate() {
        let panel = cell.margin(6, 6, 6, 6);
        let title_style = ("sans-serif", 26)
            .into_text_style(&panel)
            .pos(Pos::new(HPos::Center, VPos::Top))
            .color(&BLACK);
        let note_style = ("sans-serif", 17)
            .into_text_style(&panel)
            .pos(Pos::new(HPos::Center, VPos::Top))
            .color(&RGBColor(90, 90, 90));

        let mut chart = ChartBuilder::on(&panel).build_cartesian_2d(-0.25..1.25, -0.25..1.25)?;

        // Light background of the element.
        chart.draw_series(std::iter::once(Polygon::new(
            triangle_outline(),
            RGBColor(246, 244, 240),
        )))?;

        // Arrows on a barycentric lattice.
        let n = 8;
        let mut shafts = Vec::new();
        let mut heads = Vec::new();
        for i in 0..=n {
            for j in 0..=(n - i) {
                let point = (i as f64 / n as f64, j as f64 / n as f64);
                let q = rt.shape([point.0, point.1])[index];
                let magnitude = (q[0] * q[0] + q[1] * q[1]).sqrt();
                if magnitude < 1e-12 {
                    continue;
                }
                let color = sequential(magnitude / max_magnitude);
                let tip = (point.0 + arrow_scale * q[0], point.1 + arrow_scale * q[1]);
                let (dx, dy) = (tip.0 - point.0, tip.1 - point.1);
                let length = (dx * dx + dy * dy).sqrt();
                let (ux, uy) = (dx / length, dy / length);
                let base = (tip.0 - 0.040 * ux, tip.1 - 0.040 * uy);
                let left = (base.0 - 0.015 * uy, base.1 + 0.015 * ux);
                let right = (base.0 + 0.015 * uy, base.1 - 0.015 * ux);
                shafts.push((point, tip, color));
                heads.push((vec![tip, left, right], color));
            }
        }
        chart.draw_series(
            shafts
                .iter()
                .map(|(from, to, color)| PathElement::new(vec![*from, *to], color.stroke_width(2))),
        )?;
        chart.draw_series(
            heads
                .iter()
                .map(|(points, color)| Polygon::new(points.clone(), color.filled())),
        )?;

        // Element outline.
        chart.draw_series(std::iter::once(PathElement::new(
            triangle_outline(),
            BLACK.stroke_width(2),
        )))?;

        chart.draw_series(std::iter::once(Text::new(
            basis_label(index),
            (0.5, 1.17),
            title_style,
        )))?;
        chart.draw_series(std::iter::once(Text::new(
            format!(
                "div q{} = {:.0},  Φ{}(q{}) = 1,  colour = |q{}|",
                index + 1,
                1.0 / rt.area(),
                index + 1,
                index + 1,
                index + 1
            ),
            (0.5, 1.09),
            note_style,
        )))?;
    }

    root.present()?;
    Ok(())
}

/// One panel per (basis field, component): heat map + iso-lines + colour bar.
fn plot_components(rt: &RTTri3<2>, out: &str) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(out, (1200, 1800)).into_drawing_area();
    root.fill(&WHITE)?;
    let panels = root.split_evenly((3, 2));
    let n = 48; // triangular sub-cells per edge

    for index in 0..3 {
        for comp in 0..2 {
            let panel = panels[2 * index + comp].margin(6, 6, 6, 6);
            let title_style = ("sans-serif", 26)
                .into_text_style(&panel)
                .pos(Pos::new(HPos::Center, VPos::Top))
                .color(&BLACK);
            let tick_style = ("sans-serif", 15)
                .into_text_style(&panel)
                .color(&RGBColor(80, 80, 80));

            let mut chart = ChartBuilder::on(&panel).build_cartesian_2d(-0.2..1.2, -0.2..1.2)?;

            // Vertex values of this affine component: min/max and iso-lines.
            let values: [f64; 3] = std::array::from_fn(|k| component(rt, index, comp, TRI[k]));
            let min = values.iter().copied().fold(f64::INFINITY, f64::min);
            let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let span = (max - min).max(f64::MIN_POSITIVE);

            // Heat map on a triangular sub-grid (each cell is a small triangle
            // whose vertices lie on the lattice, so the element is tiled exactly).
            let mut cells = Vec::new();
            for i in 0..n {
                for j in 0..(n - i) {
                    let x0 = i as f64 / n as f64;
                    let y0 = j as f64 / n as f64;
                    let x1 = (i + 1) as f64 / n as f64;
                    let y1 = (j + 1) as f64 / n as f64;
                    let lower = vec![(x0, y0), (x1, y0), (x0, y1)];
                    let t = (component(
                        rt,
                        index,
                        comp,
                        ((x0 + x1 + x0) / 3.0, (y0 + y0 + y1) / 3.0),
                    ) - min)
                        / span;
                    cells.push((lower, sequential(t)));
                    if i + j + 1 < n {
                        let upper = vec![(x1, y0), (x0, y1), (x1, y1)];
                        let t = (component(
                            rt,
                            index,
                            comp,
                            ((x1 + x0 + x1) / 3.0, (y0 + y1 + y1) / 3.0),
                        ) - min)
                            / span;
                        cells.push((upper, sequential(t)));
                    }
                }
            }
            chart.draw_series(
                cells
                    .into_iter()
                    .map(|(points, color)| Polygon::new(points, color.filled())),
            )?;

            // Iso-lines of the affine component.
            for k in 1..=5 {
                let level = min + span * k as f64 / 6.0;
                if let Some((a, b)) = contour_segment(values, level) {
                    chart.draw_series(std::iter::once(PathElement::new(
                        vec![a, b],
                        RGBColor(45, 45, 45).stroke_width(1),
                    )))?;
                }
            }
            if min < 0.0 && max > 0.0 {
                if let Some((a, b)) = contour_segment(values, 0.0) {
                    chart.draw_series(std::iter::once(PathElement::new(
                        vec![a, b],
                        BLACK.stroke_width(2),
                    )))?;
                }
            }

            // Element outline.
            chart.draw_series(std::iter::once(PathElement::new(
                triangle_outline(),
                BLACK.stroke_width(2),
            )))?;

            // Horizontal colour bar under the element, with min/max labels.
            let bar = 48;
            chart.draw_series((0..bar).map(|k| {
                let t0 = k as f64 / bar as f64;
                let t1 = (k + 1) as f64 / bar as f64;
                Rectangle::new(
                    [(t0, -0.11), (t1, -0.06)],
                    sequential(0.5 * (t0 + t1)).filled(),
                )
            }))?;
            chart.draw_series(std::iter::once(PathElement::new(
                vec![
                    (0.0, -0.11),
                    (1.0, -0.11),
                    (1.0, -0.06),
                    (0.0, -0.06),
                    (0.0, -0.11),
                ],
                BLACK.stroke_width(1),
            )))?;
            chart.draw_series(std::iter::once(Text::new(
                format!("{min:.2}"),
                (0.0, -0.13),
                tick_style.clone().pos(Pos::new(HPos::Right, VPos::Top)),
            )))?;
            chart.draw_series(std::iter::once(Text::new(
                format!("{max:.2}"),
                (1.0, -0.13),
                tick_style.clone().pos(Pos::new(HPos::Left, VPos::Top)),
            )))?;

            chart.draw_series(std::iter::once(Text::new(
                component_label(index, comp),
                (0.5, 1.17),
                title_style,
            )))?;
        }
    }

    root.present()?;
    Ok(())
}
