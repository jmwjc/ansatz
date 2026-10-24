//! General 2D linear-elasticity driver driven by a TOML configuration and a
//! GMSH mesh.
//!
//! Run with:
//!
//! ```text
//! cargo run --example elasticity_cli -- examples/elasticity_cli.toml examples/msh/elasticity_cli.msh
//! ```
//!
//! The mesh and the configuration share three physical-region names, declared
//! as local constants at the top of this file:
//!
//! * `Ω`  — interior (solution domain);
//! * `Γᵍ` — essential (Dirichlet) boundary, where the prescribed displacement
//!   from `[essential]` is imposed with the penalty method;
//! * `Γᵗ` — natural (Neumann) boundary, where the prescribed traction from
//!   `[natural]` is applied.
//!
//! Body force, traction, and prescribed displacement are specified as
//! **functions of the physical coordinates `x`, `y`**: each component of a
//! `value` array is either a number or a string expression such as
//! `"sin(pi * x) * cos(pi * y)"`, written in [`fasteval`] syntax. On top of
//! fasteval's operators (`+ - * / ^`, comparisons, ...) and built-ins, the
//! bare constants `pi` and `e` as well as the functions `sqrt`, `exp`, and
//! `ln` are also accepted.
//!
//! The interior element type is chosen at run time from the GMSH element type
//! stored in the `Ω` physical group: `Tri3` (type 2) or `Quad4` (type 3).
//! The boundary groups are expected to contain `Seg2` elements (type 1).

use std::fs;

use faer::Side;
use faer::prelude::*;
use faer::sparse::SparseColMat;

use ansatz::element::quad4::Quad4;
use ansatz::element::seg2::Seg2;
use ansatz::element::tri3::Tri3;
use ansatz::element::{Element, FullIntegration};
use ansatz::io::gmsh::{FromGmsh, GmshMesh, read_gmsh};
use ansatz::operation::elasticity::{
    Lame, body_force, displacement_penalty, stiffness, traction,
};
use ansatz::sparse_vector::try_new_from_doublet;

// Physical-region naming convention for this CLI. These are intentionally
// local to the example: different problem drivers (e.g. a future `heat_cli.rs`)
// may use different physical-group names.
const INTERIOR: &str = "Ω";
const ESSENTIAL: &str = "Γᵍ";
const NATURAL: &str = "Γᵗ";

// ---------------------------------------------------------------------------
// Scalar expressions in `x` and `y`
// ---------------------------------------------------------------------------

mod expr {
    use fasteval::{Compiler, Evaler, Instruction, Parser, Slab};

    /// A compiled scalar expression in the two variables `x` and `y`.
    pub struct Expr {
        slab: Slab,
        instr: Instruction,
    }

    impl Expr {
        /// A constant expression.
        pub fn constant(value: f64) -> Expr {
            Expr {
                slab: Slab::new(),
                instr: Instruction::IConst(value),
            }
        }

        /// Parses and compiles an expression in `x` and `y`.
        ///
        /// All of fasteval's operators and built-in functions are available.
        /// In addition, `pi` and `e` may be used without parentheses and
        /// `sqrt`, `exp`, and `ln` are provided as functions.
        pub fn parse(src: &str) -> Result<Expr, fasteval::Error> {
            let mut slab = Slab::new();
            let parsed = Parser::new().parse(src, &mut slab.ps)?.from(&slab.ps);
            let instr = parsed.compile(&slab.ps, &mut slab.cs);
            Ok(Expr { slab, instr })
        }

        /// Evaluates the expression at the physical point `(x, y)`.
        pub fn eval(&self, [x, y]: [f64; 2]) -> f64 {
            let mut namespace = |name: &str, args: Vec<f64>| -> Option<f64> {
                match name {
                    "x" => Some(x),
                    "y" => Some(y),
                    "pi" => Some(std::f64::consts::PI),
                    "e" => Some(std::f64::consts::E),
                    "sqrt" => args.first().copied().map(f64::sqrt),
                    "exp" => args.first().copied().map(f64::exp),
                    "ln" => args.first().copied().map(f64::ln),
                    _ => None,
                }
            };
            self.instr
                .eval(&self.slab, &mut namespace)
                .unwrap_or_else(|err| panic!("cannot evaluate expression: {err}"))
        }
    }
}

// ---------------------------------------------------------------------------
// TOML configuration
// ---------------------------------------------------------------------------

mod config {
    use serde::Deserialize;

    use super::expr::Expr;

    /// Plane-stress / plane-strain material model selector.
    #[derive(Clone, Copy, Default, Deserialize)]
    #[serde(rename_all = "snake_case")]
    pub enum ProblemType {
        PlaneStress,
        #[default]
        PlaneStrain,
    }

    impl ProblemType {
        pub fn label(self) -> &'static str {
            match self {
                Self::PlaneStress => "plane stress",
                Self::PlaneStrain => "plane strain",
            }
        }
    }

    /// One component of a field: either a number or an expression string.
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum Component {
        Number(f64),
        Expression(String),
    }

    impl Component {
        fn compile(&self) -> Expr {
            match self {
                Self::Number(value) => Expr::constant(*value),
                Self::Expression(src) => Expr::parse(src)
                    .unwrap_or_else(|err| panic!("invalid expression `{src}`: {err}")),
            }
        }
    }

    /// A two-component field (displacement, traction, or body force) whose
    /// entries are scalar expressions in the physical coordinates `x`, `y`.
    #[derive(Deserialize)]
    #[serde(from = "RawField")]
    pub struct Field([Expr; 2]);

    impl Field {
        pub fn eval(&self, point: [f64; 2]) -> [f64; 2] {
            [self.0[0].eval(point), self.0[1].eval(point)]
        }
    }

    #[derive(Deserialize)]
    struct RawField {
        value: [Component; 2],
    }

    impl From<RawField> for Field {
        fn from(raw: RawField) -> Self {
            let [first, second] = raw.value;
            Self([first.compile(), second.compile()])
        }
    }

    /// The `[essential]` table: prescribed displacement plus penalty factor.
    #[derive(Deserialize)]
    #[serde(from = "RawEssential")]
    pub struct Essential {
        pub prescribed: Field,
        pub penalty: f64,
    }

    #[derive(Deserialize)]
    struct RawEssential {
        value: [Component; 2],
        #[serde(default = "default_penalty")]
        penalty: f64,
    }

    impl From<RawEssential> for Essential {
        fn from(raw: RawEssential) -> Self {
            Self {
                prescribed: Field::from(RawField { value: raw.value }),
                penalty: raw.penalty,
            }
        }
    }

    fn default_penalty() -> f64 {
        1e7
    }

    /// Problem configuration read from the TOML file.
    #[derive(Deserialize)]
    pub struct Config {
        #[serde(default)]
        pub problem: Problem,
        pub material: Material,
        pub body_force: Option<Field>,
        pub essential: Option<Essential>,
        pub natural: Option<Field>,
    }

    impl Config {
        pub fn parse(text: &str) -> Config {
            toml::from_str(text).unwrap_or_else(|err| panic!("invalid configuration: {err}"))
        }
    }

    #[derive(Default, Deserialize)]
    pub struct Problem {
        #[serde(rename = "type", default)]
        pub kind: ProblemType,
    }

    #[derive(Deserialize)]
    pub struct Material {
        pub young: f64,
        pub poisson: f64,
    }
}

use config::{Config, ProblemType};

// ---------------------------------------------------------------------------
// Assembly and solve
// ---------------------------------------------------------------------------

/// Returns the typed boundary elements for `name`, or `None` when the mesh has
/// no such physical group.
fn optional_elements<T: FromGmsh<2>>(mesh: &GmshMesh<2>, name: &str) -> Option<Vec<T>> {
    if mesh.physical_names.values().any(|n| n == name) {
        Some(mesh.elements(name).unwrap())
    } else {
        None
    }
}

/// Assembles and solves the plane-stress/plane-strain problem for an interior
/// element type `T` with `N` nodes and `G` full-integration points.
fn solve<T, const N: usize, const G: usize>(
    mesh: &GmshMesh<2>,
    units: Vec<T>,
    cfg: &Config,
) -> Col<f64>
where
    T: FromGmsh<2> + Element<2, 2, N> + FullIntegration<2, G>,
{
    let max_id = mesh.nodes.keys().copied().max().unwrap();
    let n_dofs = (max_id + 1) * 2;

    let (young, poisson) = (cfg.material.young, cfg.material.poisson);
    let lame = match cfg.problem.kind {
        ProblemType::PlaneStress => Lame::plane_stress(young, poisson),
        ProblemType::PlaneStrain => Lame::plane_strain(young, poisson),
    };
    let (lambda, mu) = (lame.lambda, lame.mu);

    let mut triplets = stiffness(&units, lambda, mu);
    let mut doublets = Vec::new();

    if let Some(field) = &cfg.body_force {
        doublets.extend(body_force(&units, |x| field.eval(x)));
    }

    if let (Some(elems), Some(essential)) = (
        optional_elements::<Seg2<2>>(mesh, ESSENTIAL),
        &cfg.essential,
    ) {
        let (pt, pd) =
            displacement_penalty(&elems, essential.penalty, |x| essential.prescribed.eval(x));
        triplets.extend(pt);
        doublets.extend(pd);
    }

    if let (Some(elems), Some(field)) = (optional_elements::<Seg2<2>>(mesh, NATURAL), &cfg.natural)
    {
        doublets.extend(traction(&elems, |x| field.eval(x)));
    }

    let k = SparseColMat::<usize, f64>::try_new_from_triplets(n_dofs, n_dofs, &triplets).unwrap();
    let f = try_new_from_doublet(n_dofs, &doublets).unwrap();
    k.sp_cholesky(Side::Lower).unwrap().solve(&f)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: cargo run --example elasticity_cli -- <config.toml> <mesh.msh>");
        std::process::exit(2);
    }
    let toml_path = &args[1];
    let msh_path = &args[2];

    let text = fs::read_to_string(toml_path).expect("cannot read TOML file");
    let cfg = Config::parse(&text);
    let mesh = read_gmsh::<2>(msh_path).unwrap();

    println!("== elasticity CLI ==");
    println!("  config: {toml_path}");
    println!("  mesh:   {msh_path}");
    println!("  nodes:  {}", mesh.nodes.len());
    println!("  problem: {}", cfg.problem.kind.label());
    println!(
        "  material: E = {}, nu = {}",
        cfg.material.young, cfg.material.poisson
    );

    let elem_type = mesh.element_type(INTERIOR).unwrap();
    let d = match elem_type {
        Tri3::<2>::GMSH_ELEMENT_TYPE => {
            let units: Vec<Tri3<2>> = mesh.elements(INTERIOR).unwrap();
            println!("  interior: {} Tri3 elements", units.len());
            solve::<Tri3<2>, 3, 3>(&mesh, units, &cfg)
        }
        Quad4::<2>::GMSH_ELEMENT_TYPE => {
            let units: Vec<Quad4<2>> = mesh.elements(INTERIOR).unwrap();
            println!("  interior: {} Quad4 elements", units.len());
            solve::<Quad4<2>, 4, 4>(&mesh, units, &cfg)
        }
        other => {
            panic!("unsupported interior element type {other}: expected Tri3 (2) or Quad4 (3)")
        }
    };

    let max_id = mesh.nodes.keys().copied().max().unwrap();
    let (mut max_ux, mut max_uy) = (0.0_f64, 0.0_f64);
    for id in 0..=max_id {
        max_ux = max_ux.max(d[id * 2].abs());
        max_uy = max_uy.max(d[id * 2 + 1].abs());
    }
    println!("  max |u_x|: {max_ux:.6e}");
    println!("  max |u_y|: {max_uy:.6e}");
}
