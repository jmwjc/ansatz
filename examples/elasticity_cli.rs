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
//! `"sin(pi * x) * cos(pi * y)"`. The supported syntax is documented in the
//! `expr` module below.
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
    body_force, displacement_penalty, lame_constants, plane_stress_lambda, stiffness, traction,
};
use ansatz::sparse_vector::try_new_from_doublet;

// Physical-region naming convention for this CLI. These are intentionally
// local to the example: different problem drivers (e.g. a future `heat_cli.rs`)
// may use different physical-group names.
const INTERIOR: &str = "Ω";
const ESSENTIAL: &str = "Γᵍ";
const NATURAL: &str = "Γᵗ";

// ---------------------------------------------------------------------------
// Tiny scalar expression parser for `x`, `y`
// ---------------------------------------------------------------------------

mod expr {
    use std::iter::Peekable;
    use std::str::Chars;

    /// A parsed scalar expression in the two variables `x` and `y`.
    #[derive(Clone, Debug)]
    pub struct Expr(ExprKind);

    #[derive(Clone, Debug)]
    enum ExprKind {
        Num(f64),
        Var(char),
        Neg(Box<ExprKind>),
        Add(Box<ExprKind>, Box<ExprKind>),
        Sub(Box<ExprKind>, Box<ExprKind>),
        Mul(Box<ExprKind>, Box<ExprKind>),
        Div(Box<ExprKind>, Box<ExprKind>),
        Pow(Box<ExprKind>, Box<ExprKind>),
        Func(&'static str, Box<ExprKind>),
    }

    #[derive(Clone, Debug, PartialEq)]
    enum Token {
        Num(f64),
        Ident(String),
        Plus,
        Minus,
        Star,
        Slash,
        Caret,
        LParen,
        RParen,
    }

    impl Expr {
        /// A constant expression.
        pub fn constant(value: f64) -> Expr {
            Expr(ExprKind::Num(value))
        }

        /// Parses an expression in `x` and `y`.
        ///
        /// Supported syntax: `+ - * / ^`, parentheses, the constants `pi` and
        /// `e`, and the functions `sin`, `cos`, `tan`, `exp`, `ln`, `sqrt`,
        /// `abs`, `sinh`, `cosh`.
        pub fn parse(src: &str) -> Expr {
            parse_expr(src)
        }

        /// Evaluates the expression at the physical point `(x, y)`.
        pub fn eval(&self, x: f64, y: f64) -> f64 {
            self.0.eval(x, y)
        }
    }

    impl ExprKind {
        fn eval(&self, x: f64, y: f64) -> f64 {
            match self {
                ExprKind::Num(value) => *value,
                ExprKind::Var('x') => x,
                ExprKind::Var('y') => y,
                ExprKind::Var(other) => panic!("unknown variable: {other}"),
                ExprKind::Neg(a) => -a.eval(x, y),
                ExprKind::Add(a, b) => a.eval(x, y) + b.eval(x, y),
                ExprKind::Sub(a, b) => a.eval(x, y) - b.eval(x, y),
                ExprKind::Mul(a, b) => a.eval(x, y) * b.eval(x, y),
                ExprKind::Div(a, b) => a.eval(x, y) / b.eval(x, y),
                ExprKind::Pow(a, b) => a.eval(x, y).powf(b.eval(x, y)),
                ExprKind::Func(name, a) => {
                    let value = a.eval(x, y);
                    match *name {
                        "sin" => value.sin(),
                        "cos" => value.cos(),
                        "tan" => value.tan(),
                        "exp" => value.exp(),
                        "ln" => value.ln(),
                        "sqrt" => value.sqrt(),
                        "abs" => value.abs(),
                        "sinh" => value.sinh(),
                        "cosh" => value.cosh(),
                        _ => unreachable!(),
                    }
                }
            }
        }
    }

    fn parse_expr(src: &str) -> Expr {
        let tokens = lex(src);
        assert!(!tokens.is_empty(), "empty expression");
        let mut parser = Parser { tokens, pos: 0 };
        let node = parser.parse_add();
        assert!(parser.peek().is_none(), "unexpected trailing token");
        Expr(node)
    }

    fn lex(src: &str) -> Vec<Token> {
        let mut chars = src.chars().peekable();
        let mut tokens = Vec::new();
        while let Some(&c) = chars.peek() {
            if c.is_whitespace() {
                chars.next();
                continue;
            }
            match c {
                '+' => {
                    tokens.push(Token::Plus);
                    chars.next();
                }
                '-' => {
                    tokens.push(Token::Minus);
                    chars.next();
                }
                '*' => {
                    tokens.push(Token::Star);
                    chars.next();
                }
                '/' => {
                    tokens.push(Token::Slash);
                    chars.next();
                }
                '^' => {
                    tokens.push(Token::Caret);
                    chars.next();
                }
                '(' => {
                    tokens.push(Token::LParen);
                    chars.next();
                }
                ')' => {
                    tokens.push(Token::RParen);
                    chars.next();
                }
                '0'..='9' | '.' => {
                    tokens.push(Token::Num(lex_number(&mut chars)));
                }
                'a'..='z' | 'A'..='Z' | '_' => {
                    let mut ident = String::new();
                    while let Some(&c) = chars.peek() {
                        if c.is_alphanumeric() || c == '_' {
                            ident.push(c);
                            chars.next();
                        } else {
                            break;
                        }
                    }
                    tokens.push(Token::Ident(ident));
                }
                other => panic!("unexpected character in expression: {other}"),
            }
        }
        tokens
    }

    fn lex_number(chars: &mut Peekable<Chars<'_>>) -> f64 {
        let mut text = String::new();
        while let Some(&c) = chars.peek() {
            if c.is_ascii_digit() || c == '.' {
                text.push(c);
                chars.next();
            } else {
                break;
            }
        }
        if let Some(&c) = chars.peek() {
            if c == 'e' || c == 'E' {
                text.push(c);
                chars.next();
                if let Some(&c) = chars.peek() {
                    if c == '+' || c == '-' {
                        text.push(c);
                        chars.next();
                    }
                }
                while let Some(&c) = chars.peek() {
                    if c.is_ascii_digit() {
                        text.push(c);
                        chars.next();
                    } else {
                        break;
                    }
                }
            }
        }
        text.parse().unwrap_or_else(|_| panic!("invalid number: {text}"))
    }

    struct Parser {
        tokens: Vec<Token>,
        pos: usize,
    }

    impl Parser {
        fn peek(&self) -> Option<&Token> {
            self.tokens.get(self.pos)
        }

        fn bump(&mut self) -> Token {
            let token = self.tokens[self.pos].clone();
            self.pos += 1;
            token
        }

        fn eat(&mut self, token: Token) -> bool {
            if self.peek() == Some(&token) {
                self.bump();
                true
            } else {
                false
            }
        }

        fn parse_add(&mut self) -> ExprKind {
            let mut lhs = self.parse_mul();
            while let Some(tok) = self.peek() {
                match tok {
                    Token::Plus => {
                        self.bump();
                        let rhs = self.parse_mul();
                        lhs = ExprKind::Add(Box::new(lhs), Box::new(rhs));
                    }
                    Token::Minus => {
                        self.bump();
                        let rhs = self.parse_mul();
                        lhs = ExprKind::Sub(Box::new(lhs), Box::new(rhs));
                    }
                    _ => break,
                }
            }
            lhs
        }

        fn parse_mul(&mut self) -> ExprKind {
            let mut lhs = self.parse_unary();
            while let Some(tok) = self.peek() {
                match tok {
                    Token::Star => {
                        self.bump();
                        let rhs = self.parse_unary();
                        lhs = ExprKind::Mul(Box::new(lhs), Box::new(rhs));
                    }
                    Token::Slash => {
                        self.bump();
                        let rhs = self.parse_unary();
                        lhs = ExprKind::Div(Box::new(lhs), Box::new(rhs));
                    }
                    _ => break,
                }
            }
            lhs
        }

        fn parse_unary(&mut self) -> ExprKind {
            match self.peek() {
                Some(Token::Plus) => {
                    self.bump();
                    self.parse_unary()
                }
                Some(Token::Minus) => {
                    self.bump();
                    ExprKind::Neg(Box::new(self.parse_unary()))
                }
                _ => self.parse_power(),
            }
        }

        fn parse_power(&mut self) -> ExprKind {
            let base = self.parse_primary();
            if self.eat(Token::Caret) {
                let exponent = self.parse_unary();
                ExprKind::Pow(Box::new(base), Box::new(exponent))
            } else {
                base
            }
        }

        fn parse_primary(&mut self) -> ExprKind {
            match self.bump() {
                Token::Num(value) => ExprKind::Num(value),
                Token::Ident(name) => match name.as_str() {
                    "x" => ExprKind::Var('x'),
                    "y" => ExprKind::Var('y'),
                    "pi" => ExprKind::Num(std::f64::consts::PI),
                    "e" => ExprKind::Num(std::f64::consts::E),
                    "sin" | "cos" | "tan" | "exp" | "ln" | "sqrt" | "abs" | "sinh"
                    | "cosh" => {
                        assert!(self.eat(Token::LParen), "expected '(' after {name}");
                        let arg = self.parse_add();
                        assert!(self.eat(Token::RParen), "expected ')' after {name}");
                        let func: &'static str = match name.as_str() {
                            "sin" => "sin",
                            "cos" => "cos",
                            "tan" => "tan",
                            "exp" => "exp",
                            "ln" => "ln",
                            "sqrt" => "sqrt",
                            "abs" => "abs",
                            "sinh" => "sinh",
                            "cosh" => "cosh",
                            _ => unreachable!(),
                        };
                        ExprKind::Func(func, Box::new(arg))
                    }
                    other => panic!("unknown identifier: {other}"),
                },
                Token::LParen => {
                    let expr = self.parse_add();
                    assert!(self.eat(Token::RParen), "expected ')'");
                    expr
                }
                other => panic!("unexpected token: {other:?}"),
            }
        }
    }
}

// ---------------------------------------------------------------------------
// TOML configuration
// ---------------------------------------------------------------------------

mod config {
    use toml::Value;

    use super::expr::Expr;

    /// Plane-stress / plane-strain material model selector.
    #[derive(Clone, Copy)]
    pub enum ProblemType {
        PlaneStress,
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

    /// A two-component field (displacement, traction, or body force) whose
    /// entries are scalar expressions in the physical coordinates `x`, `y`.
    #[derive(Clone)]
    pub struct Field {
        components: Vec<Expr>,
    }

    impl Field {
        pub fn eval(&self, p: [f64; 2]) -> [f64; 2] {
            [
                self.components[0].eval(p[0], p[1]),
                self.components[1].eval(p[0], p[1]),
            ]
        }
    }

    /// Problem configuration read from the TOML file.
    pub struct Config {
        pub problem_type: ProblemType,
        pub young: f64,
        pub poisson: f64,
        pub body_force: Option<Field>,
        pub essential: Option<Field>,
        pub natural: Option<Field>,
        pub penalty: f64,
    }

    fn table<'a>(root: &'a Value, name: &str) -> &'a toml::Table {
        root.get(name)
            .and_then(Value::as_table)
            .unwrap_or_else(|| panic!("missing [{name}] table"))
    }

    fn optional_table<'a>(root: &'a Value, name: &str) -> Option<&'a toml::Table> {
        root.get(name).and_then(Value::as_table)
    }

    fn float(value: &Value, what: &str) -> f64 {
        value
            .as_float()
            .unwrap_or_else(|| panic!("{what} must be a float"))
    }

    fn component(value: &Value, what: &str) -> Expr {
        if let Some(number) = value.as_float() {
            Expr::constant(number)
        } else if let Some(number) = value.as_integer() {
            Expr::constant(number as f64)
        } else if let Some(expression) = value.as_str() {
            Expr::parse(expression)
        } else {
            panic!("{what} must be a number or a string expression in x, y")
        }
    }

    fn field(value: &Value, what: &str) -> Field {
        let array = value
            .as_array()
            .unwrap_or_else(|| panic!("{what} must be an array of two components"));
        assert_eq!(array.len(), 2, "{what} must have exactly two components");
        Field {
            components: array
                .iter()
                .map(|entry| component(entry, what))
                .collect(),
        }
    }

    pub fn parse(text: &str) -> Config {
        let root: Value = toml::from_str(text).expect("invalid TOML");

        let problem = table(&root, "problem");
        let problem_type = match problem
            .get("type")
            .and_then(Value::as_str)
            .unwrap_or("plane_strain")
        {
            "plane_stress" => ProblemType::PlaneStress,
            "plane_strain" => ProblemType::PlaneStrain,
            other => panic!("unknown problem.type: {other}"),
        };

        let material = table(&root, "material");
        let young = float(
            material.get("young").expect("missing material.young"),
            "material.young",
        );
        let poisson = float(
            material.get("poisson").expect("missing material.poisson"),
            "material.poisson",
        );

        let body_force = optional_table(&root, "body_force").map(|section| {
            field(
                section.get("value").expect("missing body_force.value"),
                "body_force.value",
            )
        });

        let (essential, penalty) = match optional_table(&root, "essential") {
            Some(section) => {
                let value = field(
                    section.get("value").expect("missing essential.value"),
                    "essential.value",
                );
                let penalty = section
                    .get("penalty")
                    .and_then(Value::as_float)
                    .unwrap_or(1e7);
                (Some(value), penalty)
            }
            None => (None, 1e7),
        };

        let natural = optional_table(&root, "natural").map(|section| {
            field(
                section.get("value").expect("missing natural.value"),
                "natural.value",
            )
        });

        Config {
            problem_type,
            young,
            poisson,
            body_force,
            essential,
            natural,
            penalty,
        }
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
        Some(mesh.elements(name))
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

    let (lambda, mu) = match cfg.problem_type {
        ProblemType::PlaneStress => {
            let lambda = plane_stress_lambda(cfg.young, cfg.poisson);
            let (_, mu) = lame_constants(cfg.young, cfg.poisson);
            (lambda, mu)
        }
        ProblemType::PlaneStrain => lame_constants(cfg.young, cfg.poisson),
    };

    let mut triplets = stiffness(&units, lambda, mu);
    let mut doublets = Vec::new();

    if let Some(field) = &cfg.body_force {
        let field = field.clone();
        doublets.extend(body_force(&units, move |x| field.eval(x)));
    }

    if let (Some(elems), Some(field)) =
        (optional_elements::<Seg2<2>>(mesh, ESSENTIAL), &cfg.essential)
    {
        let field = field.clone();
        let (pt, pd) = displacement_penalty(&elems, cfg.penalty, move |x| field.eval(x));
        triplets.extend(pt);
        doublets.extend(pd);
    }

    if let (Some(elems), Some(field)) = (optional_elements::<Seg2<2>>(mesh, NATURAL), &cfg.natural) {
        let field = field.clone();
        doublets.extend(traction(&elems, move |x| field.eval(x)));
    }

    let k = SparseColMat::<usize, f64>::try_new_from_triplets(n_dofs, n_dofs, &triplets).unwrap();
    let f = try_new_from_doublet(n_dofs, &doublets);
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
    let cfg = config::parse(&text);
    let mesh = read_gmsh::<2>(msh_path);

    println!("== elasticity CLI ==");
    println!("  config: {toml_path}");
    println!("  mesh:   {msh_path}");
    println!("  nodes:  {}", mesh.nodes.len());
    println!("  problem: {}", cfg.problem_type.label());
    println!("  material: E = {}, nu = {}", cfg.young, cfg.poisson);

    let elem_type = mesh.element_type(INTERIOR);
    let d = match elem_type {
        Tri3::<2>::GMSH_ELEMENT_TYPE => {
            let units: Vec<Tri3<2>> = mesh.elements(INTERIOR);
            println!("  interior: {} Tri3 elements", units.len());
            solve::<Tri3<2>, 3, 3>(&mesh, units, &cfg)
        }
        Quad4::<2>::GMSH_ELEMENT_TYPE => {
            let units: Vec<Quad4<2>> = mesh.elements(INTERIOR);
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
