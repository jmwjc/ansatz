// Linear-elasticity CLI example: unit square plate under uniform tension.
//
// Physical groups (shared with the TOML configuration and the library):
//   Γᵍ : left edge  (essential / Dirichlet boundary)
//   Γᵗ : right edge (natural / Neumann boundary)
//   Ω  : interior   (solution domain)
//
// Generate the mesh with:
//   gmsh -2 examples/msh/elasticity_cli.geo -format msh4 -o examples/msh/elasticity_cli.msh

L = 1.0;   // side length
n = 8;     // elements per edge

Point(1) = {0.0, 0.0, 0, 1.0};
Point(2) = {L,   0.0, 0, 1.0};
Point(3) = {L,   L,   0, 1.0};
Point(4) = {0.0, L,   0, 1.0};

Line(1) = {1, 2}; // bottom edge (traction free)
Line(2) = {2, 3}; // right edge  (Γᵗ, loaded)
Line(3) = {3, 4}; // top edge    (traction free)
Line(4) = {4, 1}; // left edge   (Γᵍ, fixed)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Structured (transfinite) quadrilateral mesh.
Transfinite Curve {1, 3} = n + 1;
Transfinite Curve {2, 4} = n + 1;
Transfinite Surface(1) = {1, 2, 3, 4};
Recombine Surface(1);

Physical Curve("Γᵍ") = {4};
Physical Curve("Γᵗ") = {2};
Physical Surface("Ω") = {1};

Mesh.MshFileVersion = 4.1;
