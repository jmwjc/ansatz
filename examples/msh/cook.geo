// Cook's membrane (plane strain), structured quadrilateral mesh.
//
// Trapezoidal panel with corners
//   P1 = (0, 0),  P2 = (48, 44),  P3 = (48, 60),  P4 = (0, 44)   [mm]
// The left edge (x = 0) is clamped; the right edge (x = 48) carries a
// uniform shear traction in +y whose resultant is 1 N. The reference
// quantity is the vertical displacement at the top-right corner P3.
//
// Physical groups:
//   Γu : left edge  (essential / Dirichlet, clamped)
//   Γt : right edge (natural / Neumann, shear load)
//   Ω  : interior   (Quad4 elements)
//
// Generate a mesh with n elements per edge:
//   gmsh -2 examples/msh/cook.geo -format msh4 -o examples/msh/cook.msh
// Override n for a refinement study, e.g. sed 's/^n = .*/n = 32;/' cook.geo

n = 16;     // elements per edge

Point(1) = {0.0,  0.0, 0, 1.0};
Point(2) = {48.0, 44.0, 0, 1.0};
Point(3) = {48.0, 60.0, 0, 1.0};
Point(4) = {0.0,  44.0, 0, 1.0};

Line(1) = {1, 2}; // bottom edge
Line(2) = {2, 3}; // right edge (Γt, loaded)
Line(3) = {3, 4}; // top edge
Line(4) = {4, 1}; // left edge (Γu, clamped)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Structured transfinite quadrilateral mesh.
Transfinite Curve {1, 3} = n + 1;
Transfinite Curve {2, 4} = n + 1;
Transfinite Surface(1) = {1, 2, 3, 4};
Recombine Surface(1);

Physical Curve("Γu") = {4};
Physical Curve("Γt") = {2};
Physical Surface("Ω") = {1};

Mesh.MshFileVersion = 4.1;
