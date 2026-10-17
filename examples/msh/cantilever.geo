// Cantilever beam under a transverse tip load (plane stress).
//
// Domain: rectangle [0, L] x [-H/2, H/2], with the left end (x = 0) clamped
// by the analytical displacement field and the right end (x = L) loaded by the
// equilibrating parabolic shear traction. Top and bottom edges are traction
// free.
//
// Physical groups:
//   Γu : left end  (essential / Dirichlet boundary)
//   Γt : right end (natural / Neumann boundary)
//   Ω  : interior  (plane-stress quadrilateral elements)
//
// Generate the mesh with:
//   gmsh -2 examples/msh/cantilever.geo -format msh4 -o examples/msh/cantilever.msh

L  = 8.0;   // beam length
H  = 1.0;   // beam height
nL = 80;    // elements along the length
nH = 8;     // elements across the height (even: a node sits on the mid-line)

Point(1) = {0.0, -H/2, 0, 1.0};
Point(2) = {L,   -H/2, 0, 1.0};
Point(3) = {L,    H/2, 0, 1.0};
Point(4) = {0.0,  H/2, 0, 1.0};

Line(1) = {1, 2}; // bottom edge (traction free)
Line(2) = {2, 3}; // right edge  (Γt)
Line(3) = {3, 4}; // top edge    (traction free)
Line(4) = {4, 1}; // left edge   (Γu)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Structured (transfinite) quadrilateral mesh.
Transfinite Curve {1, 3} = nL + 1;
Transfinite Curve {2, 4} = nH + 1;
Transfinite Surface(1) = {1, 2, 3, 4};
Recombine Surface(1);

Physical Curve("Γu") = {4};
Physical Curve("Γt") = {2};
Physical Surface("Ω") = {1};

Mesh.MshFileVersion = 4.1;
