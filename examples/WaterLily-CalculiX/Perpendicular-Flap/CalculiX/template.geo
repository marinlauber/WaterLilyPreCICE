SetFactory("OpenCASCADE");
Point(1) = {0, 0, -0.5, 1.0};
Point(2) = {6.4, 0, -0.5, 1.0};
Point(3) = {6.4, 64, -0.5, 1.0};
Point(4) = {0, 64, -0.5, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};
Extrude {0, 0, 1} {
  Surface{1};
}
Transfinite Surface {1, 2, 3, 4, 5, 6};
// Physical Surface("Interface", 13) = {2, 3, 4, 5};
