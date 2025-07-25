SetFactory("OpenCASCADE");
Point(2) = {0.0, NX/32, -NX/320.0, 1.0};
Point(3) = {0.0, -NX/32, -NX/320.0, 1.0};
Point(4) = {  NX, -NX/32, -NX/320.0, 1.0};
Point(5) = {  NX, NX/32, -NX/320.0, 1.0};
Line(2) = {2, 5};
Line(3) = {5, 4};
Line(4) = {4, 3};
Line(5) = {3, 2};
Curve Loop(1) = {2, 3, 4, 5};
Plane Surface(1) = {1};
Extrude {0, 0, NX/160.0} {
    Surface{1};
}
// Physical Surface(14) = {5, 2, 4, 3};
Transfinite Surface {1, 2, 3, 4, 5, 6};
