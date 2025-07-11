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
Transfinite Curve {8, 2, 4, 12} = 48 Using Progression 1;
Transfinite Curve {6, 11, 7, 9} = 4 Using Progression 1;
// Transfinite Curve {13, 5, 10, 3} = NX/10 Using Progression 1;
