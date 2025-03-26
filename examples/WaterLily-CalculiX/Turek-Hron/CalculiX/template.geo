SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 0.5, 0, 2*Pi};
Point(2) = {0.4, 0.1, 00., 1.0};
Point(3) = {0.4, -0.1, 00., 1.0};
Point(4) = {6.1, -0.1, 00., 1.0};
Point(5) = {6.1, 0.1, 00., 1.0};
Line(2) = {2, 5};
Line(3) = {5, 4};
Line(4) = {4, 3};
Line(5) = {3, 2};
Curve Loop(1) = {1};
Curve Loop(2) = {2, 3, 4, 5};
Plane Surface(1) = {1, 2};
Curve Loop(3) = {1};
Plane Surface(2) = {3};
BooleanUnion{ Surface{2}; Delete; }{ Surface{1}; Delete; }
Extrude {0, 0, 0.1} {
    Surface{2}; 
}
// Physical Surface(15) = {4, 3, 5, 6};
Transfinite Curve {5, 11} = 32 Using Progression 1;
Transfinite Curve {9, 13, 6, 4} = 48 Using Progression 1;
Transfinite Curve {8, 7, 10, 12, 14, 3} = 2 Using Progression 1;
