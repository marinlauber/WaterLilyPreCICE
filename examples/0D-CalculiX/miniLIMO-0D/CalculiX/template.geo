SetFactory("OpenCASCADE");

R = DefineNumber[ 0.12740, Name "Parameters/R" ];
r = DefineNumber[ 0.013650, Name "Parameters/r" ];
L = DefineNumber[ 1.0, Name "Parameters/L" ];
l = DefineNumber[ 0.714370, Name "Parameters/l" ];
W = DefineNumber[ 0.92822/2, Name "Parameters/W" ];
w = DefineNumber[ 0.25935, Name "Parameters/w" ];

// outer edge of the pouch
Point(1) = {  -W,   0, 0, 1.0};
Point(2) = {  -W, L-R, 0, 1.0};
Point(3) = {-W+R, L-R, 0, 1.0};
Point(4) = {-W+R,   L, 0, 1.0};
Point(5) = {   W, L-R, 0, 1.0};
Point(6) = { W-R, L-R, 0, 1.0};
Point(7) = { W-R,   L, 0, 1.0};
Point(8) = {   W,   0, 0, 1.0};

// inner edge of the pouch
Point( 9) = {-W+R/2,  L-l, 0, 1.0};
Point(10) = {-W+R , L-l, 0, 1.0};
Point(11) = {-W+R,  L-l-R/2, 0, 1.0};
Point(12) = {-W+R/2,   L-R, 0, 1.0};
Point(13) = {-W+R, L-R/2, 0, 1.0};
Point(14) = { W-R/2,   L-R, 0, 1.0};
Point(15) = { W-R, L-R/2, 0, 1.0};
Point(16) = { W-R/2, L-l, 0, 1.0};
Point(17) = { W-R ,  L-l, 0, 1.0};
Point(18) = { W-R,  L-l-R/2, 0, 1.0};

// internal edges of the pouch
Point(19) = {-w/2,   L-R/2, 0, 1.0};
Point(20) = {-w/2, L-l-R/2, 0, 1.0};
Point(21) = { w/2,   L-R/2, 0, 1.0};
Point(22) = { w/2, L-l-R/2, 0, 1.0};

// outer seams
Line(1) = {1, 2};
Circle(2) = {2, 3, 4};
Line(3) = {4, 7};
Circle(4) = {7, 6, 5};
Line(5) = {5, 8};
Line(6) = {8, 1};

// inner seams
Line(7) = {11, 20};
Line(8) = {20, 22};
Line(9) = {22, 18};
Circle(10) = {18, 17, 16};
Line(11) = {16, 14};
Circle(12) = {14, 6, 15};
Line(13) = {15, 21};
Line(14) = {21, 19};
Line(15) = {19, 13};
Circle(16) = {13, 3, 12};
Line(17) = {12, 9};
Circle(18) = {9, 10, 11};
Line(19) = {19, 20};
Line(20) = {21, 22};

// surfaces
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Curve Loop(2) = {7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
Plane Surface(1) = {1, 2};
Reverse Surface{1}; // needed since the normals must be all in the same direction
Curve Loop(3) = {15, 16, 17, 18, 7, -19};
Plane Surface(2) = {3};
Curve Loop(4) = {14, 19, 8, -20};
Plane Surface(3) = {4};
Curve Loop(5) = {13, 20, 9, 10, 11, 12};
Plane Surface(4) = {5};

// mesh definition
Characteristic Length{ PointsOf{ Surface{1,2,3,4}; } } = 0.03;
Recombine Surface{1,2,3,4};

