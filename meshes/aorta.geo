SetFactory("OpenCASCADE");
Circle(2) = {0, 0, 0, 0.5, 0, Pi};
Point(4) = {-0.5, -1.0, 0.0, 1.0};
Point(5) = { 0.5, -1.0, 0.0, 1.0};
Point(6) = { 0.7, -1.0, 0.0, 1.0};
Point(7) = { 0.3, -1.0, 0.0, 1.0};
Point(8) = { 0.5, -1.0, -0.2, 1.0};
Point(9) = { 0.5, -1.0, 0.2, 1.0};
Line(3) = {2, 4};
Line(4) = {1, 5};
Circle(5) = {8, 5, 7};
Circle(6) = {7, 5, 9};
Circle(7) = {9, 5, 6};
Circle(8) = {6, 5, 8};
// extrude first bit
Wire(2) = {4};
Extrude { Curve{5}; Curve{6}; Curve{7}; Curve{8}; } Using Wire {2}

Wire(7) = {2};
Extrude { Curve{10}; Curve{13}; Curve{16}; Curve{19}; } Using Wire {7}

Wire(12) = {3};
Extrude { Curve{32}; Curve{24}; Curve{27}; Curve{30}; } Using Wire {12}

Characteristic Length{ PointsOf{ Surface{1,2,3,4,5,6,7,8,9,10,11,12}; } } = 0.03;
