SetFactory("OpenCASCADE");

Sphere(1) = {0, 0, 0, 0.5, -Pi/2, Pi/2, 2*Pi};
Characteristic Length{ PointsOf{ Surface{1}; } } = 0.1;
//+
Physical Surface("shell", 4) = {1};
