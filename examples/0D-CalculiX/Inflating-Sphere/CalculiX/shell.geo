SetFactory("OpenCASCADE");

Sphere(1) = {0, 0, 0, 25, -Pi/2, Pi/2, 2*Pi};
Characteristic Length{ PointsOf{ Surface{1}; } } = 5;
