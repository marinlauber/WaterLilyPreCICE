SetFactory("OpenCASCADE");

Box(1) = {-0.5, -0.5, -0.5, 1, 1, 1};
// Characteristic Length{ PointsOf{ Surface{1,2,3,4,5,6}; } } = 1;
Recombine Surface {1, 5, 4, 2, 3, 6};
//+
Transfinite Surface {1} = {3, 4, 2, 1};
//+
Transfinite Surface {6} = {3, 1, 5, 7};
//+
Transfinite Surface {2} = {7, 8, 6, 5};
//+
Transfinite Surface {5} = {8, 6, 2, 4};
//+
Transfinite Surface {4} = {3, 7, 8, 4};
//+
Transfinite Surface {3} = {5, 1, 2, 6};
