Mesh.ElementOrder = 2; // Second-order elements
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Circle(1) = {2, 1, 3};
//+
Line(2) = {3, 1};
//+
Line(3) = {1, 2};
//+
Characteristic Length {1, 2, 3} = 0.1;
