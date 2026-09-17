Mesh.ElementOrder = 2; // Second-order elements
//+
Point(1) = {1, 0, 0, 0.1};
//+
Point(2) = {0, 0, 0, 0.1};
//+
Point(3) = {0, 1, 0, 0.1};
//+
Point(4) = {-1, 0, 0, 0.1};
//+
Circle(1) = {1, 2, 3};
//+
Circle(2) = {3, 2, 4};
//+
Line(3) = {4, 2};
//+
Line(4) = {2, 1};
