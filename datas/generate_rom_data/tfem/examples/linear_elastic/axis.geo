Mesh.ElementOrder = 2; // Second-order elements
lc=0.1;
gc=1;
Point(1) = {0, 0, 0, gc};
Point(2) = {10, 0, 0, gc};
Point(3) = {10, 4, 0, gc};
Point(4) = {5, 4, 0, lc};
Point(5) = {5, 8, 0, gc};
Point(6) = {0, 8, 0, gc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
//+
Plane Surface(1) = {1};
