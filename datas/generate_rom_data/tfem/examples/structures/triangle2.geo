Mesh.ElementOrder = 2; // Second-order elements
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Characteristic Length {1, 2, 3} = 0.1;
