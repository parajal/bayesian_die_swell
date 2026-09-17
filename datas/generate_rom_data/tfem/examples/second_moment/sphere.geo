Mesh.ElementOrder = 2; // Second-order elements
//+
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 1.0, -Pi/2, Pi/2, 2*Pi};
