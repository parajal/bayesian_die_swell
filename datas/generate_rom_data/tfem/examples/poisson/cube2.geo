Mesh.ElementOrder=2;

nl = 5;
lc = 1/nl;

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc} ;
Point(3) = {1, 1, 0, lc} ;
Point(4) = {0, 1, 0, lc} ;

Line(1) = {1,2} ;
Line(2) = {3,2} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

Curve Loop(1) = {4,1,-2,3} ;

Plane Surface(1) = {1} ;

h = 1;

Extrude {0,0,h} { Surface{1}; Layers{nl}; Recombine; } 
