Mesh.ElementOrder=1;

h = 1.0;
nl = 5;
lc = h/nl;


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

Point(5) = {0, 0, h, lc};
Point(6) = {1, 0, h, lc} ;

Line(5) = {5,6} ;

Extrude {0,1,0} { Line{5}; Layers{nl}; Recombine; } 

Line(9) = {1,5} ;
Line(10) = {2,6} ;
Line(11) = {3,8} ;
Line(12) = {4,7} ;

Curve Loop(2) = {12, -7, -9, -4};

Plane Surface(10) = {2};

Curve Loop(3) = {12, 6, -11, 3};

Plane Surface(11) = {3};

Curve Loop(4) = {8, -11, 2, 10};

Plane Surface(12) = {4};

Curve Loop(5) = {9, 5, -10, -1};

Plane Surface(13) = {5};

Surface Loop(1) = {9, 13, 10, 11, 12, 1};

Volume(1) = {1};
