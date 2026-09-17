Point(1) = {0, 0, 0, 10.0};
Translate {3, 0, 0} {
  Duplicata { Point{1}; }
}
//+
Translate {3, 0, 0} {
  Duplicata { Point{2}; }
}
//+
Translate {3, 0, 0} {
  Duplicata { Point{3}; }
}
//+
Translate {3, 0, 0} {
  Duplicata { Point{4}; }
}
//+
Translate {0, 3, 0} {
  Duplicata { Point{2}; Point{3}; Point{4}; }
}
//+
Translate {0, 6, 0} {
  Duplicata { Point{5}; }
}
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 8};
//+
Line(6) = {8, 7};
//+
Line(7) = {7, 6};
//+
Line(8) = {6, 1};
//+
Line(9) = {6, 2};
//+
Line(10) = {6, 3};
//+
Line(11) = {3, 7};
//+
Line(12) = {7, 4};
//+
Line(13) = {4, 8};
//+
Line(14) = {8, 9};
//+
Line(15) = {9, 5};
//+
Physical Point("points") = {1, 2, 3, 4, 5, 9, 8, 7, 6};
//+
Physical Curve("outside") = {1, 2, 3, 4, 15, 14, 6, 7, 8};
//+
Physical Curve("inside") = {9, 10, 11, 12, 13, 5};
