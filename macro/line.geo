// Gmsh project created on Fri Aug  4 14:15:06 2023
f=1.5;
d=0.07*f;
Point(1) = {0.25*f, 0.8*f, 0, d};
//+
Point(2) = {-0.125*f, 0, 0, d};
//+
Point(3) = {0.25*f, 0, 0, d};
//+
Point(4) = {-0.125*f, 0.8*f, 0, d};
//+
Point(5) = {0,  0.8*f, 0, d/10};

Point(6) = {0, 0.133*f, 0, d/10};

Point(7) = {0, 0, 0, d/10};
//+
Line(1) = {3, 1};
//+
Line(2) = {1, 5};
//+
Line(3) = {5,6};
//+
Line(33) = {6,7};
//+
Line(4) = {7, 3};
//+
Line(5) = {7, 2};
//+
Line(6) = {2, 4};
//+
Line(7) = {4, 5};
//+
Curve Loop(1) = {1, 2, 3,33, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3,33, 5, 6, 7};
//+
Plane Surface(2) = {2};
//+
Physical Surface(0) = {1};
//+
Physical Surface(1) = {2};
//+
Physical Curve(1) = {4, 5};
//+
Physical Curve(3) = {1};
//+
Physical Curve(4) = {2, 7};
//+
Physical Curve(5) = {6};
