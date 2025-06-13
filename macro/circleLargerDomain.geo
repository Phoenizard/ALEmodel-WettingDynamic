// Gmsh project created on Tue Nov 23 15:18:19 2021
//+
Mesh.MshFileVersion = 2.2;
f = 2;
f2 = 0.5;
h = f*0.1;
hout = f*0.6;

//SetFactory("OpenCASCADE");
Point(1) = {0,0.0,0,hout/3};
Point(2) = {f2,0.0,0,h};
Point(3) = {-f2,0.0,0,h};
Circle(889) = {2,1,3};
Circle(890) = {3,1,2};
//Circle(889) = {0, 0.05, 0, 0.25, 0, 2*Pi};

Point(11) = {-f*0.5, -f*0.5, 0, h};
//+
Point(21) = {f*0.5, -f*0.5, 0, h};
//+
Point(31) = {f*0.5, f*0.5, 0, h};
//+
Point(41) = {-f*0.5, f*0.5, 0, h};
//+

Line(1) = {11, 21};
//+
Line(2) = {21, 31};
//+
Line(3) = {31, 41};
//+
Line(4) = {41, 11};

Point(441) = {-f*1, -f*1, 0, hout};
//+
Point(42) = {f*1, -f*1, 0, hout};
//+
Point(43) = {f*1, f*1, 0, hout};
//+
Point(44) = {-f*1, f*1, 0, hout};

Line(21) = {441, 42};
//+
Line(22) = {42, 43};
//+
Line(23) = {43, 44};
//+
Line(24) = {44, 441};

//+
Line(311) = {3, 1};
//+
Line(312) = {1, 2};

Line Loop(1) = {1,2,3,4};
//+
Line Loop(2) = {889,311,312};
Line Loop(3) = {-890,311,312};

Line Loop(4) = {21,22,23,24};
//+
Plane Surface(1) = {2};
Plane Surface(2) = {3};
//+
Plane Surface(0) = {1, 2,3};
//+
Plane Surface(10) = {4, 1};

Physical Line(3) = {22};
//+
Physical Line(4) = {24};
//+
Physical Line(5) = {23};
//+
Physical Line(1) = {21};
//+
Physical Surface(0) = {0,10};
//+
Physical Surface(1) = {1,2};

