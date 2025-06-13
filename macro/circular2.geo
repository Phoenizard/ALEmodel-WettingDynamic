// Gmsh project created on Tue Nov 23 15:18:19 2021
//+
Mesh.MshFileVersion = 2.2;
f = 5;
f2 = 0.5;
h = f*0.1;
hout = f*0.3;

//SetFactory("OpenCASCADE");
Point(1) = {0,0.0,0,hout/5};
Point(2) = {f2*0.25,0.0,0,h/8};
Point(3) = {-f2*0.25,0.0,0,h/8};
Circle(889) = {2,1,3};

Point(11) = {-f*0.5, 0, 0, h};
//+
Point(21) = {f*0.5, 0, 0, h};
//+
Point(31) = {f*0.5, f*0.5, 0, h};
//+
Point(41) = {-f*0.5, f*0.5, 0, h};
//+

Line(2) = {21, 31};
//+
Line(3) = {31, 41};
//+
Line(4) = {41, 11};

//+
Line(311) = {3, 1};
//+
Line(312) = {1, 2};
//+
Line(5) = {11, 3};
//+
Line(6) = {2, 21};

Line Loop(1) = {2,3,4,5,311,312,6};
//+
Line Loop(2) = {889,311,312};

//+
Plane Surface(1) = {2};
Plane Surface(2) = {1,2};

Physical Line(3) = {2};
//+
Physical Line(4) = {4};
//+
Physical Line(5) = {3};
//+
Physical Line(1) = {5,311,312,6};
//+
Physical Surface(1) = {1};
//+
Physical Surface(0) = {2};

