factor = 1; //1e-5 skalierungsfaktor
y_scale = 1.0 * factor;
width = 0.215485 * factor; //0.3563
height = 0.669058 * factor; // 0.4505
x_scale = y_scale; // 1 gives 2x1 geometry, 2 gives 4x1 ...
size = 1.5*factor;
yOffset = -0.000;

l_circ = 5*0.015 * factor;
l_box = 5* 0.1 * factor;

Point(9) = {0, height, 0, l_box}; //circle-center
Point(10) = {-width, height, 0, l_circ};
Point(11) = {width, height, 0, l_circ};
Point(12) = {0, width + height, 0, l_circ};
Point(13) = {-width,  yOffset, 0, l_circ};
Point(14) = {width, yOffset, 0, l_circ};

Circle(31) = {10, 9, 12};
Circle(32) = {12, 9, 11};

Point(441000) = {-size, -size, 0, l_box};
Point(42) = {size, -size, 0, l_box};
Point(43) = {size, size, 0, l_box};
Point(44) = {-size, size, 0, l_box};

Line(210) = {441000,42};
Line(220) = {42,43};
Line(230) = {43,44};
Line(240) = {44,441000};


Line(18) = {13, 10};
Line(19) = {11, 14};
Line(22) = {14, 13};


Line Loop(22) = {18, 31, 32, 19, 22};
Plane Surface(32) = {22};

//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi} {
  Duplicata { Point{10}; Point{11}; Point{12}; Point{13}; Point{14};  Curve{18}; Curve{19}; Curve{22}; Curve{31}; Curve{32}; Surface{32}; }
}

Line Loop(4) = {210,220,230,240};
Curve Loop(23) = {18, 31, 32, 19, 241, 245, 242, 244};
Plane Surface(333) = {4, 23};
Physical Line(1) = {210};
Physical Line(3) = {220};
Physical Line(4) = {240};
Physical Line(5) = {230};
Physical Surface(0) = {333};
Physical Surface(1) = {32,246};
Mesh.MshFileVersion = 2.2;

Field[1] = Distance;
Field[1].CurvesList = {31,32,18,19,241,242,244,245};
Field[1].Sampling = 100;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].DistMax = l_box;
Field[2].DistMin = 0.75*l_circ;
Field[2].SizeMax = 0.8*l_box;
Field[2].SizeMin = 0.75*l_circ;

Field[7] = Min;
Field[7].FieldsList = {2};
Background Field = 7;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;//+