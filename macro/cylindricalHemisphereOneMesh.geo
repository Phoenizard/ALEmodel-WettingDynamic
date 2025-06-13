Mesh.MshFileVersion = 2.2;
factor = 1; //1e-5 skalierungsfaktor
radius = 0.54 * factor;
width = 0.1 * factor;
height = 0.45 * factor;
x_scale = 2; // 1 gives 2x1 geometry, 2 gives 4x1 ...
yOffset = -0.001;

l_circ = 0.06 * factor;
l_box = 4* 0.1 * factor;

Point(1) = {-x_scale*factor, 0, 0, l_box};
Point(2) = {0, height, 0, l_box}; //circle-center
Point(22) = {0, width + height, 0, l_circ}; //circle-center
Point(3) = {- width, 0, 0, l_circ};
Point(4) = {width, 0, 0, l_circ};
Point(5) = {-width, height, 0, l_circ};
Point(55) = {width, height, 0, l_circ};
Point(6) = {x_scale * factor, 0, 0, l_box};
Point(7) = {x_scale * factor,1.094* factor, 0, l_box};
Point(8) = {-x_scale*factor, 1.094*factor, 0, l_box};

Circle(11) = {5, 2, 22};
Circle(12) = {22, 2, 55};

Point(30) = {- 1.5*width, 0, 0, l_circ};
Point(40) = {1.5*width, 0, 0, l_circ};
Point(50) = {-1.5*width, height, 0, l_circ};
Point(550) = {1.5*width, height, 0, l_circ};
Point(220) = {0, 1.5*width + height, 0, l_circ}; 

Circle(110) = {50, 2, 220};
Circle(120) = {220, 2, 550};

Line(13) = {1, 30};
Line(14) = {40, 6};
Line(15) = {6, 7};
Line(16) = {7, 8};
Line(17) = {8, 1};
Line(20) = {3, 5};
Line(21) = {55, 4};
Line(22) = {4, 3};

Line(200) = {30, 50};
Line(210) = {550, 40};
Line(211) = {30, 3};
Line(212) = {4, 40};

Physical Line(1) = {13, 22, 14, 211, 212};
//Physical Line(2) = {11, 12, 20, 21};
Physical Line(3) = {16};
Physical Line(4) = {15};
Physical Line(5) = {17};

Line Loop(21) = {13, 200, 110, 120, 210, 14, 15, 16, 17};
Line Loop(210) = {211, 20, 11, 12, 21, 212, -210, -120, -110, -200};
Line Loop(22) = {20, 11, 12, 21, 22};
Plane Surface(31) = {21};
Plane Surface(310) = {210};
Plane Surface(32) = {22};

Physical Surface(0) = {31,310};
Physical Surface(1) = {32};

