Mesh.MshFileVersion = 2.2;
factor = 1; //1e-5 skalierungsfaktor
radius = 0.54 * factor;
width = 0.1 * factor;
height = 0.45 * factor;
x_scale = 2; // 1 gives 2x1 geometry, 2 gives 4x1 ...
yOffset = -0.001;

l_circ = 0.06 * factor;
l_box = 0.4 * factor;

Point(1) = {-x_scale*factor, 0, 0, l_box};
//Point(2) = {0, height, 0, l_box}; //circle-center
//Point(22) = {0, width + height, 0, l_circ/2}; //circle-center
//Point(3) = {- width, 0, 0, l_circ/2};
//Point(4) = {width, 0, 0, l_circ/2};
//Point(5) = {-width, height, 0, l_circ/2};
//Point(55) = {width, height, 0, l_circ/2};
Point(6) = {x_scale * factor, 0, 0, l_box};
Point(7) = {x_scale * factor,1.094* factor, 0, l_box};
Point(8) = {-x_scale*factor, 1.094*factor, 0, l_box};

Point(100) = {0.0994137, 0, 0, l_circ};
Point(101) = {0.0994051, 0.0563143, 0, l_box};
Point(102) = {0.0993837, 0.11261, 0, l_box};
Point(103) = {0.0993836, 0.168854, 0, l_box};
Point(104) = {0.0995375, 0.224997, 0, l_box};
Point(105) = {0.100233, 0.280979, 0, l_box};
Point(106) = {0.10205, 0.33676, 0, l_box};
Point(107) = {0.104348, 0.392421, 0, l_box};
Point(108) = {0.101783, 0.448247, 0, l_box};
Point(109) = {0.085065, 0.497129, 0, l_box};
Point(110) = {0.0493049, 0.534306, 0, l_box};
Point(111) = {-1.67576e-06, 0.548296, 0, l_box};
Point(112) = {-0.0493021, 0.534302, 0, l_box};
Point(113) = {-0.0850662, 0.497129, 0, l_box};
Point(114) = {-0.101787, 0.448247, 0, l_box};
Point(115) = {-0.104348, 0.392424, 0, l_box};
Point(116) = {-0.102049, 0.33676, 0, l_box};
Point(117) = {-0.100233, 0.280979, 0, l_box};
Point(118) = {-0.0995371, 0.224997, 0, l_box};
Point(119) = {-0.0993832, 0.168854, 0, l_box};
Point(120) = {-0.0993832, 0.11261, 0, l_box};
Point(121) = {-0.0994047, 0.0563143, 0, l_box};
Point(122) = {-0.0994132, 0, 0, l_circ};

Line(100) = {100,101};
Line(101) = {101,102};
Line(102) = {102,103};
Line(103) = {103,104};
Line(104) = {104,105};
Line(105) = {105,106};
Line(106) = {106,107};
Line(107) = {107,108};
Line(108) = {108,109};
Line(109) = {109,110};
Line(110) = {110,111};
Line(111) = {111,112};
Line(112) = {112,113};
Line(113) = {113,114};
Line(114) = {114,115};
Line(115) = {115,116};
Line(116) = {116,117};
Line(117) = {117,118};
Line(118) = {118,119};
Line(119) = {119,120};
Line(120) = {120,121};
Line(121) = {121,122};
//Circle(11) = {5, 2, 22};
//Circle(12) = {22, 2, 55};

Line(13) = {1, 122};
Line(14) = {100, 6};
Line(15) = {6, 7};
Line(16) = {7, 8};
Line(17) = {8, 1};
Line(61) = {100,122};
//Line(20) = {3, 5};
//Line(21) = {55, 4};
//Line(22) = {4, 3};

Physical Line(1) = {122}; //13, 61, 14};
//Physical Line(2) = {11, 12, 20, 21};
Physical Line(3) = {16};
Physical Line(4) = {15,148};
Physical Line(5) = {17,138};

Line Loop(21) = {100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,-13,-17, -16, -15,-14};
Line Loop(22) = {100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,-61};
Plane Surface(31) = {21};
Plane Surface(32) = {22};

//+
Symmetry {0, 1, 0, 0} {
  Duplicata { Point{112}; Point{8}; Point{7}; Point{111}; Point{110}; Point{113}; Point{109}; Point{114}; Point{108}; Point{115}; Point{107}; Point{116}; Point{106}; Point{117}; Point{105}; Point{118}; Point{104}; Point{119}; Point{1}; Point{103}; Point{120}; Point{102}; Point{121}; Point{101}; Point{122}; Point{100}; Point{6}; Line{16}; Line{111}; Line{110}; Line{112}; Line{109}; Line{113}; Line{108}; Line{114}; Line{107}; Line{115}; Line{106}; Line{116}; Line{105}; Line{117}; Line{104}; Line{118}; Line{17}; Line{103}; Line{119}; Line{102}; Line{120}; Line{101}; Line{121}; Line{13}; Line{100}; Line{61}; Line{15}; Line{14}; }
}
//+
Line Loop(23) = {138, 13, -144, -142, -140, -137, -135, -133, -131, -129, -127, -125, -123, -124, -126, -128, -130, -132, -134, -136, -139, -141, -143, -146, 14, 148, 122};
//+
Plane Surface(33) = {23};
//+
Line Loop(24) = {144, -61, 146, 143, 141, 139, 136, 134, 132, 130, 128, 126, 124, 123, 125, 127, 129, 131, 133, 135, 137, 140, 142};
//+
Plane Surface(34) = {24};

Physical Surface(0) = {31,33};
Physical Surface(1) = {32,34};
