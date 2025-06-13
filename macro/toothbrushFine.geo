Mesh.MshFileVersion = 2.2;
width = 3.5;
widthy = 0.0;
height = 31.5;
leftright = 1.3; // distance of the first an last membrane line from left and right boundary

l_circ = 3.5;
l_box = 1;

// first rectangle
Point(1) = {0, 0, 0, l_box};
Point(2) = {-leftright*2*width, leftright*widthy, 0, l_box}; //circle-center
Point(3) = {-leftright*2*width, leftright*widthy + height, 0, l_circ};
Point(4) = {0, height, 0, l_circ};

// additional points for second rectangle
Point(5) = {width, -widthy, 0, l_box};
Point(6) = {width, height - widthy, 0, l_circ};

k = 2;
j = 0;
For i In {7:51:2}
  Point(i) = {k*width, -k*widthy, 0, l_box};
  Point(i+1) = {k*width, height - k*widthy, 0, l_circ};

  Line(i+k-1) = {i, i-2};
  Line(i+k) = {i-1, i+1};
  Line(i+k+1) = {i+1, i};
  k = k+1;
EndFor

// ninth rectangle
Point(i) = {(k+leftright)*width, -(k+leftright)*widthy, 0, l_box};
Point(i+1) = {(k+leftright)*width, height - (k+leftright)*widthy, 0, l_circ};

// upper points
Point(i+2) = {-leftright*2*width, 2*height + leftright*widthy, 0, l_circ};
Point(i+3) = {(k+leftright)*width, 2*height - (k+leftright)*widthy, 0, l_circ};


// first rectangle
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// additional lines for second rectangle
Line(5) = {5, 1};
Line(6) = {4, 6};
Line(7) = {6, 5};

// ninth rectangle
Line(i+k-1) = {i, i-2};
Line(i+k) = {i-1, i+1};
Line(i+k+1) = {i+1, i};

//outer lines
Line(i+k+2) = {3, i+2};
Line(i+k+3) = {i+2, i+3};
Line(i+k+4) = {i+3, i+1};

Physical Line(1) = {1,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77};
Physical Line(3) = {i+k+2,2};
Physical Line(4) = {i+k+4,i+k+1};
Physical Line(5) = {i+k+3};

Line Loop(1) = {1,2,3,4};
s = 2;
For l In {5:i+k+1:3}
  Line Loop(s) = {l,-l+1,l+1,l+2};
  Plane Surface(s) = {s};
  s = s+1;
EndFor
Plane Surface(1) = {1};

Line Loop(100) = {-(i+k+4),-(i+k+3),-(i+k+2),3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78};
Plane Surface(100) = {100};

Physical Surface(0) = {2,5,8,11,14,17,20,23,26};
Physical Surface(1) = {1,4,7,10,13,16,19,22,25};
Physical Surface(2) = {100,3,6,9,12,15,18,21,24};

