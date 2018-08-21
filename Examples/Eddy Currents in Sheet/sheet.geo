lc = 0.1E-3; //characteristic length
h = 1.0E-3;
w = 15E-3;

//bounding box
Point(1) = {0, 0, 0, lc};
Point(2) = {w, 0, 0, lc};
Point(3) = {w, h, 0, lc};
Point(4) = {0, h, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1,2,3,4};
Surface(1) = {1};

Physical Surface("Iron") = {1};
Physical Line("Out") = {1,2,3,4};