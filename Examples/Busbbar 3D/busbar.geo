lc = 0.07; //characteristic length

//bounding box
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1,2,3,4};

//conductor
lcc = 0.5*lc;
Point(5) = {0.4,0.4,0, lcc}; 
Point(6) = {0.6,0.4,0, lcc}; 
Point(7) = {0.6,0.6,0, lcc}; 
Point(8) = {0.4,0.6,0, lcc}; 
Line(5) = {5, 6}; 
Line(6) = {6, 7}; 
Line(7) = {7, 8}; 
Line(8) = {8, 5};
Line Loop(2) = {5,6,7,8};

//2D geometry
Plane Surface(1) = {1,2}; //bounding box minus the conductor inside
Plane Surface(2) = {2}; //conductor

//extruding geometry
ext_air[] = Extrude {0, 0, 1} { Surface{1}; };
ext_bar[] = Extrude {0, 0, 1} { Surface{2}; };

Physical Volume("Air") = ext_air[1];
Physical Volume("Bar") = ext_bar[1];

//outer surface
surf_out[] = CombinedBoundary{ Volume{ext_air[1]}; Volume{ext_bar[1]}; };
Physical Surface("Outer") = {surf_out[]};