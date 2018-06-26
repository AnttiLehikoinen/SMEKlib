//Units
cm= 1/100; //centimeters
lc = 0.2*cm;

//Dimensions
w = 6*cm;
h = w;
hy = 2*cm;
wmag = 2*cm;


//Points
Point(0) = {0,0,0,lc};
Point(1) = {w,0,0,lc};
Point(2) = {w,hy,0,lc};
Point(3) = {w-wmag,hy,0,lc};
Point(4) = {hy,hy,0,lc};
Point(5) = {hy,h-hy,0,lc};
Point(6) = {w-wmag,h-hy,0,lc};
Point(7) = {w,h-hy,0,lc};
Point(8) = {w,h,0,lc};
Point(9) = {0,h,0,lc};

//Lines
Line(0) = {0,1};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,0};

lm1 = newl; Line(lm1) = {2,7};
lm2 = newl; Line(lm2) = {3,6};

//Loops
Line Loop(0) = {0,1,2,3,4,5,6,7,8,9};
mag = newll;
Line Loop(mag) = {-2,lm1,-6,-lm2};
air = newll; Line Loop(air) = {3,4,5,-lm2};

//Physical surfaces
Plane Surface(0) = {0};
Physical Surface("Core") = {0};

Plane Surface(1) = {mag};
Physical Surface("Magnet") = (1);

Plane Surface(2) = {air};
Physical Surface("Air") = (2);

Physical Line("Out") = {0,1,lm1,7,8,9};