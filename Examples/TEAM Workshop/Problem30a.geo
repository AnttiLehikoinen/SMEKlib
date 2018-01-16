
//Units
cm= 1/100; //centimeters

//Dimensions
r1 = 2*cm;
r2 = 3*cm;
r3 = 3.2*cm;
r4 = 5.2*cm;
r5 = 5.7*cm;

//Mesh Control
mr1 = 0.3*cm;
mr2 = 0.05*cm;
mr3 = 0.08*cm;
mr4 = 0.5*cm;
mr5 = 0.2*cm;
//+

//Points r1
r1Point[]= newp; Point(newp) = {0, 0, 0, mr1};
r1Point[]+= newp; Point(newp) = {r1, 0, 0, mr1};
r1Point[]+= newp; Point(newp) = {0, r1, 0, mr1};
r1Point[]+= newp; Point(newp) = {-r1, 0, 0, mr1};
r1Point[]+= newp; Point(newp) = {0, -r1, 0, mr1};

//Lines r1
r1Line[] = newl; Circle(newl) = {r1Point[1], r1Point[0], r1Point[2]};
r1Line[]+= newl; Circle(newl) = {r1Point[2], r1Point[0], r1Point[3]};
r1Line[]+= newl; Circle(newl) = {r1Point[3], r1Point[0], r1Point[4]};
r1Line[]+= newl; Circle(newl) = {r1Point[4], r1Point[0], r1Point[1]};

//Surface r1
r1LineLoop[]= {};
r1Surface[]= {};
r1LineLoop[] += newll;  Line Loop(newll) = {r1Line[0], r1Line[1], r1Line[2], r1Line[3]};
r1Surface[]  += news;  Plane Surface(news) = {r1LineLoop[]};


//+

//Points r2
//r2Point[]= newp; Point(newp) = {0, 0, 0, mr2};
r2Point[]+= newp; Point(newp) = {r2, 0, 0, mr2};
r2Point[]+= newp; Point(newp) = {0, r2, 0, mr2};
r2Point[]+= newp; Point(newp) = {-r2, 0, 0, mr2};
r2Point[]+= newp; Point(newp) = {0, -r2, 0, mr2};

//Lines r2
r2Line[] = newl; Circle(newl) = {r2Point[0], r1Point[0], r2Point[1]};
r2Line[]+= newl; Circle(newl) = {r2Point[1], r1Point[0], r2Point[2]};
r2Line[]+= newl; Circle(newl) = {r2Point[2], r1Point[0], r2Point[3]};
r2Line[]+= newl; Circle(newl) = {r2Point[3], r1Point[0], r2Point[0]};

//Surface r2
r2LineLoop[]= {};
r2Surface[]= {};
r2LineLoop[] += newll;  Line Loop(newll) = {r2Line[0], r2Line[1], r2Line[2], r2Line[3]};
r2Surface[]  += news;  Plane Surface(news) = {r2LineLoop[],r1LineLoop[]};


//+ AirGap

//Points (r2+r3)/2
//r32Point[]= newp; Point(newp) = {0, 0, 0, mr3};
r32Point[]+= newp; Point(newp) = {(r2+r3)/2, 0, 0, mr3};
r32Point[]+= newp; Point(newp) = {0, (r2+r3)/2, 0, mr3};
r32Point[]+= newp; Point(newp) = {-(r2+r3)/2, 0, 0, mr3};
r32Point[]+= newp; Point(newp) = {0, -(r2+r3)/2, 0, mr3};

//Lines (r2+r3)/2
r32Line[] = newl; Circle(newl) = {r32Point[0], r1Point[0], r32Point[1]};
r32Line[]+= newl; Circle(newl) = {r32Point[1], r1Point[0], r32Point[2]};
r32Line[]+= newl; Circle(newl) = {r32Point[2], r1Point[0], r32Point[3]};
r32Line[]+= newl; Circle(newl) = {r32Point[3], r1Point[0], r32Point[0]};

//Surface r32
r32LineLoop[]= {};
r32Surface[]= {};
r32LineLoop[] += newll;  Line Loop(newll) = {r32Line[0], r32Line[1], r32Line[2], r32Line[3]};
r32Surface[]  += news;  Plane Surface(news) = {r32LineLoop[],r2LineLoop[]};

//Points r3
//r3Point[]= newp; Point(newp) = {0, 0, 0, mr3};
r3Point[]+= newp; Point(newp) = {r3, 0, 0, mr3};
r3Point[]+= newp; Point(newp) = {0, r3, 0, mr3};
r3Point[]+= newp; Point(newp) = {-r3, 0, 0, mr3};
r3Point[]+= newp; Point(newp) = {0, -r3, 0, mr3};

/*
//Lines r3
r3Line[] = newl; Circle(newl) = {r3Point[0], r1Point[0], r3Point[1]};
r3Line[]+= newl; Circle(newl) = {r3Point[1], r1Point[0], r3Point[2]};
r3Line[]+= newl; Circle(newl) = {r3Point[2], r1Point[0], r3Point[3]};
r3Line[]+= newl; Circle(newl) = {r3Point[3], r1Point[0], r3Point[0]};
*/

//+

//Points r4
r4Point[]+= newp; Point(newp) = {r4, 0, 0, mr4};
r4Point[]+= newp; Point(newp) = {0, r4, 0, mr4};
r4Point[]+= newp; Point(newp) = {-r4, 0, 0, mr4};
r4Point[]+= newp; Point(newp) = {0, -r4, 0, mr4};

/*
//Lines r4
r4Line[] = newl; Circle(newl) = {r4Point[0], r1Point[0], r4Point[1]};
r4Line[]+= newl; Circle(newl) = {r4Point[1], r1Point[0], r4Point[2]};
r4Line[]+= newl; Circle(newl) = {r4Point[2], r1Point[0], r4Point[3]};
r4Line[]+= newl; Circle(newl) = {r4Point[3], r1Point[0], r4Point[0]};
*/
//+

//Points r5
r5Point[]+= newp; Point(newp) = {r5, 0, 0, mr5};
r5Point[]+= newp; Point(newp) = {0, r5, 0, mr5};
r5Point[]+= newp; Point(newp) = {-r5, 0, 0, mr5};
r5Point[]+= newp; Point(newp) = {0, -r5, 0, mr5};

//Lines r5
r5Line[] = newl; Circle(newl) = {r5Point[0], r1Point[0], r5Point[1]};
r5Line[]+= newl; Circle(newl) = {r5Point[1], r1Point[0], r5Point[2]};
r5Line[]+= newl; Circle(newl) = {r5Point[2], r1Point[0], r5Point[3]};
r5Line[]+= newl; Circle(newl) = {r5Point[3], r1Point[0], r5Point[0]};



//Points coils

polePitch = Pi/3;
poles = 6;
polePairs = 2*poles;

For j In {0:poles}
coilPoint[]+= newp; Point(newp) = {r3*Cos[(Pi/4)/2+j*polePitch], r3*Sin[(Pi/4)/2+j*polePitch], 0, mr4};
coilPoint[]+= newp; Point(newp) = {r4*Cos[(Pi/4)/2+j*polePitch], r4*Sin[(Pi/4)/2+j*polePitch], 0, mr4};
coilPoint[]+= newp; Point(newp) = {r3*Cos[-(Pi/4)/2+j*polePitch], r3*Sin[-(Pi/4)/2+j*polePitch], 0, mr4};
coilPoint[]+= newp; Point(newp) = {r4*Cos[-(Pi/4)/2+j*polePitch], r4*Sin[-(Pi/4)/2+j*polePitch], 0, mr4};

coilLine[]= newl; Line(newl) = {coilPoint[0+4*j], coilPoint[1+4*j]};
coilLine[]+= newl; Line(newl) = {coilPoint[2+4*j], coilPoint[3+4*j]};

EndFor






/*
//Surface r3
r3LineLoop[]= {};
r3Surface[]= {};
r3LineLoop[] += newll;  Line Loop(newll) = {r3Line[0], r3Line[1], r3Line[2], r3Line[3]};
r3Surface[]  += news;  Plane Surface(news) = {r3LineLoop[],r32LineLoop[]};
*/


//+
Circle(37) = {14, 1, 26};
//+
Circle(38) = {26, 1, 32};
//+
Circle(39) = {32, 1, 30};
//+
Circle(40) = {30, 1, 15};
//+
Circle(41) = {15, 1, 36};
//+
Circle(42) = {36, 1, 34};
//+
Circle(43) = {34, 1, 40};
//+
Circle(44) = {40, 1, 16};
//+
Circle(45) = {16, 1, 38};
//+
Circle(46) = {38, 1, 44};
//+
Circle(47) = {44, 1, 42};
//+
Circle(48) = {42, 1, 17};
//+
Circle(49) = {17, 1, 48};
//+
Circle(50) = {48, 1, 46};
//+
Circle(51) = {46, 1, 28};
//+
Circle(52) = {28, 1, 14};
//+
Circle(53) = {18, 1, 27};
//+
Circle(54) = {27, 1, 33};
//+
Circle(55) = {33, 1, 31};
//+
Circle(56) = {31, 1, 19};
//+
Circle(57) = {19, 1, 37};
//+
Circle(58) = {37, 1, 35};
//+
Circle(59) = {35, 1, 41};
//+
Circle(60) = {41, 1, 20};
//+
Circle(61) = {20, 1, 39};
//+
Circle(62) = {39, 1, 45};
//+
Circle(63) = {45, 1, 43};
//+
Circle(64) = {43, 1, 21};
//+
Circle(65) = {21, 1, 49};
//+
Circle(66) = {49, 1, 47};
//+
Circle(67) = {47, 1, 29};
//+
Circle(68) = {29, 1, 18};
//+
Line Loop(18) = {37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52};
//+
Plane Surface(19) = {17, 18};
//+
Line Loop(19) = {53, -23, -37, -52, 24, 68};
//+
Plane Surface(20) = {19};
//+
Line Loop(20) = {23, 54, -26, -38};
//+
Plane Surface(21) = {20};
//+
Line Loop(21) = {26, 55, -25, -39};
//+
Plane Surface(22) = {21};
//+
Line Loop(22) = {40, 41, 28, -57, -56, -25};
//+
Plane Surface(23) = {22};
//+
Line Loop(23) = {42, 27, -58, -28};
//+
Plane Surface(24) = {23};
//+
Line Loop(24) = {43, 30, -59, -27};
//+
Plane Surface(25) = {24};
//+
Line Loop(25) = {44, 45, 29, -61, -60, -30};
//+
Plane Surface(26) = {25};
//+
Line Loop(26) = {46, 32, -62, -29};
//+
Plane Surface(27) = {26};
//+
Line Loop(27) = {31, -63, -32, 47};
//+
Plane Surface(28) = {27};
//+
Line Loop(28) = {34, -65, -64, -31, 48, 49};
//+
Plane Surface(29) = {28};
//+
Line Loop(29) = {66, -33, -50, 34};
//+
Plane Surface(30) = {29};
//+
Line Loop(30) = {67, -24, -51, 33};
//+
Plane Surface(31) = {30};
//+

//+
Line Loop(31) = {19, 20, 21, 22};
//+
Line Loop(32) = {55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 53, 54};
//+
Plane Surface(32) = {31, 32};
//+
Physical Surface("RotorSteel", 1) = {6};
//+
Physical Surface("RotorAluminum", 2) = {12};
//+
Physical Surface("A", 3) = {20};
//+
Physical Surface("Aminus", 4) = {26};
//+
Physical Surface("B", 5) = {28};
//+
Physical Surface("Bminus", 6) = {22};
//+
Physical Surface("C", 7) = {24};
//+
Physical Surface("Cminus", 8) = {30};
//+
Physical Surface("StatorSteel", 9) = {32};
//+
Physical Surface("AG_rotor", 10) = {18};
//+
Physical Surface("AG_stator", 11) = {19};
//+
Physical Surface("AG_coils", 12) = {21, 23, 25, 27, 29, 31};
//+
Transfinite Line {37} = 10 Using Progression 1;
//+
Transfinite Line {52, 44, 45} = 10 Using Progression 1;
//+
Transfinite Line {39, 42, 47, 50} = 20 Using Progression 1;
//+
Transfinite Line {38, 43, 46, 51} = 6 Using Progression 1;
//+
Transfinite Line {40, 41, 49, 48} = 3 Using Progression 1;
//+
Transfinite Line {14, 8, 9, 15, 16, 10, 7, 13} = 42 Using Progression 1;
//+
Transfinite Line {58, 55, 63, 66} = 20 Using Progression 1;
//+
Transfinite Line {59, 54, 67, 62} = 6 Using Progression 1;
//+
Transfinite Line {57, 56, 64, 65} = 3 Using Progression 1;
//+
Transfinite Line {53, 68, 60, 61} = 10 Using Progression 1;
//+
//+
Transfinite Line {23, 26, 25, 28, 27, 30, 29, 32, 31, 34, 33, 24} = 10 Using Progression 1;

//+
SetFactory("OpenCASCADE");
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {-0, -0, -0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {-0, -0, -0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0e, 2*Pi};
//+
Circle(69) = {0, 0, 0, 0.067, 0, 2*Pi};
//+
Line Loop(33) = {69};
//+
Surface(33) = {31, 33};
//+
Line Loop(34) = {69};
//+
Plane Surface(33) = {34};
//+
Recursive Delete {
  Surface{33}; 
}
//+
Circle(69) = {0, 0, 0, 0.067, 0, 2*Pi};
//+
Line Loop(34) = {69};
//+
Plane Surface(33) = {31, 34};
//+
Line Loop(35) = {69};
//+
Plane Surface(33) = {31, 35};
//+
Line Loop(36) = {69};
//+
Plane Surface(33) = {31, 36};
//+
Line Loop(37) = {69};
//+
Plane Surface(33) = {31, 37};
//+
Line Loop(38) = {69};
//+
Line Loop(39) = {69};
//+
Plane Surface(33) = {31, 39};
