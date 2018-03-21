function [p,t,m,FL,LL,cir,ag] = calculate_stator10(dim)

%preallocate memory
p = zeros(30,2);
t = zeros(40,3);
m = zeros(40,1);

%Help variables

line  = [dim.Sin 0.5*dim.Sin+0.5*(dim.Sin+dim.S_height/2) dim.Sin+0.5*dim.S_height dim.Sin+dim.S_height ...
        0.75*(dim.Sin+dim.S_height)+0.25*dim.Sout dim.Sout];
    
%Init the first sector

%Nodes
p(1,1) = line(1)*cos(dim.angleS(1)/2);
p(1,2) = line(1)*sin(dim.angleS(1)/2);

p(2,1) = line(1)*cos(dim.angleS(1)/2-atan(dim.S_width1/2/line(1)));
p(2,2) = line(1)*sin(dim.angleS(1)/2-atan(dim.S_width1/2/line(1)));

p(3,1) = line(2)*cos(dim.angleS(1)/2-atan(dim.S_width1/2/line(2)));
p(3,2) = line(2)*sin(dim.angleS(1)/2-atan(dim.S_width1/2/line(2)));

p(4,1) = line(3)*cos(dim.angleS(1)/2-atan(dim.S_width1/2/line(3)));
p(4,2) = line(3)*sin(dim.angleS(1)/2-atan(dim.S_width1/2/line(3)));

p(5,1) = line(4)*cos(dim.angleS(1)/2-atan(dim.S_width1/2/line(4)));
p(5,2) = line(4)*sin(dim.angleS(1)/2-atan(dim.S_width1/2/line(4)));

p(6,1) = line(4)*cos(dim.angleS(1)/2);
p(6,2) = line(4)*sin(dim.angleS(1)/2);

p(7,1) = line(4)*cos(dim.angleS(1)/2+atan(dim.S_width1/2/line(4)));
p(7,2) = line(4)*sin(dim.angleS(1)/2+atan(dim.S_width1/2/line(4)));

p(8,1) = line(3)*cos(dim.angleS(1)/2+atan(dim.S_width1/2/line(3)));
p(8,2) = line(3)*sin(dim.angleS(1)/2+atan(dim.S_width1/2/line(3)));

p(9,1) = line(2)*cos(dim.angleS(1)/2+atan(dim.S_width1/2/line(2)));
p(9,2) = line(2)*sin(dim.angleS(1)/2+atan(dim.S_width1/2/line(2)));

p(10,1) = line(1)*cos(dim.angleS(1)/2+atan(dim.S_width1/2/line(1)));
p(10,2) = line(1)*sin(dim.angleS(1)/2+atan(dim.S_width1/2/line(1)));

p(11,1) = line(2)*cos(6*dim.angleS(1)/7);
p(11,2) = line(2)*sin(6*dim.angleS(1)/7);

p(12,1) = line(2)*cos(dim.angleS(1)/2);
p(12,2) = line(2)*sin(dim.angleS(1)/2);

p(13,1) = line(2)*cos(dim.angleS(1)/7);
p(13,2) = line(2)*sin(dim.angleS(1)/7);

p(14,1) = line(3)*cos(dim.angleS(1)/2);
p(14,2) = line(3)*sin(dim.angleS(1)/2);

p(15,1) = line(5)*cos(dim.angleS(1)/2);
p(15,2) = line(5)*sin(dim.angleS(1)/2);

p(16,1) = line(1)*cos(dim.angleS(1)/7);
p(16,2) = line(1)*sin(dim.angleS(1)/7);

p(17,1) = line(1)*cos(6*dim.angleS(1)/7);
p(17,2) = line(1)*sin(6*dim.angleS(1)/7);

p(18,1) = line(1)*cos(dim.angleS(1));
p(18,2) = line(1)*sin(dim.angleS(1));

p(19,1) = line(2)*cos(dim.angleS(1));
p(19,2) = line(2)*sin(dim.angleS(1));

p(20,1) = line(3)*cos(dim.angleS(1));
p(20,2) = line(3)*sin(dim.angleS(1));

p(21,1) = line(4)*cos(dim.angleS(1));
p(21,2) = line(4)*sin(dim.angleS(1));

p(22,1) = line(5)*cos(dim.angleS(1));
p(22,2) = line(5)*sin(dim.angleS(1));

p(23,1) = line(6)*cos(dim.angleS(1));
p(23,2) = line(6)*sin(dim.angleS(1));

p(24,1) = line(6)*cos(dim.angleS(1)/2);
p(24,2) = line(6)*sin(dim.angleS(1)/2);

p(25,1) = line(6);

p(26,1) = line(5);

p(27,1) = line(4);

p(28,1) = line(3);

p(29,1) = line(2);

p(30,1) = line(1);




%Elements

t(1,:) = [2 13 3];
m(1) = 4;

t(2,:) = [2 16 13];
m(2) = 4;

t(3,:) = [16 30 13];
m(3) = 4;

t(4,:) = [30 29 13];
m(4) = 4;

t(5,:) = [3 13 4];
m(5) = 4;

t(6,:) = [13 28 4];
m(6) = 4;

t(7,:) = [13 29 28];
m(7) = 4;

t(8,:) = [4 28 5];
m(8) = 4;

t(9,:) = [28 27 5];
m(9) = 4;

t(10,:) = [27 26 5];
m(10) = 4;

t(11,:) = [5 26 15];
m(11) = 4;

t(12,:) = [6 5 15];
m(12) = 4;

t(13,:) = [15 26 24];
m(13) = 4;

t(14,:) = [26 25 24];
m(14) = 4;

t(15,:) = [24 23 22];
m(15) = 4;

t(16,:) = [24 22 15];
m(16) = 4; 

t(17,:) = [15 22 7];
m(17) = 4;

t(18,:) = [7 6 15];
m(18) = 4;

t(19,:) = [22 21 7];
m(19) = 4;

t(20,:) = [21 20 7];
m(20) = 4;

t(21,:) = [20 8 7];
m(21) = 4;

t(22,:) = [20 19 11];
m(22) = 4;

t(23,:) = [20 11 8];
m(23) = 4;

t(24,:) = [11 9 8];
m(24) = 4;

t(25,:) = [19 18 11];
m(25) = 4;

t(26,:) = [18 17 11];
m(26) = 4;

t(27,:) = [17 10 11];
m(27) = 4;

t(28,:) = [10 9 11];
m(28) = 4;

t(29,:) = [10 1 9];
m(29) = dim.SSM1;

t(30,:) = [1 12 9];
m(30) = dim.SSM1;

t(31,:) = [1 3 12];
m(31) = dim.SSM1;

t(32,:) = [1 2 3];
m(32) = dim.SSM1;

t(33,:) = [9 12 8];
m(33) = dim.SSM1;

t(34,:) = [12 14 8];
m(34) = dim.SSM1;

t(35,:) = [12 4 14];
m(35) = dim.SSM1;

t(36,:) = [12 3 4];
m(36) = dim.SSM1;

t(37,:) = [8 6 7];
m(37) = dim.SSM2;

t(38,:) = [8 14 6];
m(38) = dim.SSM2;

t(39,:) = [14 4 6];
m(39) = dim.SSM2;

t(40,:) = [4 5 6];
m(40) = dim.SSM2;


FL = [30 29 28 27 26 25];
LL = [18 19 20 21 22 23];
cir = [25 24 23];
ag = [30 16 2 1 10 17 18];















end