function [p,t,m,FL,LL,cir,ag] = calculate_stator3(dim)

%preallocate memory
p = zeros(32,2);
t = zeros(44,3);
m = zeros(44,1);


line = [dim.Sin dim.Sin/2+(dim.Sin+dim.S_height-dim.S_height3)/2 dim.Sin+dim.S_height-dim.S_height3 ...
        dim.Sin+dim.S_height-dim.S_height3/2 dim.Sin+dim.S_height 0.75*(dim.Sin+dim.S_height)+0.25*dim.Sout dim.Sout];

%Init the first layer%

p(1,1) = line(1)*cos(dim.angleS(1)/2);
p(1,2) = line(1)*sin(dim.angleS(1)/2);

p(2,1) = line(1)*cos(dim.angleS(1)/2-atan(dim.S_width1/2/line(1)));
p(2,2) = line(1)*sin(dim.angleS(1)/2-atan(dim.S_width1/2/line(1)));

p(3,1) = line(2)*cos(dim.angleS(1)/2-atan(dim.S_width3/2/line(2)));
p(3,2) = line(2)*sin(dim.angleS(1)/2-atan(dim.S_width3/2/line(2)));

p(4,1) = line(3)*cos(dim.angleS(1)/2-atan(dim.S_width2/2/line(3)));
p(4,2) = line(3)*sin(dim.angleS(1)/2-atan(dim.S_width2/2/line(3)));

p(5,1) = line(4)*cos(dim.angleS(1)/2-atan(dim.S_width2/2/line(4)));
p(5,2) = line(4)*sin(dim.angleS(1)/2-atan(dim.S_width2/2/line(4)));

p(6,1) = line(5)*cos(dim.angleS(1)/2-atan(dim.S_width2/2/line(5)));
p(6,2) = line(5)*sin(dim.angleS(1)/2-atan(dim.S_width2/2/line(5)));

p(7,1) = line(5)*cos(dim.angleS(1)/2);
p(7,2) = line(5)*sin(dim.angleS(1)/2);

p(8,1) = line(5)*cos(dim.angleS(1)/2+atan(dim.S_width2/2/line(5)));
p(8,2) = line(5)*sin(dim.angleS(1)/2+atan(dim.S_width2/2/line(5)));

p(9,1) = line(4)*cos(dim.angleS(1)/2+atan(dim.S_width2/2/line(4)));
p(9,2) = line(4)*sin(dim.angleS(1)/2+atan(dim.S_width2/2/line(4)));

p(10,1) = line(3)*cos(dim.angleS(1)/2+atan(dim.S_width2/2/line(3)));
p(10,2) = line(3)*sin(dim.angleS(1)/2+atan(dim.S_width2/2/line(3)));

p(11,1) = line(2)*cos(dim.angleS(1)/2+atan(dim.S_width3/2/line(2)));
p(11,2) = line(2)*sin(dim.angleS(1)/2+atan(dim.S_width3/2/line(2)));

p(12,1) = line(1)*cos(dim.angleS(1)/2+atan(dim.S_width1/2/line(1)));
p(12,2) = line(1)*sin(dim.angleS(1)/2+atan(dim.S_width1/2/line(1)));

p(13,1) = line(2)*cos(dim.angleS(1)/2);
p(13,2) = line(2)*sin(dim.angleS(1)/2);

p(14,1) = line(3)*cos(dim.angleS(1)/2);
p(14,2) = line(3)*sin(dim.angleS(1)/2);

p(15,1) = line(4)*cos(dim.angleS(1)/2);
p(15,2) = line(4)*sin(dim.angleS(1)/2);

p(16,1) = line(6)*cos(dim.angleS(1)/2+atan(dim.S_width2/2/line(6)));
p(16,2) = line(6)*sin(dim.angleS(1)/2+atan(dim.S_width2/2/line(6)));

p(17,1) = line(6)*cos(dim.angleS(1)/2-atan(dim.S_width2/2/line(6)));
p(17,2) = line(6)*sin(dim.angleS(1)/2-atan(dim.S_width2/2/line(6)));

p(18,1) = line(1)*cos(dim.angleS(1)/2-atan(dim.S_width3/2/line(1)));
p(18,2) = line(1)*sin(dim.angleS(1)/2-atan(dim.S_width3/2/line(1)));

p(19,1) = line(1)*cos(dim.angleS(1)/2+atan(dim.S_width3/2/line(1)));
p(19,2) = line(1)*sin(dim.angleS(1)/2+atan(dim.S_width3/2/line(1)));

p(20,1) = line(1)*cos(dim.angleS(1));
p(20,2) = line(1)*sin(dim.angleS(1));

p(21,1) = line(2)*cos(dim.angleS(1));
p(21,2) = line(2)*sin(dim.angleS(1));

p(22,1) = line(4)*cos(dim.angleS(1));
p(22,2) = line(4)*sin(dim.angleS(1));

p(23,1) = line(5)*cos(dim.angleS(1));
p(23,2) = line(5)*sin(dim.angleS(1));

p(24,1) = line(6)*cos(dim.angleS(1));
p(24,2) = line(6)*sin(dim.angleS(1));

p(25,1) = line(7)*cos(dim.angleS(1));
p(25,2) = line(7)*sin(dim.angleS(1));

p(26,1) = line(7)*cos(dim.angleS(1)/2);
p(26,2) = line(7)*sin(dim.angleS(1)/2);

p(27,1) = line(7);

p(28,1) = line(6);

p(29,1) = line(5);

p(30,1) = line(4);

p(31,1) = line(2);

p(32,1) = line(1);


t(1,:) = [2 18 3];
m(1) = dim.SM;

t(2,:) = [18 32 3];
m(2) = dim.SM;

t(3,:) = [32 31 3];
m(3) = dim.SM;

t(4,:) = [31 4 3];
m(4) = dim.SM;

t(5,:) = [31 30 4];
m(5) = dim.SM;

t(6,:) = [30 5 4];
m(6) = dim.SM;

t(7,:) = [30 29 5];
m(7) = dim.SM;

t(8,:) = [29 6 5];
m(8) = dim.SM;

t(9,:) = [29 28 17];
m(9) = dim.SM;

t(10,:) = [29 17 6];
m(10) = dim.SM;

t(11,:) = [17 7 6];
m(11) = dim.SM;

t(12,:)  = [17 16 7];
m(12)  = dim.SM;

t(13,:) = [16 8 7];
m(13) = dim.SM;

t(14,:) = [16 24 23];
m(14) = dim.SM;

t(15,:) = [16 23 8];
m(15) = dim.SM;

t(16,:) = [23 9 8];
m(16) = dim.SM;

t(17,:) = [23 22 9];
m(17) = dim.SM;

t(18,:) = [22 10 9];
m(18) = dim.SM;

t(19,:) = [22 21 10];
m(19) = dim.SM;

t(20,:) = [21 11 10];
m(20) = dim.SM;

t(21,:) = [21 20 11];
m(21) = dim.SM;

t(22,:) = [20 19 11];
m(22) = dim.SM;

t(23,:) = [19 12 11];
m(23) = dim.SM;

t(24,:) =[1 2 13];
m(24) = 9999;

t(25,:) = [2 3 13];
m(25) = 9999;

t(26,:) = [3 4 13];
m(26) = 9999;

t(27,:) = [4 14 13];
m(27) = 9999;

t(28,:) = [14 10 13];
m(28) = 9999;

t(29,:) = [10 11 13];
m(29) = 9999;

t(30,:) = [11 12 13];
m(30) = 9999;

t(31,:) = [12 1 13];
m(31) = 9999;

t(32,:) = [4 15 14];
m(32) = dim.SSM1;

t(33,:) = [4 5 15];
m(33) = dim.SSM1;

t(34,:) = [10 14 15];
m(34) = dim.SSM1;

t(35,:) = [9 10 15];
m(35) = dim.SSM1;

t(36,:) = [5 6 15];
m(36) = dim.SSM2;

t(37,:) = [6 7 15];
m(37) = dim.SSM2;

t(38,:) = [7 8 15];
m(38) = dim.SSM2;

t(39,:) = [8 9 15];
m(39) = dim.SSM2;

t(40,:) = [28 27 17];
m(40) = dim.SM;

t(41,:) = [27 26 17];
m(41) = dim.SM;

t(42,:) = [26 16 17];
m(42) = dim.SM;

t(43,:) = [26 25 16];
m(43) = dim.SM;

t(44,:) = [25 24 16];
m(44) = dim.SM;


FL = [32 31 30 29 28 27];
LL = [20 21 22 23 24 25];
cir = [27 26 25];
ag = [32 18 2 1 12 19 20];
end