function [p,t,m,FL,LL,cir,ag] = calculate_stator6(dim)

p = zeros(36,2);
t = zeros(52,3);
m = zeros(52,1);

line = [dim.Sin dim.Sin+dim.S_height1 dim.Sin+dim.S_height12 dim.Sin+dim.S_height-dim.S_height3 ...
        dim.Sin+dim.S_height-dim.S_height4 dim.Sin+dim.S_height-dim.S_height4/2 dim.Sin+dim.S_height ...
        (dim.Sin+dim.S_height)*0.75+dim.Sout*0.25 dim.Sout];


p(1,1) = line(1)*cos(dim.angleS(1)/2);
p(1,2) = line(1)*sin(dim.angleS(1)/2);

p(2,1) = line(1)*cos(dim.angleS(1)/2-atan(dim.S_width1/2/line(1)));
p(2,2) = line(1)*sin(dim.angleS(1)/2-atan(dim.S_width1/2/line(1)));

p(3,1) = line(2)*cos(dim.angleS(1)/2-atan(dim.S_width1/2/line(2)));
p(3,2) = line(2)*sin(dim.angleS(1)/2-atan(dim.S_width1/2/line(2)));
    
p(4,1) = line(3)*cos(dim.angleS(1)/2-atan(dim.S_width3/2/line(3)));
p(4,2) = line(3)*sin(dim.angleS(1)/2-atan(dim.S_width3/2/line(3)));

p(5,1) = line(4)*cos(dim.angleS(1)/2-atan(dim.S_width2/2/line(4)));
p(5,2) = line(4)*sin(dim.angleS(1)/2-atan(dim.S_width2/2/line(4)));

p(6,1) = line(5)*cos(dim.angleS(1)/2-atan(dim.S_width2/2/line(5)));
p(6,2) = line(5)*sin(dim.angleS(1)/2-atan(dim.S_width2/2/line(5)));

p(7,1) = line(6)*cos(dim.angleS(1)/2-atan(dim.S_width2/2/line(6)));
p(7,2) = line(6)*sin(dim.angleS(1)/2-atan(dim.S_width2/2/line(6)));

p(8,1) = line(7)*cos(dim.angleS(1)/2-atan(dim.S_width2/2/line(7)));
p(8,2) = line(7)*sin(dim.angleS(1)/2-atan(dim.S_width2/2/line(7)));

p(9,1) = line(7)*cos(dim.angleS(1)/2);
p(9,2) = line(7)*sin(dim.angleS(1)/2);

p(10,1) = line(7)*cos(dim.angleS(1)/2+atan(dim.S_width2/2/line(7)));
p(10,2) = line(7)*sin(dim.angleS(1)/2+atan(dim.S_width2/2/line(7)));

p(11,1) = line(6)*cos(dim.angleS(1)/2+atan(dim.S_width2/2/line(6)));
p(11,2) = line(6)*sin(dim.angleS(1)/2+atan(dim.S_width2/2/line(6)));

p(12,1) = line(5)*cos(dim.angleS(1)/2+atan(dim.S_width2/2/line(5)));
p(12,2) = line(5)*sin(dim.angleS(1)/2+atan(dim.S_width2/2/line(5)));

p(13,1) = line(4)*cos(dim.angleS(1)/2+atan(dim.S_width2/2/line(4)));
p(13,2) = line(4)*sin(dim.angleS(1)/2+atan(dim.S_width2/2/line(4)));

p(14,1) = line(3)*cos(dim.angleS(1)/2+atan(dim.S_width3/2/line(3)));
p(14,2) = line(3)*sin(dim.angleS(1)/2+atan(dim.S_width3/2/line(3)));

p(15,1) = line(2)*cos(dim.angleS(1)/2+atan(dim.S_width1/2/line(2)));
p(15,2) = line(2)*sin(dim.angleS(1)/2+atan(dim.S_width1/2/line(2)));

p(16,1) = line(1)*cos(dim.angleS(1)/2+atan(dim.S_width1/2/line(1)));
p(16,2) = line(1)*sin(dim.angleS(1)/2+atan(dim.S_width1/2/line(1)));

p(17,1) = line(2)*cos(dim.angleS(1)/2);
p(17,2) = line(2)*sin(dim.angleS(1)/2);

p(18,1) = line(4)*cos(dim.angleS(1)/2);
p(18,2) = line(4)*sin(dim.angleS(1)/2);

p(19,1) = line(5)*cos(dim.angleS(1)/2);
p(19,2) = line(5)*sin(dim.angleS(1)/2);

p(20,1) = line(8)*cos(dim.angleS(1)/2+atan(dim.S_width2/2/line(8)));
p(20,2) = line(8)*sin(dim.angleS(1)/2+atan(dim.S_width2/2/line(8)));

p(21,1) = line(8)*cos(dim.angleS(1)/2-atan(dim.S_width2/2/line(8)));
p(21,2) = line(8)*sin(dim.angleS(1)/2-atan(dim.S_width2/2/line(8)));

p(22,1) = line(1)*cos(dim.angleS(1)/2-atan(dim.S_width3/2/line(1)));
p(22,2) = line(1)*sin(dim.angleS(1)/2-atan(dim.S_width3/2/line(1)));

p(23,1) = line(1)*cos(dim.angleS(1)/2+atan(dim.S_width3/2/line(1)));
p(23,2) = line(1)*sin(dim.angleS(1)/2+atan(dim.S_width3/2/line(1)));

p(24,1) = line(1)*cos(dim.angleS(1));
p(24,2) = line(1)*sin(dim.angleS(1));

p(25,1) = line(3)*cos(dim.angleS(1));
p(25,2) = line(3)*sin(dim.angleS(1));

p(26,1) = line(6)*cos(dim.angleS(1));
p(26,2) = line(6)*sin(dim.angleS(1));

p(27,1) = line(7)*cos(dim.angleS(1));
p(27,2) = line(7)*sin(dim.angleS(1));

p(28,1) = line(8)*cos(dim.angleS(1));
p(28,2) = line(8)*sin(dim.angleS(1));

p(29,1) = line(9)*cos(dim.angleS(1));
p(29,2) = line(9)*sin(dim.angleS(1));

p(30,1) = line(9)*cos(dim.angleS(1)/2);
p(30,2) = line(9)*sin(dim.angleS(1)/2);

p(31,1) = line(9);

p(32,1) = line(8);

p(33,1) = line(7);

p(34,1) = line(6);

p(35,1) = line(3);

p(36,1) = line(1);

t(1,:) = [2 22 3];
m(1) = dim.SM;

t(2,:) = [22 4 3];
m(2) = dim.SM;

t(3,:) = [22 36 4];
m(3) = dim.SM;

t(4,:) =[36 35 4];
m(4) = dim.SM;

t(5,:) = [35 5 4];
m(5) = dim.SM;

t(6,:) =  [35 6 5];
m(6) = dim.SM;

t(7,:) = [35 34 6];
m(7) = dim.SM;

t(8,:) = [34 7 6];
m(8) = dim.SM;

t(9,:) = [34 33 7];
m(9) = dim.SM;

t(10,:) = [33 8 7];
m(10) = dim.SM;

t(11,:) = [33 32 8];
m(11) = dim.SM;

t(12,:) = [32 21 8];
m(12) = dim.SM;

t(13,:) = [21 9 8];
m(13) = dim.SM;

t(14,:) = [21 20 9];
m(14) = dim.SM;

t(15,:) = [20 10 9];
m(15) = dim.SM;

t(16,:) = [20 28 10];
m(16) = dim.SM;

t(17,:) = [28 27 10];
m(17) = dim.SM;

t(18,:) = [27 11 10];
m(18) = dim.SM;

t(19,:) = [27 26 11];
m(19) = dim.SM;

t(20,:) = [26 12 11];
m(20) = dim.SM;

t(21,:) = [26 25 12];
m(21) = dim.SM;

t(22,:) = [25 13 12];
m(22) = dim.SM;

t(23,:) = [25 14 13];
m(23) = dim.SM;

t(24,:) = [25 24 14];
m(24) = dim.SM;

t(25,:) = [24 23 14];
m(25) = dim.SM;

t(26,:) = [23 15 14];
m(26) = dim.SM;

t(27,:) = [23 16 15];
m(27) = dim.SM;

t(28,:) = [21 32 31];
m(28) = dim.SM;

t(29,:) = [21 31 30];
m(29) = dim.SM;

t(30,:) = [20 21 30];
m(30) = dim.SM;

t(31,:) = [20 30 29];
m(31) = dim.SM;

t(32,:) = [20 29 28];
m(32) = dim.SM;

t(33,:) = [1 2 17];
m(33) = 9999;

t(34,:) = [2 3 17];
m(34) = 9999;

t(35,:) = [1 17 16];
m(35) = 9999;

t(36,:) = [15 16 17];
m(36) = 9999;

t(37,:) = [3 4 5];
m(37) = 9999;

t(38,:) = [3 5 17];
m(38) = 9999;

t(39,:) = [5 18 17];
m(39) = 9999;

t(40,:) = [17 18 13];
m(40) = 9999;

t(41,:) = [17 13 15];
m(41) = 9999;

t(42,:) = [13 14 15];
m(42) = 9999;

t(43,:) = [5 6 19];
m(43) = 9999;

t(44,:) = [5 19 18];
m(44) = 9999;

t(45,:) = [18 19 13];
m(45) = 9999;

t(46,:) = [19 12 13];
m(46) = 9999;

t(47,:) = [6 7 19];
m(47) = dim.SSM1;

t(48,:) = [19 7 11];
m(48) = dim.SSM1;

t(49,:) = [19 11 12];
m(49) = dim.SSM1;

t(50,:) = [7 8 9];
m(50) = dim.SSM2;

t(51,:) = [7 9 11];
m(51) = dim.SSM2;

t(52,:) = [9 10 11];
m(52) = dim.SSM2;

FL = [36 35 34 33 32 31];
LL = [24 25 26 27 28 29];
cir = [31 30 29];
ag = [36 22 2 1 16 23 24];

end