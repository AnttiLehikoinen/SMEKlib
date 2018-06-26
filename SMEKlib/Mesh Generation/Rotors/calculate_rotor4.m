function [p,t,m,FL,LL,ag] = calculate_rotor4(dim)

%Boundaries
FL = [1 39 38 37 36 35 34];
LL = [1 23 24 25 26 27 28];
ag = [34 33 32 31 30 29 28];

%pre-allocating memory
p = zeros(39,2);
t = zeros(58,3);
m = zeros(58,1);

line = [dim.Rin dim.Rin/2+(dim.Rout-dim.R_height1-dim.R_height3)/2 0.1*dim.Rin+0.9*(dim.Rout-dim.R_height1-dim.R_height3) ...
    dim.Rout-dim.R_height1-dim.R_height3 0.8*(dim.Rout-dim.R_height1-dim.R_height3)+0.2*(dim.Rout-dim.R_height1) ...
    (dim.Rout-dim.R_height1-dim.R_height3)/2+(dim.Rout-dim.R_height1)/2 0.2*(dim.Rout-dim.R_height1-dim.R_height3)+0.8*(dim.Rout-dim.R_height1) ...
    dim.Rout-dim.R_height1 dim.Rout-dim.R_height1/2 dim.Rout];

p(2,1) = line(2)*cos(dim.angleR(1)/2);
p(2,2) = line(2)*sin(dim.angleR(1)/2);

p(3,1) = line(3)*cos(dim.angleR(1)/2);
p(3,2) = line(3)*sin(dim.angleR(1)/2);

p(4,1) = line(4)*cos(dim.angleR(1)/2);
p(4,2) = line(4)*sin(dim.angleR(1)/2);

p(5,1) = line(4)*cos(dim.angleR(1)/2-atan(dim.R_width3/2/line(4)));
p(5,2) = line(4)*sin(dim.angleR(1)/2-atan(dim.R_width3/2/line(4)));

p(6,1) = line(5)*cos(dim.angleR(1)/2-atan((0.8*dim.R_width3/2+0.2*dim.R_width2/2)/line(5)));
p(6,2) = line(5)*sin(dim.angleR(1)/2-atan((0.8*dim.R_width3/2+0.2*dim.R_width2/2)/line(5)));

p(7,1) = line(6)*cos(dim.angleR(1)/2-atan((0.5*dim.R_width3/2+0.5*dim.R_width2/2)/line(6)));
p(7,2) = line(6)*sin(dim.angleR(1)/2-atan((0.5*dim.R_width3/2+0.5*dim.R_width2/2)/line(6)));

p(8,1) = line(7)*cos(dim.angleR(1)/2-atan((0.2*dim.R_width3/2+0.8*dim.R_width2/2)/line(7)));
p(8,2) = line(7)*sin(dim.angleR(1)/2-atan((0.2*dim.R_width3/2+0.8*dim.R_width2/2)/line(7)));

p(9,1) = line(8)*cos(dim.angleR(1)/2-atan(dim.R_width2/2/line(8)));
p(9,2) = line(8)*sin(dim.angleR(1)/2-atan(dim.R_width2/2/line(8)));

p(10,1) = line(8)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(8)));
p(10,2) = line(8)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(8)));

p(11,1) = line(9)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(9)));
p(11,2) = line(9)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(9)));

p(12,1) = line(9)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(9)));
p(12,2) = line(9)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(9)));

p(13,1) = line(8)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(8)));
p(13,2) = line(8)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(8)));

p(14,1) = line(8)*cos(dim.angleR(1)/2+atan(dim.R_width2/2/line(8)));
p(14,2) = line(8)*sin(dim.angleR(1)/2+atan(dim.R_width2/2/line(8)));

p(15,1) = line(7)*cos(dim.angleR(1)/2+atan((0.2*dim.R_width3/2+0.8*dim.R_width2/2)/line(7)));
p(15,2) = line(7)*sin(dim.angleR(1)/2+atan((0.2*dim.R_width3/2+0.8*dim.R_width2/2)/line(7)));

p(16,1) = line(6)*cos(dim.angleR(1)/2+atan((0.5*dim.R_width3/2+0.5*dim.R_width2/2)/line(6)));
p(16,2) = line(6)*sin(dim.angleR(1)/2+atan((0.5*dim.R_width3/2+0.5*dim.R_width2/2)/line(6)));

p(17,1) = line(5)*cos(dim.angleR(1)/2+atan((0.8*dim.R_width3/2+0.2*dim.R_width2/2)/line(5)));
p(17,2) = line(5)*sin(dim.angleR(1)/2+atan((0.8*dim.R_width3/2+0.2*dim.R_width2/2)/line(5)));

p(18,1) = line(4)*cos(dim.angleR(1)/2+atan(dim.R_width3/2/line(4)));
p(18,2) = line(4)*sin(dim.angleR(1)/2+atan(dim.R_width3/2/line(4)));

p(19,1) = line(6)*cos(dim.angleR(1)/2);
p(19,2) = line(6)*sin(dim.angleR(1)/2);

p(20,1) = line(7)*cos(dim.angleR(1)/2);
p(20,2) = line(7)*sin(dim.angleR(1)/2);

p(21,1) = line(8)*cos(dim.angleR(1)/2);
p(21,2) = line(8)*sin(dim.angleR(1)/2);

p(22,1) = line(9)*cos(dim.angleR(1)/2);
p(22,2) = line(9)*sin(dim.angleR(1)/2);

p(23,1) = line(1)*cos(dim.angleR(1));
p(23,2) = line(1)*sin(dim.angleR(1));

p(24,1) = line(3)*cos(dim.angleR(1));
p(24,2) = line(3)*sin(dim.angleR(1));

p(25,1) = line(5)*cos(dim.angleR(1));
p(25,2) = line(5)*sin(dim.angleR(1));

p(26,1) = line(6)*cos(dim.angleR(1));
p(26,2) = line(6)*sin(dim.angleR(1));

p(27,1) = line(8)*cos(dim.angleR(1));
p(27,2) = line(8)*sin(dim.angleR(1));

p(28,1) = line(10)*cos(dim.angleR(1));
p(28,2) = line(10)*sin(dim.angleR(1));

p(29,1) = line(10)*cos(dim.angleR(1)/2+atan(dim.R_width2/2/line(10)));
p(29,2) = line(10)*sin(dim.angleR(1)/2+atan(dim.R_width2/2/line(10)));

p(30,1) = line(10)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(10)));
p(30,2) = line(10)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(10)));

p(31,1) = line(10)*cos(dim.angleR(1)/2);
p(31,2) = line(10)*sin(dim.angleR(1)/2);

p(32,1) = line(10)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(10)));
p(32,2) = line(10)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(10)));

p(33,1) = line(10)*cos(dim.angleR(1)/2-atan(dim.R_width2/2/line(10)));
p(33,2) = line(10)*sin(dim.angleR(1)/2-atan(dim.R_width2/2/line(10)));

p(34,1) = line(10);

p(35,1) = line(8);

p(36,1) = line(6);

p(37,1) = line(5);

p(38,1) = line(3);

p(39,1) = line(1);

t(1,:) = [1 39 23];
m(1) = dim.SFM;

t(2,:) = [39 2 23];
m(2) = dim.RM;

t(3,:) = [39 38 2];
m(3) = dim.RM;

t(4,:) = [38 3 2];
m(4) = dim.RM;

t(5,:) = [2 24 23];
m(5) = dim.RM;

t(6,:) = [2 3 24];
m(6) = dim.RM;

t(7,:) = [3 5 4];
m(7) = dim.RM;

t(8,:) = [38 5 3];
m(8) = dim.RM;

t(9,:) = [38 37 5];
m(9) = dim.RM;

t(10,:) = [37 6 5];
m(10) = dim.RM;

t(11,:) = [37 36 6];
m(11) = dim.RM;

t(12,:) = [36 7 6];
m(12) = dim.RM;

t(13,:) = [36 8 7];
m(13) = dim.RM;

t(14,:) = [36 35 8];
m(14) = dim.RM;

t(15,:) = [35 9 8];
m(15) = dim.RM;

t(16,:) = [35 34 33];
m(16) = dim.RM;

t(17,:) = [35 33 9];
m(17) = dim.RM;

t(18,:) = [33 32 11];
m(18) = dim.RM;

t(19,:) = [33 11 9];
m(19) = dim.RM;

t(20,:) = [11 10 9];
m(20) = dim.RM;

t(21,:) = [14 13 12];
m(21) = dim.RM;

t(22,:) = [14 12 29];
m(22) = dim.RM;

t(23,:) = [12 30 29];
m(23) = dim.RM;

t(24,:) = [29 28 27];
m(24) = dim.RM;

t(25,:) = [29 27 14];
m(25) = dim.RM;

t(26,:) = [27 26 15];
m(26) = dim.RM;

t(27,:) = [27 15 14];
m(27) = dim.RM;

t(28,:) = [26 16 15];
m(28) = dim.RM;

t(29,:) = [26 17 16];
m(29) = dim.RM;

t(30,:) = [26 25 17];
m(30) = dim.RM;

t(31,:) = [25 18 17];
m(31) = dim.RM;

t(32,:) = [25 24 18];
m(32) = dim.RM;

t(33,:) = [24 3 18];
m(33) = dim.RM;

t(34,:) = [3 4 18];
m(34) = dim.RM;

t(35,:) = [18 4 17];
m(35) = 999;

t(36,:) = [4 5 6];
m(36) = 999;
  
t(37,:) = [4 6 17];
m(37) = 999;

t(38,:) = [6 7 19];
m(38) = 999;

t(39,:) = [6 19 17];
m(39) = 999;

t(40,:) = [19 16 17];
m(40) = 999;

t(41,:) = [7 8 20];
m(41) = 999;

t(42,:) = [7 20 19];
m(42) = 999;

t(43,:) = [19 20 16];
m(43) = 999;

t(44,:) = [20 15 16];
m(44) = 999;

t(45,:) = [8 9 10];
m(45) = 999;

t(46,:) = [8 10 20];
m(46) = 999;

t(47,:) = [10 21 20];
m(47) = 999;

t(48,:) = [20 21 13];
m(48) = 999;

t(49,:) = [20 13 15];
m(49) = 999;

t(50,:) = [13 14 15];
m(50) = 999;

t(51,:) = [10 11 22];
m(51) = 9999;

t(52,:) = [11 32 22];
m(52) = 9999;

t(53,:) = [32 31 22];
m(53) = 9999;

t(54,:) = [22 21 10];
m(54) = 9999;

t(55,:) = [31 30 22];
m(55) = 9999;

t(56,:) = [30 12 22];
m(56) = 9999;

t(57,:) = [12 13 22];
m(57) = 9999;

t(58,:) = [13 21 22];
m(58) = 9999;

end