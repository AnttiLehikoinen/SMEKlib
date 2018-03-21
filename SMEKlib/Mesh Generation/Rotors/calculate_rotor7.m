function [p,t,m,FL,LL,ag] = calculate_rotor7(dim)

%Boundries
FL = [1 41 40 39 38 37 36];
LL = [1 25 26 27 28 29 30];
ag = [36 35 34 33 32 31 30];

%prealloacte memory
p = zeros(41,2);
t = zeros(62,3);
m = zeros(62,1);


%help variables
line = [dim.Rin 0.5*dim.Rin+0.5*(dim.Rout-dim.R_height) 0.3*dim.Rin+0.7*(dim.Rout-dim.R_height) 0.1*dim.Rin+0.9*(dim.Rout-dim.R_height) ...
        dim.Rout-dim.R_height dim.Rout-dim.R_height+dim.R_width3/2 0.7*(dim.Rout-dim.R_height+dim.R_width3/2)+0.3*(dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height2) ...
        dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height2 0.3*(dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height2)+0.7*(dim.Rout-dim.R_height1-dim.R_width1/2) ...
        dim.Rout-dim.R_height1-dim.R_width1/2 dim.Rout-dim.R_height1 dim.Rout];
c1 = sqrt(line(6)^2+(dim.R_width3/2)^2-2*line(6)*(dim.R_width3/2)*cosd(45));
alfa1 = asin((dim.R_width3/2*sind(45))/c1);
c2 = sqrt(line(10)^2+(dim.R_width1/2)^2-2*line(10)*(dim.R_width1/2)*cosd(135));
alfa2 = asin((dim.R_width1/2*sind(135))/c2);

p(2,1) = line(2)*cos(dim.angleR(1)/2);
p(2,2) = line(2)*sin(dim.angleR(1)/2);

p(3,1) = line(4)*cos(dim.angleR(1)/2);
p(3,2) = line(4)*sin(dim.angleR(1)/2);

p(4,1) = line(5)*cos(dim.angleR(1)/2);
p(4,2) = line(5)*sin(dim.angleR(1)/2);

p(5,1) = c1*cos(dim.angleR(1)/2-alfa1);
p(5,2) = c1*sin(dim.angleR(1)/2-alfa1);

p(6,1) = line(6)*cos(dim.angleR(1)/2-atan(dim.R_width3/2/line(6)));
p(6,2) = line(6)*sin(dim.angleR(1)/2-atan(dim.R_width3/2/line(6)));

p(7,1) = line(7)*cos(dim.angleR(1)/2-atan((0.7*dim.R_width3/2+0.3*dim.R_width2/2)/line(7)));
p(7,2) = line(7)*sin(dim.angleR(1)/2-atan((0.7*dim.R_width3/2+0.3*dim.R_width2/2)/line(7)));

p(8,1) = line(8)*cos(dim.angleR(1)/2-atan(dim.R_width2/2/line(8)));
p(8,2) = line(8)*sin(dim.angleR(1)/2-atan(dim.R_width2/2/line(8)));

p(9,1) = line(9)*cos(dim.angleR(1)/2-atan((0.3*dim.R_width2/2+0.7*dim.R_width1/2)/line(9)));
p(9,2) = line(9)*sin(dim.angleR(1)/2-atan((0.3*dim.R_width2/2+0.7*dim.R_width1/2)/line(9)));

p(10,1) = line(10)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(10)));
p(10,2) = line(10)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(10)));

p(11,1) = c2*cos(dim.angleR(1)/2-alfa2);
p(11,2) = c2*sin(dim.angleR(1)/2-alfa2);

p(12,1) = line(11)*cos(dim.angleR(1)/2);
p(12,2) = line(11)*sin(dim.angleR(1)/2);

p(13,1) = c2*cos(dim.angleR(1)/2+alfa2);
p(13,2) = c2*sin(dim.angleR(1)/2+alfa2);

p(14,1) = line(10)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(10)));
p(14,2) = line(10)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(10)));

p(15,1) = line(9)*cos(dim.angleR(1)/2+atan((0.3*dim.R_width2/2+0.7*dim.R_width1/2)/line(9)));
p(15,2) = line(9)*sin(dim.angleR(1)/2+atan((0.3*dim.R_width2/2+0.7*dim.R_width1/2)/line(9)));

p(16,1) = line(8)*cos(dim.angleR(1)/2+atan(dim.R_width2/2/line(8)));
p(16,2) = line(8)*sin(dim.angleR(1)/2+atan(dim.R_width2/2/line(8)));

p(17,1) = line(7)*cos(dim.angleR(1)/2+atan((0.7*dim.R_width3/2+0.3*dim.R_width2/2)/line(7)));
p(17,2) = line(7)*sin(dim.angleR(1)/2+atan((0.7*dim.R_width3/2+0.3*dim.R_width2/2)/line(7)));

p(18,1) = line(6)*cos(dim.angleR(1)/2+atan(dim.R_width3/2/line(6)));
p(18,2) = line(6)*sin(dim.angleR(1)/2+atan(dim.R_width3/2/line(6)));

p(19,1) = c1*cos(dim.angleR(1)/2+alfa1);
p(19,2) = c1*sin(dim.angleR(1)/2+alfa1);

p(20,1) = line(6)*cos(dim.angleR(1)/2);
p(20,2) = line(6)*sin(dim.angleR(1)/2);

p(21,1) = line(8)*cos(dim.angleR(1)/2);
p(21,2) = line(8)*sin(dim.angleR(1)/2);

p(22,1) = line(10)*cos(dim.angleR(1)/2);
p(22,2) = line(10)*sin(dim.angleR(1)/2);

p(23,1) = line(10)*cos(3*dim.angleR(1)/4);
p(23,2) = line(10)*sin(3*dim.angleR(1)/4);

p(24,1) = line(10)*cos(dim.angleR(1)/4);
p(24,2) = line(10)*sin(dim.angleR(1)/4);

p(25,1) = line(1)*cos(dim.angleR(1));
p(25,2) = line(1)*sin(dim.angleR(1));

p(26,1) = line(3)*cos(dim.angleR(1));
p(26,2) = line(3)*sin(dim.angleR(1));

p(27,1) = line(6)*cos(dim.angleR(1));
p(27,2) = line(6)*sin(dim.angleR(1));

p(28,1) = line(8)*cos(dim.angleR(1));
p(28,2) = line(8)*sin(dim.angleR(1));

p(29,1) = line(10)*cos(dim.angleR(1));
p(29,2) = line(10)*sin(dim.angleR(1));

p(30,1) = line(12)*cos(dim.angleR(1));
p(30,2) = line(12)*sin(dim.angleR(1));

p(31,1) = line(12)*cos(3*dim.angleR(1)/4);
p(31,2) = line(12)*sin(3*dim.angleR(1)/4);

p(32,1) = line(12)*cos(dim.angleR(1)/2+alfa2);
p(32,2) = line(12)*sin(dim.angleR(1)/2+alfa2);

p(33,1) = line(12)*cos(dim.angleR(1)/2);
p(33,2) = line(12)*sin(dim.angleR(1)/2);

p(34,1) = line(12)*cos(dim.angleR(1)/2-alfa2);
p(34,2) = line(12)*sin(dim.angleR(1)/2-alfa2);

p(35,1) = line(12)*cos(dim.angleR(1)/4);
p(35,2) = line(12)*sin(dim.angleR(1)/4);

p(36,1) = line(12);

p(37,1) = line(10);

p(38,1) = line(8);

p(39,1) = line(6);

p(40,1) = line(3);

p(41,1) = line(1);


t(1,:) = [1 41 25];
m(1)= dim.SFM;

t(2,:) = [41 2 25];
m(2) = dim.RM;

t(3,:) = [41 40 2];
m(3) = dim.RM;

t(4,:) = [2 26 25];
m(4) = dim.RM;

t(5,:) = [2 40 26];
m(5) = dim.RM;

t(6,:) = [40 3 26];
m(6) = dim.RM;

t(7,:) = [40 5 3];
m(7) = dim.RM;

t(8,:) = [3 5 4];
m(8) = dim.RM;

t(9,:) = [40 39 5];
m(9) = dim.RM;

t(10,:) = [39 6 5];
m(10) = dim.RM;

t(11,:) = [39 7 6];
m(11) = dim.RM;

t(12,:) = [39 38 7];
m(12) = dim.RM;

t(13,:) = [38 8 7];
m(13) = dim.RM;

t(14,:) = [38 37 24];
m(14) = dim.RM;

t(15,:) = [38 24 8];
m(15) = dim.RM;

t(16,:) = [8 24 9];
m(16) = dim.RM;

t(17,:) = [24 10 9];
m(17) = dim.RM;

t(18,:) = [37 36 35];
m(18) = dim.RM;

t(19,:) = [37 35 24];
m(19) = dim.RM;

t(20,:) = [24 35 11];
m(20) = dim.RM;

t(21,:) = [24 11 10];
m(21) = dim.RM;

t(22,:) = [11 35 34];
m(22) = dim.RM;

t(23,:) = [12 11 34];
m(23) = dim.RM;

t(24,:) = [12 34 33];
m(24) = dim.RM;

t(25,:) = [12 33 32];
m(25) = dim.RM;

t(26,:) = [13 12 32];
m(26) = dim.RM;

t(27,:) = [13 32 31];
m(27) = dim.RM;

t(28,:) = [13 31 23];
m(28) = dim.RM;

t(29,:) = [14 13 23];
m(29) = dim.RM;

t(30,:) = [23 31 29];
m(30) = dim.RM;

t(31,:) = [31 30 29];
m(31) = dim.RM;

t(32,:) = [23 29 28];
m(32) = dim.RM;

t(33,:) = [23 28 16];
m(33) = dim.RM;

t(34,:) = [23 16 15];
m(34) = dim.RM;

t(35,:) = [23 15 14];
m(35) = dim.RM;

t(36,:) = [28 17 16];
m(36) = dim.RM;

t(37,:) = [28 27 17];
m(37) = dim.RM;

t(38,:) = [27 18 17];
m(38) = dim.RM;

t(39,:) = [27 19 18];
m(39) = dim.RM;

t(40,:) = [27 26 19];
m(40) = dim.RM;

t(41,:) = [26 3 19];
m(41) = dim.RM;

t(42,:) = [3 4 19];
m(42) = dim.RM;

t(43,:) = [4 5 20];
m(43) = 9999;

t(44,:) = [5 6 20];
m(44) = 9999;

t(45,:) = [4 20 19];
m(45) = 9999;

t(46,:) = [18 19 20];
m(46) = 9999;

t(47,:) = [6 7 20];
m(47) = 9999;

t(48,:) = [7 17 20];%
m(48) = 9999;

t(49,:) = [17 18 20];
m(49) = 9999;

t(50,:) = [7 8 21];
m(50) = 9999;

t(51,:) = [7 21 17];
m(51) = 9999;

t(52,:) = [21 16 17];
m(52) = 9999;

t(53,:) = [8 9 21];
m(53) = 9999;

t(54,:) = [9 15 21];
m(54) = 9999;

t(55,:) = [21 15 16];
m(55) = 9999;

t(56,:) = [9 22 15];
m(56) = 9999;

t(57,:) = [9 10 22];
m(57) = 9999;

t(58,:) = [22 14 15];
m(58) = 9999;

t(59,:) = [10 11 22];
m(59) = 9999;

t(60,:) = [11 12 22];
m(60) = 9999;

t(61,:) =[12 13 22];
m(61) = 9999;

t(62,:) = [13 14 22];
m(62) = 9999;

end