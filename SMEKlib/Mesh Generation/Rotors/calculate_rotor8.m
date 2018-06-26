function [p,t,m,FL,LL,ag] = calculate_rotor8(dim)

%Boundaries
FL = [1 53 52 51 50 49 48 47];
LL = [1 25 26 27 28 29 30 31];
ag = [47 46 41 43 35 32 31];

%prealloacte memory
p = zeros(53,2);
t = zeros(84,3);
m = zeros(84,1);

%help variables
line = [dim.Rin 0.5*dim.Rin+0.5*(dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height4-dim.R_height5) 0.3*dim.Rin+0.7*(dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height4-dim.R_height5) 0.1*dim.Rin+0.9*(dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height4-dim.R_height5) ...
        dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height4-dim.R_height5 dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height4-dim.R_height5+dim.R_width4/2 0.7*(dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height4-dim.R_height5+dim.R_width4/2)+0.3*(dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height4-dim.R_width2/2-dim.R_height2) ...
        dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height4-dim.R_width2/2-dim.R_height2 0.3*(dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height4-dim.R_width2/2-dim.R_height2)+0.7*(dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height4-dim.R_width2/2) ...
        dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height4-dim.R_width2/2 dim.Rout-dim.R_height1-dim.R_width1/2-dim.R_height4 dim.Rout-dim.R_height1-dim.R_width1/2 dim.Rout-dim.R_height1 dim.Rout];

c1 = sqrt(line(6)^2+(dim.R_width4/2)^2-2*line(6)*(dim.R_width4/2)*cosd(45));
alfa1 = asin((dim.R_width4/2*sind(45))/c1);
c2 = sqrt(line(10)^2+(dim.R_width2/2)^2-2*line(10)*(dim.R_width2/2)*cosd(135));
alfa2 = asin((dim.R_width2/2*sind(135))/c2);
c3 = sqrt(line(13)^2+(dim.R_width1/2)^2-2*line(6)*(dim.R_width1/2)*cosd(45));
alfa3 = asin((dim.R_width1/2*sind(45))/c2);

p(2,1) = line(2)*cos(dim.angleR(1)/2);
p(2,2) = line(2)*sin(dim.angleR(1)/2);

p(3,1) = line(4)*cos(dim.angleR(1)/2);
p(3,2) = line(4)*sin(dim.angleR(1)/2);

p(4,1) = line(5)*cos(dim.angleR(1)/2);
p(4,2) = line(5)*sin(dim.angleR(1)/2);

p(5,1) = c1*cos(dim.angleR(1)/2-alfa1);
p(5,2) = c1*sin(dim.angleR(1)/2-alfa1);

p(6,1) = line(6)*cos(dim.angleR(1)/2-atan(dim.R_width4/2/line(6)));
p(6,2) = line(6)*sin(dim.angleR(1)/2-atan(dim.R_width4/2/line(6)));

p(7,1) = line(7)*cos(dim.angleR(1)/2-atan((0.7*dim.R_width4/2+0.3*dim.R_width3/2)/line(7)));
p(7,2) = line(7)*sin(dim.angleR(1)/2-atan((0.7*dim.R_width4/2+0.3*dim.R_width3/2)/line(7)));

p(8,1) = line(8)*cos(dim.angleR(1)/2-atan(dim.R_width3/2/line(8)));
p(8,2) = line(8)*sin(dim.angleR(1)/2-atan(dim.R_width3/2/line(8)));

p(9,1) = line(9)*cos(dim.angleR(1)/2-atan((0.3*dim.R_width3/2+0.7*dim.R_width2/2)/line(9)));
p(9,2) = line(9)*sin(dim.angleR(1)/2-atan((0.3*dim.R_width3/2+0.7*dim.R_width2/2)/line(9)));

p(10,1) = line(10)*cos(dim.angleR(1)/2-atan(dim.R_width2/2/line(10)));
p(10,2) = line(10)*sin(dim.angleR(1)/2-atan(dim.R_width2/2/line(10)));

p(11,1) = c2*cos(dim.angleR(1)/2-alfa2);
p(11,2) = c2*sin(dim.angleR(1)/2-alfa2);

p(12,1) = line(11)*cos(dim.angleR(1)/2);
p(12,2) = line(11)*sin(dim.angleR(1)/2);

p(13,1) = c2*cos(dim.angleR(1)/2+alfa2);
p(13,2) = c2*sin(dim.angleR(1)/2+alfa2);

p(14,1) = line(10)*cos(dim.angleR(1)/2+atan(dim.R_width2/2/line(10)));
p(14,2) = line(10)*sin(dim.angleR(1)/2+atan(dim.R_width2/2/line(10)));

p(15,1) = line(9)*cos(dim.angleR(1)/2+atan((0.3*dim.R_width3/2+0.7*dim.R_width2/2)/line(9)));
p(15,2) = line(9)*sin(dim.angleR(1)/2+atan((0.3*dim.R_width3/2+0.7*dim.R_width2/2)/line(9)));

p(16,1) = line(8)*cos(dim.angleR(1)/2+atan(dim.R_width3/2/line(8)));
p(16,2) = line(8)*sin(dim.angleR(1)/2+atan(dim.R_width3/2/line(8)));

p(17,1) = line(7)*cos(dim.angleR(1)/2+atan((0.7*dim.R_width4/2+0.3*dim.R_width3/2)/line(7)));
p(17,2) = line(7)*sin(dim.angleR(1)/2+atan((0.7*dim.R_width4/2+0.3*dim.R_width3/2)/line(7)));

p(18,1) = line(6)*cos(dim.angleR(1)/2+atan(dim.R_width4/2/line(6)));
p(18,2) = line(6)*sin(dim.angleR(1)/2+atan(dim.R_width4/2/line(6)));

p(19,1) = c1*cos(dim.angleR(1)/2+alfa1);
p(19,2) = c1*sin(dim.angleR(1)/2+alfa1);

p(20,1) = line(6)*cos(dim.angleR(1)/2);
p(20,2) = line(6)*sin(dim.angleR(1)/2);

p(21,1) = line(8)*cos(dim.angleR(1)/2);
p(21,2) = line(8)*sin(dim.angleR(1)/2);

p(22,1) = line(10)*cos(dim.angleR(1)/2);
p(22,2) = line(10)*sin(dim.angleR(1)/2);

p(23,1) = line(10)*cos(3*dim.angleR(1)/4);%3/4
p(23,2) = line(10)*sin(3*dim.angleR(1)/4);

p(24,1) = line(10)*cos(dim.angleR(1)/4);%1/4
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

p(30,1) = line(13)*cos(dim.angleR(1));
p(30,2) = line(13)*sin(dim.angleR(1));

p(31,1) = line(14)*cos(dim.angleR(1));
p(31,2) = line(14)*sin(dim.angleR(1));

p(32,1) = line(14)*cos(3/4*dim.angleR(1));
p(32,2) = line(14)*sin(3/4*dim.angleR(1));

p(33,1) = line(13)*cos(3/4*dim.angleR(1));
p(33,2) = line(13)*sin(3/4*dim.angleR(1));

p(34,1) = (line(11)/2+line(12)/2)*cos(dim.angleR(1)/2+atan((dim.R_width2/4+dim.R_width1/4)/line(11)));
p(34,2) = (line(11)/2+line(12)/2)*sin(dim.angleR(1)/2+atan((dim.R_width2/4+dim.R_width1/4)/line(11)));

p(35,1) = line(14)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(14)));
p(35,2) = line(14)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(14)));

p(36,1) = line(13)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(13)));
p(36,2) = line(13)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(13)));

p(37,1) = c3*cos(dim.angleR(1)/2+alfa3);
p(37,2) = c3*sin(dim.angleR(1)/2+alfa3);

p(38,1) = line(12)*cos(dim.angleR(1)/2);
p(38,2) = line(12)*sin(dim.angleR(1)/2);

p(39,1) = c3*cos(dim.angleR(1)/2-alfa3);
p(39,2) = c3*sin(dim.angleR(1)/2-alfa3);

p(40,1) = line(13)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(13)));
p(40,2) = line(13)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(13)));

p(41,1) = line(14)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(14)));
p(41,2) = line(14)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(14)));

p(42,1) = line(13)*cos(dim.angleR(1)/2);
p(42,2) = line(13)*sin(dim.angleR(1)/2);

p(43,1) = line(14)*cos(dim.angleR(1)/2);
p(43,2) = line(14)*sin(dim.angleR(1)/2);

p(44,1) = (line(11)/2+line(12)/2)*cos(dim.angleR(1)/2-atan((dim.R_width2/4+dim.R_width1/4)/line(11)));
p(44,2) = (line(11)/2+line(12)/2)*sin(dim.angleR(1)/2-atan((dim.R_width2/4+dim.R_width1/4)/line(11)));

p(45,1) = line(13)*cos(1/4*dim.angleR(1));
p(45,2) = line(13)*sin(1/4*dim.angleR(1));

p(46,1) = line(14)*cos(1/4*dim.angleR(1));
p(46,2) = line(14)*sin(1/4*dim.angleR(1));

p(47,1) = line(14);

p(48,1) = line(13);

p(49,1) = line(10);

p(50,1) = line(8);

p(51,1) = line(6);

p(52,1) = line(3);

p(53,1) = line(1);


t(1,:) = [1 53 25];
m(1) = dim.SFM;

t(2,:) = [53 2 25];
m(2) = dim.RM;

t(3,:) = [53 52 2];
m(3) = dim.RM;

t(4,:) = [2 26 25];
m(4) = dim.RM;

t(5,:) = [2 52 26];
m(5) = dim.RM;

t(6,:) = [52 3 26];
m(6) = dim.RM;

t(7,:) = [52 5 3];
m(7) = dim.RM;

t(8,:) = [5 4 3];
m(8) = dim.RM;

t(9,:) = [52 51 5];
m(9) = dim.RM;

t(10,:) = [51 6 5];
m(10) = dim.RM;

t(11,:) = [51 7 6];
m(11) = dim.RM;

t(12,:) = [51 50 7];
m(12) = dim.RM;

t(13,:) = [50 8 7];
m(13) = dim.RM;

t(14,:) = [50 49 24];
m(14) = dim.RM;

t(15,:) = [50 24 8];
m(15) = dim.RM;

t(16,:) = [24 9 8];
m(16) = dim.RM;

t(17,:) = [24 10 9];
m(17) = dim.RM;

t(18,:) = [49 48 45];
m(18) = dim.RM;

t(19,:) = [48 47 45];
m(19) = dim.RM;

t(20,:) = [47 46 45];
m(20) = dim.RM;

t(21,:) = [49 45 24];
m(21) = dim.RM;

t(22,:) = [46 41 45];
m(22) = dim.RM;

t(23,:) = [45 41 40];
m(23) = dim.RM;

t(24,:) = [45 40 44];
m(24) = dim.RM;

t(25,:) = [45 44 24];
m(25) = dim.RM;

t(26,:) = [44 40 39];
m(26) = dim.RM;

t(27,:) = [44 39 11];
m(27) = dim.RM;

t(28,:) = [24 44 11];
m(28) = dim.RM;

t(29,:) = [24 11 10];
m(29) = dim.RM;

t(30,:) = [39 12 11];
m(30) = dim.RM;

t(31,:) = [39 38 12];
m(31) = dim.RM;

t(32,:) = [38 37 12];
m(32) = dim.RM;

t(33,:) = [37 13 12];
m(33) = dim.RM;

t(34,:) = [23 14 13];
m(34) = dim.RM;

t(35,:) = [23 13 34];
m(35) = dim.RM;

t(36,:) = [13 37 34];
m(36) = dim.RM;

t(37,:) = [34 37 36];
m(37) = dim.RM;

t(38,:) = [34 36 33];
m(38) = dim.RM;

t(39,:) = [33 36 35];
m(39) = dim.RM;

t(40,:) = [33 35 32];
m(40) = dim.RM;

t(41,:) = [23 34 33];
m(41) = dim.RM;

t(42,:) = [33 32 31];
m(42) = dim.RM;

t(43,:) = [31 30 33];
m(43) = dim.RM;

t(44,:) = [33 30 29];
m(44) = dim.RM;

t(45,:) = [33 29 23];
m(45) = dim.RM;

t(46,:) = [15 14 23];
m(46) = dim.RM;

t(47,:) = [16 15 23];
m(47) = dim.RM;

t(48,:) = [23 29 28];
m(48) = dim.RM;

t(49,:) = [23 28 16];
m(49) = dim.RM;

t(50,:) = [28 17 16];
m(50) = dim.RM;

t(51,:) = [28 27 17];
m(51) = dim.RM;

t(52,:) = [27 18 17];
m(52) = dim.RM;

t(53,:) = [27 19 18];
m(53) = dim.RM;

t(54,:) = [27 26 19];
m(54) = dim.RM;

t(55,:) = [26 3 19];
m(55) = dim.RM;

t(56,:) = [3 4 19];
m(56) = dim.RM;

t(57,:) = [4 5 20];
m(57) = 9999;

t(58,:) = [5 6 20];
m(58) = 9999;

t(59,:) = [4 20 19];
m(59) = 9999;

t(60,:) = [20 18 19];
m(60) = 9999;

t(61,:) = [6 7 20];
m(61) = 9999;

t(62,:) = [7 17 20];
m(62) = 9999;

t(63,:) = [17 18 20];
m(63) = 9999;

t(64,:) = [7 8 21];
m(64) = 9999;

t(65,:) = [7 21 17];
m(65) = 9999;

t(66,:) = [21 16 17];
m(66) = 9999;

t(67,:) = [8 9 21];
m(67) = 9999;

t(68,:) = [9 15 21];
m(68) = 9999;

t(69,:) = [21 15 16];
m(69) = 9999;

t(70,:) = [9 10 22];
m(70) = 9999;

t(71,:) = [9 22 15];
m(71) = 9999;

t(72,:) = [22 14 15];
m(72) = 9999;

t(73,:) = [10 11 22];
m(73) = 9999;

t(74,:) = [11 12 22];
m(74) = 9999;

t(75,:) = [12 13 22];
m(75) = 9999;

t(76,:) = [13 14 22];
m(76) = 9999;

t(77,:) = [38 39 42];
m(77) = 999;

t(78,:) = [39 40 42];
m(78) = 999;

t(79,:) = [40 41 42];
m(79) = 999;

t(80,:) = [41 43 42];
m(80) = 999;

t(81,:) = [43 35 42];
m(81) = 999;

t(82,:) = [35 36 42];
m(82) = 999;

t(83,:) = [36 37 42];
m(83) = 999;

t(84,:) = [37 38 42];
m(84) = 999;

end