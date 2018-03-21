function [p,t,m,FL,LL,ag] = calculate_rotor22(dim)

%Boundaries
FL = [1 53 52 51 50 49 48 47 46];
LL = [1 34 35 36 37 38 39 40 41];
ag = [46 45 44 43 42 41];

%prealloacte memory
p = zeros(53,2);
t = zeros(83,3);
m = zeros(83,1);

line = [dim.Rin 0.5*dim.Rin+0.5*(dim.Rout-dim.R_height) 0.25*dim.Rin+0.75*(dim.Rout-dim.R_height) ...
        0.1*dim.Rin+0.9*(dim.Rout-dim.R_height) 0.05*dim.Rin+0.95*(dim.Rout-dim.R_height) dim.Rout-dim.R_height ...
        dim.Rout-dim.R_height+dim.R_width2/2 dim.Rout-dim.R_height+dim.R_width2 0.5*(dim.Rout-dim.R_height1)+0.5*(dim.Rout-dim.R_height+dim.R_width2) ...
        dim.Rout-dim.R_height1 dim.Rout];


c1 = sqrt(line(7)^2+(dim.R_width2/2)^2-2*line(7)*(dim.R_width2/2)*cosd(30));
alfa1 = asin((dim.R_width2/2*sind(30))/c1);
c2 = sqrt(line(7)^2+(dim.R_width2/2)^2-2*line(7)*(dim.R_width2/2)*cosd(60));
alfa2 = asin((dim.R_width2/2*sind(60))/c2);
c3 = sqrt(line(7)^2+(dim.R_width2/2)^2-2*line(7)*(dim.R_width2/2)*cosd(120));
alfa3 = asin((dim.R_width2/2*sind(120))/c3);

c4 = sqrt(line(7)^2+(dim.R_width3/2)^2-2*line(7)*(dim.R_width3/2)*cosd(30));
alfa4 = asin((dim.R_width3/2*sind(30))/c3);
c5 = sqrt(line(7)^2+(dim.R_width3/2)^2-2*line(7)*(dim.R_width3/2)*cosd(150));
alfa5 = asin((dim.R_width3/2*sind(150))/c3);


p(2,1) = line(2)*cos(dim.angleR(1)/2);
p(2,2) = line(2)*sin(dim.angleR(1)/2);

p(3,1) = line(4)*cos(dim.angleR(1)/2);
p(3,2) = line(4)*sin(dim.angleR(1)/2);
    
p(4,1) = line(6)*cos(dim.angleR(1)/2);
p(4,2) = line(6)*sin(dim.angleR(1)/2);
    
p(5,1) = c1*cos(dim.angleR(1)/2-alfa1);
p(5,2) = c1*sin(dim.angleR(1)/2-alfa1); 

p(6,1) = c2*cos(dim.angleR(1)/2-alfa2);
p(6,2) = c2*sin(dim.angleR(1)/2-alfa2); 

p(7,1) = line(7)*cos(dim.angleR(1)/2-atan(dim.R_width2/2/line(7)));
p(7,2) = line(7)*sin(dim.angleR(1)/2-atan(dim.R_width2/2/line(7)));

p(8,1) = c3*cos(dim.angleR(1)/2-alfa3);
p(8,2) = c3*sin(dim.angleR(1)/2-alfa3); 

p(9,1) = line(8)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(8)));
p(9,2) = line(8)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(8)));

p(10,1) = line(8)*cos(dim.angleR(1)/2);
p(10,2) = line(8)*sin(dim.angleR(1)/2);
    
p(11,1) = line(8)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(8)));
p(11,2) = line(8)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(8)));

p(12,1) = c3*cos(dim.angleR(1)/2+alfa3);
p(12,2) = c3*sin(dim.angleR(1)/2+alfa3); 

p(13,1) = line(7)*cos(dim.angleR(1)/2+atan(dim.R_width2/2/line(7)));
p(13,2) = line(7)*sin(dim.angleR(1)/2+atan(dim.R_width2/2/line(7)));

p(14,1) = c2*cos(dim.angleR(1)/2+alfa2);
p(14,2) = c2*sin(dim.angleR(1)/2+alfa2); 

p(15,1) = c1*cos(dim.angleR(1)/2+alfa1);
p(15,2) = c1*sin(dim.angleR(1)/2+alfa1); 

p(16,1) = c4*cos(dim.angleR(1)/2+alfa4);
p(16,2) = c4*sin(dim.angleR(1)/2+alfa4); 

p(17,1) = c4*cos(dim.angleR(1)/2-alfa4);
p(17,2) = c4*sin(dim.angleR(1)/2-alfa4); 

p(18,1) = line(7)*cos(dim.angleR(1)/2-atan(dim.R_width3/2/line(7)));
p(18,2) = line(7)*sin(dim.angleR(1)/2-atan(dim.R_width3/2/line(7)));

p(19,1) = c5*cos(dim.angleR(1)/2-alfa5);
p(19,2) = c5*sin(dim.angleR(1)/2-alfa5); 

p(20,1) = c5*cos(dim.angleR(1)/2+alfa5);
p(20,2) = c5*sin(dim.angleR(1)/2+alfa5); 

p(21,1) = line(7)*cos(dim.angleR(1)/2+atan(dim.R_width3/2/line(7)));
p(21,2) = line(7)*sin(dim.angleR(1)/2+atan(dim.R_width3/2/line(7)));

p(22,1) = line(7)*cos(dim.angleR(1)/2);
p(22,2) = line(7)*sin(dim.angleR(1)/2);

p(23,1) = line(5)*cos(3/4*dim.angleR(1));
p(23,2) = line(5)*sin(3/4*dim.angleR(1));

p(24,1) = line(5)*cos(1/4*dim.angleR(1));
p(24,2) = line(5)*sin(1/4*dim.angleR(1));

p(25,1) = line(9)*cos(3/4*dim.angleR(1));
p(25,2) = line(9)*sin(3/4*dim.angleR(1));

p(26,1) = line(9)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(9)));
p(26,2) = line(9)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(9)));

p(27,1) = line(9)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(9)));
p(27,2) = line(9)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(9)));

p(28,1) = line(9)*cos(1/4*dim.angleR(1));
p(28,2) = line(9)*sin(1/4*dim.angleR(1));

p(29,1) = line(10)*cos(3/4*dim.angleR(1));
p(29,2) = line(10)*sin(3/4*dim.angleR(1));

p(30,1) = line(10)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(10)));
p(30,2) = line(10)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(10)));

p(31,1) = line(10)*cos(dim.angleR(1)/2);
p(31,2) = line(10)*sin(dim.angleR(1)/2);

p(32,1) = line(10)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(10)));
p(32,2) = line(10)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(10)));

p(33,1) = line(10)*cos(1/4*dim.angleR(1));
p(33,2) = line(10)*sin(1/4*dim.angleR(1));

p(34,1) = line(1)*cos(dim.angleR(1));
p(34,2) = line(1)*sin(dim.angleR(1));

p(35,1) = line(3)*cos(dim.angleR(1));
p(35,2) = line(3)*sin(dim.angleR(1));

p(36,1) = line(6)*cos(dim.angleR(1));
p(36,2) = line(6)*sin(dim.angleR(1));

p(37,1) = line(7)*cos(dim.angleR(1));
p(37,2) = line(7)*sin(dim.angleR(1));

p(38,1) = line(8)*cos(dim.angleR(1));
p(38,2) = line(8)*sin(dim.angleR(1));

p(39,1) = line(9)*cos(dim.angleR(1));
p(39,2) = line(9)*sin(dim.angleR(1));

p(40,1) = line(10)*cos(dim.angleR(1));
p(40,2) = line(10)*sin(dim.angleR(1));

p(41,1) = line(11)*cos(dim.angleR(1));
p(41,2) = line(11)*sin(dim.angleR(1));

p(42,1) = line(11)*cos(3/4*dim.angleR(1));
p(42,2) = line(11)*sin(3/4*dim.angleR(1));

p(43,1) = line(11)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(11)));
p(43,2) = line(11)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(11)));

p(44,1) = line(11)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(11)));
p(44,2) = line(11)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(11)));

p(45,1) = line(11)*cos(1/4*dim.angleR(1));
p(45,2) = line(11)*sin(1/4*dim.angleR(1));

p(46,1) = line(11);

p(47,1) = line(10);

p(48,1) = line(9);

p(49,1) = line(8);

p(50,1) = line(7);

p(51,1) = line(6);

p(52,1) = line(3);

p(53,1) = line(1);


t(1,:) = [1 53 34];
m(1) = 1;

t(2,:) = [53 2 34];
m(2) = 4;

t(3,:) = [53 52 2];
m(3) = 4;

t(4,:) = [34 2 35];
m(4) = 4;

t(5,:) =[2 52 35];
m(5) = 4;

t(6,:) = [52 3 35];
m(6) = 4;

t(7,:) = [52 24 3];
m(7) = 4;

t(8,:) =[52 51 24];
m(8) = 4;

t(9,:) = [24 4 3];
m(9) = 4;

t(10,:) = [24 5 4];
m(10) = 4;

t(11,:) = [51 5 24];
m(11) = 4;

t(12,:) =[51 6 5];
m(12) = 4;

t(13,:) = [51 50 6];
m(13) = 4;

t(14,:) = [50 7 6];
m(14) = 4;

t(15,:) = [50 49 7];
m(15) = 4;

t(16,:) = [49 8 7];
m(16) = 4;

t(17,:) = [49 28 8];
m(17) = 4;

t(18,:) = [49 48 28];
m(18) = 4;

t(19,:) = [28 9 8];
m(19) = 4;

t(20,:) = [28 27 9];
m(20) = 4;

t(21,:) = [48 33 28];
m(21) = 4;

t(22,:) = [48 47 33];
m(22) = 4;

t(23,:) = [47 46 33];
m(23) = 4;

t(24,:) = [46 45 33];
m(24) = 4;

t(25,:) = [33 45 44];
m(25) = 4;

t(26,:) = [44 32 33];
m(26) = 4;

t(27,:) = [32 33 27];
m(27) = 4;

t(28,:) = [33 27 28];
m(28) = 4;

t(29,:) = [43 42 29];
m(29) = 4;

t(30,:) = [43 29 30];
m(30) = 4;

t(31,:) = [26 30 29];
m(31) = 4;

t(32,:) = [29 25 26];
m(32) = 4;

t(33,:) = [42 41 29];
m(33) = 4;

t(34,:) = [41 40 29];
m(34) = 4;

t(35,:) =[40 39 29];
m(35) = 4;

t(36,:) = [39 25 29];
m(36) = 4;

t(37,:) = [25 11 26];
m(37) = 4;

t(38,:) = [25 12 11];
m(38) = 4;

t(39,:) = [25 39 38];
m(39) = 4;

t(40,:) = [25 38 12];
m(40) = 4;

t(41,:) = [38 13 12];
m(41) = 4;

t(42,:) = [38 37 13];
m(42) = 4;

t(43,:) = [37 14 13];
m(43) = 4;

t(44,:) = [37 36 14];
m(44) = 4;

t(45,:) = [36 15 14];
m(45) = 4;

t(46,:) =[36 23 15];
m(46) = 4;

t(47,:) =[23 4 15];
m(47) = 4;

t(48,:) = [23 3 4];
m(48) = 4;

t(49,:) = [36 35 23];
m(49) = 4;

t(50,:) = [16 17 22];
m(50) = 999;% inner circle 50-55

t(51,:) = [17 18 22];
m(51) = 999;%

t(52,:) = [18 19 22];
m(52) = 999;%

t(53,:) = [19 20 22];
m(53) = 999;%

t(54,:) = [20 21 22];
m(54) = 999;%

t(55,:) = [16 22 21];
m(55) = 999;%

%outer ring 56-73
t(56,:) = [4 17 16];
m(56) = 9999;

t(57,:) = [4 5 17];
m(57) = 9999;

t(58,:) = [5 6 17];
m(58) = 9999;

t(59,:) = [6 18 17];
m(59) = 9999;

t(60,:) = [6 7 18];
m(60) = 9999;

t(61,:) = [7 8 18];
m(61) = 9999;

t(62,:) = [8 19 18];
m(62) = 9999;

t(63,:) = [8 9 19];
m(63) = 9999;

t(64,:) = [9 10 19];
m(64) = 9999;

t(65,:) = [10 20 19];
m(65) = 9999;

t(66,:) = [10 11 20];
m(66) = 9999;

t(67,:) = [11 12 20];
m(67) = 9999;

t(68,:) = [12 21 20];
m(68) = 9999;

t(69,:) =[12 13 21];
m(69) = 9999;

t(70,:) =[13 14 21];
m(70) = 9999;

t(71,:) = [14 16 21];
m(71) = 9999;

t(72,:) = [14 15 16];
m(72) = 9999;

t(73,:) = [4 16 15];
m(73) = 9999;

%gap 74-end
t(74,:) = [26 11 10];
m(74) = 999;

t(75,:) = [10 27 26];
m(75) = 999;

t(76,:) = [27 10 9];
m(76) = 999;

t(77,:) = [31 26 27];
m(77) = 999;

t(78,:) = [27 32 31];
m(78) = 999;

t(79,:) = [32 44 31];
m(79) = 4;

t(80,:) = [44 43 31];
m(80) = 4;

t(81,:) = [43 30 31];
m(81) = 4;

t(82,:) = [31 30 26];
m(82) = 999;

t(83,:) = [35 3 23];
m(83) = 4;

end