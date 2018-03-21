function [p,t,m,FL,LL,ag] = calculate_rotor3(dim)

%Boundaries
FL = [1 38 37 36 35 34 33 32];
LL = [1 20 21 22 23 24 25 26];
ag = [32 31 30 29 28 27 26];

%allocating memory
p = zeros(38,2);
t = zeros(55,3);
m = zeros(55,1);

%help variables
L = dim.Rout-dim.R_height;
line = [dim.Rin dim.Rin/2+L/2 dim.Rin/4+L/4+L/2 L dim.Rout-dim.R_height+dim.R_width2/2 L/4+dim.Rout/2-dim.R_width2/2 ...
        dim.Rout-dim.R_height+dim.R_width2 dim.Rout/2+(dim.Rout-dim.R_height+dim.R_width2)/2 dim.Rout];

c = sqrt(line(5)^2+(dim.R_width2/2)^2-2*((line(5)+line(5))/2)*(dim.R_width2/2)*cosd(135));  
alfa = asin((dim.R_width2/2*sind(135))/c);

c1 = sqrt(line(5)^2+(dim.R_width2/2)^2-2*((line(5)/2+line(4)/2))*(dim.R_width2/2)*cosd(45));  
alfa1 = asin((dim.R_width2/2*sind(45))/c1);

p(2,1) = line(3)*cos(dim.angleR(1)/2);
p(2,2) = line(3)*sin(dim.angleR(1)/2);

p(3,1) = line(4)*cos(dim.angleR(1)/2);
p(3,2) = line(4)*sin(dim.angleR(1)/2);

p(4,1) = c1*cos(dim.angleR(1)/2-alfa1);
p(4,2) = c1*sin(dim.angleR(1)/2-alfa1);

p(5,1) = sqrt(line(5)^2+(dim.R_width2/2)^2)*cos(dim.angleR(1)/2-atan(dim.R_width2/2/line(5)));
p(5,2) = sqrt(line(5)^2+(dim.R_width2/2)^2)*sin(dim.angleR(1)/2-atan(dim.R_width2/2/line(5)));

p(6,1) = c*cos(dim.angleR(1)/2-alfa);
p(6,2) = c*sin(dim.angleR(1)/2-alfa);

p(7,1) = line(7)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(7)));
p(7,2) = line(7)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(7)));

p(8,1) = line(7)*cos(dim.angleR(1)/2);
p(8,2) = line(7)*sin(dim.angleR(1)/2);

p(9,1) = line(7)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(7)));
p(9,2) = line(7)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(7)));

p(10,1) = c*cos(dim.angleR(1)/2+alfa);
p(10,2) = c*sin(dim.angleR(1)/2+alfa);

p(11,1) = sqrt(line(5)^2+(dim.R_width2/2)^2)*cos(dim.angleR(1)/2+atan(dim.R_width2/2/line(5)));
p(11,2) = sqrt(line(5)^2+(dim.R_width2/2)^2)*sin(dim.angleR(1)/2+atan(dim.R_width2/2/line(5)));

p(12,1) = c1*cos(dim.angleR(1)/2+alfa1);
p(12,2) = c1*sin(dim.angleR(1)/2+alfa1);

p(13,1) = line(5)*cos(dim.angleR(1)/2);
p(13,2) = line(5)*sin(dim.angleR(1)/2);

p(14,1) = c*cos(dim.angleR(1)/2);
p(14,2) = c*sin(dim.angleR(1)/2);

p(15,1) = line(8)*cos(dim.angleR(1)/2+dim.angleR(1)/4);
p(15,2) = line(8)*sin(dim.angleR(1)/2+dim.angleR(1)/4);

p(16,1) = line(8)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(8)));
p(16,2) = line(8)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(8)));

p(17,1) = line(8)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(8)));
p(17,2) = line(8)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(8)));

p(18,1) = line(8)*cos(dim.angleR(1)/2-dim.angleR(1)/4);
p(18,2) = line(8)*sin(dim.angleR(1)/2-dim.angleR(1)/4);

p(19,1) = (line(1)/2+line(2)/2)*cos(dim.angleR(1)/2);
p(19,2) = (line(1)/2+line(2)/2)*sin(dim.angleR(1)/2);

p(20,1) = line(1)*cos(dim.angleR(1));
p(20,2) = line(1)*sin(dim.angleR(1)); 

p(21,1) = line(2)*cos(dim.angleR(1));
p(21,2) = line(2)*sin(dim.angleR(1)); 

p(22,1) = line(3)*cos(dim.angleR(1));
p(22,2) = line(3)*sin(dim.angleR(1)); 

p(23,1) = line(5)*cos(dim.angleR(1));
p(23,2) = line(5)*sin(dim.angleR(1)); 

p(24,1) = line(7)*cos(dim.angleR(1));
p(24,2) = line(7)*sin(dim.angleR(1)); 

p(25,1) = line(8)*cos(dim.angleR(1));
p(25,2) = line(8)*sin(dim.angleR(1)); 

p(26,1) = line(9)*cos(dim.angleR(1));
p(26,2) = line(9)*sin(dim.angleR(1)); 

p(27,1) = line(9)*cos(dim.angleR(1)-dim.angleR(1)/4);
p(27,2) = line(9)*sin(dim.angleR(1)-dim.angleR(1)/4);

p(28,1) = line(9)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(9)));
p(28,2) = line(9)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(9)));

p(29,1) = line(9)*cos(dim.angleR(1)/2);
p(29,2) = line(9)*sin(dim.angleR(1)/2);

p(30,1) = line(9)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(9)));
p(30,2) = line(9)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(9)));

p(31,1) = line(9)*cos(dim.angleR(1)/4);
p(31,2) = line(9)*sin(dim.angleR(1)/4);

p(32,1) = line(9);

p(33,1) = line(8);

p(34,1) = line(7);

p(35,1) = line(5);

p(36,1) = line(3);

p(37,1) = line(2);

p(38,1) = line(1);

t(1,:) = [1 38 20];
m(1) = dim.SFM;

t(2,:) = [38 19 20];
m(2) = dim.RM;

t(3,:) = [38 37 19];
m(3) = dim.RM;

t(4,:) = [19 21 20];
m(4) = dim.RM;

t(5,:) = [37 21 19];
m(5) = dim.RM;

t(6,:) = [37 36 2];
m(6) = dim.RM;

t(7,:) = [37 2 21];
m(7) = dim.RM;

t(8,:) = [2 22 21];
m(8) = dim.RM;

t(9,:) = [36 35 4];
m(9) = dim.RM;

t(10,:) = [36 4 2];
m(10) = dim.RM;

t(11,:) = [4 3 2];
m(11) = dim.RM;

t(12,:) = [2 3 12];
m(12) = dim.RM;

t(13,:) = [2 12 22];
m(13) = dim.RM;

t(14,:) = [12 23 22];
m(14) = dim.RM;

t(15,:) = [35 5 4];
m(15) = dim.RM;

t(16,:) = [35 34 5];
m(16) = dim.RM;

t(17,:) = [34 6 5];
m(17) = dim.RM;

t(18,:) = [34 18 6];
m(18) = dim.RM;

t(19,:) = [34 33 18];
m(19) = dim.RM;

t(20,:) = [33 31 18];
m(20) = dim.RM;

t(21,:) = [33 32 31];
m(21) = dim.RM;

t(22,:) = [18 7 6];
m(22) = dim.RM;

t(23,:) = [18 17 7];
m(23) = dim.RM;

t(24,:) = [18 31 17];
m(24) = dim.RM;

t(25,:) = [31 30 17];
m(25) = dim.RM;

t(26,:) = [16 28 27];
m(26) = dim.RM;

t(27,:) = [16 27 15];
m(27) = dim.RM;

t(28,:) = [16 27 15];
m(28) = dim.RM;

t(29,:) = [9 16 15];
m(29) = dim.RM;

t(30,:) = [10 9 15];
m(30) = dim.RM;

t(31,:) = [27 26 25];
m(31) = dim.RM;

t(32,:) = [15 27 25];
m(32) = dim.RM;

t(33,:) = [15 25 24];
m(33) = dim.RM;

t(34,:) = [10 15 24];
m(34) = dim.RM;

t(35,:) = [11 10 24];
m(35) = dim.RM;

t(36,:) = [11 24 23];
m(36) = dim.RM;

t(37,:) = [12 11 23];
m(37) = dim.RM;

t(38,:) = [17 30 29];
m(38) = dim.RO;%

t(39,:) = [17 29 16];
m(39) = dim.RO;%

t(40,:) = [16 29 28];
m(40) = dim.RO;%

t(41,:) = [7 17 8];
m(41) = dim.RO;%

t(42,:) = [17 16 8];
m(42) = dim.RO;%

t(43,:) = [8 16 9];
m(43) = dim.RO;%

t(44,:) = [3 4 13];
m(44) = 9999;

t(45,:) = [4 5 13];
m(45) = 9999;

t(46,:) = [5 6 13];
m(46) = 9999;

t(47,:) = [6 14 13];
m(47) = 9999;

t(48,:) = [6 7 14];
m(48) = 9999;

t(49,:) = [7 8 14];
m(49) = 9999;

t(50,:) = [8 9 14];
m(50) = 9999;

t(51,:) = [9 10 14];
m(51) = 9999;

t(52,:) = [14 10 13];
m(52) = 9999;

t(53,:)= [10 11 13];
m(53) = 9999;

t(54,:) = [11 12 13];
m(54) = 9999;

t(55,:) = [12 3 13];
m(55) = 9999;

end