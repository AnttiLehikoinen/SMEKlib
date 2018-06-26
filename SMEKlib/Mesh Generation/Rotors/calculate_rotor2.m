function [p,t,m,FL,LL,ag] = calculate_rotor2(dim)

%Boundaries
FL = [1 42 41 40 39 38 37];
LL = [1 24 25 26 27 28 29];
ag = [37 36 35 34 33 32 31 30 29];

%prealloacte memory
p = zeros(42,2);
t = zeros(62,3);
m = zeros(62,1);

%help variables
L = dim.R_height1+dim.R_width2/2+dim.R_height3+dim.R_width3/2;    %length of slot
line = [(dim.Rin+dim.Rout-L)/2 dim.Rout-L dim.Rout-(dim.R_height1+dim.R_width2/2+dim.R_height3) dim.Rout-(dim.R_height1+dim.R_width2/2+dim.R_height3/2) dim.Rout-(dim.R_height1+dim.R_width2/2) dim.Rout-dim.R_height1 dim.Rout-dim.R_height1/2];
c = line(5)/2+line(6)/2;
alfa = atan((dim.R_width2/4+dim.R_width1/4)/(line(5)/2+line(6)/2));

%init the first sector nodes and elements


%Nodes
p(2,1) = line(1)*cos(dim.angleR(1)/2);
p(2,2) = line(1)*sin(dim.angleR(1)/2);

p(3,1) = (0.2*line(1)+0.8*line(2))*cos(dim.angleR(1)/2);
p(3,2) = (0.2*line(1)+0.8*line(2))*sin(dim.angleR(1)/2);

p(4,1) = line(2)*cos(dim.angleR(1)/2);
p(4,2) = line(2)*sin(dim.angleR(1)/2);

p(5,1) = sqrt(line(3)^2+(dim.R_width3/2)^2)*cos(dim.angleR(1)/2-atan((dim.R_width3/2)/line(3)));
p(5,2) = sqrt(line(3)^2+(dim.R_width3/2)^2)*sin(dim.angleR(1)/2-atan((dim.R_width3/2)/line(3)));

p(6,1) = sqrt(line(4)^2+((dim.R_width3/2+dim.R_width2/2)/2)^2)*cos(dim.angleR(1)/2-atan(((dim.R_width3/2+dim.R_width2/2)/2)/line(4)));
p(6,2) = sqrt(line(4)^2+((dim.R_width3/2+dim.R_width2/2)/2)^2)*sin(dim.angleR(1)/2-atan(((dim.R_width3/2+dim.R_width2/2)/2)/line(4)));

p(7,1) = sqrt(line(5)^2+(dim.R_width2/2)^2)*cos(dim.angleR(1)/2-atan((dim.R_width2/2)/line(5)));
p(7,2) = sqrt(line(5)^2+(dim.R_width2/2)^2)*sin(dim.angleR(1)/2-atan((dim.R_width2/2)/line(5)));

p(8,1) = c*cos(dim.angleR(1)/2-alfa);
p(8,2) = c*sin(dim.angleR(1)/2-alfa);

p(9,1) = sqrt(line(6)^2+(dim.R_width1/2)^2)*cos(dim.angleR(1)/2-atan((dim.R_width1/2)/line(6)));
p(9,2) = sqrt(line(6)^2+(dim.R_width1/2)^2)*sin(dim.angleR(1)/2-atan((dim.R_width1/2)/line(6)));

p(10,1) = line(7)*cos(dim.angleR(1)/2);
p(10,2) = line(7)*sin(dim.angleR(1)/2);

p(11,1) = sqrt(line(6)^2+(dim.R_width1/2)^2)*cos(dim.angleR(1)/2+atan((dim.R_width1/2)/line(6)));
p(11,2) = sqrt(line(6)^2+(dim.R_width1/2)^2)*sin(dim.angleR(1)/2+atan((dim.R_width1/2)/line(6)));

p(12,1) = c*cos(dim.angleR(1)/2+alfa);
p(12,2) = c*sin(dim.angleR(1)/2+alfa);

p(13,1) = sqrt(line(5)^2+(dim.R_width2/2)^2)*cos(dim.angleR(1)/2+atan((dim.R_width2/2)/line(5)));
p(13,2) = sqrt(line(5)^2+(dim.R_width2/2)^2)*sin(dim.angleR(1)/2+atan((dim.R_width2/2)/line(5)));

p(14,1) = sqrt(line(4)^2+((dim.R_width3/2+dim.R_width2/2)/2)^2)*cos(dim.angleR(1)/2+atan(((dim.R_width3/2+dim.R_width2/2)/2)/line(4)));
p(14,2) = sqrt(line(4)^2+((dim.R_width3/2+dim.R_width2/2)/2)^2)*sin(dim.angleR(1)/2+atan(((dim.R_width3/2+dim.R_width2/2)/2)/line(4)));

p(15,1) = sqrt(line(3)^2+(dim.R_width3/2)^2)*cos(dim.angleR(1)/2+atan((dim.R_width3/2)/line(3)));
p(15,2) = sqrt(line(3)^2+(dim.R_width3/2)^2)*sin(dim.angleR(1)/2+atan((dim.R_width3/2)/line(3)));

p(16,1) = line(3)*cos(dim.angleR(1)/2);
p(16,2) = line(3)*sin(dim.angleR(1)/2);

p(17,1) = line(5)*cos(dim.angleR(1)/2);
p(17,2) = line(5)*sin(dim.angleR(1)/2);

p(18,1) = sqrt(line(6)^2+(dim.R_width1/2)^2)*cos(5*dim.angleR(1)/6);
p(18,2) = sqrt(line(6)^2+(dim.R_width1/2)^2)*sin(5*dim.angleR(1)/6);

p(19,1) = sqrt(line(6)^2+(dim.R_width1/2)^2)*cos(dim.angleR(1)/2+atan((dim.R_width2/2)/line(5)));
p(19,2) = sqrt(line(6)^2+(dim.R_width1/2)^2)*sin(dim.angleR(1)/2+atan((dim.R_width2/2)/line(5)));

p(20,1) = sqrt(line(6)^2+(dim.R_width1/2)^2)*cos(dim.angleR(1)/2-atan((dim.R_width2/2)/line(5)));
p(20,2) = sqrt(line(6)^2+(dim.R_width1/2)^2)*sin(dim.angleR(1)/2-atan((dim.R_width2/2)/line(5)));

p(21,1) = sqrt(line(6)^2+(dim.R_width1/2)^2)*cos(dim.angleR(1)/6);
p(21,2) = sqrt(line(6)^2+(dim.R_width1/2)^2)*sin(dim.angleR(1)/6);

p(22,1) = line(2)*cos(5*dim.angleR(1)/6);
p(22,2) = line(2)*sin(5*dim.angleR(1)/6);

p(23,1) = line(2)*cos(dim.angleR(1)/6);
p(23,2) = line(2)*sin(dim.angleR(1)/6);

p(24,1) = dim.Rin*cos(dim.angleR(1));
p(24,2) = dim.Rin*sin(dim.angleR(1));

p(25,1) = (line(1)+line(2))/2*cos(dim.angleR(1));                       
p(25,2) = (line(1)+line(2))/2*sin(dim.angleR(1));               

p(26,1) = line(3)*cos(dim.angleR(1));
p(26,2) = line(3)*sin(dim.angleR(1));

p(27,1) = line(5)*cos(dim.angleR(1));
p(27,2) = line(5)*sin(dim.angleR(1));

p(28,1) = sqrt(line(6)^2+(dim.R_width1/2)^2)*cos(dim.angleR(1));
p(28,2) = sqrt(line(6)^2+(dim.R_width1/2)^2)*sin(dim.angleR(1));

p(29,1) = dim.Rout*cos(dim.angleR(1));
p(29,2) = dim.Rout*sin(dim.angleR(1));

p(30,1) = dim.Rout*cos(5*dim.angleR(1)/6);
p(30,2) = dim.Rout*sin(5*dim.angleR(1)/6);

p(31,1) = dim.Rout*cos(dim.angleR(1)/2+atan((dim.R_width2/2)/line(5)));
p(31,2) = dim.Rout*sin(dim.angleR(1)/2+atan((dim.R_width2/2)/line(5)));

p(32,1) = sqrt(dim.Rout^2+(dim.R_width1/2)^2)*cos(dim.angleR(1)/2+atan((dim.R_width1/2)/dim.Rout));
p(32,2) = sqrt(dim.Rout^2+(dim.R_width1/2)^2)*sin(dim.angleR(1)/2+atan((dim.R_width1/2)/dim.Rout));

p(33,1) = dim.Rout*cos(dim.angleR(1)/2);
p(33,2) = dim.Rout*sin(dim.angleR(1)/2);

p(34,1) = sqrt(dim.Rout^2+(dim.R_width1/2)^2)*cos(dim.angleR(1)/2-atan((dim.R_width1/2)/dim.Rout));
p(34,2) = sqrt(dim.Rout^2+(dim.R_width1/2)^2)*sin(dim.angleR(1)/2-atan((dim.R_width1/2)/dim.Rout));

p(35,1) = dim.Rout*cos(dim.angleR(1)/2-atan((dim.R_width2/2)/line(5)));
p(35,2) = dim.Rout*sin(dim.angleR(1)/2-atan((dim.R_width2/2)/line(5)));

p(36,1) = dim.Rout*cos(dim.angleR(1)/6);
p(36,2) = dim.Rout*sin(dim.angleR(1)/6);

p(37,1) = dim.Rout;

p(38,1) = sqrt(line(6)^2+(dim.R_width1/2)^2);

p(39,1) = line(5);

p(40,1) = line(3);

p(41,1) = (line(1)+line(2))/2;

p(42,1) = dim.Rin;


%Elements

t(1,:) = [1 42 24];
m(1) = dim.SFM;

t(2,:) = [42 2 24];
m(2) = dim.RM;

t(3,:) = [42 41 2];
m(3) = dim.RM;

t(4,:) = [2 25 24];
m(4) = dim.RM;

t(5,:) = [41 3 2];
m(5) = dim.RM;

t(6,:) = [2 3 25];
m(6) = dim.RM;

t(7,:) = [41 40 23];
m(7) = dim.RM;

t(8,:) = [41 23 3];
m(8) = dim.RM;

t(9,:) = [23 4 3];
m(9) = dim.RM;

t(10,:) = [23 5 4];
m(10) = dim.RM;

t(11,:) = [40 5 23];
m(11) = dim.RM;

t(12,:) = [3 4 22];
m(12) = dim.RM;

t(13,:) = [4 15 22];
m(13) = dim.RM;

t(14,:) = [3 22 25];
m(14) = dim.RM;

t(15,:) = [22 26 25];
m(15) = dim.RM;

t(16,:) = [15 26 22];
m(16) = dim.RM;

t(17,:) = [40 6 5];
m(17) = dim.RM;

t(18,:) = [40 39 6];
m(18) = dim.RM;

t(19,:) = [39 7 6];
m(19) = dim.RM;

t(20,:) = [39 38 21];
m(20) = dim.RM;

t(21,:) = [39 21 7];
m(21) = dim.RM;

t(22,:) = [38 37 21];
m(22) = dim.RM;

t(23,:) = [37 36 21];
m(23) = dim.RM;

t(24,:) = [36 35 21];
m(24) = dim.RM;

t(25,:) = [35 20 21];
m(25) = dim.RM;

t(26,:) = [21 20 7];
m(26) = dim.RM;

t(27,:) = [20 8 7];
m(27) = dim.RM;

t(28,:) = [20 9 8];
m(28) = dim.RM;

t(29,:) = [20 35 9];
m(29) = dim.RM;

t(30,:) = [35 34 9];
m(30) = dim.RM;

t(31,:) = [32 31 11];
m(31) = dim.RM;

t(32,:) = [31 19 11];
m(32) = dim.RM;

t(33,:) = [19 12 11];
m(33) = dim.RM;

t(34,:) = [19 13 12];
m(34) = dim.RM;

t(35,:) = [31 30 18];
m(35) = dim.RM;

t(36,:) = [31 18 19];
m(36) = dim.RM;

t(37,:) = [30 29 18];
m(37) = dim.RM;

t(38,:) = [29 28 18];
m(38) = dim.RM;

t(39,:) = [13 19 18];
m(39) = dim.RM;

t(40,:) = [13 18 27];
m(40) = dim.RM;

t(41,:) = [28 27 18];
m(41) = dim.RM;

t(42,:) = [27 14 13];
m(42) = dim.RM;

t(43,:) = [27 26 14];
m(43) = dim.RM;

t(44,:) = [26 15 14];
m(44) = dim.RM;

t(45,:) = [4 5 16];
m(45) = 9999;

t(46,:) = [4 16 15];
m(46) = 9999;

t(47,:) = [5 6 16];
m(47) = 9999;

t(48,:) = [6 14 16];
m(48) = 9999;

t(49,:) = [14 15 16];
m(49) = 9999;

t(50,:) = [6 7 17];
m(50) = 9999;

t(51,:) = [6 17 14];
m(51) = 9999;

t(52,:) = [17 13 14];
m(52) = 9999;

t(53,:) = [7 8 17];
m(53) = 9999;

t(54,:) = [8 9 17];
m(54) = 9999;

t(55,:) = [9 11 17];
m(55) = 9999;

t(56,:) = [11 12 17];
m(56) = 9999;

t(57,:) = [12 13 17];
m(57) = 9999;

t(58,:) = [9 10 11];
m(58) = dim.RO;

t(59,:) = [9 34 10];
m(59) = dim.RO;

t(60,:) = [34 33 10];
m(60) = dim.RO;

t(61,:) = [33 32 10];
m(61) = dim.RO;

t(62,:) = [32 11 10];
m(62) = dim.RO;

end