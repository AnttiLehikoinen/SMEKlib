function [p,t,m,FL,LL,ag] = calculate_rotor6(dim)

FL = [1 41 40 39 38 37 36];
LL = [1 25 26 27 28 29 30];
ag = [36 35 34 33 32 31 30];

%prealloacte memory
p = zeros(41,2);
t = zeros(62,3);    
m = zeros(62,1);     


line = [dim.Rin dim.Rin/2+(dim.Rout-dim.R_height)/2 dim.Rout-dim.R_height ...
        dim.Rout-dim.R_height+dim.R_height3/2 dim.Rout-dim.R_height2 dim.Rout-dim.R_height1-dim.R_width2/2 ...
        dim.Rout-dim.R_height1 dim.Rout];
    
c = sqrt(line(6)^2+(dim.R_width2/2)^2-2*((line(6)+line(7))/2)*(dim.R_width2/2)*cosd(120));  
alfa = asin((dim.R_width2/2*sind(120))/c);

%Init the first sector

%Nodes

p(2,1) = line(2)*cos(dim.angleR(1)/2);
p(2,2) = line(2)*sin(dim.angleR(1)/2);

p(3,1) = (line(2)/4 +3*line(3)/4)*cos(dim.angleR(1)/2);
p(3,2) = (line(2)/4 +3*line(3)/4)*sin(dim.angleR(1)/2);

p(4,1) = line(3)*cos(dim.angleR(1)/2-atan(dim.R_width4/2/line(3)));
p(4,2) = line(3)*sin(dim.angleR(1)/2-atan(dim.R_width4/2/line(3)));

p(5,1) = line(4)*cos(dim.angleR(1)/2-atan((dim.R_width4/4+dim.R_width3/4)/line(4)));
p(5,2) = line(4)*sin(dim.angleR(1)/2-atan((dim.R_width4/4+dim.R_width3/4)/line(4)));

p(6,1) = line(5)*cos(dim.angleR(1)/2-atan(dim.R_width3/2/line(5)));
p(6,2) = line(5)*sin(dim.angleR(1)/2-atan(dim.R_width3/2/line(5)));

p(7,1) = line(5)*cos(dim.angleR(1)/2-atan(dim.R_width2/2/line(5)));
p(7,2) = line(5)*sin(dim.angleR(1)/2-atan(dim.R_width2/2/line(5)));

p(8,1) = line(6)*cos(dim.angleR(1)/2-atan(dim.R_width2/2/line(6)));
p(8,2) = line(6)*sin(dim.angleR(1)/2-atan(dim.R_width2/2/line(6)));

p(9,1) = c*cos(dim.angleR(1)/2-alfa);
p(9,2) = c*sin(dim.angleR(1)/2-alfa);

p(10,1) = line(7)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(7)));
p(10,2) = line(7)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(7)));

p(11,1) = line(7)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(7)));
p(11,2) = line(7)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(7)));

p(12,1) = c*cos(dim.angleR(1)/2+alfa);
p(12,2) = c*sin(dim.angleR(1)/2+alfa);

p(13,1) = line(6)*cos(dim.angleR(1)/2+atan(dim.R_width2/2/line(6)));
p(13,2) = line(6)*sin(dim.angleR(1)/2+atan(dim.R_width2/2/line(6)));

p(14,1) = line(5)*cos(dim.angleR(1)/2+atan(dim.R_width2/2/line(5)));
p(14,2) = line(5)*sin(dim.angleR(1)/2+atan(dim.R_width2/2/line(5)));

p(15,1) = line(5)*cos(dim.angleR(1)/2+atan(dim.R_width3/2/line(5)));
p(15,2) = line(5)*sin(dim.angleR(1)/2+atan(dim.R_width3/2/line(5)));

p(16,1) = line(4)*cos(dim.angleR(1)/2+atan((dim.R_width4/4+dim.R_width3/4)/line(4)));
p(16,2) = line(4)*sin(dim.angleR(1)/2+atan((dim.R_width4/4+dim.R_width3/4)/line(4)));

p(17,1) = line(3)*cos(dim.angleR(1)/2+atan(dim.R_width4/2/line(3)));
p(17,2) = line(3)*sin(dim.angleR(1)/2+atan(dim.R_width4/2/line(3)));

p(18,1) = line(4)*cos(dim.angleR(1)/2);
p(18,2) = line(4)*sin(dim.angleR(1)/2);

p(19,1) = line(5)*cos(dim.angleR(1)/2);
p(19,2) = line(5)*sin(dim.angleR(1)/2);

p(20,1) = line(6)*cos(dim.angleR(1)/2);
p(20,2) = line(6)*sin(dim.angleR(1)/2);

p(21,1) = line(6)*cos(dim.angleR(1)/2+atan(dim.R_width3/2/line(6)));
p(21,2) = line(6)*sin(dim.angleR(1)/2+atan(dim.R_width3/2/line(6)));

p(22,1) = line(7)*cos(dim.angleR(1)/2+atan(dim.R_width3/2/line(7)));
p(22,2) = line(7)*sin(dim.angleR(1)/2+atan(dim.R_width3/2/line(7)));

p(23,1) = line(6)*cos(dim.angleR(1)/2-atan(dim.R_width3/2/line(6)));
p(23,2) = line(6)*sin(dim.angleR(1)/2-atan(dim.R_width3/2/line(6)));

p(24,1) = line(7)*cos(dim.angleR(1)/2-atan(dim.R_width3/2/line(7)));
p(24,2) = line(7)*sin(dim.angleR(1)/2-atan(dim.R_width3/2/line(7)));

p(25,1) = line(1)*cos(dim.angleR(1));
p(25,2) = line(1)*sin(dim.angleR(1));

p(26,1) = (line(2)/2+line(3)/2)*cos(dim.angleR(1));
p(26,2) = (line(2)/2+line(3)/2)*sin(dim.angleR(1));

p(27,1) = line(3)*cos(dim.angleR(1));
p(27,2) = line(3)*sin(dim.angleR(1));

p(28,1) = line(5)*cos(dim.angleR(1));
p(28,2) = line(5)*sin(dim.angleR(1));

p(29,1) = line(7)*cos(dim.angleR(1));
p(29,2) = line(7)*sin(dim.angleR(1));

p(30,1) = line(8)*cos(dim.angleR(1));
p(30,2) = line(8)*sin(dim.angleR(1));

p(31,1) = line(8)*cos(dim.angleR(1)/2+atan(dim.R_width3/2/line(8)));
p(31,2) = line(8)*sin(dim.angleR(1)/2+atan(dim.R_width3/2/line(8)));

p(32,1) = line(8)*cos(dim.angleR(1)/2+atan(dim.R_width1/2/line(8)));
p(32,2) = line(8)*sin(dim.angleR(1)/2+atan(dim.R_width1/2/line(8)));

p(33,1) = line(8)*cos(dim.angleR(1)/2);
p(33,2) = line(8)*sin(dim.angleR(1)/2);

p(34,1) = line(8)*cos(dim.angleR(1)/2-atan(dim.R_width1/2/line(8)));
p(34,2) = line(8)*sin(dim.angleR(1)/2-atan(dim.R_width1/2/line(8)));

p(35,1) = line(8)*cos(dim.angleR(1)/2-atan(dim.R_width3/2/line(8)));
p(35,2) = line(8)*sin(dim.angleR(1)/2-atan(dim.R_width3/2/line(8)));

p(36,1) = line(8);

p(37,1) = line(7);

p(38,1) = line(5);

p(39,1) = line(3);

p(40,1) = line(2)/2+line(3)/2;

p(41,1) = line(1);


%Elements

t(1,:) = [1 41 25];
m(1) = dim.SFM;

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

t(7,:) = [40 39 3];
m(7) = dim.RM;

t(8,:) = [39 4 3];
m(8) = dim.RM;

t(9,:) = [39 5 4];
m(9) = dim.RM;

t(10,:) = [39 38 5];
m(10) = dim.RM;

t(11,:) = [38 6 5];
m(11) = dim.RM;

t(12,:) = [38 23 6];
m(12) = dim.RM;

t(13,:) = [38 37 23];
m(13) = dim.RM;

t(14,:) = [37 24 23];
m(14) = dim.RM;

t(15,:) = [37 36 24];
m(15) = dim.RM;

t(16,:) = [36 35 24];
m(16) = dim.RM;

t(17,:) = [6 23 7];
m(17) = dim.RM;

t(18,:) = [23 8 7 ];
m(18) = dim.RM;

t(19,:) = [23 9 8];
m(19) = dim.RM;

t(20,:) = [23 24 9];
m(20) = dim.RM;

t(21,:) = [24 10 9];
m(21) = dim.RM;

t(22,:)  = [24 34 10];
m(22) = dim.RM;

t(23,:) = [24 35 34];
m(23) = dim.RM;

t(24,:) = [32 31 22];
m(24) = dim.RM;

t(25,:)  = [32 22 11];
m(25) = dim.RM;

t(26,:) = [22 12 11];
m(26) = dim.RM;

t(27,:) = [22 21 12];
m(27) = dim.RM;

t(28,:) = [21 13 12];
m(28) = dim.RM;

t(29,:) = [21 14 13];
m(29) = dim.RM;

t(30,:) = [21 15 14];
m(30) = dim.RM;

t(31,:) = [31 30 22];
m(31) = dim.RM;

t(32,:) = [30 29 22];
m(32) = dim.RM;

t(33,:) = [29 21 22];
m(33) = dim.RM;

t(34,:) = [29 28 21];
m(34) = dim.RM;

t(35,:) = [28 15 21];
m(35) = dim.RM;

t(36,:) = [28 16 15];
m(36) = dim.RM;

t(37,:) = [28 27 16];
m(37) = dim.RM;

t(38,:) = [27 17 16];
m(38) = dim.RM;

t(39,:) = [27 3 17];
m(39) = dim.RM;

t(40,:) = [27 26 3];
m(40) = dim.RM;

t(41,:) = [3 4 17];
m(41) = dim.RM;

t(42,:) = [4 5 18];
m(42) = 9999;

t(43,:) = [4 18 17];
m(43) = 9999;

t(44,:) = [16 17 18];
m(44) = 9999;

t(45,:) = [5 6 7];
m(45) = 9999;

t(46,:) = [5 7 18];
m(46) = 9999;

t(47,:) = [7 19 18];
m(47) = 9999;

t(48,:) = [18 19 14];
m(48) = 9999;

t(49,:) = [18 14 16];
m(49) = 9999;

t(50,:) = [14 15 16];
m(50) = 9999;

t(51,:) = [7 8 19];
m(51) = 9999;

t(52,:) =[8 20 19];
m(52) = 9999;

t(53,:) = [19 20 13];
m(53) = 9999;

t(54,:) = [19 13 14];
m(54) = 9999;

t(55,:) = [8 9 20];
m(55) = 9999;

t(56,:) = [9 10 20];
m(56) = 9999;

t(57,:) = [10 11 20];
m(57) = 9999;

t(58,:) = [11 12 20];
m(58) = 9999;

t(59,:) = [12 13 20];
m(59) = 9999;

t(60,:) = [10 34 33];
m(60) = dim.RO;

t(61,:) = [10 33 11];
m(61) = dim.RO;

t(62,:) = [33 32 11];
m(62) =  dim.RO;

end