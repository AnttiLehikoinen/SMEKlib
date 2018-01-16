function msh = calculate_rotor5(msh,dim)

%preallocation of memory
p = zeros(42,2);
t = zeros(62,3);
m = zeros(62,1);


%help variables
line = [dim.Rin dim.Rout-dim.hc_r2-dim.hc_r-dim.wc2_r/2 dim.Rout-dim.hc_r2-dim.hc_r ...
        dim.Rout-dim.hc_r2-dim.hc_r/2 dim.Rout-dim.hc_r2 dim.Rout-dim.hc_r2+dim.wc3_r/2 ...
        dim.Rout-dim.htt_r-dim.wc1_r dim.Rout-dim.htt_r-dim.wc1_r/2 dim.Rout-dim.htt_r dim.Rout];
c = sqrt(line(8)^2+(dim.wc1_r/2)^2-2*((line(8)+line(9))/2)*(dim.wc1_r/2)*cosd(135));  
alfa = asin((dim.wc1_r/2*sind(135))/c);


%Nodes

p(2,1) = (line(1)+line(2))/2*cos(dim.angleR(1)/2);
p(2,2) = (line(1)+line(2))/2*sin(dim.angleR(1)/2);

p(3,1) = line(2)*cos(dim.angleR(1)/2);
p(3,2) = line(2)*sin(dim.angleR(1)/2);

p(4,1) = sqrt(line(3)^2+(dim.wc2_r/2)^2)*cos(dim.angleR(1)/2-atan((dim.wc2_r/2)/line(3)));
p(4,2) = sqrt(line(3)^2+(dim.wc2_r/2)^2)*sin(dim.angleR(1)/2-atan((dim.wc2_r/2)/line(3)));

p(5,1) = sqrt(line(4)^2+(dim.wc2_r/4+dim.wc3_r/4)^2)*cos(dim.angleR(1)/2-atan((dim.wc2_r/4+dim.wc3_r/4)/line(4)));
p(5,2) = sqrt(line(4)^2+(dim.wc2_r/4+dim.wc3_r/4)^2)*sin(dim.angleR(1)/2-atan((dim.wc2_r/4+dim.wc3_r/4)/line(4)));

p(6,1) = sqrt(line(5)^2+(dim.wc3_r/2)^2)*cos(dim.angleR(1)/2-atan((dim.wc3_r/2)/line(5)));
p(6,2) = sqrt(line(5)^2+(dim.wc3_r/2)^2)*sin(dim.angleR(1)/2-atan((dim.wc3_r/2)/line(5)));

p(7,1) = sqrt(line(6)^2+(dim.wc2_r/2)^2)*cos(dim.angleR(1)/2-atan((dim.wc2_r/2)/line(6)));
p(7,2) = sqrt(line(6)^2+(dim.wc2_r/2)^2)*sin(dim.angleR(1)/2-atan((dim.wc2_r/2)/line(6)));

p(8,1) = sqrt(line(7)^2+(dim.wc2_r/2)^2)*cos(dim.angleR(1)/2-atan((dim.wc2_r/2)/line(7)));
p(8,2) = sqrt(line(7)^2+(dim.wc2_r/2)^2)*sin(dim.angleR(1)/2-atan((dim.wc2_r/2)/line(7)));

p(9,1) = sqrt(line(8)^2+(dim.wc1_r/2)^2)*cos(dim.angleR(1)/2-atan((dim.wc1_r/2)/line(8)));
p(9,2) = sqrt(line(8)^2+(dim.wc1_r/2)^2)*sin(dim.angleR(1)/2-atan((dim.wc1_r/2)/line(8)));

p(10,1) = c*cos(dim.angleR(1)/2-alfa);
p(10,2) = c*sin(dim.angleR(1)/2-alfa);

p(11,1) = line(9)*cos(dim.angleR(1)/2);
p(11,2) = line(9)*sin(dim.angleR(1)/2);

p(12,1) = c*cos(dim.angleR(1)/2+alfa);
p(12,2) = c*sin(dim.angleR(1)/2+alfa);

p(13,1) = sqrt(line(8)^2+(dim.wc1_r/2)^2)*cos(dim.angleR(1)/2+atan((dim.wc1_r/2)/line(8)));
p(13,2) = sqrt(line(8)^2+(dim.wc1_r/2)^2)*sin(dim.angleR(1)/2+atan((dim.wc1_r/2)/line(8)));

p(14,1) = sqrt(line(7)^2+(dim.wc2_r/2)^2)*cos(dim.angleR(1)/2+atan((dim.wc2_r/2)/line(7)));
p(14,2) = sqrt(line(7)^2+(dim.wc2_r/2)^2)*sin(dim.angleR(1)/2+atan((dim.wc2_r/2)/line(7)));

p(15,1) = sqrt(line(6)^2+(dim.wc2_r/2)^2)*cos(dim.angleR(1)/2+atan((dim.wc2_r/2)/line(6)));
p(15,2) = sqrt(line(6)^2+(dim.wc2_r/2)^2)*sin(dim.angleR(1)/2+atan((dim.wc2_r/2)/line(6)));

p(16,1) = sqrt(line(5)^2+(dim.wc3_r/2)^2)*cos(dim.angleR(1)/2+atan((dim.wc3_r/2)/line(5)));
p(16,2) = sqrt(line(5)^2+(dim.wc3_r/2)^2)*sin(dim.angleR(1)/2+atan((dim.wc3_r/2)/line(5)));

p(17,1) = sqrt(line(4)^2+(dim.wc2_r/4+dim.wc3_r/4)^2)*cos(dim.angleR(1)/2+atan((dim.wc2_r/4+dim.wc3_r/4)/line(4)));
p(17,2) = sqrt(line(4)^2+(dim.wc2_r/4+dim.wc3_r/4)^2)*sin(dim.angleR(1)/2+atan((dim.wc2_r/4+dim.wc3_r/4)/line(4)));

p(18,1) = sqrt(line(3)^2+(dim.wc2_r/2)^2)*cos(dim.angleR(1)/2+atan((dim.wc2_r/2)/line(3)));
p(18,2) = sqrt(line(3)^2+(dim.wc2_r/2)^2)*sin(dim.angleR(1)/2+atan((dim.wc2_r/2)/line(3)));

p(19,1) = line(3)*cos(dim.angleR(1)/2);
p(19,2) = line(3)*sin(dim.angleR(1)/2);

p(20,1) = line(5)*cos(dim.angleR(1)/2);
p(20,2) = line(5)*sin(dim.angleR(1)/2);

p(21,1) = line(8)*cos(dim.angleR(1)/2);
p(21,2) = line(8)*sin(dim.angleR(1)/2);

p(22,1) = line(9)*cos(4*dim.angleR(1)/5);
p(22,2) = line(9)*sin(4*dim.angleR(1)/5);

p(23,1) = line(9)*cos(dim.angleR(1)/5);
p(23,2) = line(9)*sin(dim.angleR(1)/5);

p(24,1) = line(1)*cos(dim.angleR(1));
p(24,2) = line(1)*sin(dim.angleR(1));

p(25,1) = ((line(1)+line(2))/2+line(2))/2*cos(dim.angleR(1));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%((line(1)+line(2))/2+line(2))/2
p(25,2) = ((line(1)+line(2))/2+line(2))/2*sin(dim.angleR(1));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p(26,1) = line(3)*cos(dim.angleR(1));
p(26,2) = line(3)*sin(dim.angleR(1));

p(27,1) = line(5)*cos(dim.angleR(1));
p(27,2) = line(5)*sin(dim.angleR(1));

p(28,1) = line(8)*cos(dim.angleR(1));
p(28,2) = line(8)*sin(dim.angleR(1));

p(29,1) = line(9)*cos(dim.angleR(1));
p(29,2) = line(9)*sin(dim.angleR(1));

p(30,1) = line(10)*cos(dim.angleR(1));
p(30,2) = line(10)*sin(dim.angleR(1));

p(31,1) = line(10)*cos(4*dim.angleR(1)/5);
p(31,2) = line(10)*sin(4*dim.angleR(1)/5);

p(32,1) = sqrt(line(10)^2+(dim.wc1_r/2)^2)*cos(dim.angleR(1)/2+atan((dim.wc1_r/2)/line(10)));
p(32,2) = sqrt(line(10)^2+(dim.wc1_r/2)^2)*sin(dim.angleR(1)/2+atan((dim.wc1_r/2)/line(10)));

p(33,1) = line(10)*cos(dim.angleR(1)/2);
p(33,2) = line(10)*sin(dim.angleR(1)/2);

p(34,1) = sqrt(line(10)^2+(dim.wc1_r/2)^2)*cos(dim.angleR(1)/2-atan((dim.wc1_r/2)/line(10)));
p(34,2) = sqrt(line(10)^2+(dim.wc1_r/2)^2)*sin(dim.angleR(1)/2-atan((dim.wc1_r/2)/line(10)));

p(35,1) = line(10)*cos(dim.angleR(1)/5);
p(35,2) = line(10)*sin(dim.angleR(1)/5);

p(36,1) = line(10);

p(37,1) = line(9);

p(38,1) = line(8);

p(39,1) = line(5);

p(40,1) = line(3);

p(41,1) = ((line(1)+line(2))/2+line(2))/2;%%%%%%%%%%%%%%%%%%%((line(1)+line(2))/2+line(2))/2

p(42,1) = line(1);



%Elements

t(1,:) = [1 42 24];
m(1) = dim.SFM;

t(2,:) = [42 2 24];
m(2) = dim.RM;

t(3,:) = [42 41 2];
m(3) = dim.RM;

t(4,:) = [2 25 24];
m(4) = dim.RM;

t(5,:) = [41 25 2];
m(5) = dim.RM;

t(6,:) = [41 3 25];
m(6) = dim.RM;

t(7,:) = [41 40 3];
m(7) = dim.RM;

t(8,:) = [3 26 25];
m(8) = dim.RM;

t(9,:) = [40 4 3];
m(9) = dim.RM;

t(10,:) = [40 5 4];
m(10) = dim.RM;

t(11,:) = [40 39 5];
m(11) = dim.RM; 

t(12,:) = [39 6 5];
m(12) = dim.RM;

t(13,: ) = [39 7 6];
m(13) = dim.RM;

t(14,:) = [39 38 7];
m(14) = dim.RM;

t(15,:) = [38 8 7];
m(15) = dim.RM;

t(16,:) = [38 9 8];
m(16) = dim.RM;

t(17,:) = [38 23 9];
m(17) = dim.RM;

t(18,:) = [38 37 23];
m(18) = dim.RM;

t(19,:) = [37 36 23];
m(19) = dim.RM;

t(20,:) = [36 35 23];
m(20) = dim.RM;

t(21,:) = [35 34 23];
m(21) = dim.RM;

t(22,:) = [34 10 23];
m(22) = dim.RM;

t(23,:) = [23 10 9];
m(23) = dim.RM;

t(24,:) = [34 11 10];
m(24) = dim.RM;

t(25,:) = [34 33 11];
m(25) = dim.RM;

t(26,:) = [11 33 32];
m(26) = dim.RM;

t(27,:) = [32 12 11];
m(27) = dim.RM;

t(28,:) = [32 31 22];
m(28) = dim.RM;

t(29,:) = [32 22 12];
m(29) = dim.RM;

t(30,:) = [22 13 12];
m(30) = dim.RM;

t(31,:) = [31 30 22];
m(31) = dim.RM;

t(32,:) = [30 29 22];
m(32) = dim.RM;

t(33,:) = [29 28 22];
m(33) = dim.RM;

t(34,:) = [28 13 22];
m(34) = dim.RM;

t(35,:) = [28 14 13];
m(35) = dim.RM;

t(36,:) = [28 15 14];
m(36) = dim.RM;

t(37,:) = [28 27 15];
m(37) = dim.RM;

t(38,:) = [27 16 15];
m(38) = dim.RM;

t(39,:) = [27 17 16];
m(39) = dim.RM;

t(40,:) = [27 26 17];
m(40) = dim.RM;

t(41,:) = [26 18 17];
m(41) = dim.RM;

t(42,:) = [26 3 18];
m(42) = dim.RM;

t(43,:) =[19 3 4];
m(43) = 9999;

t(44,:) = [3 19 18];
m(44) = 9999;

t(45,:) = [4 5 19];
m(45) = 9999;

t(46,:) = [5 17 19];
m(46) = 9999;

t(47,:) = [17 18 19];
m(47) = 9999;

t(48,:) = [5 6 20];
m(48) = 9999;

t(49,:) = [5 20 17];
m(49) = 9999;

t(50,:) = [16 17 20];
m(50) = 9999;

t(51,:) = [6 7 20];
m(51) = 9999;

t(52,:) = [7 15 20];
m(52) = 9999;

t(53,:) = [20 15 16];
m(53) = 9999;

t(54,:) = [7 8 14];
m(54) = 9999;

t(55,:) = [7 14 15];
m(55) = 9999;

t(56,:) = [8 21 14];
m(56) = 9999;

t(57,:) = [8 9 21];
m(57) = 9999;

t(58,:) = [9 10 21];
m(58) = 9999;

t(59,:) = [10 11 21];
m(59) = 9999;

t(60,:) = [11 12 21];
m(60) = 9999;

t(61,:) = [12 13 21];
m(61) = 9999;

t(62,:) = [13 14 21];
m(62) = 9999;

%additional quality improvements
p(37,:) = 0.5*p(36,:) + 0.5*p(38,:);
p(29,:) = 0.5*p(30,:) + 0.5*p(28,:);

msh.t = t;
msh.matel = m;
msh.p = p;
msh.index_p = size(p,1);

end
