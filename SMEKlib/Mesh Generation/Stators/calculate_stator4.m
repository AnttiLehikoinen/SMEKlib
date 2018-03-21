function [p,t,m,FL,LL,cir,ag] = calculate_stator4(dim)
    

%preallocate memory
p = zeros(32,2);
t = zeros(44,3);
m = zeros(44,1);

%help variables
line = [dim.Sin dim.Sin+dim.S_height1 dim.Sin+dim.S_height-dim.S_width3/2-dim.S_height3 ...
        dim.Sin+dim.S_height-dim.S_width3/2 dim.Sin+dim.S_height (dim.Sin+dim.S_height+dim.Sout)/2 dim.Sout];
c = sqrt(line(4)^2+(dim.S_width3/2)^2-2*(dim.S_width3/2)*line(4)*cosd(135));
alfa = asin((dim.S_width3/2*sind(135))/c);

%Init the first sector

%Nodes
p(1,1) = line(1)*cos(dim.angleS(1)/2);
p(1,2) = line(1)*sin(dim.angleS(1)/2);

p(2,1) = line(1)*cos(dim.angleS(1)/2-asin(dim.S_width1/2/line(1)));
p(2,2) = line(1)*sin(dim.angleS(1)/2-asin(dim.S_width1/2/line(1)));

p(3,1) = line(2)*cos(dim.angleS(1)/2-asin(dim.S_width1/2/line(2)));
p(3,2) = line(2)*sin(dim.angleS(1)/2-asin(dim.S_width1/2/line(2)));

p(4,1) = line(3)*cos(dim.angleS(1)/2-atan((dim.S_width2/2)/line(3)));
p(4,2) = line(3)*sin(dim.angleS(1)/2-atan((dim.S_width2/2)/line(3)));

p(5,1) = line(4)*cos(dim.angleS(1)/2-atan(dim.S_width3/2/line(4)));
p(5,2) = line(4)*sin(dim.angleS(1)/2-atan(dim.S_width3/2/line(4)));

p(6,1) = c*cos(dim.angleS(1)/2-alfa);
p(6,2) = c*sin(dim.angleS(1)/2-alfa);

p(7,1) = line(5)*cos(dim.angleS(1)/2);
p(7,2) = line(5)*sin(dim.angleS(1)/2);

p(8,1) = c*cos(dim.angleS(1)/2+alfa);
p(8,2) = c*sin(dim.angleS(1)/2+alfa);

p(9,1) = line(4)*cos(dim.angleS(1)/2+atan(dim.S_width3/2/line(4)));
p(9,2) = line(4)*sin(dim.angleS(1)/2+atan(dim.S_width3/2/line(4)));

p(10,1) = line(3)*cos(dim.angleS(1)/2+atan((dim.S_width2/2)/line(3)));
p(10,2) = line(3)*sin(dim.angleS(1)/2+atan((dim.S_width2/2)/line(3)));

p(11,1) = line(2)*cos(dim.angleS(1)/2+asin(dim.S_width1/2/line(2)));
p(11,2) = line(2)*sin(dim.angleS(1)/2+asin(dim.S_width1/2/line(2)));

p(12,1) = line(1)*cos(dim.angleS(1)/2+asin(dim.S_width1/2/line(1)));
p(12,2) = line(1)*sin(dim.angleS(1)/2+asin(dim.S_width1/2/line(1)));

p(13,1) = line(2)*cos(dim.angleS(1)/2);
p(13,2) = line(2)*sin(dim.angleS(1)/2);

p(14,1) = line(4)*cos(dim.angleS(1)/2);
p(14,2) = line(4)*sin(dim.angleS(1)/2);

p(15,1) = ((dim.Sin+dim.S_height)*0.75+0.25*dim.Sout)*cos(dim.angleS(1)/2+dim.angleS(1)/4);
p(15,2) = ((dim.Sin+dim.S_height)*0.75+0.25*dim.Sout)*sin(dim.angleS(1)/2+dim.angleS(1)/4);

p(16,1) = ((dim.Sin+dim.S_height)*0.75+0.25*dim.Sout)*cos(dim.angleS(1)/2);
p(16,2) = ((dim.Sin+dim.S_height)*0.75+0.25*dim.Sout)*sin(dim.angleS(1)/2);

p(17,1) = ((dim.Sin+dim.S_height)*0.75+0.25*dim.Sout)*cos(dim.angleS(1)/2-dim.angleS(1)/4);
p(17,2) = ((dim.Sin+dim.S_height)*0.75+0.25*dim.Sout)*sin(dim.angleS(1)/2-dim.angleS(1)/4);

p(18,1) = line(1)*cos(dim.angleS(1)/6);
p(18,2) = line(1)*sin(dim.angleS(1)/6);

p(19,1) = line(1)*cos(dim.angleS(1)/2-atan(dim.S_width2/2/line(1)));
p(19,2) = line(1)*sin(dim.angleS(1)/2-atan(dim.S_width2/2/line(1)));

p(20,1) = line(1)*cos(dim.angleS(1)/2+(dim.S_width2/2/line(1)));
p(20,2) = line(1)*sin(dim.angleS(1)/2+atan(dim.S_width2/2/line(1)));

p(21,1) = line(1)*cos(5*dim.angleS(1)/6);
p(21,2) = line(1)*sin(5*dim.angleS(1)/6);

p(22,1) = line(1)*cos(dim.angleS(1));
p(22,2) = line(1)*sin(dim.angleS(1));

p(23,1) = line(3)*cos(dim.angleS(1));
p(23,2) = line(3)*sin(dim.angleS(1));

p(24,1) = line(4)*cos(dim.angleS(1));
p(24,2) = line(4)*sin(dim.angleS(1));

p(25,1) = ((dim.Sin+dim.S_height)*0.75+0.25*dim.Sout)*cos(dim.angleS(1));
p(25,2) = ((dim.Sin+dim.S_height)*0.75+0.25*dim.Sout)*sin(dim.angleS(1));

p(26,1) = line(7)*cos(dim.angleS(1));
p(26,2) = line(7)*sin(dim.angleS(1));

p(27,1) = line(7)*cos(dim.angleS(1)/2);
p(27,2) = line(7)*sin(dim.angleS(1)/2);

p(28,1) = line(7);

p(29,1) = ((dim.Sin+dim.S_height)*0.75+0.25*dim.Sout);

p(30,1) = line(4);

p(31,1) = line(3);

p(32,1) = line(1);


%Elements
    

t(1,:) = [2 19 3];
m(1) = dim.SM;

t(2,:) = [19 4 3];
m(2) = dim.SM;

t(3,:) = [19 18 4];
m(3) = dim.SM;

t(4,:) = [18 31 4];
m(4) = dim.SM;

t(5,:) = [18 32 31];
m(5) = dim.SM;

t(6,:) = [31 5 4];
m(6) = dim.SM;

t(7,:) = [31 30 5];
m(7) = dim.SM;

t(8,:) = [30 6 5];
m(8) = dim.SM;

t(9,:) = [30 29 6];
m(9) = dim.SM;

t(10,:) = [6 29 17];
m(10) = dim.SM;

t(11,:) = [17 7 6];
m(11) = dim.SM;

t(12,:) = [7 17 16];
m(12) = dim.SM;

t(13,:) = [29 28 17];
m(13) = dim.SM;

t(14,:) = [28 27 17];
m(14) = dim.SM;

t(15,:) = [27 16 17];
m(15) = dim.SM;

t(16,:) = [27 15 16];
m(16) = dim.SM;

t(17,:) = [27 26 15];
m(17) = dim.SM;

t(18,:) = [26 25 15];
m(18) = dim.SM;

t(19,:) = [7 16 15];
m(19) = dim.SM;

t(20,:) = [15 8 7];
m(20) = dim.SM;

t(21,:) = [25 8 15];
m(21) = dim.SM;

t(22,:) =[25 24 8];
m(22) = dim.SM;

t(23,:) = [24 9 8];
m(23) = dim.SM;

t(24,:) = [24 23 9];
m(24) = dim.SM;

t(25,:) = [23 10 9];
m(25) = dim.SM;

t(26,:) = [23 22 21];
m(26) = dim.SM;

t(27,:) = [23 21 10];
m(27) = dim.SM;

t(28,:) = [21 20 10];
m(28) = dim.SM;

t(29,:) = [20 11 10];
m(29) = dim.SM;

t(30,:) = [20 12 11];
m(30) = dim.SM;

t(31,:) = [1 2 3];
m(31) = 9999;

t(32,:) = [1 3 13];
m(32) = 9999;

t(33,:) = [1 13 11];
m(33) = 9999;

t(34,:) = [1 11 12];
m(34) = 9999;

t(35,:) = [13 3 4];
m(35) = 9999;

t(36,:) = [13 4 10];
m(36) = 9999;

t(37,:) = [10 11 13];
m(37) = 9999;

t(38,:) = [4 5 14];
m(38) = dim.SSM1;

t(39,:) = [4 14 10];
m(39) = dim.SSM1;

t(40,:) = [14 9 10];
m(40) = dim.SSM1;

t(41,:) = [5 6 14];
m(41) = dim.SSM1;

t(42,:) = [6 7 14];
m(42) = dim.SSM1;

t(43,:) = [7 8 14];
m(43) = dim.SSM1;

t(44,:) = [8 9 14];
m(44) = dim.SSM1;
  

FL = [32 31 30 29 28];
LL = [22 23 24 25 26];
cir = [28 27 26];
ag = [32 18 19 2 1 12 20 21 22];
end
    