function [p,t,m,FL,LL,cir,ag] = calculate_stator8(dim)
    
p = zeros(29,2);
t = zeros(44,3);
m = zeros(44,1);

    %help variables
    line = [dim.Sin dim.Sin+dim.S_height1 dim.Sin+dim.S_height-dim.S_width3/2-dim.S_height3 ...
            dim.Sin+dim.S_height-dim.S_width3/2 dim.Sin+dim.S_height ((dim.Sin+dim.S_height)*0.75+0.25*dim.Sout) dim.Sout];
    c = sqrt(line(4)^2+(dim.S_width3/2)^2-2*(dim.S_width3/2)*line(4)*cosd(135));
    alfa = asin((dim.S_width3/2*sind(135))/c);
    
    % init the first sector
    p(1,1) = line(1)*cos(dim.angleS(1)/2);
    p(1,2) = line(1)*sin(dim.angleS(1)/2);
    
    p(2,1) = line(1)*cos(dim.angleS(1)/2-asin(dim.S_width1/2/line(1)));
    p(2,2) = line(1)*sin(dim.angleS(1)/2-asin(dim.S_width1/2/line(1)));
    
    p(3,1) = line(2)*cos(dim.angleS(1)/2-asin(dim.S_width1/2/line(2)));
    p(3,2) = line(2)*sin(dim.angleS(1)/2-asin(dim.S_width1/2/line(2)));
    
    p(4,1) = line(3)*cos(dim.angleS(1)/2-atan((dim.S_width2/2)/line(3)));
    p(4,2) = line(3)*sin(dim.angleS(1)/2-atan((dim.S_width2/2)/line(3)));
    
    p(5,1) = (line(3)+dim.S_height4*(line(4)-line(3)))*cos(dim.angleS(1)/2-atan(((1-dim.S_height4)*dim.S_width2/2+dim.S_height4*dim.S_width3/2)/(line(3)+dim.S_height4*(line(4)-line(3)))));
    p(5,2) = (line(3)+dim.S_height4*(line(4)-line(3)))*sin(dim.angleS(1)/2-atan(((1-dim.S_height4)*dim.S_width2/2+dim.S_height4*dim.S_width3/2)/(line(3)+dim.S_height4*(line(4)-line(3)))));
    
    p(6,1) = line(4)*cos(dim.angleS(1)/2-atan(dim.S_width3/2/line(4)));
    p(6,2) = line(4)*sin(dim.angleS(1)/2-atan(dim.S_width3/2/line(4)));
    
    p(7,1) = c*cos(dim.angleS(1)/2-alfa);
    p(7,2) = c*sin(dim.angleS(1)/2-alfa);
    
    p(8,1) = line(5)*cos(dim.angleS(1)/2);
    p(8,2) = line(5)*sin(dim.angleS(1)/2);
    
    p(9,1) = c*cos(dim.angleS(1)/2+alfa);
    p(9,2) = c*sin(dim.angleS(1)/2+alfa);
    
    p(10,1) = line(4)*cos(dim.angleS(1)/2+atan(dim.S_width3/2/line(4)));
    p(10,2) = line(4)*sin(dim.angleS(1)/2+atan(dim.S_width3/2/line(4)));
    
    p(11,1) = (line(3)+dim.S_height4*(line(4)-line(3)))*cos(dim.angleS(1)/2+atan(((1-dim.S_height4)*dim.S_width2/2+dim.S_height4*dim.S_width3/2)/(line(3)+dim.S_height4*(line(4)-line(3)))));
    p(11,2) = (line(3)+dim.S_height4*(line(4)-line(3)))*sin(dim.angleS(1)/2+atan(((1-dim.S_height4)*dim.S_width2/2+dim.S_height4*dim.S_width3/2)/(line(3)+dim.S_height4*(line(4)-line(3)))));
    
    p(12,1) = line(3)*cos(dim.angleS(1)/2+atan((dim.S_width2/2)/line(3)));
    p(12,2) = line(3)*sin(dim.angleS(1)/2+atan((dim.S_width2/2)/line(3)));
    
    p(13,1) = line(2)*cos(dim.angleS(1)/2+asin(dim.S_width1/2/line(2)));
    p(13,2) = line(2)*sin(dim.angleS(1)/2+asin(dim.S_width1/2/line(2)));
    
    p(14,1) = line(1)*cos(dim.angleS(1)/2+asin(dim.S_width1/2/line(1)));
    p(14,2) = line(1)*sin(dim.angleS(1)/2+asin(dim.S_width1/2/line(1)));
    
    p(15,1) = line(2)*cos(dim.angleS(1)/2);
    p(15,2) = line(2)*sin(dim.angleS(1)/2);
    
    p(16,1) = line(3)*cos(dim.angleS(1)/2);
    p(16,2) = line(3)*sin(dim.angleS(1)/2);
    
    p(17,1) = line(4)*cos(dim.angleS(1)/2);
    p(17,2) = line(4)*sin(dim.angleS(1)/2);
    
    p(18,1) = line(6)*cos(dim.angleS(1)/2+dim.angleS(1)/4);
    p(18,2) = line(6)*sin(dim.angleS(1)/2+dim.angleS(1)/4);
    
    p(19,1) = line(6)*cos(dim.angleS(1)/2);
    p(19,2) = line(6)*sin(dim.angleS(1)/2);
    
    p(20,1) = line(6)*cos(dim.angleS(1)/2-dim.angleS(1)/4);
    p(20,2) = line(6)*sin(dim.angleS(1)/2-dim.angleS(1)/4);
    
    p(21,1) = line(1)*cos(dim.angleS(1));
    p(21,2) = line(1)*sin(dim.angleS(1));
    
    p(22,1) = line(4)*cos(dim.angleS(1));
    p(22,2) = line(4)*sin(dim.angleS(1));
    
    p(23,1) = line(6)*cos(dim.angleS(1));
    p(23,2) = line(6)*sin(dim.angleS(1));
    
    p(24,1) = line(7)*cos(dim.angleS(1));
    p(24,2) = line(7)*sin(dim.angleS(1));
    
    p(25,1) = line(7)*cos(dim.angleS(1)/2);
    p(25,2) = line(7)*sin(dim.angleS(1)/2);
    
    p(26,1) = line(7);
    
    p(27,1) = line(6);
    
    p(28,1) = line(4);
    
    p(29,1) = line(1);
    
    
    t(1,:) = [29 3 2];
m(1) = dim.SM;

t(2,:) = [29 4 3];
m(2) = dim.SM;

t(3,:) = [29 5 4];
m(3) = dim.SM;

t(4,:) = [29 28 5];
m(4) = dim.SM;

t(5,:) = [28 6 5];
m(5) = dim.SM;

t(6,:) = [2 3 1];
m(6) = 9999;

t(7,:) = [3 15 1];
m(7) = 9999;

t(8,:) = [1 15 13];
m(8) = 9999;

t(9,:) = [1 13 14];
m(9) = 9999;

t(10,:) = [3 4 15];
m(10) = 9999;

t(11,:) = [4 16 15];
m(11) = 9999;

t(12,:) = [15 12 13];
m(12) = 9999;

t(13,:) = [15 16 12];
m(13) = 9999;

t(14,:) = [4 5 16];
m(14) = dim.SSM1;

t(15,:) = [5 11 16];
m(15) = dim.SSM1;

t(16,:) = [16 11 12];
m(16) = dim.SSM1;

t(17,:) = [5 6 17];

t(18,:) = [5 17 11];

t(19,:) = [17 10 11];

t(20,:) = [6 7 17];

t(21,:) = [7 8 17];

t(22,:) = [8 9 17];

t(23,:) = [9 10 17];

t(24,:) = [14 13 21];
m(24) = dim.SM;

t(25,:) = [13 12 21];
m(25) = dim.SM;

t(26,:) = [12 11 21];
m(26) = dim.SM;

t(27,:) = [11 22 21];
m(27) = dim.SM;

t(28,:) = [11 10 22];
m(28) = dim.SM;

t(29,:) = [28 7 6];
m(29) = dim.SM;

t(30,:) = [28 27 20];
m(30) = dim.SM;

t(31,:) = [28 20 7];
m(31) = dim.SM;

t(32,:) = [7 20 8];
m(32) = dim.SM;

t(33,:) = [8 20 19];
m(33) = dim.SM;

t(34,:) = [8 19 18];
m(34) = dim.SM;

t(35,:) = [8 18 9];
m(35) = dim.SM;

t(36,:) = [22 10 9];
m(36) = dim.SM;

t(37,:) = [22 9 18];
m(37) = dim.SM;

t(38,:) = [22 18 23];
m(38) = dim.SM;

t(39,:) = [27 26 20];
m(39) = dim.SM;

t(40,:) = [26 25 20];
m(40) = dim.SM;

t(41,:) = [20 25 19];
m(41) =4;

t(42,:) = [19 25 18];
m(42) = dim.SM;

t(43,:) = [18 25 24];
m(43) = dim.SM;

t(44,:) = [18 24 23];
m(44) = dim.SM;
    
      
FL = [29 28 27 26];
LL = [21 22 23 24];
cir = [26 25 24];
ag = [29 2 1 14 21];
end