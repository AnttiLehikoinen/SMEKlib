function [p,t,m,FL,LL,cir,ag] = calculate_stator5(dim)

%preallocate memory
p = zeros(32,2);
t = zeros(44,3);
m = zeros(44,1);

L = dim.Sin/2 + (dim.Sin+dim.S_height-dim.S_height3)/2;
line = [dim.Sin L dim.Sin+dim.S_height-dim.S_height3 dim.Sin+dim.S_height-dim.S_height3/2 ...
        dim.Sin+dim.S_height (dim.Sin+dim.S_height)/2+dim.Sout/2 dim.Sout];
    
    
%Init first layer%

    p(1,1) = line(1)*cos(dim.angleS(1)/2);
    p(1,2) = line(1)*sin(dim.angleS(1)/2);
    
    p(2,1) = line(1)*cos(dim.angleS(1)/2-atan(dim.S_width1/2/line(1)));
    p(2,2) = line(1)*sin(dim.angleS(1)/2-atan(dim.S_width1/2/line(1)));

    p(3,1) = line(2)*cos(dim.angleS(1)/2-atan(dim.S_width1/2/line(2)));
    p(3,2) = line(2)*sin(dim.angleS(1)/2-atan(dim.S_width1/2/line(2)));
    
    p(4,1) = line(3)*cos(dim.angleS(1)/2-atan(dim.S_width1/2/line(3)));
    p(4,2) = line(3)*sin(dim.angleS(1)/2-atan(dim.S_width1/2/line(3)));
    
    p(5,1) = line(4)*cos(dim.angleS(1)/2-atan(dim.S_width1/2/line(4)));
    p(5,2) = line(4)*sin(dim.angleS(1)/2-atan(dim.S_width1/2/line(4)));
    
    p(6,1) = line(5)*cos(dim.angleS(1)/2-atan(dim.S_width1/2/line(5)));
    p(6,2) = line(5)*sin(dim.angleS(1)/2-atan(dim.S_width1/2/line(5)));
    
    p(7,1) = line(5)*cos(dim.angleS(1)/2);
    p(7,2) = line(5)*sin(dim.angleS(1)/2);
    
    p(8,1) = line(5)*cos(dim.angleS(1)/2+atan(dim.S_width1/2/line(5)));
    p(8,2) = line(5)*sin(dim.angleS(1)/2+atan(dim.S_width1/2/line(5)));

    p(9,1) = line(4)*cos(dim.angleS(1)/2+atan(dim.S_width1/2/line(4)));
    p(9,2) = line(4)*sin(dim.angleS(1)/2+atan(dim.S_width1/2/line(4)));
    
    p(10,1) = line(3)*cos(dim.angleS(1)/2+atan(dim.S_width1/2/line(3)));
    p(10,2) = line(3)*sin(dim.angleS(1)/2+atan(dim.S_width1/2/line(3)));
    
    p(11,1) = line(2)*cos(dim.angleS(1)/2+atan(dim.S_width1/2/line(2)));
    p(11,2) = line(2)*sin(dim.angleS(1)/2+atan(dim.S_width1/2/line(2)));
    
    p(12,1) = line(1)*cos(dim.angleS(1)/2+atan(dim.S_width1/2/line(1)));
    p(12,2) = line(1)*sin(dim.angleS(1)/2+atan(dim.S_width1/2/line(1)));
    
    p(13,1) = line(2)*cos(dim.angleS(1)/2);
    p(13,2) = line(2)*sin(dim.angleS(1)/2);
    
    p(14,1) = line(3)*cos(dim.angleS(1)/2);
    p(14,2) = line(3)*sin(dim.angleS(1)/2);
    
    p(15,1) = line(4)*cos(dim.angleS(1)/2);
    p(15,2) = line(4)*sin(dim.angleS(1)/2);
    
    p(16,1) = (line(5)*0.75+0.25*dim.Sout)*cos(dim.angleS(1)/2+dim.angleS(1)/4);
    p(16,2) = (line(5)*0.75+0.25*dim.Sout)*sin(dim.angleS(1)/2+dim.angleS(1)/4); %line(6)
    
    p(17,1) = (line(5)*0.75+0.25*dim.Sout)*cos(dim.angleS(1)/2-dim.angleS(1)/4);
    p(17,2) = (line(5)*0.75+0.25*dim.Sout)*sin(dim.angleS(1)/2-dim.angleS(1)/4);
    
    p(18,1) = line(1)*cos(dim.angleS(1));
    p(18,2) = line(1)*sin(dim.angleS(1));
    
    p(19,1) = line(2)*cos(dim.angleS(1));
    p(19,2) = line(2)*sin(dim.angleS(1));
    
    p(20,1) = line(3)*cos(dim.angleS(1));
    p(20,2) = line(3)*sin(dim.angleS(1));
    
    p(21,1) = line(4)*cos(dim.angleS(1));
    p(21,2) = line(4)*sin(dim.angleS(1));
    
    p(22,1) = line(5)*cos(dim.angleS(1));
    p(22,2) = line(5)*sin(dim.angleS(1));
    
    p(23,1) = (line(5)*0.75+0.25*dim.Sout)*cos(dim.angleS(1));
    p(23,2) = (line(5)*0.75+0.25*dim.Sout)*sin(dim.angleS(1));
    
    p(24,1) = line(7)*cos(dim.angleS(1));
    p(24,2) = line(7)*sin(dim.angleS(1));

    p(25,1) = line(7)*cos(dim.angleS(1)/2);
    p(25,2) = line(7)*sin(dim.angleS(1)/2);
    
    p(26,1) = line(7);
    
    p(27,1) = (line(5)*0.75+0.25*dim.Sout);
    
    p(28,1) = line(5); 
   
    p(29,1) = line(4);
    
    p(30,1) = line(3);
    
    p(31,1) = line(2); 
    
    p(32,1) = line(1);
    
    
    t(1,:) = [32 31 2];
m(1) = dim.SM;

t(2,:) = [31 3 2];
m(2) = dim.SM;

t(3,:) = [31 4 3];
m(3) = dim.SM;

t(4,:) = [31 30 4];
m(4) = dim.SM;

t(5,:) = [30 29 5];
m(5) = dim.SM;

t(6,:) = [30 5 4];
m(6) = dim.SM;

t(7,:) = [29 28 5];
m(7) = dim.SM;

t(8,:) = [28 6 5];
m(8) = dim.SM;

t(9,:) = [28 27 17];
m(9) = dim.SM;

t(10,:) = [28 17 6];
m(10) = dim.SM;

t(11,:) = [6 17 7];
m(11) = dim.SM;

t(12,:) = [17 16 7];
m(12) = dim.SM;

t(13,:) = [7 16 8];
m(13) = dim.SM;

t(14,:) = [16 23 22];
m(14) = dim.SM;

t(15,:) = [8 16 22];
m(15) = dim.SM;

t(16,:) = [8 22 9];
m(16) = dim.SM;

t(17,:) = [22 21 9];
m(17) = dim.SM;

t(18,:) = [9 21 20];
m(18)  = dim.SM;

t(19,:) = [9 20 10];
m(19) = dim.SM;

t(20,:) = [10 20 19];
m(20) = dim.SM;

t(21,:) = [10 19 11];
m(21) = dim.SM;

t(22,:) = [19 11 12];
m(22) = dim.SM;

t(23,:) = [19 18 12];
m(23) = dim.SM;

t(24,:) = [23 16 24];
m(24) = dim.SM;

t(25,:) = [16 25 24]; 
m(25) = dim.SM;

t(26,:) = [16 17 25];
m(26) = dim.SM;

t(27,:) = [17 26 25];
m(27) = dim.SM;

t(28,:) = [17 27 26];
m(28) = dim.SM;

t(29,:) = [12 1 11];
m(29) = 9999;

t(30,:) = [1 13 11];
m(30) = 9999;

t(31,:) = [1 2 3];
m(31) = 9999;

t(32,:) = [1 3 13];
m(32) = 9999;

t(33,:) = [10 11 14];
m(33) = 9999;

t(34,:) = [11 13 14];
m(34) = 9999;

t(35,:) = [14 13 3];
m(35) = 9999;

t(36,:) = [3 4 14];
m(36) = 9999;

t(37,:) = [14 4 15];
m(37) = dim.SSM1;

t(38,:) = [4 5 15];
m(38) = dim.SSM1;

t(39,:) = [10 14 15];
m(39) = dim.SSM1;

t(40,:) = [10 15 9];
m(40) = dim.SSM1;

t(41,:) = [8 9 15];
m(41) = dim.SSM2;

t(42,:) = [7 8 15];
m(42) = dim.SSM2;

t(43,:) = [6 7 15];
m(43) = dim.SSM2;

t(44,:) = [5 6 15];
m(44) = dim.SSM2;
    
FL = [32 31 30 29 28 27 26];
LL = [18 19 20 21 22 23 24];
cir = [26 25 24];
ag = [32 2 1 12 18];

end