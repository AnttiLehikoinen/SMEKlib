function [p,t,m,FL,LL,cir,ag] = calculate_stator2(dim)
    
%preallocate memory
p = zeros(44,2);
t = zeros(66,3);
m = zeros(66,1);

    line = [dim.Sin dim.Sin+dim.S_height1 dim.Sin+dim.S_height-dim.S_height3-dim.S_height4 ...
            dim.Sin+dim.S_height-dim.S_height4 dim.Sin+dim.S_height ((dim.Sin+dim.S_height)*0.75+0.25*dim.Sout) dim.Sout];
       
    c =0.4*line(4)+0.6*line(5);
    alfa = atan((0.7*dim.S_width3/2+0.3*dim.S_width2/2)/c);

      
    p(1,1) = line(1)*cos(dim.angleS(1)/2);
    p(1,2) = line(1)*sin(dim.angleS(1)/2);
    
    p(2,1) = line(1)*cos(dim.angleS(1)/2-asin(dim.S_width1/2/line(1)));
    p(2,2) = line(1)*sin(dim.angleS(1)/2-asin(dim.S_width1/2/line(1)));
    
    p(3,1) = line(2)*cos(dim.angleS(1)/2-asin(dim.S_width1/2/line(2)));
    p(3,2) = line(2)*sin(dim.angleS(1)/2-asin(dim.S_width1/2/line(2)));
    
    p(4,1) = line(3)*cos(dim.angleS(1)/2-atan((dim.S_width2/2)/line(3)));
    p(4,2) = line(3)*sin(dim.angleS(1)/2-atan((dim.S_width2/2)/line(3)));
    
    p(5,1) = (line(4)+ line(3))/2*cos(dim.angleS(1)/2-atan((dim.S_width2+dim.S_width3)/(line(4)+ line(3))/2));
    p(5,2) = (line(4)+ line(3))/2*sin(dim.angleS(1)/2-atan((dim.S_width2+dim.S_width3)/(line(4)+ line(3))/2));
    
    p(6,1) = line(4)*cos(dim.angleS(1)/2-atan(dim.S_width3/2/line(4)));
    p(6,2) = line(4)*sin(dim.angleS(1)/2-atan(dim.S_width3/2/line(4)));
    
    p(7,1) = c*cos(dim.angleS(1)/2-alfa);%%%%%%%%%%%%%%%%%%%%%%%%%%c alfa
    p(7,2) = c*sin(dim.angleS(1)/2-alfa);
    
    p(8,1) = line(5)*cos(dim.angleS(1)/2-atan(dim.S_width2/2/line(5)));
    p(8,2) = line(5)*sin(dim.angleS(1)/2-atan(dim.S_width2/2/line(5)));
    
    p(9,1) = line(5)*cos(dim.angleS(1)/2);
    p(9,2) = line(5)*sin(dim.angleS(1)/2);
    
    p(10,1) = line(5)*cos(dim.angleS(1)/2+atan(dim.S_width2/2/line(5)));
    p(10,2) = line(5)*sin(dim.angleS(1)/2+atan(dim.S_width2/2/line(5)));
    
    p(11,1) = c*cos(dim.angleS(1)/2+alfa);%%%%%%%%%%%%%%%%%%%%%%%%%%c alfa
    p(11,2) = c*sin(dim.angleS(1)/2+alfa);
    
    p(12,1) = line(4)*cos(dim.angleS(1)/2+atan(dim.S_width3/2/line(4)));
    p(12,2) = line(4)*sin(dim.angleS(1)/2+atan(dim.S_width3/2/line(4)));
    
    p(13,1) = (line(4)+ line(3))/2*cos(dim.angleS(1)/2+atan((dim.S_width2+dim.S_width3)/(line(4)+ line(3))/2));
    p(13,2) = (line(4)+ line(3))/2*sin(dim.angleS(1)/2+atan((dim.S_width2+dim.S_width3)/(line(4)+ line(3))/2));
    
    p(14,1) = line(3)*cos(dim.angleS(1)/2+atan((dim.S_width2/2)/line(3)));
    p(14,2) = line(3)*sin(dim.angleS(1)/2+atan((dim.S_width2/2)/line(3)));
    
    p(15,1) = line(2)*cos(dim.angleS(1)/2+asin(dim.S_width1/2/line(2)));
    p(15,2) = line(2)*sin(dim.angleS(1)/2+asin(dim.S_width1/2/line(2)));
    
    p(16,1) = line(1)*cos(dim.angleS(1)/2+asin(dim.S_width1/2/line(1)));
    p(16,2) = line(1)*sin(dim.angleS(1)/2+asin(dim.S_width1/2/line(1)));
    
    p(17,1) = line(2)*cos(dim.angleS(1)/2);
    p(17,2) = line(2)*sin(dim.angleS(1)/2);
    
    p(18,1) = line(3)*cos(dim.angleS(1)/2);
    p(18,2) = line(3)*sin(dim.angleS(1)/2);
    
    p(19,1) = (line(4)+ line(3))/2*cos(dim.angleS(1)/2);
    p(19,2) = (line(4)+line(3))/2*sin(dim.angleS(1)/2);
    
    p(20,1) = line(4)*cos(dim.angleS(1)/2+atan(dim.S_width1/2/line(4)));
    p(20,2) = line(4)*sin(dim.angleS(1)/2+atan(dim.S_width1/2/line(4)));
    
    p(21,1) = line(4)*cos(dim.angleS(1)/2);
    p(21,2) = line(4)*sin(dim.angleS(1)/2);
    
    p(22,1) = line(4)*cos(dim.angleS(1)/2-atan(dim.S_width1/2/line(4)));
    p(22,2) = line(4)*sin(dim.angleS(1)/2-atan(dim.S_width1/2/line(4)));
    
    p(23,1) = line(6)*cos(3/4*dim.angleS(1));
    p(23,2) = line(6)*sin(3/4*dim.angleS(1));
    
    p(24,1) = line(6)*cos(dim.angleS(1)/2);
    p(24,2) = line(6)*sin(dim.angleS(1)/2);
    
    p(25,1) = line(6)*cos(1/4*dim.angleS(1));
    p(25,2) = line(6)*sin(1/4*dim.angleS(1));
    
    p(26,1) =line(3)*cos(4*dim.angleS(1)/5);
    p(26,2) =line(3)*sin(4*dim.angleS(1)/5);
    
    p(27,1) =line(3)*cos(dim.angleS(1)/5);
    p(27,2) =line(3)*sin(dim.angleS(1)/5);
    
    p(28,1) =line(1)*cos(dim.angleS(1)/5);
    p(28,2) =line(1)*sin(dim.angleS(1)/5);
    
    p(29,1) =line(1)*cos(dim.angleS(1)/2-atan(dim.S_width2/2/line(1)));
    p(29,2) =line(1)*sin(dim.angleS(1)/2-atan(dim.S_width2/2/line(1)));
    
    p(30,1) =line(1)*cos(dim.angleS(1)/2+atan(dim.S_width2/2/line(1)));
    p(30,2) =line(1)*sin(dim.angleS(1)/2+atan(dim.S_width2/2/line(1)));
    
    p(31,1) =line(1)*cos(4*dim.angleS(1)/5);
    p(31,2) =line(1)*sin(4*dim.angleS(1)/5);
    
    p(32,1) = line(1)*cos(dim.angleS(1));
    p(32,2) = line(1)*sin(dim.angleS(1));
    
    p(33,1) = line(3)*cos(dim.angleS(1));
    p(33,2) = line(3)*sin(dim.angleS(1));
    
    p(34,1) = (line(3)+line(4))/2*cos(dim.angleS(1));
    p(34,2) = (line(3)+line(4))/2*sin(dim.angleS(1));
    
    p(35,1) = line(5)*cos(dim.angleS(1));
    p(35,2) = line(5)*sin(dim.angleS(1));
    
    p(36,1) = line(6)*cos(dim.angleS(1));
    p(36,2) = line(6)*sin(dim.angleS(1));
    
    p(37,1) = line(7)*cos(dim.angleS(1));
    p(37,2) = line(7)*sin(dim.angleS(1));
    
    p(38,1) = line(7)*cos(dim.angleS(1)/2);
    p(38,2) = line(7)*sin(dim.angleS(1)/2);
    
    p(39,1) = line(7);
    
    p(40,1) = line(6);
    
    p(41,1) = line(5);
    
    p(42,1) = (line(3)+line(4))/2;
    
    p(43,1) = line(3);
    
    p(44,1) = line(1);
    
    
    
    
    
    
    %Elements
    
    t(1,:) = [2 29 3];
    m(1) = dim.SM;
    
    t(2,:) = [29 3 4];
    m(2) = dim.SM;
    
    t(3,:) = [29 27 4];
    m(3) = dim.SM;
    
    t(4,:) = [29 28 27];
    m(4) = dim.SM;
    
    t(5,:) = [28 44 27];
    m(5) = dim.SM;
    
    t(6,:) = [44 43 27];
    m(6) = dim.SM;
    
    t(7,:) = [4 27 5];
    m(7) = dim.SM;
    
    t(8,:) = [27 42 5];
    m(8) = dim.SM;
    
    t(9,:) = [27 43 42];
    m(9) = dim.SM;
    
    t(10,:) = [42 6 5];
    m(10) = dim.SM;
    
    t(11,:) = [42 41 6];
    m(11) = dim.SM;
    
    t(12,:) = [41 7 6];
    m(12) = dim.SM;
    
    t(13,:) = [41 40 7];
    m(13) = dim.SM;
    
    t(14,:) = [40 25 7];
    m(14) = dim.SM;
    
    t(15,:) = [25 8 7];
    m(15) = dim.SM;
    
    t(16,:) = [25 24 8];
    m(16) = dim.SM;
    
    t(17,:) = [24 9 8];
    m(17) = dim.SM;
    
    t(18,:) = [40 39 25];
    m(18) = dim.SM;
    
    t(19,:) = [39 38 25];
    m(19) = dim.SM;
    
    t(20,:) = [38 24 25];
    m(20) = dim.SM;
    
    t(21,:) = [38 23 24];
    m(21) = dim.SM;
    
    t(22,:) = [38 37 23];
    m(22) = dim.SM;
    
    t(23,:) = [37 36 23];
    m(23) = dim.SM;
    
    t(24,:) = [36 11 23];
    m(24) = dim.SM;
    
    t(25,:) = [23 11 10];
    m(25) = dim.SM;
    
    t(26,:) = [10 24 23];
    m(26) = dim.SM;
    
    t(27,:) = [10 9 24];
    m(27) = dim.SM;
    
    t(28,:) = [36 35 11];
    m(28) = dim.SM;
    
    t(29,:) = [35 12 11];
    m(29) = dim.SM;
    
    t(30,:) = [35 34 12];
    m(30) = dim.SM;
    
    t(31,:) = [34 13 12];
    m(31) = dim.SM;
    
    t(32,:) = [34 26 13];
    m(32) = dim.SM;
    
    t(33,:) = [26 14 13];
    m(33) = dim.SM;
    
    t(34,:) = [34 33 26];
    m(34) = dim.SM;
    
    t(35,:) = [33 32 26];
    m(35) = dim.SM;
    
    t(36,:) = [32 31 26];
    m(36) = dim.SM;
    
    t(37,:) = [31 30 26];
    m(37) = dim.SM;
    
    t(38,:) = [30 14 26];
    m(38) = dim.SM;
    
    t(39,:) = [30 15 14];
    m(39) = dim.SM;
    
    t(40,:) = [30 16 15];
    m(40) = dim.SM;
    
    t(41,:) = [1 2 3];
    m(41) = 9999;
    
    t(42,:) = [1 3 17];
    m(42) = 9999;
    
    t(43,:) = [1 17 15];
    m(43) = 9999;
        
    t(44,:) = [1 15 16];
    m(44) = 9999;
    
    t(45,:) = [3 4 18];
    m(45) = 9999;
    
    t(46,:) = [3 18 17];
    m(46) = 9999;
    
    t(47,:) = [17 18 15];
    m(47) = 9999;
    
    t(48,:) = [14 15 18];
    m(48) = 9999;
    
    t(49,:) = [18 4 5];
    m(49) = dim.SSM1;
    
    t(50,:) = [19 18 5];
    m(50) = dim.SSM1;
    
    t(51,:) = [13 18 19];
    m(51) = dim.SSM1;
    
    t(52,:) = [13 14 18];
    m(52) = dim.SSM1;
    
    t(53,:) = [5 6 22];
    m(53) = dim.SSM2;
    
    t(54,:) = [5 22 21];
    m(54) = dim.SSM2;
    
    t(55,:) = [5 21 19];
    m(55) = dim.SSM2;
    
    t(56,:) = [19 21 13];
    m(56) = dim.SSM2;
    
    t(57,:) = [21 20 13];
    m(57) = dim.SSM2; 
    
    t(58,:) = [20 12 13];
    m(58) = dim.SSM2;
    
    t(59,:) = [6 7 22];
    m(59) = dim.SSM2;
    
    t(60,:) = [7 8 22];
    m(60) = dim.SSM2; 
    
    t(61,:) = [8 21 22];
    m(61) = dim.SSM2;
    
    t(62,:) = [8 9 21];
    m(62) = dim.SSM2;
    
    t(63,:) = [9 10 21];
    m(63) = dim.SSM2;
    
    t(64,:) =[10 20 21];
    m(64) = dim.SSM2;
    
    t(65,:) = [10 11 20];
    m(65) = dim.SSM2;
    
    t(66,:) = [11 12 20];
    m(66) = dim.SSM2;


    
FL = [44 43 42 41 40 39];
LL = [32 33 34 35 36 37];
cir = [39 38 37];
ag = [44 28 29 2 1 16 30 31 32];
end    
    
    
    