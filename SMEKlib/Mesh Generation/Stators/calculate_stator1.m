function [p,t,m,FL,LL,cir,ag] = calculate_stator1(dim)
%preallocate memory
p = zeros(37,2);
t = zeros(52,3);
m = zeros(52,1);

line = [dim.Sin dim.Sin+dim.S_height1 dim.Sin+dim.S_height-dim.S_height3 ...
        dim.Sin+dim.S_height-dim.S_height3/2 dim.Sin+dim.S_height ((dim.Sin+dim.S_height)*0.75+0.25*dim.Sout) dim.Sout];
       
        
    p(1,1) = line(1)*cos(dim.angleS(1)/2);
    p(1,2) = line(1)*sin(dim.angleS(1)/2);
    
    p(2,1) = line(1)*cos(dim.angleS(1)/2-asin(dim.S_width1/2/line(1)));
    p(2,2) = line(1)*sin(dim.angleS(1)/2-asin(dim.S_width1/2/line(1)));
    
    p(3,1) = line(2)*cos(dim.angleS(1)/2-asin(dim.S_width1/2/line(2)));
    p(3,2) = line(2)*sin(dim.angleS(1)/2-asin(dim.S_width1/2/line(2)));
    
    p(4,1) = line(3)*cos(dim.angleS(1)/2-atan((dim.S_width2/2)/line(3)));
    p(4,2) = line(3)*sin(dim.angleS(1)/2-atan((dim.S_width2/2)/line(3)));
    %
    %p(5,1) 0.5*... == dim.S_height4 to choose 
    p(5,1) = (0.5*line(3)+0.5*line(5))*cos(dim.angleS(1)/2-atan((dim.S_width3/4+dim.S_width2/4)/(0.5*line(3)+0.5*line(5))));
    p(5,2) = (0.5*line(3)+0.5*line(5))*sin(dim.angleS(1)/2-atan((dim.S_width3/4+dim.S_width2/4)/(0.5*line(3)+0.5*line(5))));
    
    p(6,1) = line(5)*cos(dim.angleS(1)/2-atan(dim.S_width3/2/line(5)));
    p(6,2) = line(5)*sin(dim.angleS(1)/2-atan(dim.S_width3/2/line(5)));
    
    p(7,1) = line(5)*cos(dim.angleS(1)/2);
    p(7,2) = line(5)*sin(dim.angleS(1)/2);
    
    p(8,1) = line(5)*cos(dim.angleS(1)/2+atan(dim.S_width3/2/line(5)));
    p(8,2) = line(5)*sin(dim.angleS(1)/2+atan(dim.S_width3/2/line(5)));
    %
    %0.5 or dim.S_height4
    p(9,1) = (0.5*line(3)+0.5*line(5))*cos(dim.angleS(1)/2+atan((dim.S_width3/4+dim.S_width2/4)/(0.5*line(3)+0.5*line(5))));
    p(9,2) = (0.5*line(3)+0.5*line(5))*sin(dim.angleS(1)/2+atan((dim.S_width3/4+dim.S_width2/4)/(0.5*line(3)+0.5*line(5))));
    
    p(10,1) = line(3)*cos(dim.angleS(1)/2+atan((dim.S_width2/2)/line(3)));
    p(10,2) = line(3)*sin(dim.angleS(1)/2+atan((dim.S_width2/2)/line(3)));
    
    p(11,1) = line(2)*cos(dim.angleS(1)/2+asin(dim.S_width1/2/line(2)));
    p(11,2) = line(2)*sin(dim.angleS(1)/2+asin(dim.S_width1/2/line(2)));
    
    p(12,1) = line(1)*cos(dim.angleS(1)/2+asin(dim.S_width1/2/line(1)));
    p(12,2) = line(1)*sin(dim.angleS(1)/2+asin(dim.S_width1/2/line(1)));
    
    p(13,1) = line(2)*cos(dim.angleS(1)/2);
    p(13,2) = line(2)*sin(dim.angleS(1)/2);
    
    p(14,1) = line(3)*cos(dim.angleS(1)/2);
    p(14,2) = line(3)*sin(dim.angleS(1)/2);
    
    %0.5 or dim.S_height4
    p(15,1) = (line(3)+0.5*(line(5)-line(3)))*cos(dim.angleS(1)/2);
    p(15,2) = (line(3)+0.5*(line(5)-line(3)))*sin(dim.angleS(1)/2);
    
    p(16,1) = line(6)*cos(dim.angleS(1)/2+dim.angleS(1)/4);
    p(16,2) = line(6)*sin(dim.angleS(1)/2+dim.angleS(1)/4);
    
    p(17,1) = line(6)*cos(dim.angleS(1)/2);
    p(17,2) = line(6)*sin(dim.angleS(1)/2);
    
    p(18,1) = line(6)*cos(dim.angleS(1)/2-dim.angleS(1)/4);
    p(18,2) = line(6)*sin(dim.angleS(1)/2-dim.angleS(1)/4);
    
    p(19,1) = line(3)*cos(4*dim.angleS(1)/5);
    p(19,2) = line(3)*sin(4*dim.angleS(1)/5);
    
    p(20,1) = line(3)*cos(dim.angleS(1)/5);
    p(20,2) = line(3)*sin(dim.angleS(1)/5);
    
    p(21,1) = line(1)*cos(dim.angleS(1)/5);
    p(21,2) = line(1)*sin(dim.angleS(1)/5);
    
    p(22,1) = line(1)*cos(dim.angleS(1)/2-atan(dim.S_width2/2/line(1)));
    p(22,2) = line(1)*sin(dim.angleS(1)/2-atan(dim.S_width2/2/line(1)));
    
    p(23,1) = line(1)*cos(dim.angleS(1)/2+atan(dim.S_width2/2/line(1)));
    p(23,2) = line(1)*sin(dim.angleS(1)/2+atan(dim.S_width2/2/line(1)));
    
    p(24,1) = line(1)*cos(4*dim.angleS(1)/5);
    p(24,2) = line(1)*sin(4*dim.angleS(1)/5);
    
    p(25,1) = line(1)*cos(dim.angleS(1));
    p(25,2) = line(1)*sin(dim.angleS(1));
    
    p(26,1) = line(3)*cos(dim.angleS(1));
    p(26,2) = line(3)*sin(dim.angleS(1));
    
    p(27,1) = line(4)*cos(dim.angleS(1));
    p(27,2) = line(4)*sin(dim.angleS(1));
    
    p(28,1) = line(5)*cos(dim.angleS(1));
    p(28,2) = line(5)*sin(dim.angleS(1));
    
    p(29,1) = line(6)*cos(dim.angleS(1));
    p(29,2) = line(6)*sin(dim.angleS(1));
    
    p(30,1) = line(7)*cos(dim.angleS(1));
    p(30,2) = line(7)*sin(dim.angleS(1));
    
    p(31,1) = line(7)*cos(dim.angleS(1)/2);
    p(31,2) = line(7)*sin(dim.angleS(1)/2);
    
    p(32,1) = line(7);
    
    p(33,1) = line(6);
    
    p(34,1) = line(5);
    
    p(35,1) = line(4);
    
    p(36,1) = line(3);
    
    p(37,1) = line(1);
    
    
    %Elements
    t(1,:) = [2 3 22];
    m(1) = dim.SM;
    
    t(2,:) = [22 4 3];
    m(2)  = dim.SM;
    
    t(3,:) = [22 21 20];
    m(3) = dim.SM;
    
    t(4,:) = [20 4 22];
    m(4) = dim.SM;
    
    t(5,:) = [21 37 20];
    m(5) = dim.SM;
    
    t(6,:) = [37 36 20];
    m(6) = dim.SM;
    
    t(7,:) = [20 5 4];
    m(7) = dim.SM;
    
    t(8,:) = [20 35 5];
    m(8) = dim.SM;
    
    t(9,:) = [36 35 20];
    m(9) = dim.SM;
    
    t(10,:) = [35 6 5];
    m(10) = dim.SM;
    
    t(11,:) = [35 34 6];
    m(11) = dim.SM;
    
    t(12,:) = [34 18 6];
    m(12) = dim.SM;
    
    t(13,:) = [34 33 18];
    m(13) = dim.SM;
    
    t(14,:) = [18 7 6];
    m(14) = dim.SM;
    
    t(15,:) = [18 17 7];
    m(15) = dim.SM;
    
    t(16,:) = [33 32 18];
    m(16) = dim.SM;
    
    t(17,:) = [32 31 18];
    m(17) = dim.SM;
    
    t(18,:) = [31 17 18];
    m(18) = dim.SM;
    
    t(19,:) = [31 16 17];
    m(19) = dim.SM;
    
    t(20,:) = [31 30 16];
    m(20) = dim.SM;
    
    t(21,:) = [30 29 16];
    m(21) = dim.SM;
    
    t(22,:) = [7 17 16];
    m(22) = dim.SM;
     
    t(23,:) = [7 16 8];
    m(23) = dim.SM;
    
    t(24,:) = [16 29 28];
    m(24) = dim.SM;
    
    t(25,:) = [16 28 8];
    m(25) = dim.SM;
    
    t(26,:) = [28 27 8];
    m(26) = dim.SM;
    
    t(27,:) = [27 9 8];
    m(27) = dim.SM;
    
    t(28,:) = [27 26 19];
    m(28) = dim.SM;
    
    t(29,:) = [27 19 9];
    m(29) = dim.SM;
    
    t(30,:) = [19 10 9];
    m(30) = dim.SM;
    
    t(31,:) = [26 25 19];
    m(31) = dim.SM;
    
    t(32,:) = [25 24 19];
    m(32) = dim.SM;
    
    t(33,:) = [24 23 19];
    m(33) = dim.SM;
    
    t(34,:) = [23 10 19];
    m(34) = dim.SM;
    
    t(35,:) = [23 11 10];
    m(35) = dim.SM;
    
    t(36,:) = [23 12 11];
    m(36) = dim.SM;
    
    t(37,:) = [1 2 3];
    m(37) = 9999;
    
    t(38,:) = [1 3 13];
    m(38) = 9999;
    
    t(39,:) = [1 13 11];
    m(39) = 9999;
    
    t(40,:) = [1 11 12];
    m(40) = 9999;
    
    t(41,:) = [13 3 4];
    m(41) = 9999;
    
    t(42,:) = [13 4 14];
    m(42) = 9999;
    
    t(43,:) = [13 14 10];
    m(43) = 9999;
    
    t(44,:) = [13 10 11];
    m(44) = 9999;
    
    t(45,:) = [14 4 5];
    m(45) = dim.SSM1;
    
    t(46,:) = [14 5 15];
    m(46) = dim.SSM1;
    
    t(47,:) = [14 15 9];
    m(47) = dim.SSM1;
    
    t(48,:) = [14 9 10];
    m(48) = dim.SSM1;
    
    t(49,:) = [5 6 7];
    m(49) = dim.SSM2;
    
    t(50,:) = [15 5 7];
    m(50) = dim.SSM2;
    
    t(51,:) = [15 7 9];
    m(51) = dim.SSM2;
    
    t(52,:) = [7 8 9];
    m(52) = dim.SSM2;

    
    FL = [37 36 35 34 33 32];
    LL = [25 26 27 28 29 30];
    cir = [32 31 30];
    ag = [37 21 22 2 1 12 28 24 25];
    
end