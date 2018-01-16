%[data, dim] = parameter_example() 
%return

dim = struct();
data = struct();

%core materials
dim.SM = 2; 
dim.RM = 2;
dim.SWM = 3; %FIXME

dim.SSM1 = 0;
dim.SSM2 = 0;
dim.RSM = 0;


dim.Sout = 310.0E-03  / 2;
dim.Sin = 200.0E-03 / 2;
dim.Rout = dim.Sin - 0.8e-3;
dim.Rin = 70e-3 / 2;

dim.num = 4; dim.num2 = 4;
dim.Qs = 48; dim.Qr = 40;

dim.angleR =(2*pi)/dim.Qr:(2*pi)/dim.Qr:2*pi;
dim.angleS = (2*pi)/dim.Qs:(2*pi)/dim.Qs:2*pi;


%setting stator dimensions
dim.S = 4;
dim.hs_s = 23.90E-03;
dim.htt_s = 1e-3;
dim.hc_s = 17.5e-3;
dim.wso_s = 3.5e-3;
dim.wc1_c = 6.5e-3;
dim.wc2_s = 8.8e-3;

%setting rotor dimensions
dim.R = 5;
dim.htt_r = 0.7e-3;
dim.hc_r = 17.8e-3;
dim.hc_r2 = 16.1e-3;
dim.wc1_r = 6e-3;
dim.wc2_r = 2.5e-3;
dim.wc3_r = 5.85e-3;

%rotor slot opening
dim.RO = 2;
dim.SFM = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meshing

data = rotor5(data,dim); 
data.rotel = 1:data.index_t;
data = stator4(data,dim);
data.SC = data.SC';
data.RC = data.RC';

%