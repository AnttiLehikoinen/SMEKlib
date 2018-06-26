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
dim.S_height = 23.90E-03;
dim.S_height1 = 1e-3;
dim.S_height3 = 17.5e-3;
dim.S_width1 = 3.5e-3;
dim.S_width2 = 6.5e-3;
dim.S_width3 = 8.8e-3;

%setting rotor dimensions
dim.R = 5;
dim.R_height1 = 0.7e-3; %tooth tip height
dim.R_height3 = 17.8e-3; %distance from bottom bar top segment center to bottom bar bottom segment center
dim.R_height4 = 16.1e-3; %distance from airgap to bottom bar top segment center
dim.R_width3 = 6e-3; %top bar top segment diameter
dim.R_width4 = 2.5e-3; %area between top bar and bottom bar width
dim.R_width5 = 5.85e-3; %bottom bar top segment diameter
dim.R_width6 = 2.95e-3; %bottom bar bottom segment diameter

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