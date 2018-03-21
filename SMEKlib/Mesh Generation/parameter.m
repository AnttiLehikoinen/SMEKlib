%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameter-file of the Mesh Generator
%Ver. 0.85
%Copyright (c) 2017-2018 Timo Davidsson / Aalto University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [msh,dim] = parameter()

%%%DIMENSIONS%%%
    %%%%%Add wanted SLOT-INDEX %%%%%%
    R = 1;                              %Index for the shape of rotor slots
    S = 1;                              %Index for the shape of the stator slots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %MACHINE DIMENSIONS%
    stator_in = 200.000E-03;            %Inner diameter of the stator core
    stator_out = 310.000E-03;           %Outer diameter of the stator core
    rotor_in = 75.000E-03;              %Inner diameter of the rotor core
    rotor_out = 198.400E-03;            %Outer diameter of the rotor core
    Qs = 48;                            %Number of stator slots
    Qr = 40;                            %Number of rotor slots                    

    
    %%ROTOR SLOT DIMENSIONS%%
    R_height = 27.000E-03;     
    R_height1 = 3.000E-03;
    R_height2 = 10.000E-03;
    R_height3 = 25.000E-03;
    R_height4 = 1.000E-03;
    R_height5 = 27.000E-03;       
    R_width1 = 2.500E-03;        
    R_width2 = 5.000E-03;      
    R_width3 = 3.500E-03;        
    R_width4 = 2.500E-03;        
    R_width5 = 0.00;
    R_width6 = 0.00;
    R_height0 = 0.00;
    beta = 0.00;

    
    %%STATOR SLOT DIMENSIONS%%
    S_height = 25.500E-03;    
    S_height1 = 1.000E-03;
    S_height12 = 2.000E-03;
    S_height2 = 0.000;
    S_height3 = 22.500E-03;   
    S_height4 = 1.500E-03;
    S_height5 = 0.00;
    S_width1 = 3.500E-03;
    S_width2 = 5.500E-03;    
    S_width3 = 8.500E-03;   
    S_width4 = 0.00;
    S_width5 = 0.00;
    
       
  
%%%MATERIALS%%%                     

    SM = 4;                         %Material index for the stator core
    SWM = 0;                        %Material index for the stator slot wedges
    SSM1 = 0;                       %Material index for statorslot (inside), SSM1 is used if the slot consists only from 1 material
    SSM2 = 0;                       %Material index for statorslot (outside), so usually this can be left unedited

    RM = 4;                         %Material index for the rotor core
    RSM = 0;                        %Material index for the rotorslot
    RO = 0;                         %Material index for opening of rotorslot

    SFM = 1;                        %Material index for the shaft


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%DO NOT EDIT DATA BELOW%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mesh-generation parameter 
  % num = 1;    num2 = 1;       %calculate dimensions for whole machine
    num = 4;    num2 = 4;       %calculate dimensions for quater of the machine. (Use this for simulating!)
  % num = 2;    num2 = 2;       %calculate dimensions for 180-degree symmetry sector machines. (Use this for simulating!) !!!Remember to change symmetrysectors from setup file!!!!
  % num = Qr;   num2 = Qs;      %calculate dimensions for a single sector of the machine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Angle data for mesh%
angleR =(2*pi)/Qr:(2*pi)/Qr:2*pi;
angleS = (2*pi)/Qs:(2*pi)/Qs:2*pi;

%%%Creating mesh structure%%%
%p = nodal coordinates
%t = elements
%matel = element material-list
%rotel = rotor element-list
%FL = first layer coordinates in p, LL = last layer... (periodic boundary data)
%RC = rotor slot elements, SC = ...

msh = struct('p',[0 0],'t',[0 0 0],'matel',0,'rotel',0,'FL',0,'LL',0,'RC',0,'SC',0,'index_p',0,'index_t',0);

%%%Assigning dimension data%%%
dim = struct('R', R, 'S', S, 'Rin', rotor_in/2, 'Rout', rotor_out/2, 'Sin', stator_in/2, 'Sout', stator_out/2,...   
    'Qs', Qs, 'Qr', Qr, 'angleR', angleR, 'angleS', angleS, 'num', num, 'num2', num2,...
    'SM', SM, 'SWM', SWM, 'SSM1', SSM1, 'SSM2', SSM2, 'RM', RM, 'RSM', RSM, 'RO', RO, 'SFM', SFM,...
    'S_width1', S_width1, 'S_width2', S_width2, 'S_width3', S_width3, 'S_width4', S_width4, 'S_width5', S_width5,...
    'S_height', S_height, 'S_height1', S_height1, 'S_height12', S_height12, 'S_height2', S_height2, 'S_height3', S_height3, 'S_height4', S_height4, 'S_height5', S_height5,...
    'R_width1', R_width1, 'R_width2', R_width2, 'R_width3', R_width3, 'R_width4', R_width4, 'R_width5', R_width5, 'R_width6', R_width6, ...
    'R_height', R_height, 'R_height1', R_height1, 'R_height2', R_height2, 'R_height3', R_height3, 'R_height4', R_height4, 'R_height5', R_height5, 'R_height0', R_height0, 'beta', beta); 
    

end
