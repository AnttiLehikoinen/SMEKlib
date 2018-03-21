function msh = surf_perma_mgt(msh,dim)
%Creates mesh for surface-mounted permanent magnet machine
%v.0.4
%Copyright (c) 2017 Timo Davidsson / Aalto University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim.alfa = 45;  

%Calculate nodes and elements
msh = SPM_new2(msh,dim); 

%Finalize mesh
msh.index_p = size(msh.p,1);
msh.index_t = size(msh.t,1);
end