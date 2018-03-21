%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main function of the Mesh Generator
%v0.9
%Copyright (c) 2017 Timo Davidsson / Aalto University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%init structures
[msh,dim] = parameter();  
tic 

%Calculate mesh according to wanted slot-types
switch dim.R
    case 1
        msh = rotor1(msh,dim);  
    case 2
        msh = rotor2(msh,dim);  
    case 3
        msh = rotor3(msh,dim);  
    case 4
        msh = rotor4(msh,dim);  
    case 5
        msh = rotor5(msh,dim);  
    case 6
        msh = rotor6(msh,dim);
    case 7
        msh = rotor7(msh,dim);  
    case 8
        msh = rotor8(msh,dim);  
    case 9
        msh = rotor9(msh,dim);
    case 20
        msh = rotor20(msh,dim);    
    case 21
        msh = rotor21(msh,dim);   
    case 22
        msh = rotor22(msh,dim);
    case 99
        msh = surf_perma_mgt(msh,dim);  %Surface-mounted permanentmagnet machine 
    otherwise
        error('Check your rotorslot index!');
end

msh.rotel = 1:msh.index_t;

switch dim.S
    case 1
        msh = stator1(msh,dim); 
    case 2
        msh = stator2(msh,dim); 
    case 3
        msh = stator3(msh,dim); 
    case 4
        msh = stator4(msh,dim); 
    case 5
        msh = stator5(msh,dim); 
    case 6
        msh = stator6(msh,dim); 
    case 7
        msh = stator7(msh,dim); 
    case 8
        msh = stator8(msh,dim); 
    case 9
        msh = stator9(msh,dim); 
    case 10
        msh = stator10(msh,dim);
    otherwise
        error('Check your statorslot index!');
end
toc


%%%%Plotting stuff for testing purposes%%%%
%figure(1)
%plot(msh.p(:,1),msh.p(:,2),'.');   %testing nodes
%axis equal
%grid on
%figure(2)
%triplot(msh.t,msh.p(:,1),msh.p(:,2));  %testing triangulation
%axis equal

%%%%Functions used for datapoints documentation%%%%
%hold on;
%number_nodes(msh.p');
%number_elements(msh.p', msh.t');
%axis off
%xlim([0 0.1]);
%ylim([-0.005 0.021]);

%%%%Just changing name of the struct for SMEKlib%%%%
data = msh;

