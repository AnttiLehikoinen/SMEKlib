function msh = SPM_new2(msh,dim)
%This function creates mesh for surface mounted permanent magnet machine
%Ver. 0.4
%Work in progress
%Copyright (c) 2017-2018 Timo Davidsson / Aalto University

%%%%New testing%%%
angle1 = linspace(0,pi/2,3);
angle2 = linspace(0,pi/2,5);
angle3 = linspace(0,pi/2,8);
angle4 = linspace(0,pi/2,9);
angle5 = linspace(0,pi/2,11);
angle6 = linspace(0,pi/2,44);
angle7 = linspace(0,pi/2,66);
angle8 = linspace(0,pi/2,90);
angle9 = linspace(0,pi/2,90);

%%%%%%%%%%%%%%%%%%
angle61 = find(angle6 < pi/2/90*(dim.alfa-dim.beta/2));
angle62 = find(angle6 > pi/2/90*(dim.alfa+dim.beta/2));
angle71 = find(angle7 < pi/2/90*(dim.alfa-dim.beta/2));
angle72 = find(angle7 > pi/2/90*(dim.alfa+dim.beta/2));
angle81 = find(angle8 < pi/2/90*(dim.alfa-dim.beta/2));
angle82 = find(angle8 > pi/2/90*(dim.alfa+dim.beta/2));
angle91 = find(angle9 < pi/2/90*(dim.alfa-dim.beta/2));
angle92 = find(angle9 > pi/2/90*(dim.alfa+dim.beta/2));

a61 = angle61;
a62 = angle62;
a71 = angle71;
a72 = angle72;
a81 = angle81;
a82 = angle82;
a91 = angle91;
a92 = angle92;

angle61 = angle6(1:angle61(end));
angle62 = angle6(angle62(1):end);
angle63 = angle6(a61(end)+2:a62(1)-2);
angle71 = angle7(1:angle71(end));
angle72 = angle7(angle72(1):end);
angle73 = angle7(a71(end)+2:a72(1)-2);
angle81 = angle8(1:angle81(end));
angle82 = angle8(angle82(1):end);
angle83 = angle8(a81(end)+2:a82(1)-2);
angle91 = angle9(1:angle91(end));
angle92 = angle9(angle92(1):end);
angle93 = angle9(a91(end)+2:a92(1)-2);

if angle63(1) ~= angle93(1)
    angle63 = [angle93(1) angle63 angle93(end)];
    
end
if angle73(1) ~= angle93(1)
    if angle93(1) > angle73(1)
       angle73(1) = angle93(1);
       angle73(end) = angle93(end);
else
    angle73 = [angle93(1) angle73 angle93(end)];
    end
end

%Assigning radii data
r = dim.Rout-dim.Rin;
r1 = dim.Rin/2;
r2 = dim.Rin;
r3 = dim.Rin+0.25*r;
r4 = dim.Rin+0.5*r;
r5 = dim.Rin+0.75*r;
r6 = dim.Rout;
r7 = dim.Rout+0.45*dim.R_height0; 
r8 = dim.Rout+0.85*dim.R_height0; 
r9 = dim.Rout+dim.R_height0;

%%Making nodes
a1 = r1*[cos(angle1); sin(angle1)]';
a2 = r2*[cos(angle2); sin(angle2)]';
a3 = r3*[cos(angle3); sin(angle3)]';
a4 = r4*[cos(angle4); sin(angle4)]';
a5 = r5*[cos(angle5); sin(angle5)]';
a61 = r6*[cos(angle61); sin(angle61)]';
a71 = r7*[cos(angle71); sin(angle71)]';
a81 = r8*[cos(angle81); sin(angle81)]';
a91 = r9*[cos(angle91); sin(angle91)]';
a62 = r6*[cos(angle62); sin(angle62)]';
a72 = r7*[cos(angle72); sin(angle72)]';
a82 = r8*[cos(angle82); sin(angle82)]';
a92 = r9*[cos(angle92); sin(angle92)]';

a63 = r6*[cos(angle63); sin(angle63)]';
a73 = r7*[cos(angle73); sin(angle73)]';
a83 = r8*[cos(angle83); sin(angle83)]';
a93 = r9*[cos(angle93); sin(angle93)]';

%Assembling nodal coordinates
p = [0 0 ;a1;a2;a3;a4;a5;a61;a71;a81;a91;a62;a72;a82;a92;a63(1,:);a73(1,:);a83(1,:);a93(1,:);...
    a63(end,:);a73(end,:);a83(end,:);a93(end,:)];
koko1 = size(p,1) +1;
p = [p;a63(2:end-1,:);a73(2:end-1,:);a83(2:end-1,:);a93(2:end-1,:)];
koko2 = size(p,1);

%Creating triangulation
dt = delaunayTriangulation(p(:,1),p(:,2));
tr = triangulation(dt.ConnectivityList,dt.Points);
msh.dt = dt;
msh.p = tr.Points;
msh.t = tr.ConnectivityList;
msh.matel = zeros(size(msh.t,1),1);

%Find elements for shaft, core and magnet
[row,~] = find(msh.t <=4);         %Yeah, just want the row information 
shaft = unique(row);                            %Elements that are part of the shaft
[row,~] = find(msh.t >= 10 & msh.t <= 37);    
core = unique(row);                             %Elements that are part of the core       
[row,~] = find(msh.t >= koko1 & msh.t <= koko2);
magnet = unique(row);

%Make elements list         %NOTE TO SELF, CHANGE THE PARAMETERS TO dim
msh.matel(magnet) = 9; 
msh.matel(shaft) = dim.SFM;
msh.matel(core) = dim.RM;


%Make air gap vector
n = [a91;a92;a93];                    %Nodes on the outer line   
[~,row,~] = intersect(msh.p,n,'rows');
msh.n_ag_r = row;

end