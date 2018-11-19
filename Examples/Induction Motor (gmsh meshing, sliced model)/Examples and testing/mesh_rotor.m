%Meshing rotor.

tol_ag = 1e-3;
tol_core1 = 2e-3;
tol_core = 2.3e-3;

gw = gwrap(gmsh_path); %gwrap object

wt = 2.8e-3; %tooth width
ht = 10e-3;
rto = dim.Rout - 2.8e-3; %this influences the iron bridge thickness on top of rotor slots
rti = rto - ht;

Xto = [rto; wt/2];
Xti = [rti; wt/2];

Xout = [dim.Rout;0];
Xin = [dim.Rin;0];

alpha_slot = 2*pi/dim.Qr; %slot pitch angle
Dt = [cos(alpha_slot/2) -sin(alpha_slot/2);sin(alpha_slot/2) cos(alpha_slot/2)]; %rotation matrix
D = Dt*[1 0;0 -1]*Dt'; %mirror matrix over the slot center line (radial)

%adding rotor slot
%%{
gw.addPcws('arc', Dt*[rto;0], Xto, alpha_slot/2, tol_core1, ...
    'arc', Dt*[rto;0], alpha_slot/2, D*Xto, tol_core1, ...
    'line', D*Xto, D*Xti, tol_core1, ...
    'arc', Dt*[rti;0], D*Xti, pi+alpha_slot/2, tol_core1, ...
    'arc', Dt*[rti;0], pi+alpha_slot/2, Xti, tol_core1, ...
    'line', Xti, Xto, tol_core1, ...
    'Slot');
%}

%adding bounding box
gw.addPcws('arc', [0;0], Xout, D*Xout, tol_ag, 'linename', 'nag_r', ...
    'line',  D*Xout, D*Xin, tol_core, 'linename', 'n_ccl', ...
    'arc', [0;0], D*Xin, Xin, tol_core, ...
    'line', Xin, Xout, tol_core, 'linename', 'n_cl',...
    'RotorCore', 'Slot');

%adding shaft
gw.addPcws('arc', [0;0], Xin, D*Xin, tol_core, ...
    'line', D*Xin, [0;0], tol_core, 'linename', 'n_ccl', ...
    'line', [0;0], Xin, tol_core, 'linename', 'n_cl', ...
    'Shaft');

figure(2); clf; hold on; axis equal;
gw.plotSurface('Slot', 'k.-');
gw.plotSurface('RotorCore', 'r.-');
gw.plotSurface('Shaft', 'g.-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meshing

%meshing
gw.removeDuplicates(); 
gw.writeFile();
gw.mesh();
[p_r, t_r, Surfaces_r] = gw.loadMesh();
Ne_r_orig = size(t_r,2);

figure(3); clf; hold on; box on; axis equal;
triplot(t_r', p_r(1,:), p_r(2,:), 'b');

n_ag_r = sortSegmentEdges(p_r, Surfaces_r.get('nag_r'));
n_cl_r = sortRadialEdges(p_r, Surfaces_r.get('n_cl'));
n_ccl_r = sortRadialEdges(p_r, Surfaces_r.get('n_ccl'));
n_dir_r = intersect(n_cl_r, n_ccl_r);

plot(p_r(1, n_cl_r), p_r(2, n_cl_r), 'ro-');
plot(p_r(1, n_ccl_r), p_r(2, n_ccl_r), 'go-');
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replicating

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replicating

%assigning material
m_r = zeros(1, size(t_r,2));
m_r( Surfaces_r.get('RotorCore') ) = dim.RM;
m_r( Surfaces_r.get('Shaft') ) = dim.sM;

[p_r, t_r, n_ccl_r, n_ag_r, n_dir_r] = ...
    replicate_sector_fixed(p_r, t_r, dim.Qr/symmetrySectors, alpha_slot, n_cl_r, n_ccl_r, n_ag_r, n_dir_r);

%n_cl_r = setdiff(n_cl_r, n_dir_r);
%n_ccl_r = setdiff(n_ccl_r, n_dir_r);

%replicating and finalizing material data
m_r = repmat(m_r, 1, dim.Qr/symmetrySectors);
RC = cell(1, dim.Qr/symmetrySectors);
for k = 1:size(RC,2);
    RC{k} = Surfaces_r.get('Slot') + (k-1)*Ne_r_orig;
end

figure(4); clf; hold on; box on; axis equal;
triplot(t_r', p_r(1,:), p_r(2,:), 'b');
plot(p_r(1,n_ccl_r), p_r(2,n_ccl_r), 'go-');
plot(p_r(1,n_cl_r), p_r(2,n_cl_r), 'ro-');
plot(p_r(1,n_ag_r), p_r(2,n_ag_r), 'cv-');