%Meshing a stator.
% 
% A single slot pitch is first meshed with gmsh, using the gwrap interface.
% The full geometry (for a symmetry sector) is then obtained by replicating
% the single-slot mesh.
%
% (c) 2018 Antti Lehikoinen / Smeklab Ltd

%dimensions
symmetrySectors = 8;
dim.SM = 2; %stator core material index
dim.RM = 2; %rotor core material index
dim.sM = 1; %shaft material
dim.Sout = 150e-3 / 2; %outer radius
dim.Sin = 104e-3 / 2; %inner radius
dim.Rout = dim.Sin - 0.5e-3;
dim.Rin = 20e-3;

dim.Qs = 48; %number of stator slots
dim.Qr = 36;

%(rectangular) slot dimensions
dim.S_height = 12e-3; %airgap to slot bottom distance
h_tt = 1.7e-3; %tooth tip height

dim.S_height1 = h_tt; %slot wedge height

dim.S_height3 = dim.S_height - h_tt; %conductor area height
dim.S_width1 = 1.2e-3; %slot opening width

%dim.S_width2 = 2.7e-3; %slot width (opening side)
%dim.S_width3 = 4.2e-3; %slot width (bottom side)
w_tooth = 4e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating gwrap geometry

gw = gwrap(gmsh_path); %gwrap object

%tolerances = max distances between points
tol_ag = 1.1e-3; %airgap
tol_core1 = 1e-3; %slot opening area
tol_core2 = 2e-3; %slot sides
tol_core3 = 10e-3; %stator yoke
tol_out = 15e-3; %stator back-iron

alpha_slot = 2*pi/dim.Qs; %slot pitch angle
%rotation matrix
Dt = [cos(alpha_slot/2) -sin(alpha_slot/2);sin(alpha_slot/2) cos(alpha_slot/2)];
D = Dt*[1 0;0 -1]*Dt'; %mirror matrix over the slot center line (radial)
%Dt = D'*Dt;


%angular pitch of slot opening
aso = alpha_slot/2 - dim.S_width1/dim.Sin/2;

%defining control points for slot
Xso = Dt*[dim.Sin; -dim.S_width1/2]; %slot opening corner (airgap side)
Xso = Xso*dim.Sin/norm(Xso);
Xsoi = Dt*[dim.Sin+h_tt/2; -dim.S_width1/2]; %slot opening corner (inner side)

%conductor area, airgap side
Xco = [dim.Sin+dim.S_height-dim.S_height3; w_tooth/2];

%slot bottom
Xcb = [dim.Sin+dim.S_height; w_tooth/2]; %corner point
Xcb1 = Xcb-[1e-3;0]; %adding "fillet"
Xcb2 = Xcb+[0;1e-3];

r_lr = 0.5; %layer height ratio, FIX
Xcc = r_lr*Xco + (1-r_lr)*Xcb; %corner point between layers


%adding surfaces
gw.addPcws('line', Xso, Xsoi, tol_core1, ...
    'line', Xsoi, Xco, tol_core1, ...
    'line', Xco, D*Xco, tol_core1, ...
    'line', D*Xco, D*Xsoi, tol_core1, ...
    'line', D*Xsoi, D*Xso, tol_core1, ...
    'arc', [0;0], D*Xso, Xso, tol_ag, 'linename', 'n_ag_s', ......
    'Wedge');

%conductor layers
gw.addPcws('line', Xco, Xcc, tol_core2, ...
    'line', Xcc, D*Xcc, tol_core2, ...
    'line', D*Xcc, D*Xco, tol_core2, ...
    'line', D*Xco, Xco, tol_core1, ...
    'Layer1');

gw.addPcws('line', Xcc, Xcb1, tol_core2, ...
    'line', Xcb1, Xcb2, tol_core2, ...
    'line', Xcb2, D*Xcb2, tol_core2, ...
    'line', D*Xcb2, D*Xcb1, tol_core2, ...
    'line', D*Xcb1, D*Xcc, tol_core2, ...
    'line', D*Xcc, Xcc, tol_core2, ...
    'Layer2');

%bounding box
Xin = [dim.Sin; 0];
Xout = [dim.Sout; 0];
Xmid = [dim.Sin+dim.S_height; 0];
gw.addPcws('line', Xin, Xmid, -5, 'linename', 'n_cl', ...
    'line', Xmid, Xout, -3, 'linename', 'n_cl', ...
    'arc', [0;0], Xout, D*Xout, tol_out, 'linename', 'n_dir', ...
    'line', D*Xout, D*Xmid, -3, 'linename', 'n_ccl', ...
    'line', D*Xmid, D*Xin, -5, 'linename', 'n_ccl', ...
    'arc', [0;0], D*Xin, D*Xso, tol_ag, 'linename', 'nag_s', ...
    'arc', [0;0], D*Xso, Xso, tol_ag, 'linename', 'nag_s', ...
    'arc', [0;0], Xso, Xin, tol_ag, 'linename', 'nag_s', ...
    'Iron', 'Wedge', 'Layer1', 'Layer2');

figure(2); clf; hold on; axis equal;
gw.plotSurface('Wedge', 'ro-');
gw.plotSurface('Layer1', 'gv-');
gw.plotSurface('Layer2', 'bo-');
gw.plotSurface('Iron', 'co-');

plot( Xsoi(1), Xsoi(2), 'ro', 'MarkerSize', 10);

%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meshing


%meshing
gw.removeDuplicates(); 
gw.writeFile();
gw.mesh();
[p_s, t_s, Surfaces_s] = gw.loadMesh();
Ne_s_orig = size(t_s,2);

figure(3); clf; hold on; box on; axis equal;
%triplot(t_s', p_s(1,:), p_s(2,:), 'b');
triplot(t_s(:, Surfaces_s.get('Wedge'))', p_s(1,:), p_s(2,:), 'c')
triplot(t_s(:, Surfaces_s.get('Layer1'))', p_s(1,:), p_s(2,:), 'b')
triplot(t_s(:, Surfaces_s.get('Layer2'))', p_s(1,:), p_s(2,:), 'r')
triplot(t_s(:, Surfaces_s.get('Iron'))', p_s(1,:), p_s(2,:), 'k')
%triplot(t_r(:, Surfaces_r.get('Magnet3'))', p_r(1,:), p_r(2,:), 'r')
%triplot(t_r(:, Surfaces_r.get('Magnet5'))', p_r(1,:), p_r(2,:), 'r')

n_ag_s = sortSegmentEdges(p_s, Surfaces_s.get('nag_s'));
n_cl_s = sortRadialEdges(p_s, Surfaces_s.get('n_cl'));
n_ccl_s = sortRadialEdges(p_s, Surfaces_s.get('n_ccl'));
n_dir_s = sortSegmentEdges(p_s, Surfaces_s.get('n_dir'));

plot(p_s(1,n_ag_s), p_s(2,n_ag_s), 'ko-');
plot(p_s(1,n_cl_s), p_s(2,n_cl_s), 'ro-');
plot(p_s(1,n_ccl_s), p_s(2,n_ccl_s), 'go-');

plot( gw.p(1, gw.p(3,:)>0), gw.p(2, gw.p(3,:)>0), 'ro');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replicating

%assigning material
m_s = zeros(1, size(t_s,2));
m_s( Surfaces_s.get('Iron') ) = dim.SM;

[p_s, t_s, n_ccl_s, n_ag_s, n_dir_s] = ...
    replicate_sector_fixed(p_s, t_s, dim.Qs/symmetrySectors, alpha_slot, n_cl_s, n_ccl_s, n_ag_s, n_dir_s);

n_cl_s = setdiff(n_cl_s, n_dir_s);
n_ccl_s = setdiff(n_ccl_s, n_dir_s);

%replicating and finalizing material data
m_s = repmat(m_s, 1, dim.Qs/symmetrySectors);
SC = cell(1, 2*dim.Qs/symmetrySectors);
for k = 1:(dim.Qs/symmetrySectors)
    SC{(k-1)*2+1} = Surfaces_s.get('Layer1') + (k-1)*Ne_s_orig;
    SC{(k-1)*2+2} = Surfaces_s.get('Layer2') + (k-1)*Ne_s_orig;
end

figure(4); clf; hold on; box on; axis equal;
triplot(t_s', p_s(1,:), p_s(2,:), 'b');
plot(p_s(1,n_ccl_s), p_s(2,n_ccl_s), 'go-');
plot(p_s(1,n_cl_s), p_s(2,n_cl_s), 'ro-');

plot(p_s(1,n_ag_s), p_s(2,n_ag_s), 'cv-');
%plot(pnew(1,n_dir), pnew(2,n_dir), 'ko');

%checking
if any( abs( sum(p_s(:,n_cl_s).^2,1) - sum(p_s(:,n_ccl_s).^2,1) ) > 1e-4 )
    error('Node ordering, again');
end

dat = struct('p', p_s', 't', t_s', 'matel', m_s', 'LL', n_ccl_s', ...
    'FL', n_cl_s', 'cir', n_dir_s, 'n_ag_s', n_ag_s');
dat.SC = SC;