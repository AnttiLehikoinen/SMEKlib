%mesh_rotor Meshing example rotor
% (FCMSMEK rotor shape 5)

PLOTTING_ON = true; %detailed plotting, for debugging

gmsh_path = dim.gmsh_path;


%extracting dimensions
Rout = dim.Rout; %outer diameter
Rin = dim.Rin; %inner diameter

Qr = dim.Qr;

%dependent dimensions
aslot = 2*pi/Qr;

tol_shaft = 5e-3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting geometry
O = [0;0]; %origin
gw = gwrap(gmsh_path);

%adding slot
D12 = [cos(aslot/2) -sin(aslot/2); sin(aslot/2) cos(aslot/2)]; %rotation over half-slot pitch
D1 = [cos(aslot) -sin(aslot); sin(aslot) cos(aslot)]; %rotation over slot pitch
D = D12*[1 0;0 -1]*D12'; %mirror over slot center line


%adding bar...
rbout = dim.w_bar_out/2;

Xbout1 = D12*[Rout - dim.h_bridge_r - rbout - sqrt( rbout^2 - dim.w_bar_mid^2/4 ); -dim.w_bar_mid/2];
Xboutc = D12*[Rout - dim.h_bridge_r - rbout; 0];
Xbout2 = D12*[Rout - dim.h_bridge_r; 0];

Xbmidc = D12*[Rout - dim.h_bridge_r - dim.h_bar_out; 0];
Xbmid1 = D12*[Rout - dim.h_bridge_r - dim.h_bar_out; -rbout];
Xbmid2 = D12*[Rout - dim.h_bridge_r - dim.h_bar_out + sqrt( rbout^2 - dim.w_bar_mid^2/4 ); -dim.w_bar_mid/2];

rbin = dim.w_bar_mid/2;
Xbinc = D12*[Rout - dim.h_bridge_r - dim.h_bar_out - dim.h_bar_in; 0];
Xbin2 = D12*[Rout - dim.h_bridge_r - dim.h_bar_out - dim.h_bar_in; -rbin];
Xbin2(2) = Xbmid1(2);
Xbin1 = D12*[Rout - dim.h_bridge_r - dim.h_bar_out - dim.h_bar_in - norm(Xbin2-Xbinc); 0];

%adding bar
tol_bar_out = 1.5e-3;
tol_bar_mid = 2.0e-3;
tol_bar_in = 2.0e-3;
gw.addPcws('arc', Xboutc, Xbout1, Xbout2, tol_bar_out, ...
    'arc', Xboutc, Xbout2, D*Xbout1, tol_bar_out, ...
    'line', D*Xbout1, Xbout1, tol_bar_out, ...
    'Bar_out');
gw.addPcws('line', Xbout1, D*Xbout1, tol_bar_out, ...
    'line', D*Xbout1, D*Xbmid2, tol_bar_mid, ...
    'line', D*Xbmid2, Xbmid2, tol_bar_mid, ...
    'line', Xbmid2, Xbout1, tol_bar_mid, ...
    'Bar_mid');

gw.addPcws('arc', Xbmidc, Xbmid1, Xbmid2, tol_bar_in, ...
    'line',  Xbmid2, D*Xbmid2, tol_bar_mid, ...
    'arc', Xbmidc, D*Xbmid2, D*Xbmid1, tol_bar_in, ...
    'line', D*Xbmid1, D*Xbin2, tol_bar_in, ...
    'arc', Xbinc, D*Xbin2, Xbin1, tol_bar_in, ...
    'arc', Xbinc, Xbin1, Xbin2, tol_bar_in, ...
    'line', Xbin2, Xbmid1, tol_bar_in, ...
    'Bar_in');


%bounding box, core part
Xin = [Rin; 0];
Xout = [Rout; 0];
Xmid = [Rout - dim.h_bar_out; 0];
Nper1 = ceil( (Xmid(1)-Xin(1)) / (tol_core) ) + 1;
Nper2 = ceil( (Xout(1)-Xmid(1)) / (0.4*tol_core) ) + 1;

gw.addPcws('line', Xin, Xmid, -Nper1, 'linename', 'n_cw', ...
    'line', Xmid, Xout, -Nper2, 'linename', 'n_cw', ...
    'arc', [0;0], Xout, D*Xout, tol_ag, 'linename', 'n_ag_r', ...
    'line', D*Xout, D*Xmid, -Nper2, 'linename', 'n_ccw', ...
    'line', D*Xmid, D*Xin, -Nper1, 'linename', 'n_ccw', ...
    'arc', [0;0], D*Xin, Xin, tol_core, ...
    'Core', 'Bar_out', 'Bar_mid', 'Bar_in');

%shaft
Nper = ceil( Rin/tol_shaft) + 1;
gw.addPcws('line', O, Xin, -Nper, 'linename', 'n_cw', ...
    'arc', [0;0], Xin, D*Xin, tol_core, ...
    'line', D*Xin, O, -Nper, 'linename', 'n_ccw', ...
    'Shaft');

if PLOTTING_ON  
    figure(2); clf; hold on; box on; axis equal;
    %msh_fill(msh, rotorConductors{1}, 'y');
    
    gw.plotSurface('Bar_out', 'b.-');
    gw.plotSurface('Bar_mid', 'k.-');
    gw.plotSurface('Bar_in', 'k.-');
    %gw.plotSurface('Layer2', 'rx-');    
    %gw.plotSurface('Wedge', 'mx-');    
    gw.plotSurface('Core', 'ro-');
    
    pointplot(Xbout1, 'Xbout1', 'ro');
    pointplot(Xboutc, 'Xboutc', 'ro');
    pointplot(Xbout2, 'Xbout2', 'ro');
    
    pointplot(Xbmid1, 'Xbmid1', 'ro');
    pointplot(Xbmidc, 'Xbmidc', 'ro');
    pointplot(Xbmid2, 'Xbmid2', 'ro');
    
    pointplot(Xbinc, 'Xbinc', 'ro');
    pointplot(Xbin1, 'Xbin1', 'ro');
    pointplot(Xbin2, 'Xbin2', 'ro');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meshing
fname = 'rotor_example.geo';

%meshing
gw.removeDuplicates(); 
gw.writeFile(fname);
gw.mesh(fname);
[p_r, t_r, Surfaces_r] = gw.loadMesh(fname);
Ne_r_orig = size(t_r,2);

%deleting temporary files
delete(strrep(fname, 'geo', 'msh'));
delete(fname);

n_ag_r = sortSegmentEdges(p_r, Surfaces_r.get('n_ag_r'));
n_cw_r = sortRadialEdges(p_r, Surfaces_r.get('n_cw'));
n_ccw_r = sortRadialEdges(p_r, Surfaces_r.get('n_ccw'));
n_dir_r = intersect(n_cw_r, n_ccw_r);

n_ccw_r = setdiff(n_ccw_r, n_dir_r, 'stable');
n_cw_r = setdiff(n_cw_r, n_dir_r, 'stable');


if PLOTTING_ON
    figure(3); clf; hold on; box on; axis equal;
    triplot(t_r', p_r(1,:), p_r(2,:), 'b');
    triplot(t_r(:, Surfaces_r.get('Bar_out'))', p_r(1,:), p_r(2,:), 'r')
    %triplot(t_r(:, Surfaces_r.get('Layer1'))', p_r(1,:), p_r(2,:), 'b')
    %triplot(t_r(:, Surfaces_r.get('Layer2'))', p_r(1,:), p_r(2,:), 'r')
    triplot(t_r(:, Surfaces_r.get('Core'))', p_r(1,:), p_r(2,:), 'k')
    %triplot(t_r(:, Surfaces_r.get('Magnet3'))', p_r(1,:), p_r(2,:), 'r')
    %triplot(t_r(:, Surfaces_r.get('Magnet5'))', p_r(1,:), p_r(2,:), 'r')


    plot(p_r(1,n_ag_r), p_r(2,n_ag_r), 'bo-');
    plot(p_r(1,n_cw_r), p_r(2,n_cw_r), 'ro-');
    plot(p_r(1,n_ccw_r), p_r(2,n_ccw_r), 'go-');
    plot(p_r(1,n_dir_r), p_r(2,n_dir_r), 'ko-');

    %plot( gw.p(1, gw.p(3,:)>0), gw.p(2, gw.p(3,:)>0), 'ro');
end

%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replicating

%assigning material
m_r_temp = zeros(1, size(t_r,2));
m_r_temp( Surfaces_r.get('Core') ) = dim.RM;
m_r_temp( Surfaces_r.get('Shaft') ) = 1;

[p_r, t_r, n_ccw_r, n_ag_r, n_dir_r] = ...
    replicate_sector_fixed(p_r, t_r, Qr/dim.symmetrySectors, aslot, n_cw_r, n_ccw_r, n_ag_r, n_dir_r);

%replicating and finalizing material data
m_r_temp = repmat(m_r_temp, 1, Qr/dim.symmetrySectors);

%setting conductors
RC = cell(1, Qr/dim.symmetrySectors);
for k = 1:(Qr/dim.symmetrySectors)
    
    RC{k} = (k-1)*Ne_r_orig + [ Surfaces_r.get('Bar_out') Surfaces_r.get('Bar_mid') Surfaces_r.get('Bar_in')];
end


if PLOTTING_ON
    figure(4); clf; hold on; box on; axis equal;
    triplot(t_r', p_r(1,:), p_r(2,:), 'b');
    plot(p_r(1,n_ccw_r), p_r(2,n_ccw_r), 'go-');
    plot(p_r(1,n_cw_r), p_r(2,n_cw_r), 'ro-');

    plot(p_r(1,n_ag_r), p_r(2,n_ag_r), 'cv-');
    plot(p_r(1,n_dir_r), p_r(2,n_dir_r), 'ko');
end

datr = struct();
datr.t_r = t_r;
datr.p_r = p_r;
datr.n_dir_r = n_dir_r;
datr.n_ccw_r = n_ccw_r;
datr.n_cw_r = n_cw_r;
datr.m_r = m_r_temp;
datr.n_ag_r = n_ag_r;
datr.RC = RC;