%mesh_stator Meshing example stator, with semicircle slot bottom
% (FCSMEK slot shape 4)


PLOTTING_ON = true; %detailed plotting, for debugging

gmsh_path = dim.gmsh_path;

%extracting dimensions
Sout = dim.Sout; %outer diameter
Sin = dim.Sin; %inner diameter

Qs = dim.Qs;

%dependent dimensions
aslot = 2*pi/Qs;

%parsing tolerances
if isfield(dim, 'tol_slot')
    tol_slot = dim.tol_slot_s;
else
    tol_slot = dim.hslot_s / 10;
end
if isfield(dim, 'tol_core')
    tol_core = dim.tol_core_s;
else
    tol_core = abs( dim.Sout - dim.Sin ) / 10;
end
if isfield(dim, 'tol_ag')
    tol_ag = dim.tol_ag_s;
else
    tol_ag = dim.delta*0.9;
end
if isfield(dim, 'tol_opening')
    tol_opening = dim.tol_opening_s;
else
    tol_opening = 2*tol_ag;
end
if isfield(dim, 'tol_out')
    tol_out = dim.tol_out;
else
    tol_out = 5*tol_core;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting geometry
O = [0;0]; %origin
gw = gwrap(gmsh_path);

%adding slot
D12 = [cos(aslot/2) -sin(aslot/2); sin(aslot/2) cos(aslot/2)]; %rotation over half-slot pitch
D1 = [cos(aslot) -sin(aslot); sin(aslot) cos(aslot)]; %rotation over slot pitch
%D12 = D;
D = D12*[1 0;0 -1]*D12'; %mirror over slot center line


Xco = D12*[Sin + dim.htt_s; -dim.ws_o/2];
Xcb = D12*[Sin + dim.hslot_s-dim.ws_b/2; -dim.ws_b/2]; %slot bottom corner arc, start
Xcb_mid = D12*[Sin + dim.hslot_s; 0];
Xcb_center = D12*[Sin + dim.hslot_s-dim.ws_b/2; 0];


%computing layer height ratio
ht = norm( 0.5*(D*Xcb + Xcb) - 0.5*(D*Xco + Xco) ); %height of trapezoidal part
wb = dim.ws_b; %width of trap. area, bottom-side
wo = dim.ws_o; %width of trap. area, opening-side
Af = pi*(wb/2)^2 / 2; %half-circle area
Ao = @(r)( r*ht * (wo + wo+(wb-wo)*r)/2 ); %opening-side trapezoid area
Ab = @(r)( (1-r)*ht * (wb + wb+(wo-wb)*(1-r))/2 + Af );
r_lr = fzero(@(r)( Ao(r) - Ab(r) ), 0.5);

Xcc = (1-r_lr)*Xco + r_lr*Xcb; %corner point between layers


%adding layers
gw.addPcws('line', Xco, Xcc, tol_slot, ...
'line', Xcc, D*Xcc, tol_slot, ...
'line', D*Xcc, D*Xco, tol_slot, ...
'line', D*Xco, Xco, tol_opening, ...
'Layer1');


gw.addPcws('line', Xcc, Xcb, tol_slot, ...
    'arc', Xcb_center, Xcb, Xcb_mid, tol_slot, ...
    'arc', Xcb_center, Xcb_mid, D*Xcb, tol_slot, ...
    'line', D*Xcb, D*Xcc, tol_slot, ...
    'line', D*Xcc, Xcc, tol_slot, ...
    'Layer2');

%adding slot opening
Xso = D12*[Sin; -dim.w_slotOpening_s/2]; Xso = Xso / norm(Xso)*Sin;
Xsoi = D12*[Sin+dim.htt_taper_s; -dim.w_slotOpening_s/2];
gw.addPcws('line', Xso, Xsoi, tol_opening, ...
    'line', Xsoi, Xco, tol_opening, ...
    'line', Xco, D*Xco, tol_opening, ...
    'line', D*Xco, D*Xsoi, tol_opening, ...
    'line', D*Xsoi, D*Xso, tol_opening, ...
    'arc', [0;0], D*Xso, Xso, tol_ag, 'linename', 'n_ag_s', ...
    'Wedge');

%bounding box
Xin = [Sin; 0];
Xout = [Sout; 0];
Xmid = [Sin + dim.hslot_s; 0];
Nper1 = ceil( (Xmid(1)-Xin(1)) / (0.6*tol_core) ) + 1;
Nper2 = ceil( (Xout(1)-Xmid(1)) / tol_core ) + 1;

gw.addPcws('line', Xin, Xmid, -Nper1, 'linename', 'n_cw', ...
    'line', Xmid, Xout, -Nper2, 'linename', 'n_cw', ...
    'arc', [0;0], Xout, D*Xout, tol_out, 'linename', 'n_dir', ...
    'line', D*Xout, D*Xmid, -Nper2, 'linename', 'n_ccw', ...
    'line', D*Xmid, D*Xin, -Nper1, 'linename', 'n_ccw', ...
    'arc', [0;0], D*Xin, D*Xso, tol_ag, 'linename', 'n_ag_s', ...
    'arc', [0;0], D*Xso, Xso, tol_ag, 'linename', 'n_ag_s', ...
    'arc', [0;0], Xso, Xin, tol_ag, 'linename', 'n_ag_s', ...
    'Core', 'Wedge', 'Layer1', 'Layer2');

if PLOTTING_ON
    figure(2); clf; hold on; box on; axis equal;
    gw.plotSurface('Layer1', 'bo-');
    gw.plotSurface('Layer2', 'rx-');    
    gw.plotSurface('Wedge', 'mx-');    
    gw.plotSurface('Core', 'ro-');
    
    %plot(Xcb1(1), Xcb1(2), 'ko');
    %plot(Xcbc(1), Xcbc(2), 'kv');
    %plot(Xcb2(1), Xcb2(2), 'ko');
end

%return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meshing
fname = ['stator_example.geo'];

%meshing
gw.removeDuplicates(); 
gw.writeFile(fname);
gw.mesh(fname);
[p_s, t_s, Surfaces_s] = gw.loadMesh(fname);
Ne_s_orig = size(t_s,2);

%deleting temporary files
delete(strrep(fname, 'geo', 'msh'));
delete(fname);

n_ag_s = sortSegmentEdges(p_s, Surfaces_s.get('n_ag_s'));
n_cw_s = sortRadialEdges(p_s, Surfaces_s.get('n_cw'));
n_ccw_s = sortRadialEdges(p_s, Surfaces_s.get('n_ccw'));
n_dir_s = sortSegmentEdges(p_s, Surfaces_s.get('n_dir'));


if PLOTTING_ON
    figure(3); clf; hold on; box on; axis equal;
    %triplot(t_s', p_s(1,:), p_s(2,:), 'b');
    triplot(t_s(:, Surfaces_s.get('Wedge'))', p_s(1,:), p_s(2,:), 'c')
    triplot(t_s(:, Surfaces_s.get('Layer1'))', p_s(1,:), p_s(2,:), 'b')
    triplot(t_s(:, Surfaces_s.get('Layer2'))', p_s(1,:), p_s(2,:), 'r')
    triplot(t_s(:, Surfaces_s.get('Core'))', p_s(1,:), p_s(2,:), 'k')
    %triplot(t_r(:, Surfaces_r.get('Magnet3'))', p_r(1,:), p_r(2,:), 'r')
    %triplot(t_r(:, Surfaces_r.get('Magnet5'))', p_r(1,:), p_r(2,:), 'r')


    plot(p_s(1,n_ag_s), p_s(2,n_ag_s), 'bo-');
    plot(p_s(1,n_cw_s), p_s(2,n_cw_s), 'ro-');
    plot(p_s(1,n_ccw_s), p_s(2,n_ccw_s), 'go-');
    plot(p_s(1,n_dir_s), p_s(2,n_dir_s), 'ko-');

    %plot( gw.p(1, gw.p(3,:)>0), gw.p(2, gw.p(3,:)>0), 'ro');
end

%verifying layer areas
%%{
N = Nodal2D(Operators.I);
m = SimpleMesh(p_s, t_s);
Mslot = MatrixConstructor().assemble_vector(N, 1, 1, Surfaces_s.get('Layer1'), m).finalize();
Acopper = full(sum( Mslot, 1));
Mslot = MatrixConstructor().assemble_vector(N, 1, 1, Surfaces_s.get('Layer2'), m).finalize();
Acopper2 = full(sum( Mslot, 1));
%%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replicating

%assigning material
m_s_temp = zeros(1, size(t_s,2));
m_s_temp( Surfaces_s.get('Core') ) = dim.SM;

[p_s, t_s, n_ccw_s, n_ag_s, n_dir_s] = ...
    replicate_sector_fixed(p_s, t_s, Qs/dim.symmetrySectors, aslot, n_cw_s, n_ccw_s, n_ag_s, n_dir_s);

%replicating and finalizing material data
m_s_temp = repmat(m_s_temp, 1, Qs/dim.symmetrySectors);

%setting conductors
if isfield(dim, 'N_layers') && dim.N_layers == 2
    SC = cell(2, Qs/dim.symmetrySectors);
    for k = 1:(Qs/dim.symmetrySectors)
        SC{(k-1)*2+1} = Surfaces_s.get('Layer1') + (k-1)*Ne_s_orig;
        SC{(k-1)*2+2} = Surfaces_s.get('Layer2') + (k-1)*Ne_s_orig;
    end
else
    SC = cell(1, Qs/dim.symmetrySectors);
    for k = 1:(Qs/dim.symmetrySectors)
        SC{k} = (k-1)*Ne_s_orig + [Surfaces_s.get('Layer1') Surfaces_s.get('Layer2')];
    end
end


if PLOTTING_ON
    figure(4); clf; hold on; box on; axis equal;
    triplot(t_s', p_s(1,:), p_s(2,:), 'b');
    plot(p_s(1,n_ccw_s), p_s(2,n_ccw_s), 'go-');
    plot(p_s(1,n_cw_s), p_s(2,n_cw_s), 'ro-');

    plot(p_s(1,n_ag_s), p_s(2,n_ag_s), 'cv-');
    plot(p_s(1,n_dir_s), p_s(2,n_dir_s), 'ko');
end

dats = struct();
dats.t_s = t_s;
dats.p_s = p_s;
dats.n_dir_s = n_dir_s;
dats.n_ccw_s = n_ccw_s;
dats.n_cw_s = n_cw_s;
dats.m_s = m_s_temp;
dats.n_ag_s = n_ag_s;
dats.SC = SC;
dats.Aslot = [Acopper Acopper2];