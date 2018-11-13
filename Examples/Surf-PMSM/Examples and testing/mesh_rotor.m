%Meshing a simple surf-PMSM rotor pole with segmented magnets.
%
% (c) 2018 Antti Lehikoinen / Smeklab Ltd

mu0 = pi*4e-7;

%magnet parameters
Br = 1.25;
mur_PM = 1.05;

%dimensions
delta = 2e-3; %airgap
hpm = 3.2e-3; %magnet height
hyr = 12e-3; %yoke thickness
p = 4; %number of pole-pairs
N_segments = 3; %number of segments
alpha_pm = 0.8; %magnet pitch

%dependent variables
r_magnet = dim.Sin - delta;
r_rotorCore = r_magnet - hpm;
r_shaft = r_rotorCore - hyr;
dim.Rout = r_magnet;

gamma_pole = 2*pi / 8; %pole pitch
gamma_mag = gamma_pole*alpha_pm; %magnet pitch


tol_core = 3e-3;
tol_air = 1.1e-3;
tol_mag = 0.25e-2;

O = [0;0]; %origin

gw = gwrap(gmsh_path);

%setting permanent magnets and shaft
alpha_temp = (gamma_pole - gamma_mag)/2;
mag_angles = linspace(alpha_temp, gamma_pole-alpha_temp, N_segments+1);
arg_corelines = cell(1, 5*N_segments);
arg_maglines = cell(1, 5*N_segments);
for k = 1:N_segments
    Xi1 = r_rotorCore*[cos(mag_angles(k)); sin(mag_angles(k))];
    Xi2 = r_rotorCore*[cos(mag_angles(k+1)); sin(mag_angles(k+1))];
    Xo1 = r_magnet*[cos(mag_angles(k)); sin(mag_angles(k))];
    Xo2 = r_magnet*[cos(mag_angles(k+1)); sin(mag_angles(k+1))];
    
    [arg_corelines{(k-1)*5 + (1:5)}] = deal('arc', O, Xi1, Xi2, tol_mag);
    [arg_maglines{(k-1)*5 + (1:5)}] = deal('arc', O, Xo1, Xo2, tol_air);
    
    gw.addPcws('arc', O, mag_angles(k), mag_angles(k+1), r_rotorCore, tol_mag, ...
        'line', Xi2, Xo2, tol_mag, ...
        'arc', O, Xo2, Xo1, tol_air, 'linename', 'nag_r', ...
        'line', Xo1, Xi1, tol_mag, ['Magnet' num2str(k)]);    
end

%adding shaft
gw.addPcws('line', O, [r_shaft;0], tol_core, 'linename', 'n_cl', ...
    'arc', O, 0, gamma_pole, r_shaft, tol_core, ...
    'line', r_shaft*[cos(gamma_pole);sin(gamma_pole)], O, tol_core, 'linename', 'n_ccl', ...
    'Shaft');

%adding rotor core
%%{
gw.addPcws('line', [r_shaft;0], [r_rotorCore;0], tol_core, 'linename', 'n_cl', ...
    'arc', O, 0, mag_angles(1), r_rotorCore, tol_core, ...
    arg_corelines{:}, ...
    'arc', O, mag_angles(end), gamma_pole, r_rotorCore, tol_core, ...
    'line', r_rotorCore*[cos(gamma_pole);sin(gamma_pole)], r_shaft*[cos(gamma_pole);sin(gamma_pole)], tol_core, 'linename', 'n_ccl', ...
    'arc', O, gamma_pole, 0, r_shaft, tol_core, ...
    'RotorCore');


%adding air on rotor side
str_magnets = strcat({'Magnet'}, num2str((1:N_segments)'));
D = [cos(gamma_pole) -sin(gamma_pole);sin(gamma_pole) cos(gamma_pole)];
gw.addPcws('arc', O, 0, mag_angles(1), r_rotorCore, tol_core, ...
    'line', r_rotorCore*[cos(mag_angles(1)); sin(mag_angles(1))], r_magnet*[cos(mag_angles(1)); sin(mag_angles(1))], tol_mag, ...
    'arc', O, mag_angles(1), 0, r_magnet, tol_air, 'linename', 'nag_r', ...
    'line', [r_magnet;0], [r_rotorCore;0], tol_core, 'linename', 'n_cl', ...
    'Air');
gw.addPcws('arc', O, gamma_pole, mag_angles(end), r_rotorCore, tol_core, ...
    'line', r_rotorCore*[cos(mag_angles(end)); sin(mag_angles(end))], r_magnet*[cos(mag_angles(end)); sin(mag_angles(end))], tol_mag, ...
    'arc', O, mag_angles(end), gamma_pole, r_magnet, tol_air, 'linename', 'nag_r', ...
    'line', D*[r_magnet;0], D*[r_rotorCore;0], tol_core, 'linename', 'n_ccl', ...
    'Air');


figure(1); clf; hold on; axis equal;
gw.plotSurface('Shaft', 'ko-');
gw.plotSurface('RotorCore', 'gx-');
gw.plotSurface('Magnet1', 'ro-');
gw.plotSurface('Magnet2', 'ro-');

gw.plotSurface('Air', 'mx-');

%meshing
gw.removeDuplicates();
gw.writeFile();
gw.mesh();
[p_r, t_r, Surfaces_r] = gw.loadMesh();

figure(2); clf; hold on; box on; axis equal;
triplot(t_r', p_r(1,:), p_r(2,:), 'b');
triplot(t_r(:, Surfaces_r.get('Shaft'))', p_r(1,:), p_r(2,:), 'k')
triplot(t_r(:, Surfaces_r.get('RotorCore'))', p_r(1,:), p_r(2,:), 'Color', [1 1 1]*0.2)
triplot(t_r(:, Surfaces_r.get('Magnet1'))', p_r(1,:), p_r(2,:), 'r')
triplot(t_r(:, Surfaces_r.get('Magnet3'))', p_r(1,:), p_r(2,:), 'r')

%getting boundary nodes of different kinds
n_ag_r = sortSegmentEdges(p_r, Surfaces_r.get('nag_r'));
n_cl_r = sortRadialEdges(p_r, Surfaces_r.get('n_cl'));
n_ccl_r = sortRadialEdges(p_r, Surfaces_r.get('n_ccl'));

n_dir_r = intersect(n_cl_r, n_ccl_r);
n_cl_r = setdiff(n_cl_r, n_dir_r, 'stable');
n_ccl_r = setdiff(n_ccl_r, n_dir_r, 'stable');

plot(p_r(1,n_ag_r), p_r(2,n_ag_r), 'ko-');
plot(p_r(1,n_cl_r), p_r(2,n_cl_r), 'ro-');
plot(p_r(1,n_ccl_r), p_r(2,n_ccl_r), 'go-');

%material indices for rotor
m_r = zeros(1, size(t_r,2));
m_r( Surfaces_r.get('Shaft') ) = 0;
m_r( Surfaces_r.get('RotorCore') ) = dim.RM;
PMs = cell(2, N_segments);
for k = 1:N_segments
    m_r( Surfaces_r.get(['Magnet' num2str(k)]) ) = 6;
    
    angle_pm = 0.5*(mag_angles(k) + mag_angles(k+1));
    
    PMs{2, k} = Surfaces_r.get(['Magnet' num2str(k)]);
    PMs{1, k} = -repmat(Br/(mu0*mur_PM)*[cos(angle_pm); sin(angle_pm)], 1, size(PMs{2,k},2));
end

%material data for rotor
sigma_shaft = 1 / (0.23e-6);
sigma_rotorCore = 1 / (0.23e-6);

%checking
if any( abs( sum(p_r(:,n_ccl_r).^2,1) - sum(p_r(:,n_cl_r).^2,1) ) > 1e-4 )
    error('Node ordering, again');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ne_r_e = size(t_r,2);

% replicating
Nrep = 2*pi/gamma_pole;
[p_r, t_r, n_ccl_r, n_ag_r, n_dir_r] = ...
    replicate_sector_fixed(p_r, t_r, Nrep, gamma_pole, n_cl_r, n_ccl_r, n_ag_r, n_dir_r);

n_cl_r = setdiff(n_cl_r, n_dir_r, 'stable');
n_ccl_r = setdiff(n_ccl_r, n_dir_r, 'stable');
n_dir_r = [];

figure(4); clf; hold on;
triplot(t_r', p_r(1,:), p_r(2,:), 'k');
plot(p_r(1,n_ccl_r), p_r(2,n_ccl_r), 'bv-');
plot(p_r(1,n_cl_r), p_r(2,n_cl_r), 'ro-');

plot(p_r(1,n_ag_r), p_r(2,n_ag_r), 'cv-');
plot(p_r(1,n_dir_r), p_r(2,n_dir_r), 'ko');

%return

%replicating and finalizing material data
m_r = repmat(m_r, 1, Nrep);

%magnets
mu0 = pi*4e-7;

PMs = cell(2, 1*Nrep);

for kr = 1:Nrep
    m_r( Surfaces_r.get('Magnet1') ) = 6;
    PMs{2, kr} = Surfaces_r.get('Magnet1') + (kr-1)*Ne_r_e;
    angle_pm = (gamma_pole/2 + (kr-1)*gamma_pole) + mod(kr+1,2)*pi;
    PMs{1, kr} = -repmat(Br/(mu0*mur_PM)*[cos(angle_pm); sin(angle_pm)], 1, size(PMs{2,1},2));

end

figure(4);
km = 4;
triplot(t_r(:, PMs{2, km})', p_r(1,:), p_r(2,:), 'b');

sum(p_r(:,n_ccl_r).^2,1)
sum(p_r(:,n_cl_r).^2,1)
