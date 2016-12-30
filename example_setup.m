%example script.
% Initializes the pre-defined geometry, and runs both time-harmonic and
% time-stepping simulations.
%
% Note that the machine has parallel paths outside the symmetry sector, so
% the coefficient 2 appears a lot :)
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

addpath(genpath('SMEKlib'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting example data

%loading the following data:
% p = nodal coordinates (2x874)
% t = elements (3x1340)
% t_ag = air-gap elements (3x352)
% rotel = list of elements in the rotor
% matel = element materials (0 = air/copper, 1 = shaft, 2 = steel sheet)
% statorConductors_all = elements belonging to stator conductors
% rotorConductors_all = self-evident
geometry;

%setting some dimensions:
% D_so = outer diameter of stator
% p = number of pole-pairs
% a = number of parallel paths
% Qs, Qr = numbers of stator and rotor slot
% fillingFactor = stator winding filling factor
% leff = effective length
% sigma_rotor, sigma_stator = stator and rotor winding conductivities
% A_ring = end-ring cross-sectional area (radial cut)
% D_ring = effective diameter of end-ring
% slip = whaddya think?
% N_series = number of turns in stator
dims = struct('D_so', 310e-3, ...
    'p', 2, 'a', 2, 'Qs', 48, 'Qr', 40, ...
    'fillingFactor', 0.4, 'leff', 0.246, ...
    'sigma_rotor', 25.8E+06 *(230.0 + 80)/(230.0 + 80), 'sigma_stator', 45e6, 'A_ring', 520.00E-06 , 'D_ring', 154.00E-03, ...
    'slip', 0.03, 'N_series', 12);
dims.q = dims.Qs/(3*2*dims.p); %number of slots per pole and phase

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

msh = inittri(p, t);
msh.rotel = rotel;
msh.matel = matel;

%assigning some data
msh.symmetrySectors = 4;
msh.periodicityCoeff = get_periodicityFactor(msh.symmetrySectors, dims.p);
Np = size(msh.p, 2);
Ne = size(msh.t, 2);

% plotting mesh etc
figure(1); clf; hold on; box on; 
axis(dims.D_so*[0 0.5 0 0.5]); daspect([1 1 1]); axis square;

%air (and conductors, but they are overdrawn later)
air = find( msh.matel == 3 );
msh_fill(msh, air, 'r', 'EdgeColor', 'w');

%rotor bars
msh_fill(msh, rotorConductors_all, 'y');

%assigning rotor conductors slot-by-slot
temp_angles = linspace(pi/dims.Qr, 2*pi/msh.symmetrySectors*(1-1/dims.Qr), dims.Qr/msh.symmetrySectors);
rotorConductors = binElementsWithAngle(rotorConductors_all, msh, temp_angles);
msh_fill(msh, rotorConductors{5}, 'r'); 

%stator conductors
msh_fill(msh, statorConductors_all, 'b'); 

%assigning stator conductors slot-by-slot
temp_angles = linspace(pi/dims.Qs, 2*pi/msh.symmetrySectors*(1-1/dims.Qs), dims.Qs/msh.symmetrySectors);
statorConductors = binElementsWithAngle(statorConductors_all, msh, temp_angles);
msh_fill(msh, statorConductors{11}, 'g');

%iron parts
iron = find( msh.matel == 2 );
msh_fill(msh, iron, [1 1 1]*0.5);

%shaft
shaft = find( msh.matel == 1 );
msh_fill(msh, shaft, [1 0.5 0.2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determining winding configuration etc.

%stator winding
W = windingConfiguration_1(dims.q, dims.p);

%plotting winding
figure(2); clf; box on; drawWindingConfiguration(W, dims); axis equal;

%loop matrix for stator
Ls = statorConnectionMatrix(W, 1, 1);
Ls = Ls(1:(dims.Qs/msh.symmetrySectors),:); 
Ls = Ls(:, sum(abs(Ls),1)>0) * dims.N_series;

%loop matrix for rotor
r_er = 1/dims.Qr*pi*dims.D_ring/(dims.sigma_rotor*dims.A_ring); %end-ring segment resistance
[Lr, Zr] = rotorConnectionMatrix(dims.Qr, msh.symmetrySectors, dims.p, 0, r_er);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% material data
%setting up reluctivity interpolation tables and -function
mu0 = pi*4e-7;
nu_struct = initialize_reluctivityStruct_interp1(msh, true, [0 10;0 10/mu0]');
nu_fun = @(B)( calculate_reluctivity(B, nu_struct) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting boundary conditions
%
% Redundant boundary variables (Dirichlet, periodic) are eliminated by
% writing
% a_ll = P*a_free
% where a_free contains the free variables, and P is a Np x Nfree matrix.
% Hence, the final problem to solve will be
% P'*S*P*a_free = P'*f_free

% finding Dirichlet boundary
tol = 1e-4;
n_dir = [find( abs(sum(msh.p.^2,1) - 0.25*dims.D_so^2).^0.5 < tol ) ...
    find( abs(sum(msh.p.^2,1) ).^0.5 < tol )];

%setting periodic boundary data in an ugly fashion
np_master = find( (msh.p(2,:) < tol) & (msh.p(1,:) > tol) );
[~, IA] = sort( msh.p(1, np_master) ); np_master = np_master(IA);
np_slave = find( (msh.p(1,:) < tol) );
[~, IA] = sort( msh.p(2, np_slave) ); np_slave = np_slave(IA);

np_master = setdiff(np_master, n_dir); np_slave = setdiff(np_slave, n_dir);

%setting interpolation matrix from free to all variables
P_data = {[np_slave; np_master; msh.periodicityCoeff*ones(1, numel(np_master))], ...
    [n_dir; zeros(2, numel(n_dir))]};
P = assemble_TotalMasterSlaveMatrix(Np, P_data, []);

%plotting Dirichlet and periodic boundary nodes
figure(1);
msh_plot(msh, np_master, 'ro-');
msh_plot(msh, np_slave, 'go-');
msh_plot(msh, n_dir, 'ko');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting air-gap data for rotor movement
msh.bandData = initializeBandData(msh, [], msh.rotel, t_ag);

figure(1);
msh_fill(struct('t', msh.bandData.t_ag, 'p', msh.bandData.p_ag_virt), [], 'c', 'EdgeColor', 'w'); drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulating
example_timeHarmonicSimulation; drawnow;
example_timeSteppingSimulation; drawnow;