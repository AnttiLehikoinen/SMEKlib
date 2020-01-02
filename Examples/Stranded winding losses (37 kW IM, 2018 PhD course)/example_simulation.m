%A reference simulation.
%
% Copyright (c) 2018 Antti Lehikoinen / Aalto University

addpath(genpath('SMEKlib'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting example data

mesh_example2; %meshing

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
dims = struct('D_so', dim.Sout*2, ...
    'p', 2, 'a', 2, 'Qs', dim.Qs, 'Qr', dim.Qr, ...
    'fillingFactor', 0.4, 'leff', 0.249, ...
    'sigma_rotor', 35.5E+06 *(230.0 + 20)/(230.0 + 80), 'sigma_stator', 45e6, 'A_ring', 520.00E-06 , 'D_ring', 154.00E-03, ...
    'slip', 0.015, 'N_series', 12, 'l_halfCoil', 0.520);
dims.q = dims.Qs/(3*2*dims.p); %number of slots per pole and phase

dims.type_statorWinding = defs.stranded;
dims.type_rotorWinding = defs.cage;
dims.connection_stator = defs.delta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

msh_ref = MachineMesh(data.p', data.t', data.matel');
msh_ref.namedElements.add('rotel', data.rotel);
msh_ref.setSymmetrySectors(4, dims);
msh_ref.setMachineBoundaryNodes(dims);


%assigning some data
Np = size(msh_ref.p, 2);
Ne = size(msh_ref.t, 2);

% plotting mesh etc
figure(1); clf; hold on; box on; 
axis(dims.D_so*[0 0.5 0 0.5]); daspect([1 1 1]); axis square;

%air (and conductors, but they are overdrawn later)
air = find( msh_ref.matel == 3 );
msh_fill(msh_ref, air, 'r', 'EdgeColor', 'w');

%assigning conductors
rotorConductors = cell(1, size(data.RC, 1));
for k = 1:numel(rotorConductors)
    rotorConductors{k} = data.RC(k, :);
end
msh_ref.namedElements.add('rotorConductors', rotorConductors);
statorConductors = cell(1, size(data.SC, 1));
for k = 1:numel(statorConductors)
    statorConductors{k} = data.SC(k, :);
end
msh_ref.namedElements.add('statorConductors', statorConductors);

msh_fill(msh_ref, cell2mat(statorConductors), 'b');
msh_fill(msh_ref, statorConductors{11}, 'g');

msh_fill(msh_ref, cell2mat(rotorConductors), 'y');
msh_fill(msh_ref, rotorConductors{5}, 'r'); 

%iron parts
iron = find( msh_ref.matel == 2 );
msh_fill(msh_ref, iron, [1 1 1]*0.5);

%shaft
shaft = find( msh_ref.matel == 1 );
msh_fill(msh_ref, shaft, [1 0.5 0.2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

msh_ref.namedNodes.add('n_ag_s', data.n_ag_s');
msh_ref.namedNodes.add('n_ag_r', data.n_ag_r');
msh_ref.generateMovingBand();

%plotting air-gap triangulation
%[tag, p_ag_virt, ~] = msh_ref.bandData.t_ag(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting boundary conditions
np_master = msh_ref.namedNodes.get('Periodic_master');
np_slave = msh_ref.namedNodes.get('Periodic_slave');
n_dir = msh_ref.namedNodes.get('Dirichlet');
figure(1);
msh_plot(msh_ref, np_master, 'ro-');
msh_plot(msh_ref, np_slave, 'go-');
msh_plot(msh_ref, n_dir, 'ko');

n_ag_s = msh_ref.namedNodes.get('n_ag_s'); n_ag_r = msh_ref.namedNodes.get('n_ag_r');
msh_plot(msh_ref, n_ag_s, 'bo-');
msh_plot(msh_ref, n_ag_r, 'ro-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sim_refulating
sim_ref = MachineSimulation(msh_ref, dims);
pars_ref = SimulationParameters('U', 400/sqrt(3), 'slip', 1.5e-2, 'N_periods', 1, 'N_stepsPerPeriod', 200);

sim_ref.run_harmonic(pars_ref);

figure(11); clf; hold on; box on;
A = sim_ref.results.Xh(1:sim_ref.Np, 1);
drawFluxDensity(msh_ref, A, 'LineStyle', 'none'); colormap('jet'); colorbar; caxis([0 2]);
drawFluxLines(msh_ref, A, 16, 'k');
axis(dims.D_so/2*[-1 1 0 1]); box on; axis tight; daspect([1 1 1]);

sim_ref.init(pars_ref);
tic
sim_ref.run_timestepping(pars_ref);
toc

figure(12); clf; hold on;
plot(pars_ref.ts, 2*sim_ref.Is');

figure(13); clf;
sim_ref.fluxplot(25, pars_ref);

%saving variables
Ish_ref = sim_ref.Ish;
Is_ref = sim_ref.Is;