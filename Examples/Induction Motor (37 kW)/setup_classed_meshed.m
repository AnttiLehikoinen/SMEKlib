%example script.
% Initializes the pre-defined geometry, and runs both time-harmonic and
% time-stepping simulations.
%
% Note that the machine has parallel paths outside the symmetry sector, so
% the coefficient 2 appears a lot :)
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

addpath(genpath('..\..\SMEKlib'));

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

mshc = MachineMesh(data.p', data.t', data.matel');
mshc.namedElements.add('rotel', data.rotel);
mshc.setSymmetrySectors(4, dims);
mshc.setMachineBoundaryNodes(dims);


%assigning some data
Np = size(mshc.p, 2);
Ne = size(mshc.t, 2);

% plotting mesh etc
figure(1); clf; hold on; box on; 
axis(dims.D_so*[0 0.5 0 0.5]); daspect([1 1 1]); axis square;

%air (and conductors, but they are overdrawn later)
air = find( mshc.matel == 3 );
msh_fill(mshc, air, 'r', 'EdgeColor', 'w');

%assigning conductors
rotorConductors = cell(1, size(data.RC, 1));
for k = 1:numel(rotorConductors)
    rotorConductors{k} = data.RC(k, :);
end
mshc.namedElements.add('rotorConductors', rotorConductors);
statorConductors = cell(1, size(data.SC, 1));
for k = 1:numel(statorConductors)
    statorConductors{k} = data.SC(k, :);
end
mshc.namedElements.add('statorConductors', statorConductors);

msh_fill(mshc, cell2mat(statorConductors), 'b');
msh_fill(mshc, statorConductors{11}, 'g');

msh_fill(mshc, cell2mat(rotorConductors), 'y');
msh_fill(mshc, rotorConductors{5}, 'r'); 

%iron parts
iron = find( mshc.matel == 2 );
msh_fill(mshc, iron, [1 1 1]*0.5);

%shaft
shaft = find( mshc.matel == 1 );
msh_fill(mshc, shaft, [1 0.5 0.2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mshc.namedNodes.add('n_ag_s', data.n_ag_s');
mshc.namedNodes.add('n_ag_r', data.n_ag_r');
mshc.generateMovingBand();

%plotting air-gap triangulation
%[tag, p_ag_virt, ~] = mshc.bandData.t_ag(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting boundary conditions
np_master = mshc.namedNodes.get('Periodic_master');
np_slave = mshc.namedNodes.get('Periodic_slave');
n_dir = mshc.namedNodes.get('Dirichlet');
figure(1);
msh_plot(mshc, np_master, 'ro-');
msh_plot(mshc, np_slave, 'go-');
msh_plot(mshc, n_dir, 'ko');

n_ag_s = mshc.namedNodes.get('n_ag_s'); n_ag_r = mshc.namedNodes.get('n_ag_r');
msh_plot(mshc, n_ag_s, 'bo-');
msh_plot(mshc, n_ag_r, 'ro-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulating
sim = MachineSimulation(mshc, dims);
pars = SimulationParameters('U', 400/sqrt(3), 'slip', 1.5e-2, 'N_periods', 1, 'N_stepsPerPeriod', 200);

sim.run_harmonic(pars);

figure(11); clf; hold on; box on;
A = sim.results.Xh(1:sim.Np, 1);
drawFluxDensity(mshc, A, 'LineStyle', 'none'); colormap('jet'); colorbar; caxis([0 2]);
drawFluxLines(mshc, A, 16, 'k');
axis(dims.D_so/2*[-1 1 0 1]); box on; axis tight; daspect([1 1 1]);

sim.init(pars);
tic
sim.run_timestepping(pars);
toc

figure(12); clf; hold on;
plot(pars.ts, 2*sim.Is');
%plot(pars.ts, 2*Xsamples(indI, :)', 'b--'); %comparison to earlier example

figure(13); clf;
sim.fluxplot(25, pars);

Tem = sim_compute_torque(sim, pars, 'stepping')
figure(7); clf; hold on; box on; grid on;
plot(pars.ts, Tem);