%example script.
% Meshes an example geometry with gmsh, and manual-runs a simulation.
%
% Copyright (c) 2018 Antti Lehikoinen / Smeklab

%meshing; see the script files for details
mesh_stator; 
mesh_rotor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beginning to parse models

%some dimensions
dimsc = struct('D_so', dim.Sout*2, ...
    'p', p, 'a', 1, 'Qs', dim.Qs, 'Qr', dim.Qr, ...
    'fillingFactor', 0.34, 'Nc_slot', 1, ...
    'leff', 120e-3, ...
    'sigma_rotor', 32e6*1e-6, 'sigma_stator', 45e6, ...
    'N_conductorsInSlot', 1, 'N_series', 16, 'N_layers', 2);
dimsc.q = dimsc.Qs/(3*2*dimsc.p); %number of slots per pole and phase

dimsc.type_statorWinding = defs.stranded;
dimsc.connection_stator = defs.delta;

dimsc.c = 2; %chording by two slots

%setting up mesh object
Ne_r = size(t_r, 2);
Np_r = size(p_r, 2);
mshc = MachineMesh([p_r p_s], [t_r Np_r+t_s], [m_r m_s]);
mshc.setSymmetrySectors(symmetrySectors, dimsc);

%setting special nodes
mshc.namedNodes.set( 'Dirichlet', [n_dir_r Np_r+n_dir_s] );
mshc.namedNodes.set( 'n_ag_s', Np_r+n_ag_s );
mshc.namedNodes.set( 'n_ag_r', n_ag_r );
mshc.namedNodes.set( 'Periodic_slave', [n_ccl_r Np_r+n_ccl_s] );
mshc.namedNodes.set( 'Periodic_master', [n_cl_r Np_r+n_cl_s] );

%setting special elements
mshc.namedElements.set('rotel', 1:Ne_r); %rotor elements
statorConductors = cellfun(@(x)(x+Ne_r), SC, 'UniformOutput', false);
mshc.namedElements.set('statorConductors', statorConductors);
mshc.namedElements.set('PMs', PMs);

%setting PMs as solid conductors, shorted at their by end-plate
%Note the huuuuuge conductivity sigma_rotor btw :)
%mshc.namedElements.set('rotorConductors', {PMs{2,:}});
dimsc.type_rotorWinding = defs.none;
%dimsc.type_rotorWinding = defs.user_defined;
%dimsc.Lr = [1 0;-1 -1;0 1]; %rotor loop matrix
%dimsc.Zew_r = zeros(2); %rotor end-plate impedance matrix


%generating airgap triangulation
mshc.generateMovingBand();

%plotting mesh
figure(1); clf; hold on; box on; axis equal;
msh_triplot(mshc, [], 'b');
msh_fill(mshc, find(mshc.matel == dim.SM), [1 1 1]*0.7);
msh_fill(mshc, find(mshc.matel == dim.RM), [1 1 1]*0.7);

msh_plot(mshc, mshc.namedNodes.get('Periodic_master'), 'go-');
msh_plot(mshc, mshc.namedNodes.get('Periodic_slave'), 'ro-');
msh_plot(mshc, mshc.namedNodes.get('Dirichlet'), 'ko');
msh_plot(mshc, mshc.namedNodes.get('n_ag_s'), 'cd-');
msh_plot(mshc, mshc.namedNodes.get('n_ag_r'), 'md-');

statorConductors = mshc.namedElements.get('statorConductors');
rotorConductors = mshc.namedElements.get('rotorConductors');
msh_fill(mshc, statorConductors{6}, 'b');
msh_fill(mshc, [PMs{2,:}], 'r');
%msh_fill(mshc, rotorConductors{1}, 'm');

%plotting airgap triangulation
AGT = mshc.bandData;
[~, pag, tag] = mshc.bandData.t_ag(0);
triplot(tag', mshc.p(1,:), mshc.p(2,:), 'm');
drawnow;

detF = abs(mappingDeterminant(mshc.getMappingMatrix()));
[~, ind] = min(detF);
msh_triplot(mshc, ind, 'r', 'linewidth', 2);

%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running simulation
%%{

%determining supply angles
phi_rotor = pi/8; %angle of rotor axis
%phi_stator = pi/(2*dimsc.p) + pi/dimsc.Qs*dimsc.q*(1/2 + 1/2*(1 - 0.5*dimsc.c)); %angle of stator winding axis

angle_layer_1 = pi/(2*dimsc.p) + pi/dimsc.Qs*dimsc.q - 2*pi/dimsc.Qs*dimsc.c; %chorded layer
angle_layer_2 = pi/(2*dimsc.p) + pi/dimsc.Qs*dimsc.q; %non-chorded layer
phi_stator = 0.5*angle_layer_1 + 0.5*angle_layer_2;
phi_bias = -phi_rotor + phi_stator;

sim = MachineSimulation(mshc, dimsc);
pars = SimulationParameters('U', 525, 'f', 200, 'slip', 0, 'N_periods', 2, 'N_stepsPerPeriod', 600, 'isDC', true, ...
    'phi0', phi_bias-pi/2 - pi/180*11);

%uncomment for elegant no-load simulation
%sim.matrices.Zew_s = eye(3)*1e9;

%harmonic simulation
sim.run_harmonic(pars);

%flux plot from harmonic analysis
figure(2); clf; hold on; box on;
sim.fluxplot(-1, pars);
drawnow
%return

%running time-stepping single-slice simulation
sim.init(pars); %initial conditions
tic
sim.run_timestepping(pars);
toc
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESSING

%plotting currents
figure(3); clf; hold on;
I_rated = sim.Is;
plot(pars.ts, I_rated');
plot(pars.ts, sum(I_rated,1), 'k');

%flux plot
figure(4); clf;
sim.fluxplot(2, pars);

%torque plot
figure(5); clf; hold on; box on;
T_rated = sim_compute_torque(sim, pars, 'stepping');
plot(pars.ts, T_rated);
wm = 2*pi*pars.f/dimsc.p;
Pave = mean(T_rated)*wm;
Sin = pars.U * sum(rms(sim.Is,2));

%current density in slots
Alayer = dimsc.leff / (dimsc.sigma_stator*full(mean(diag(sim.matrices.DRs)))) / dimsc.fillingFactor;
Jrms = dimsc.N_series*rms(sim.Is,2) / Alayer *1e-6

%back-emf
ts = pars.ts; dt = ts(2)-ts(1);
Phi_rated = sim.matrices.Ls' * dimsc.leff*sim.matrices.Cs' * sim.results.Xt(1:sim.Np,:) * mshc.symmetrySectors/sim.dims.a;
E_rated = diff(Phi_rated, [], 2); 
E_rated = [E_rated(:,1) E_rated]/dt;
w = 2*pi*pars.f; phi0 = pars.phi0; %phi0 = pi/180 * 60;
Uplot = pars.U * sim.dims.a * sqrt(2)* ...
    [cos(w*ts-phi0); cos(w*ts - 2*pi/3-phi0); cos(w*ts - 4*pi/3-phi0)];
%Uloop = pars.U/sim.msh.symmetrySectors * sim.dims.a * sqrt(2)*[cos(w*ts-phi0); cos(w*ts - pi/3-phi0)];
%Uplot = [1 -1 0;0 1 -1;-1 0 1]*Uloop;

figure(6); clf; hold on; box on;
plot(ts, E_rated(1, :), 'b');
plot(ts, E_rated(2, :), 'r');
plot(ts, E_rated(3, :), 'k');
plot(ts, Uplot(1, :), 'b--');
plot(ts, Uplot(2, :), 'r--');
plot(ts, Uplot(3, :), 'k--');

%Idq_plot;
save('data_rated.mat', 'I_rated', 'T_rated', 'E_rated', 'Phi_rated', 'Babs_all', 'Bcirc_all', 'Brad_all');