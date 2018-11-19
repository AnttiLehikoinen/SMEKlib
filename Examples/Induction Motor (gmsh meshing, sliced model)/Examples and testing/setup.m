%example script.
% Meshes an example geometry with gmsh, and manual-runs a simulation.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

%meshing; see the script files for details
mesh_stator; 
mesh_rotor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beginning to parse models

%some dimensions
dimsc = struct('D_so', dim.Sout*2, ...
    'p', 2, 'a', 1, 'Qs', dim.Qs, 'Qr', dim.Qr, ...
    'fillingFactor', 0.34, 'Nc_slot', 1, ...
    'leff', 120e-3, ...
    'sigma_rotor', 32e6, 'sigma_stator', 45e6, ...
    'N_conductorsInSlot', 1, 'N_series', 42);
dimsc.q = dimsc.Qs/(3*2*dimsc.p); %number of slots per pole and phase

dimsc.type_statorWinding = defs.stranded;
dimsc.type_rotorWinding = defs.cage;
dimsc.connection_stator = defs.delta;

%for computing end-ring resistnace
dimsc.D_ring = (dim.Rout - ht/2)*2;
dimsc.A_ring = ht*ht;

%setting up mesh object
Ne_r = size(t_r, 2);
Np_r = size(p_r, 2);
mshc = MachineMesh([p_r p_s], [t_r Np_r+t_s], [m_r m_s]);
mshc.setSymmetrySectors(symmetrySectors, dimsc);

%setting special nodes
mshc.namedNodes.set( 'Dirichlet', Np_r+n_dir_s );
mshc.namedNodes.set( 'n_ag_s', Np_r+n_ag_s );
mshc.namedNodes.set( 'n_ag_r', n_ag_r );
mshc.namedNodes.set( 'Periodic_slave', [n_ccl_r Np_r+n_ccl_s] );
mshc.namedNodes.set( 'Periodic_master', [n_cl_r Np_r+n_cl_s] );

%setting special elements
mshc.namedElements.set('rotel', 1:Ne_r); %rotor elements
statorConductors = cellfun(@(x)(x+Ne_r), SC, 'UniformOutput', false);
mshc.namedElements.set('statorConductors', statorConductors);
mshc.namedElements.set('rotorConductors', RC);

%use 2nd-order elements?
%The method transforms the mesh elements into non-curved second-order
%elements. For the named nodes (airgap bnd, periodic, Dirichlet) to be
%updated correctly, the following criteria must be met:
%   - The method is only called after the nodes have been set with e.g.
%       msh.namedNodes.set('n_ag_s', stator_ag_nodes);
%   - The named nodes are ordered either radially (periodic boundaries) or
%       circumferentially (airgap nodes, Dirichlet nodes excluding possible
%       center node).
%   - Air gap mesh generation (msh.generateMovingBand()) is only called
%   AFTER msh.2ndOrder().
mshc.to2ndOrder();

%generating airgap triangulation
mshc.generateMovingBand();

%plotting mesh
figure(1); clf; hold on; box on; axis equal;
msh_triplot(mshc, [], 'b');
msh_fill(mshc, find(mshc.matel == dim.SM), [1 1 1]*0.7);
msh_fill(mshc, find(mshc.matel == dim.RM), [1 1 1]*0.7);
msh_fill(mshc, find(mshc.matel == dim.sM), [1 1 1]*0.85);

msh_plot(mshc, mshc.namedNodes.get('Periodic_master'), 'go-');
msh_plot(mshc, mshc.namedNodes.get('Periodic_slave'), 'ro-');
msh_plot(mshc, mshc.namedNodes.get('Dirichlet'), 'ko');
msh_plot(mshc, mshc.namedNodes.get('n_ag_s'), 'cd-');
msh_plot(mshc, mshc.namedNodes.get('n_ag_r'), 'md-');

statorConductors = mshc.namedElements.get('statorConductors');
msh_fill(mshc, statorConductors{12}, 'b');
rotorConductors = mshc.namedElements.get('rotorConductors');
msh_fill(mshc, rotorConductors{9}, 'r');

%plotting airgap triangulation
AGT = mshc.bandData;
[~, pag, tag] = mshc.bandData.t_ag(0);
triplot(tag(1:3,:)', mshc.p(1,:), mshc.p(2,:), 'm');
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running simulation
sim = MachineSimulation(mshc, dimsc);
pars = SimulationParameters('U', 400/sqrt(3), 'slip', 1.5e-2, 'N_periods', 1, 'N_stepsPerPeriod', 400);

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

%plotting currents
figure(3); clf; hold on;
plot(pars.ts, sim.Is');

%flux plot
figure(4); clf;
sim.fluxplot(25, pars);

%torque plot
figure(5); clf; hold on; box on;
T = sim_compute_torque(sim, pars, 'stepping');
plot(pars.ts, T);