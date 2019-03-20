%example script.
% Meshes an example geometry with gmsh, and manual-runs a simulation.
%
% Copyright (c) 2018 Antti Lehikoinen / Smeklab

%meshing; see the script files for details
%mesh_stator; 
%mesh_rotor;

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
dimsc.supply_type = defs.current_supply;

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
%msh_fill(mshc, [PMs{2,:}], 'r');
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

%determining supply angles
phi_rotor = pi/8; %angle of rotor axis
phi_stator = pi/(2*dimsc.p) + pi/dimsc.Qs*dimsc.q*(1/2 + 1/2*(1 - 0.5*dimsc.c)); %angle of stator winding axis
phi_bias = -phi_rotor + phi_stator;

angle_layer_1 = pi/(2*dimsc.p) + pi/dimsc.Qs*dimsc.q - 2*pi/dimsc.Qs*dimsc.c; %chorded layer
angle_layer_2 = pi/(2*dimsc.p) + pi/dimsc.Qs*dimsc.q; %non-chorded layer
phi_stator = 0.5*angle_layer_1 + 0.5*angle_layer_2;
phi_bias = (phi_rotor - phi_stator)*dimsc.p;

angles = linspace(0, 2*pi, 100);
Imain = sqrt(2)*6*[cos(phi_bias+angles); cos(phi_bias-2*pi/3+angles); cos(phi_bias-4*pi/3+angles)];
Imain = Imain + randn(size(Imain));

msch.matel(mshc.rotel) = 0;
sim = MachineSimulation(mshc, dimsc);
pars = SimulationParameters('Is', Imain, 'slip', 1, 'rotorAngle', angles/dimsc.p);

%harmonic simulation
sim.run_static(pars);

%flux plot from harmonic analysis
figure(2); clf; hold on; box on;
sim.fluxplot(1, pars, 'static');
drawnow
%return