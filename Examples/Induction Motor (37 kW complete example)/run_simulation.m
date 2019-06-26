%running simulation, postprocessing

%%{

sim = MachineSimulation(mshc, dim);
pars = SimulationParameters('U', 400, 'slip', 1.25e-2, 'N_periods', 1.05, 'N_stepsPerPeriod', 400);

%return
sim.run_harmonic(pars);

figure(5); clf; hold on; axis equal tight;
sim.fluxplot(-1, pars);
drawnow;

T = sim_compute_torque(sim, pars, 'harmonic')

sim.init(pars);
%}


sim.run_timestepping(pars);

%figure(6); clf; hold on; axis equal tight;
%sim.fluxplot(100, pars);

T = sim_compute_torque(sim, pars, 'stepping');
figure(7); clf; hold on; box on; grid on;


%computing and plotting iron losses
[Ptot, Physt, Peddy, Pexcess] = sim_IronLosses(sim, pars, true);

figure(8); caxis([0 30]);
figure(9); caxis([0 100]);