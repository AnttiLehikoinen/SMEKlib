%running simulation, postprocessing

%%{
dim.N_series = 6;
dim.l_halfCoil = 0.520;
dim.a = 2;
dim.N_layers = 1;
dim.A_ring = 520.00E-06;
dim.D_ring = 154.00E-03;

dim.connection_stator = defs.star;

dim.Lew = 0.1e-3;

sim = MachineSimulation(mshc, dim);
pars = SimulationParameters('U', 400, 'slip', 1.25e-2, 'N_periods', 200/400, 'N_stepsPerPeriod', 400, 'maxIter', 50);

%return
sim.run_harmonic(pars);

figure(5); clf; hold on; axis equal tight;
sim.fluxplot(-1, pars);

T = sim_compute_torque(sim, pars, 'harmonic')

sim.init(pars);
%}

%{
sim_initialConditions_temp(sim, pars)

figure(6); clf; hold on; box on; axis equal tight;
drawFluxDensity(mshc, sim.results.X0, 'LineStyle', 'none'); 
colormap('jet'); colorbar; caxis([0 2]);
drawFluxLines(mshc, sim.results.X0, 16, 'k');
%}
%return


sim.run_timestepping(pars);

%figure(6); clf; hold on; axis equal tight;
%sim.fluxplot(100, pars);

T = sim_compute_torque(sim, pars, 'stepping')
figure(7); clf; hold on; box on; grid on;
plot(pars.ts, T);
