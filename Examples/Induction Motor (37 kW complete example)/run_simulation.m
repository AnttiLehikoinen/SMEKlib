%running simulation, postprocessing

%{

sim = MachineSimulation(mshc, dim);
pars = SimulationParameters('U', 400, 'slip', 1.25e-2, 'N_periods', 1.05, 'N_stepsPerPeriod', 400);

%return
sim.run_harmonic(pars);

figure(5); clf; hold on; axis equal tight;
sim.fluxplot(-1, pars);
drawnow;

T = sim_compute_torque(sim, pars, 'harmonic')

sim.init(pars);


sim.run_timestepping(pars);

%figure(6); clf; hold on; axis equal tight;
%sim.fluxplot(100, pars);

%}

figure(6); clf; hold on; box on; grid on;
plot(pars.ts, dim.a*sim.Is);
xlabel('Time (s)');
ylabel('Phase current (A)');

T = sim_compute_torque(sim, pars, 'stepping');
figure(7); clf; hold on; box on; grid on;
plot(pars.ts, T);
xlabel('Time (s)');
ylabel('Torque (Nm)');

%computing copper losses
[Pmean, P_bar, J] = sim_compute_CageLosses(sim, pars, false, 350:400);
Pcage = sum(Pmean)

Iring = sim.results.Xt( (1:sim.results.Ni_r) + sim.Np + sim.results.Nu_s+sim.results.Nu_r+sim.results.Ni_s, :);
Pring = sum( mean( sim.matrices.Zew_r * Iring.^2, 2) ) * mshc.symmetrySectors

Pcu_stator = sum( mean(real(sim.matrices.Zew_s) * (sim.Is).^2, 2) )  * mshc.symmetrySectors


%computing and plotting iron losses
[Ptot, Physt, Peddy, Pexcess] = sim_IronLosses(sim, pars, true);