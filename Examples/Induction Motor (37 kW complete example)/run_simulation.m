%run_simulation.
%   - Runs harmonic analysis
%   - computes initial conditions for time-stepping
%   - performs time-stepping with circuit equations and sinusoidal voltage
%   - post-processes torque etc

sim = MachineSimulation(mshc, dim); %simulation object

%simulation parameters
pars = SimulationParameters('U', 400, 'slip', 1.25e-2, 'N_periods', 2, 'N_stepsPerPeriod', 200);

sim.run_harmonic(pars); %harmonic analysis

%plotting flux from harmonic analysis
figure(5); clf; hold on; axis equal tight;
sim.fluxplot(-1, pars);
drawnow;

%torque from harmonic analysis
Th = sim_compute_torque(sim, pars, 'harmonic')


sim.init(pars); %initial conditions for time-stepping
sim.run_timestepping(pars); %time-stepping

%plotting flux from time-stepping
figure(6); clf; hold on; axis equal tight;
sim.fluxplot(100, pars);

%plotting phase currents
figure(7); clf; hold on; box on; grid on;
plot(pars.ts, dim.a*sim.Is);
xlabel('Time (s)');
ylabel('Phase current (A)');

%plotting torque
T = sim_compute_torque(sim, pars, 'stepping');
figure(8); clf; hold on; box on; grid on;
plot(pars.ts, T);
xlabel('Time (s)');
ylabel('Torque (Nm)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loss computation

%bar losses
[Pmean, P_bar, J] = sim_compute_CageLosses(sim, pars, true, 1:200);
Pcage = sum(Pmean)

%end-ring losses
Iring = sim.results.Xt( (1:sim.results.Ni_r) + sim.Np + sim.results.Nu_s+sim.results.Nu_r+sim.results.Ni_s, :);
Pring = sum( mean( sim.matrices.Zew_r * Iring.^2, 2) ) * mshc.symmetrySectors

%stator copper losses
Pcu_stator = sum( mean(real(sim.matrices.Zew_s) * (sim.Is).^2, 2) )  * mshc.symmetrySectors

%computing and plotting iron losses
[PFe_tot, Physt, Peddy, Pexcess] = sim_IronLosses(sim, pars, true);

%computing input power
Is = dim.a * sim.Is;
Irms = mean( Is.^2, 2).^0.5;
Sin = sum( pars.U / sqrt(3) * Irms ) %apparent input power

%active input power
phi0 = -pi/6;
Pin = mean( sqrt(2)*pars.U/sqrt(3)*( cos(2*pi*pars.f*pars.ts - phi0).*Is(1, :) + ...
    cos(2*pi*pars.f*pars.ts-2*pi/3 - phi0).*Is(2,:) + cos(2*pi*pars.f*pars.ts-4*pi/3 - phi0).*Is(3,:) ) )

%shaft power
Pout = mean(T(200:end)) * 2*pi*pars.f * (1-pars.slip) / dim.p

%efficiency (%)
eta = 100 * (Pout - PFe_tot) / Pin

%power factor
pf = Pin / Sin
