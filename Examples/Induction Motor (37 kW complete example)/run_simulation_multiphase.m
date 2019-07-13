%run_simulation Run simulations for a multiphase machine.
%   - Runs harmonic analysis
%   - computes initial conditions for time-stepping
%   - performs time-stepping with circuit equations and sinusoidal voltage
%   - post-processes torque etc

sim = MachineSimulation(mshc, dim); %simulation object

%simulation parameters
pars = SimulationParameters('U', Uh, 'slip', 1.25e-2, 'N_periods', 2, 'N_stepsPerPeriod', 200);

sim.run_harmonic(pars); %harmonic analysis

%plotting flux from harmonic analysis
figure(5); clf; hold on; axis equal tight;
sim.fluxplot(-1, pars);
drawnow;

%torque from harmonic analysis
Thm = sim_compute_torque(sim, pars, 'harmonic');

%setting voltage for time-stepping
pars.U = Ut;

sim.init(pars); %initial conditions for time-stepping
sim.run_timestepping(pars); %time-stepping

%plotting flux from time-stepping
figure(6); clf; hold on; axis equal tight;
sim.fluxplot(100, pars);

%plotting phase currents
figure(7); clf; hold on; box on; grid on;
plot(pars.ts, sim.Is);
xlabel('Time (s)');
ylabel('Phase current (A)');

%plotting torque
Tm = sim_compute_torque(sim, pars, 'stepping');
figure(8); clf; hold on; box on; grid on;
plot(pars.ts, Tm);
xlabel('Time (s)');
ylabel('Torque (Nm)');


%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loss computation

%bar losses
[Pmean, P_bar, J] = sim_compute_CageLosses(sim, pars, false, 1:200);
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
Sin = sum( Uampl/sqrt(2) * Irms )*mshc.symmetrySectors %apparent input power

%active input power
Usamples = Ut( pars.ts );
Pin = mean( sum(Usamples(:,200:end).*Is(:,200:end),1) ) * mshc.symmetrySectors

%shaft power
Pout = mean(Tm(200:end)) * 2*pi*pars.f * (1-pars.slip) / dim.p

%efficiency (%)
eta = 100 * (Pout - PFe_tot) / Pin

%power factor
pf = Pin / Sin
