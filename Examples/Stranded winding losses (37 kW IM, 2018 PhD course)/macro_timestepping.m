%Running time-stepping analysis.
%
% Copyright (c) 2018 Antti Lehikoinen / Aalto University

pars = SimulationParameters('U', 400/sqrt(3), 'slip', 1.5e-2, 'N_periods', 1, 'N_stepsPerPeriod', 100);

%initializing time-stepping:
%sim = sim_computeSlaveDomainSolution_timestepping(sim, pars); %computing the impulse response function
%sim = sim_initialConditions_DDM(sim, pars); %initial conditions for time-stepping analysis


%single-pulse PWM supply
%{
U1 = sqrt(2)*230/2; %amplidute of fundamental (division by 2 has to be there due to symmetry sectors and parallel paths)
Ubus = U1; %DC-bus voltage
pars.U = @(t)( PWM_1pulse(2*pi*pars.f*t-[0;2*pi/3;4*pi/3], U1, Ubus) );
%}

%sine-triangle PWM supply
%%{
U1 = sqrt(2)*230/2; %amplidute of fundamental (division by 2 has to be there due to symmetry sectors and parallel paths)
Ubus = 400; %DC-bus voltage
fs = 1000; %switching frequency
pars.U = @(t)( ThreePhasePWM(U1, pars.f, t, Ubus, fs) );
%}

sim = sim_runtimeSteppingSimulation_DDM(sim, pars);

%plotting currents
Is = sim.Is;

figure(7); clf; hold on;
plot(pars.ts, Is);