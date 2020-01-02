%Running time-harmonic analysis.
%
% Copyright (c) 2018 Antti Lehikoinen / Aalto University

%simulation parameters
pars = SimulationParameters('U', 400/sqrt(3), 'slip', 1.5e-2, 'N_periods', 1, 'N_stepsPerPeriod', 200);

%uncomment this to change the number of turns
%{
dims.N_series = 12 / 2;
pars.U = 230 / 2;
%}

%uncomment this to change the winding transposition
%{
dims.transpositionArgs = {4, 1};
%}

sim = MachineSimulation(mshc, dims);
sim.dims.W = sim.matrices.W;

sim = sim_computeSlaveDomainSolution_harmonic(sim, pars); %elimination operation
sim = sim_runtimeHarmonicSimulation_DDM(sim, pars); %actual simulation

figure(8); clf;
Ish = sim.Ish;
compass(real(Ish), imag(Ish));

%INSERT YOUR kcc CALCULATION HERE:
disp('The circulating current factor is:')

figure(9); clf;
subplot(1,2,1);
spy(sim.matrices.Ls(1:dims.Nc_slot,:));
title('Winding matrices for slots 1 and 2 (phase A)')
ylabel('Strand number');
xlabel('Slot 1');
subplot(1,2,2);
spy(sim.matrices.Ls(dims.Nc_slot+(1:dims.Nc_slot),:));
xlabel('Slot 2');