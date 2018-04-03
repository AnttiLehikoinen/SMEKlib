function sim = sim_computeSlaveDomainSolution_timestepping(sim, pars)
%sim_computeSlaveDomainSolution_timestepping Compute impulse-response
%solution for the slave domain.
%
% (c) 2018 Antti Lehikoinen / Aalto University

msh = sim.msh.misc.msh_slave;
nd = sim.msh.misc.nd_slave;
conductors = sim.msh.misc.conductors_slave;
Nu = numel(conductors);
Qs_sector = sim.dims.Qs / sim.msh.symmetrySectors;

Np = size(msh.p,2);

%slave-domain matrices
S = sim.matrices.S_slave;
M = sim.matrices.M_slave;
C = sim.matrices.C_slave;
DR = sim.matrices.DR_slave;


end