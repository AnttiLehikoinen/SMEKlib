function sim = sim_computeSlaveDomainSolution_harmonic(sim, pars)
%sim_computeSlaveDomainSolution_harmonic slave-domain solution in frequency
%domain.
%
% (c) 2018 Antti Lehikoinen / Aalto University

msh = sim.msh.misc.msh_slave;
nd = sim.msh.misc.nd_slave;
conductors = sim.msh.misc.conductors_slave;
Nu = numel(conductors);
Qs_sector = sim.dims.Qs / sim.msh.symmetrySectors;

Np = size(msh.p,2);
nu0 = 1/(pi*4e-7);

S = MatrixConstructor(Nodal2D(Operators.grad), Nodal2D(Operators.grad), nu0, [], msh).finalize();
M = MatrixConstructor(Nodal2D(Operators.I), Nodal2D(Operators.I), sim.dims.sigma_stator, horzcat(conductors{:}), msh, []).finalize(Np,Np);

Cc = MatrixConstructor();
for kc = 1:Nu
    Cc.assemble_vector(Nodal2D(Operators.I), kc, 1, conductors{kc}, msh);
end
C = Cc.finalize(Np, Nu);

%conductor areas
cA_slave = sum(C*speye(Nu, Nu), 1);
DR = sparsediag( sim.dims.leff ./ cA_slave );
DR_ew = kron(speye(Qs_sector, Qs_sector),...
    sparsediag( (sim.dims.l_halfCoil - sim.dims.leff) ./ cA_slave ));
sim.matrices.Zew_s = sim.matrices.Ls'*DR_ew*sim.matrices.Ls;

P_D2s = sim.msh.misc.P_D2s;

w = 2*pi*pars.f;

%actual computation
nfree = setdiff(1:Np, nd); ind_freeVars = [nfree (Np+1):(Np+Nu)]; 
ND = size(P_D2s, 2);

Q = [S(nfree, nfree)+1i*w*M(nfree, nfree) -1/sim.dims.leff*C(nfree,:);
    DR*1i*w*transpose(C(nfree,:)) -speye(Nu, Nu)];
Qbnd_slave = [S(nfree, nd)+1i*w*M(nfree, nd);
    DR*1i*w*transpose(C(nd,:))];

Ximp_H = zeros(Np+Nu, ND + Nu);
Ximp_H(ind_freeVars, :) = Q \ [...
    -Qbnd_slave*P_D2s [sparse(numel(nfree), Nu); -DR]];
Ximp_H(nd, 1:ND) = P_D2s;

sim.misc.Ximp_H = Ximp_H;


%assigning reduced matrices
X_AA = Ximp_H(1:Np, 1:ND); X_uA = Ximp_H((Np+1):end, 1:ND);
X_AI = Ximp_H(1:Np, (ND+1):end); X_uI = Ximp_H((Np+1):end, (ND+1):end);

sim.misc.R_AA = X_AA' * ( (S + 1i*w*M)*X_AA - 1/sim.dims.leff*C*X_uA );    
sim.misc.R_AI = X_AA'*( (S + 1i*w*M)*X_AI - 1/sim.dims.leff*C*X_uI );
sim.misc.R_UA = X_uA; 
sim.misc.R_UI = X_uI;

end
