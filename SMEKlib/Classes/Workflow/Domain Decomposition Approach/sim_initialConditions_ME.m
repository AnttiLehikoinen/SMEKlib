function sim = sim_initialConditions_ME(sim, pars, varargin)
%sim_initialConditions_ME initial conditions for macro-element simulation.
% 
% (c) 2018 Antti Lehikoinen / Aalto University

f = pars.f;
w = 2*pi*f;

if ~isempty(pars.slip)
    slips = pars.slip;
else
    slips = sim.dims.slip;
end
kslip = 1;
slip = slips(kslip);

Jc = JacobianConstructor(sim.msh, Nodal2D(Operators.curl), Nodal2D(Operators.curl), false);

[Stot, Mtot] = get_circuitMatrices(sim, slip);

%numbers of variables
Ntot = size(sim.results.Xh, 1) / 2;
Nui = Ntot - size(sim.matrices.P,1);
Np = sim.Np;

%free variables: only those nodal potentials not coupled to the reduced
%domains
[~,J] = find(sim.msh.misc.P_m2D);
P = blkdiag(sim.matrices.P, eye(Nui));
[~, Jcond] = find(P(setdiff(1:Np, J),:));

PT = P(:,unique(Jcond));

%setting reluctivity function
if ~numel(varargin)
    nu_fun = sim.nu_fun;
else
    nu_fun = varargin{1};
end


%effect of derivatives
Finit = -(  w*Mtot*sim.results.Xh(Ntot + (1:Ntot), kslip) )  + [sim.matrices.F; zeros(Nui,1)];

X0 = sim.results.Xh(1:Ntot, kslip);
Sag = sim.msh.get_AGmatrix(0, Ntot);
for kiter = 1:15
    [J, res] = Jc.eval(X0, nu_fun);
    
    Jtot = PT'*( J + Sag + Stot )*PT;
    res_tot = PT'*( (Sag+Stot)*X0 - Finit + res );
    
    %checking convergence
    disp( norm(res_tot) / norm(Finit) )
    if (norm(res_tot) / norm(Finit)) < 1e-6
        break;
    end
    
    %Newton step
    dX = - Jtot \ res_tot;
    
    X0 = X0 + PT*dX;
end

sim.results.X0 = X0;
sim.results.res_prev = (  w*Mtot*sim.results.Xh(Ntot + (1:Ntot), kslip) )  - 0*[sim.matrices.F; zeros(Nui,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial conditions for slave domain
Qs_sector = sim.dims.Qs / sim.msh.symmetrySectors;

%indices
nd = sim.msh.misc.nd_slave; %Dirichlet nodes
Np_slave = size(sim.msh.misc.msh_slave.p,2);
Nu_slave = numel(sim.msh.misc.conductors_slave);

%slave-domain matrices
S = sim.matrices.S_slave;
M = sim.matrices.M_slave;
C = sim.matrices.C_slave;
DR = sim.matrices.DR_slave;
P_m2D = sim.msh.misc.P_m2D;
P_D2s = sim.msh.misc.P_D2s;
ND = size(P_D2s, 2);

%conductor nodes
[Itemp, ~] = find(C); nc = toRow(unique(Itemp));
np_free = setdiff(1:Np_slave, [nd nc]);
nfree = [np_free (Np_slave+1):(Np_slave+Nu_slave)]; 

%assembling final matrices
S_slave = [S -1/sim.dims.leff*C;
    sparse(Nu_slave, Np_slave) -speye(Nu_slave, Nu_slave)];
M_slave = [M sparse(Np_slave, Nu_slave);
    DR*transpose(C) sparse(Nu_slave, Nu_slave)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%slave-domain solution from harmonic analysis
indA = 1:sim.Np;
indI = (sim.Np + sim.results.Nu_r) + (1:sim.results.Ni_s);
Xh = sim.results.Xh;
Xh = reshape(Xh, [], 2);

D0 = P_m2D*Xh(indA,:);
D0 = reshape(D0(:,1) + 1i*D0(:,2), [], Qs_sector);
I0 = sim.matrices.Ls*Xh(indI,:);
I0 = reshape(I0(:,1) + 1i*I0(:,2), [], Qs_sector);

Xslave0 = sim.misc.Ximp_H * [D0; I0];

sim.results.Xslave0 = reshape(real(Xslave0), 1, []);

end
