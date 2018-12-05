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

[Stot, Mtot] = get_circuitMatrices_2(sim, slip);

%numbers of variables
Ntot = size(sim.results.Xh, 1) / 2;
Nui = Ntot - size(sim.matrices.P,1);
Np = sim.Np;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slave-domain contribution
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
%np_free = setdiff(1:Np_slave, [nd nc]);
np_free = setdiff(1:Np_slave, nd);
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

L_s = sim.matrices.Ls;
L_s = [L_s zeros(size(L_s,1), sim.results.Ni_s-size(L_s,2))]; %for dynamic current supply case

Xh = sim.results.Xh;
Xh = reshape(Xh, [], 2);

D0 = P_m2D*Xh(indA,:);
D0 = reshape(D0(:,1) + 1i*D0(:,2), [], Qs_sector);
I0 = L_s*Xh(indI,:);
I0 = reshape(I0(:,1) + 1i*I0(:,2), [], Qs_sector);

Xslave0 = sim.misc.Ximp_H * [D0; I0];

xs = w*M_slave*imag(Xslave0);

%{
nd = setdiff( 1:size(S_slave,1), nfree );
Fbnd = S_slave(nfree,nd)*Xslave0(nd,:);
Fi = [zeros(numel(np_free), size(I0,2)); DR*real(I0)];

Xslave0(nfree,:) =  S_slave(nfree, nfree) \ ( -Fbnd - Fi + xs(nfree,:) );
%}
hslave = zeros(size(S_slave,1), Qs_sector);
hslave(nfree, :) = S_slave(nfree, nfree) \ xs(nfree,:);

%contribution of slave-domain time-derivative to master-domain load vector
FD = zeros(Ntot, 1);
FD(indA) = FD(indA) - P_m2D'*reshape(P_D2s'*S_slave(nd, :)*hslave, [], 1);

%"impulse" solutions
Nu = sim.results.Nu_s + sim.results.Nu_r; 
Ni_r = sim.results.Ni_r; Ni_s = sim.results.Ni_s;
XX = zeros(Np_slave+Nu_slave, ND + Nu_slave);
XX(nfree, :) = S_slave(nfree,nfree) \ [...
    -S_slave(nfree, nd)*P_D2s [sparse(numel(np_free), Nu_slave); -DR]];
XX(nd, 1:ND) = P_D2s;
Hvar = [P_D2s'*S_slave(nd, :)*XX(:,:);
    XX((Np_slave+1):end,:)];
tempcell_SDD = cell(1, Qs_sector); [tempcell_SDD{:}] = deal( Hvar(1:ND, 1:ND) );
tempcell_SDI = cell(1, Qs_sector); [tempcell_SDI{:}] = deal( Hvar(1:ND, (ND+1):end) );
tempcell_ID = cell(1, Qs_sector); [tempcell_ID{:}] = deal( Hvar((ND+1):end, 1:ND) );
tempcell_II = cell(1, Qs_sector); [tempcell_II{:}] = deal( Hvar((ND+1):end, (ND+1):end) );
Q_SDD = P_m2D'*blkdiag(tempcell_SDD{:})*P_m2D; clear tempcell_SDD;
Q_SDI = P_m2D'*blkdiag(tempcell_SDI{:})*L_s; clear tempcell_SDI;
Q_ID = L_s'*blkdiag(tempcell_ID{:})*P_m2D; clear tempcell_ID;
Q_II = L_s'*blkdiag(tempcell_II{:})*L_s; clear tempcell_II;
Qred = [Q_SDD sparse(sim.Np, Nu) Q_SDI sparse(sim.Np, Ni_r);
    sparse(Nu, sim.Np + Nu + Ni_r + Ni_s);
    Q_ID sparse(Ni_s,Nu) Q_II];
Stot = Stot + Qred;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%free variables: only those nodal potentials not coupled to the reduced
%domains
%[~,J] = find(sim.msh.misc.P_m2D);
%P = blkdiag(sim.matrices.P, eye(Nui));
%[~, Jcond] = find(P(setdiff(1:Np, J),:));
%PT = P(:,unique(Jcond));

%scratch that, free variables as normally
Np = sim.Np;
P = blkdiag(sim.matrices.P, eye(Nui));
[~, Jcond] = find(P(1:Np,:));

PT = P(:,unique(Jcond));

%setting reluctivity function
if ~numel(varargin)
    nu_fun = sim.nu_fun;
else
    nu_fun = varargin{1};
end


%effect of derivatives
Finit = (  w*Mtot*sim.results.Xh(Ntot + (1:Ntot), kslip) )+ FD + [sim.matrices.F; zeros(Nui,1)];
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
%sim.results.res_prev = (  -w*Mtot*sim.results.Xh(Ntot + (1:Ntot), kslip) )  - [sim.matrices.F; zeros(Nui,1)];
sim.results.res_prev = res + 0*FD;


%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial conditions for slave domain

Xslave0 = XX * [reshape(P_m2D*X0(indA), [], Qs_sector);
    reshape(L_s*X0(indI), [], Qs_sector)] + hslave;

sim.results.Xslave0 = reshape(Xslave0, [], 1);

end
