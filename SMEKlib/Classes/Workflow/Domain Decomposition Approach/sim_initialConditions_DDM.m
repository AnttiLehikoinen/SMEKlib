function sim = sim_initialConditions_DDM(sim, pars, varargin)
%sim_initialConditions_DDM
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
Finit = -(  w*Mtot*sim.results.Xh(Ntot + (1:Ntot), kslip) );

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
sim.results.res_prev = Finit;

end
