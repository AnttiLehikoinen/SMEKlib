function sim = sim_initialConditions_CN(sim, pars, varargin)
%sim_initialConditions_CN initial conditions for time-stepping.
%
% (c) 2017 Antti Lehikoinen / Aalto University

f = pars.f;
w = 2*pi*f;

if ~isempty(pars.slip)
    slips = pars.slip;
else
    slips = sim.dims.slip;
end

Jc = JacobianConstructor(sim.msh, Nodal2D(Operators.curl), Nodal2D(Operators.curl), false);

Ntot = size(sim.results.Xh, 1) / 2;
Nui = Ntot - size(sim.matrices.P,1);

kslip = 1;

slip = slips(kslip);
[Stot, Mtot] = get_circuitMatrices(sim, slip);

%{
%finding indices to free non-damped variables (non-conductive regions,
%mostly)
[Itemp, Jtemp] = find(Mtot);
Inc = setdiff(1:Ntot, [Itemp; Jtemp]);

P = blkdiag(sim.matrices.P, eye(Nui));
[~, Jcond] = find(P(Inc,:));

PT = P(:,unique(Jcond));
%}
%%{
Np = sim.Np;
P = blkdiag(sim.matrices.P, eye(Nui));
[~, Jcond] = find(P(1:Np,:));

PT = P(:,unique(Jcond));
%}

%{
Np = sim.Np;
PT = blkdiag(sim.matrices.P, eye(Nui));
%voltage function
phi0 = 0;
U = abs(pars.U) / sim.msh.symmetrySectors * sim.dims.a * sqrt(2);
if sim.dims.connection_stator == defs.delta
    Ufun = @(t)(  U*[cos(w*t-phi0); cos(w*t - 2*pi/3-phi0); cos(w*t - 4*pi/3-phi0)] );
else
    Ufun = @(t)(  U*[cos(w*t-phi0); cos(w*t - pi/3-phi0)] );
end
Finit = Finit + [zeros(Ntot-sim.results.Ni_s, 1); Ufun(0)];
%}

%setting reluctivity function
if ~numel(varargin)
    nu_fun = sim.nu_fun;
else
    nu_fun = varargin{1};
end


%effect of derivatives
Finit = -(  w*Mtot*sim.results.Xh(Ntot + (1:Ntot), kslip) ) + [sim.matrices.F; zeros(Nui,1)];

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

end