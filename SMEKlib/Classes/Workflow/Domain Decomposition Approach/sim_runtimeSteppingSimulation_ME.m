function sim = sim_runtimeSteppingSimulation_ME(sim, pars, varargin)
%sim_runtimeHarmonicSimulation_DDM domain-reduced time-harmonic simulation.
%
% (c) 2018 Antti Lehikoinen / Aalto University

%TODO
% [x] CN initial conditions
% [ ] remove quick fixes from e.g. voltage
% [ ] remove unnecessary data from impulse response computation
% [ ] switch to harmonic reluctivity? (maybe not required yet

% adjusted CN for stability
alpha2 = 2; %weight for implicit (k+1) step; 1 for CN, 2 for BE
alpha1 = 2 - alpha2;

%basic setup
f = pars.f;
w = 2*pi*f;

if ~isempty(pars.slip)
    slip = pars.slip;
else
    slip = sim.dims.slip;
end

%voltage-function
L_s = sim.matrices.Ls;
if isa(pars.U, 'function_handle')
    Nin =  nargin(pars.U);
    %N_phases = numel(pars.U(0));
    N_phases = 3; %quick fix
    N_inParallel = size(L_s,2) / N_phases;
    Mphase =  kron(eye(N_phases), ones(N_inParallel,1));
    Ufun = @(t, varargin)(Mphase*pars.U(t, varargin{:}) );
else
    U = pars.U / sim.msh.symmetrySectors * sim.dims.a * sqrt(2);
    Nin = 1;
    %U = 400/2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %voltage function
    phi0 = 0;
    if sim.dims.connection_stator == defs.delta
        Ufun = @(t)(  kron(eye(3), ones(size(L_s,2)/3,1))*U*[cos(w*t-phi0); cos(w*t - 2*pi/3-phi0); cos(w*t - 4*pi/3-phi0)] );
    else
        Ufun = @(t)(  U*[cos(w*t-phi0); cos(w*t - pi/3-phi0)] );
    end
end

%timestamps etc
wm = w/sim.dims.p * (1-slip);
tsamples = pars.ts;
Nsamples = numel(tsamples);
dt = tsamples(2) - tsamples(1);

%numbers of variables
Nu = sim.results.Nu_r;
Ni_r = sim.results.Ni_r;
Ni_s = sim.results.Ni_s;
indI = (sim.Np + Nu) + (1:Ni_s);
indA = 1:sim.Np;

%circuit and other matrices
[Sc, Mtot] = get_circuitMatrices(sim);
Mtot = Mtot / dt;
Jc = JacobianConstructor(sim.msh, Nodal2D(Operators.curl), Nodal2D(Operators.curl), false);
Ntot = size(Sc,1);
Nui = Ntot - sim.Np;
PT = blkdiag(sim.matrices.P, speye(Nui));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reduced-basis coupling blocks
Qs_sector = sim.dims.Qs / sim.msh.symmetrySectors;

%indices
nd = sim.msh.misc.nd_slave; %Dirichlet nodes
Np_slave = size(sim.msh.misc.msh_slave.p,2);
Nu_slave = numel(sim.msh.misc.conductors_slave);
np_free = setdiff(1:Np_slave, nd);
nfree = [np_free (Np_slave+1):(Np_slave+Nu_slave)]; 

%slave-domain matrices
S = sim.matrices.S_slave;
M = sim.matrices.M_slave;
C = sim.matrices.C_slave;
DR = sim.matrices.DR_slave;
P_m2D = sim.msh.misc.P_m2D;
P_D2s = sim.msh.misc.P_D2s;
ND = size(P_D2s, 2);

Nh = numel(nd);
L = sparse(1:Nh, nd, ones(1, Nh), Nh, Np_slave);

%assembling final matrices
S_slave = [S -1/sim.dims.leff*C;
    sparse(Nu_slave, Np_slave) -speye(Nu_slave, Nu_slave)];
M_slave = [M sparse(Np_slave, Nu_slave);
    DR*transpose(C) sparse(Nu_slave, Nu_slave)];
Q_slave = S_slave + 2/(dt*alpha2)*M_slave;

%"impulse" solutions
XX = zeros(Np_slave+Nu_slave, ND + Nu_slave);
XX(nfree, :) = Q_slave(nfree,nfree) \ [...
    -Q_slave(nfree, nd)*P_D2s [sparse(numel(np_free), Nu_slave); -DR]];
XX(nd, 1:ND) = P_D2s;

sim.misc.temp = XX;

Hvar = [P_D2s'*Q_slave(nd, 1:Np_slave)*XX(1:Np_slave,:);
    XX((Np_slave+1):end,:)];

sim.misc.Hvar = Hvar;

%dependencies on the newest time-step
tempcell_SDD = cell(1, Qs_sector); [tempcell_SDD{:}] = deal( Hvar(1:ND, 1:ND) );
tempcell_SDI = cell(1, Qs_sector); [tempcell_SDI{:}] = deal( Hvar(1:ND, (ND+1):end) );
tempcell_ID = cell(1, Qs_sector); [tempcell_ID{:}] = deal( Hvar((ND+1):end, 1:ND) );
tempcell_II = cell(1, Qs_sector); [tempcell_II{:}] = deal( Hvar((ND+1):end, (ND+1):end) );
Q_SDD = P_m2D'*blkdiag(tempcell_SDD{:})*P_m2D; clear tempcell_SDD;
Q_SDI = P_m2D'*blkdiag(tempcell_SDI{:})*L_s; clear tempcell_SDI;
Q_ID = L_s'*blkdiag(tempcell_ID{:})*P_m2D; clear tempcell_ID;
Q_II = L_s'*blkdiag(tempcell_II{:})*L_s; clear tempcell_II;

Xslave = zeros(Qs_sector*(Np_slave+Nu_slave), Nsamples);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sc = Sc + [Q_SDD sparse(sim.Np, Nu) Q_SDI sparse(sim.Np, Ni_r);
    sparse(Nu, sim.Np + Nu + Ni_r + Ni_s);
    Q_ID sparse(Ni_s,Nu) Q_II];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulating
Xsamples = zeros(Ntot, Nsamples);

%setting initial condition
if isfield(sim.results, 'X0') && ~isempty( sim.results.X0 )
    Xsamples(:, 1) = 0*sim.results.X0;
    %Xsamples(:, 1) = sim.results.Xh(1:Ntot);
else
    Xsamples(:, 1) = 0*sim.results.Xh(1:Ntot);
    %error('Initial conditions not computed.')
end

%initializing previous residual term
%res_prev = sim.results.res_prev;
res_prev = zeros(size(Xsamples,1),1);


for kt = 2:Nsamples
    disp(['Time step ' num2str(kt) '...']);
    
    S_ag = sim.msh.get_AGmatrix(wm*tsamples(kt), Ntot);
    Qconst = S_ag + Sc + (2/alpha2)*Mtot;
    
    %load function
    if Nin == 1
        Ustep = Ufun(tsamples(kt));
    elseif Nin == 3
        Ustep = Ufun(Mphase'*Xsamples(indI, kt-1), wm*tsamples(kt), tsamples(kt));
    end
    FL = (2/alpha2)*Mtot*Xsamples(:,kt-1) + [sim.matrices.F; zeros(Nu,1); Ustep(1:sim.results.Ni_s)];
    
    % updating the decomposed-domain contribution to F
    FD = zeros(Ntot, 1);
    xslave_prev = (2/(dt*alpha2))*M_slave*reshape(Xslave(:,kt-1), [], Qs_sector);
    hslave = Q_slave(nfree,nfree) \ xslave_prev(nfree,:);
    FD(indA) = FD(indA) - P_m2D'*reshape(P_D2s'*Q_slave(nd, nfree)*hslave(1:numel(nfree),:), [], 1);
    FD(indI) = FD(indI) - L_s' * reshape(hslave((numel(np_free)+1):end,:), [], 1);
    
    if kt == 3
        sim.misc.Fa = P_m2D'*reshape(P_D2s'*Q_slave(nd, nfree)*hslave(1:numel(nfree),:), [], 1);
        sim.misc.Fi = L_s' * reshape(hslave((numel(np_free)+1):end,:), [], 1);
    end
    
    FL = FL + FD;
    
    Xsamples(:,kt) = Xsamples(:,kt-1); %initial condition for NR
    for kiter = 1:50
        [J, res] = Jc.eval(Xsamples(:,kt), sim.nu_fun); 
        
        res_tot = PT'*(res + Qconst*Xsamples(:,kt) - FL - (alpha1/alpha2)*res_prev);
        
        %sim.misc.J = J;
        %sim.misc.res = res;
        
        resNorm = norm(res_tot) / norm(FL);
        disp(['    Newton step ' num2str(kiter) ', relative residual ' num2str(resNorm) '.']);
        %disp(['    Newton step ' num2str(kiter) ', relative residual ' num2str(resNorm) '.']); fflush(stdout); uncomment in Octave
        if resNorm < 1e-6
            break;
        end
        
        Jtot = PT' * (J + Qconst) * PT;
        
        Xsamples(:,kt) = Xsamples(:,kt) - PT*( Jtot\res_tot );
    end
    
    %updating prev-residual term
    res_prev = -res - (S_ag + Sc)*Xsamples(:,kt) + [sim.matrices.F; zeros(Nu,1); Ustep(1:sim.results.Ni_s)] ...
        + FD;
    
    %updating slave-domain solution    
    xslav = XX*[reshape(P_m2D*Xsamples(indA, kt),[],Qs_sector); reshape(L_s*Xsamples(indI, kt),[],Qs_sector)];
    xslav(nfree,:) = xslav(nfree,:) + hslave;
    Xslave(:,kt) = reshape(xslav, [], 1);
    
    %plotting currents
    %%{
    Mphase_plot = kron(eye(N_phases), ones(N_inParallel,1))';
    Is = Xsamples(indI(:), 1:kt);
    Iphase = Mphase_plot*Is;
    figure(13); clf; hold on;
    plot( Iphase(1:3,:)' ); 
    plot( Is(1:N_inParallel,:)', 'b--');
    drawnow;
    %}
    
    %{
    figure(13); clf;
    Xtemp = reshape(Xslave(:,kt), [], Qs_sector);
    %msh_trimesh( sim.msh.misc.msh_slave, Xtemp(1:Np_slave, 3), []);
    msh_trimesh(sim.msh, Xsamples(:,kt), []);
    figure(14); clf;
    plot(Xtemp(:,3))
    %}
end

sim.results.Xt = Xsamples;


end