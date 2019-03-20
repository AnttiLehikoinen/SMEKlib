function sim = sim_runTimeSteppingSimulation_CN(sim, pars)
%sim_runTimeSteppingSimulation_CN Crank-Nicolson time-stepping
%
% (c) 2017 Antti Lehikoinen / Aalto University

f = pars.f;
w = 2*pi*f;

if ~isempty(pars.slip)
    slip = pars.slip;
else
    slip = sim.dims.slip;
end
angle0 = pars.rotorAngle;

%voltage-function
Nin = 1;
if isa(pars.U, 'function_handle')
    Ufun = pars.U;
    Nin = nargin(Ufun);
else
    U = pars.U / sim.msh.symmetrySectors * sim.dims.a * sqrt(2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %voltage function
    phi0 = pars.phi0;
    if sim.dims.connection_stator == defs.delta
        Ufun = @(t)(  U*[cos(w*t-phi0); cos(w*t - 2*pi/3-phi0); cos(w*t - 4*pi/3-phi0)] );
    else
        Ufun = @(t)(  U*[cos(w*t-phi0); cos(w*t - pi/3-phi0)] );
    end
end

%{
dt = 1/f / pars.N_stepsPerPeriod; %time-step length
Nsamples = ceil( (pars.N_periods/f) / dt );
tsamples = (0:(Nsamples-1))*dt;
%}

wm = w/sim.dims.p * (1-slip);
tsamples = pars.ts;
Nsamples = numel(tsamples);
dt = tsamples(2) - tsamples(1);


%[Sc, Mtot] = get_circuitMatrices(sim);
[Sc, Mtot] = get_circuitMatrices_2(sim);
Mtot = Mtot/dt;

Ntot = size(Mtot, 1);
Nui = Ntot - size(sim.matrices.P,1);
Nu = sim.results.Nu_s + sim.results.Nu_r;
PT = blkdiag(sim.matrices.P, speye(Nui));
indI = (1:sim.results.Ni_s) + sim.Np+Nu;

%Jacobian constructor object
Jc = JacobianConstructor(sim.msh, Nodal2D(Operators.curl), Nodal2D(Operators.curl), true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulating
Xsamples = zeros(size(Sc,1), Nsamples);

%setting initial condition
if isfield(sim.results, 'X0') && ~isempty( sim.results.X0 )
    Xsamples(:, 1) = sim.results.X0;
else
    error('Initial conditions not computed.')
end

% adjusted CN for stability
%alpha2 = 1.1; %weight for implicit (k+1) step; 1 for CN, 2 for BE
alpha2 = pars.alpha2;
%alpha2 = 2
alpha1 = 2 - alpha2;


%initializing previous residual term
[~, res_prev] = Jc.eval(Xsamples(:, 1), sim.nu_fun);
if Nin == 1
    Ustep = Ufun(0);
elseif Nin == 3
    Ustep = Ufun(Xsamples(indI, 1), 0, 0);
end
if size(Ustep,1) == sim.results.Ni_s
    Ustep = [Ustep; zeros(sim.results.Ni_r,size(Ustep,2))];
end
res_prev = -res_prev - (sim.msh.get_AGmatrix(angle0, Ntot) + Sc)*Xsamples(:,1) + [sim.matrices.F; zeros(Nu, 1); Ustep];

for kt = 2:Nsamples
    disp(['Time step ' num2str(kt) '...']);
    
    %S_ag = get_MovingBandMatrix(wm*tsamples(kt), sim.msh, Ntot);
    S_ag = sim.msh.get_AGmatrix(angle0 + wm*tsamples(kt), Ntot);
    
    Qconst = S_ag + Sc + (2/alpha2)*Mtot;
    
    if Nin == 1
        Ustep = Ufun(tsamples(kt));
    elseif Nin == 3
        Ustep = Ufun(Xsamples(indI, kt-1), wm*tsamples(kt), tsamples(kt));
    end
    if size(Ustep,1) == sim.results.Ni_s
        Ustep = [Ustep; zeros(sim.results.Ni_r,size(Ustep,2))];
    end
        
    FL = (2/alpha2)*Mtot*Xsamples(:,kt-1) + [sim.matrices.F; zeros(Nu, 1); Ustep];
    
    Xsamples(:,kt) = Xsamples(:,kt-1); %initial condition for NR
    for kiter = 1:pars.maxIter
        %[J, res] = assemble_Jacobian(sim.nu_fun, Xsamples(:,kt), [], sim.msh);
        [J, res] = Jc.eval(Xsamples(:,kt), sim.nu_fun); 
        
        res_tot = PT'*(res + Qconst*Xsamples(:,kt) - FL - (alpha1/alpha2)*res_prev);
        
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
    res_prev = -res - (S_ag + Sc)*Xsamples(:,kt) + ...
        [sim.matrices.F; zeros(Nu, 1); Ustep];
end

sim.results.Xt = Xsamples;

end