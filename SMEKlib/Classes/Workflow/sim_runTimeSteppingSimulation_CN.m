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

%voltage-function
if isa(pars.U, 'function_handle')
    Ufun = pars.U;
else
    U = pars.U / sim.msh.symmetrySectors * sim.dims.a * sqrt(2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %voltage function
    phi0 = 0;
    if sim.dims.connection_stator == defs.delta
        Ufun = @(t)(  U*[cos(w*t-phi0); cos(w*t - 2*pi/3-phi0); cos(w*t - 4*pi/3-phi0)] );
    else
        Ufun = @(t)(  U*[cos(w*t-phi0); cos(w*t - pi/3-phi0)] );
    end
end

dt = 1/f / pars.N_stepsPerPeriod; %time-step length
wm = w/sim.dims.p * (1-slip);

Nsamples = ceil( (pars.N_periods/f) / dt );
tsamples = (0:(Nsamples-1))*dt;


[Sc, Mtot] = get_circuitMatrices(sim);
Mtot = Mtot/dt;

Ntot = size(Mtot, 1);
Nui = Ntot - size(sim.matrices.P,1);
PT = blkdiag(sim.matrices.P, speye(Nui));

%Jacobian constructor object
Jc = JacobianConstructor(sim.msh, Nodal2D(Operators.curl), Nodal2D(Operators.curl), false);

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
alpha2 = 1.1; %weight for implicit (k+1) step; 1 for CN, 2 for BE
alpha1 = 2 - alpha2;


%initializing previous residual term
[~, res_prev] = Jc.eval(Xsamples(:, 1), sim.nu_fun);
res_prev = -res_prev - (sim.msh.get_AGmatrix(0, Ntot) + Sc)*Xsamples(:,1) + [zeros(Ntot-sim.results.Ni_s, 1); Ufun(0)];

for kt = 2:Nsamples
    disp(['Time step ' num2str(kt) '...']);
    
    %S_ag = get_MovingBandMatrix(wm*tsamples(kt), sim.msh, Ntot);
    S_ag = sim.msh.get_AGmatrix(wm*tsamples(kt), Ntot);
    
    Qconst = S_ag + Sc + (2/alpha2)*Mtot;
    
    Ustep = Ufun(tsamples(kt));
    FL = (2/alpha2)*Mtot*Xsamples(:,kt-1) + [zeros(Ntot-sim.results.Ni_s, 1); Ustep(1:sim.results.Ni_s)];
    
    Xsamples(:,kt) = Xsamples(:,kt-1); %initial condition for NR
    for kiter = 1:50
        %[J, res] = assemble_Jacobian(sim.nu_fun, Xsamples(:,kt), [], sim.msh);
        [J, res] = Jc.eval(Xsamples(:,kt), sim.nu_fun); 
        
        res_tot = PT'*(res + Qconst*Xsamples(:,kt) - FL - (alpha1/alpha2)*res_prev);
        
        resNorm = norm(res_tot) / norm(FL);
        disp([9 'Newton step ' num2str(kiter) ', relative residual ' num2str(resNorm) '.']);
        %disp(['    Newton step ' num2str(kiter) ', relative residual ' num2str(resNorm) '.']); fflush(stdout); uncomment in Octave
        if resNorm < 1e-6
            break;
        end
        
        Jtot = PT' * (J + Qconst) * PT;
        
        Xsamples(:,kt) = Xsamples(:,kt) - PT*( Jtot\res_tot );
    end
    
    %updating prev-residual term
    res_prev = -res - (S_ag + Sc)*Xsamples(:,kt) + ...
        [zeros(Ntot-sim.results.Ni_s, 1); Ustep(1:sim.results.Ni_s)];
end

sim.results.Xt = Xsamples;

end