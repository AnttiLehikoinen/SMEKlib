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

%time-step
tsamples = pars.ts;
dt_base = tsamples(2) - tsamples(1);

%setting source function (unnecessarily complex for discrete impulse
%response, but whatevs
Nsteps = 500; %steps for impulse response
Npeak = 3;
d_Tpeak = dt_base;

%T_decay = 1/f * 1.5;
T_decay = dt_base*Nsteps;
N_decay = ceil( (T_decay-d_Tpeak)/d_Tpeak ) + 1;
ts_peak = linspace(-d_Tpeak, d_Tpeak, Npeak); f_peak = (0.5 + 0.5*sawtooth(pi*ts_peak/d_Tpeak-pi, 0.5));
ts_imp = [ts_peak linspace(ts_peak(end) + T_decay/N_decay, ts_peak(end)+T_decay, N_decay)];
fs_imp = [f_peak zeros(1, N_decay)];

X_imp = zeros(N_imp, (Np_slave+Nu_slave)*(N_Dir + Nu_slave));
Xprev = zeros(Np_slave + Nu_slave, N_Dir + Nu_slave); Xnew = Xprev;

%stepping
for kt = 2:N_imp
    %setting source
    if kt > Npeak
        %zero excitation
        FL = Mprev_slave*Xprev(ind_freeVars_slave,:);
        Xnew(ind_freeVars_slave, :) = Qfree_slave \ FL;
    else
        %changing boundary/current
        Xnew(nd_slave, 1:N_Dir) = P_c2s*fs_imp(kt);

        bndTerm = -(Sbnd_slave + Mbnd_slave) * Xnew(nd_slave, :) + ...
            Mbnd_slave * Xprev(nd_slave, :);

        Xnew(ind_freeVars_slave, :) = Qfree_slave \ ...
            (Mprev_slave * Xprev(ind_freeVars_slave,:) + bndTerm + ...
                blkdiag(sparse(Np_free_slave, N_Dir), -DR_slave*fs_imp(kt)) );
    end

    %reshaping results to column-major
    X_imp(kt,:) = reshape(Xnew, 1, []);

    if kt == (Npeak+1)
        % maxes for plotting
        %temp_maxes_plot = max(abs(Xnew), [], 1);
        %axis([1 N_imp 0 1]);
    elseif kt == ceil(Npeak/2)
        % maxes for checking convergence
        temp_maxes_comp = max(abs(Xnew), [], 1);            
    end
    if kt >= (Npeak+1)
        %plot(kt, max(abs(Xnew), [], 1)./temp_maxes_plot, 'k.'); pause(0.05);
        %disp( max( max(abs(Xnew), [], 1)./temp_maxes_comp ) )

        if (max( max(abs(Xnew), [], 1)./temp_maxes_comp)) < 1e-14
            break;
        end
    end

    Xprev = Xnew;
end
if (kt == N_imp)
    warning('Impulse response did not decay sufficiently.');
end

sim.misc.X_imp_T = X_imp(2:kt,:);

%assigning the constant reduced matrices
Ximp_temp = reshape(X_imp(1,:), Np_slave + Nu_slave, []);
sim.misc.X_AA = Ximp_temp(1:Np_slave, 1:N_Dir); 
sim.misc.X_uA = Ximp_temp((Np_slave+1):end, 1:N_Dir);
sim.misc.X_AI = Ximp_temp(1:Np_slave, (N_Dir+1):end); 
sim.misc.X_uI = Ximp_temp((Np_slave+1):end, (N_Dir+1):end);

sim.misc.RS_AA_T = X_AA' * ( (S_slave + 1/dt*M_slave)*X_AA - 1/dims.leff*CF_slave*X_uA );    
sim.misc.RS_AI_T = X_AA'*( (S_slave + 1/dt*M_slave)*X_AI - 1/dims.leff*CF_slave*X_uI );
sim.misc.R_UA_T = X_uA; 
sim.misc.R_UI_T = X_uI;

end