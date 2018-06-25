function sim = sim_computeSlaveDomainSolution_timestepping(sim, pars)
%sim_computeSlaveDomainSolution_timestepping Compute impulse-response
%solution for the slave domain.
%
% (c) 2018 Antti Lehikoinen / Aalto University

msh = sim.msh.misc.msh_slave;
nd = sim.msh.misc.nd_slave; %Dirichlet nodes
Nu = numel(sim.msh.misc.conductors_slave);

Np = size(msh.p,2);
nfree = setdiff(1:Np, nd);
Np_free = numel(nfree);

ind_freeVars = [nfree (Np+1):(Np+Nu)]; 

dims = sim.dims;

%slave-domain matrices
S = sim.matrices.S_slave;
M = sim.matrices.M_slave;
C = sim.matrices.C_slave;
DR = sim.matrices.DR_slave;
P_D2s = sim.msh.misc.P_D2s;
P_m2D = sim.msh.misc.P_m2D;
ND = size(P_D2s, 2);

%for extracting Lagrange multipliers from the solution
Nh = numel(nd); %number of Lagrance multipliers
L = sparse(1:Nh, nd, ones(1, Nh), Nh, Np);

%time-step
tsamples = pars.ts;
dt_base = tsamples(2) - tsamples(1);

%assembling final matrices
Qfree_slave = [S(nfree, nfree)+1/dt_base*M(nfree, nfree) -1/dims.leff*C(nfree,:);
    1/dt_base*DR*transpose(C(nfree,:)) -speye(Nu, Nu)];
Mprev_slave = 1/dt_base*[M(nfree, nfree) sparse(Np_free, Nu);
    DR*transpose(C(nfree,:)) sparse(Nu, Nu)];
Sbnd = [S(nfree, nd);
    sparse(Nu, numel(nd))];
Mbnd = 1/dt_base*[M(nfree, nd);
    DR*transpose(C(nd,:))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% impulse response computation

%setting source function (unnecessarily complex for discrete impulse
%response, but whatevs
Nsteps = 50; %steps for impulse response
Npeak = 3;
d_Tpeak = dt_base;

%T_decay = 1/f * 1.5;
T_decay = dt_base*Nsteps;
N_decay = ceil( (T_decay-d_Tpeak)/d_Tpeak ) + 1;
ts_peak = linspace(-d_Tpeak, d_Tpeak, Npeak); f_peak = (0.5 + 0.5*sawtooth(pi*ts_peak/d_Tpeak-pi, 0.5));
ts_imp = [ts_peak linspace(ts_peak(end) + T_decay/N_decay, ts_peak(end)+T_decay, N_decay)];
fs_imp = [f_peak zeros(1, N_decay)];

N_imp = numel(fs_imp);

X_imp = zeros(N_imp, (Np+Nu)*(ND + Nu));
H_imp = zeros(N_imp, (ND+Nu)*(ND + Nu));
Xprev = zeros(Np + Nu, ND + Nu); Xnew = Xprev;

H_imp2 = zeros(N_imp, (Nh+Nu)*(ND + Nu));

%stepping
for kt = 2:N_imp
    %setting source
    if kt > Npeak
        %zero excitation
        FL = Mprev_slave*Xprev(ind_freeVars,:);
        Xnew(ind_freeVars, :) = Qfree_slave \ FL;
    else
        %changing boundary/current
        Xnew(nd, 1:ND) = P_D2s*fs_imp(kt);

        bndTerm = -(Sbnd + Mbnd) * Xnew(nd, :) + ...
            Mbnd * Xprev(nd, :);

        Xnew(ind_freeVars, :) = Qfree_slave \ ...
            (Mprev_slave * Xprev(ind_freeVars,:) + bndTerm + ...
                blkdiag(sparse(Np_free, ND), -DR*fs_imp(kt)) );
    end
    
    
    if kt == 2
        sim.misc.temp = Xnew;
    end
    
    %extracting Lagrance multipliers and voltages
    h_imp_temp = [-P_D2s'*L*( S*Xnew(1:Np,:) + M*(Xnew(1:Np,:)-Xprev(1:Np,:))/dt_base ...
        - 1/dims.leff*C*Xnew((Np+1):end,:) );
        Xnew((Np+1):end,:)];
    H_imp(kt,:) = reshape(h_imp_temp, 1, []);
    
    h_imp_temp = [-L*( S*Xnew(1:Np,:) + M*(Xnew(1:Np,:)-Xprev(1:Np,:))/dt_base ...
        - 1/dims.leff*C*Xnew((Np+1):end,:) );
        Xnew((Np+1):end,:)];
    H_imp2(kt,:) = reshape(h_imp_temp, 1, []);

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

H_imp = H_imp(2:kt,:);
H_var = reshape(H_imp(1,:), ND + Nu, []);
sim.misc.H_imp = H_imp;
sim.misc.H_var = H_var;

sim.misc.H_imp2 = H_imp2(2:kt,:);
sim.misc.X_imp_T = X_imp(2:kt,:);

%assigning the constant reduced matrices
Ximp_temp = reshape(X_imp(1,:), Np + Nu, []);
sim.misc.X_AA = Ximp_temp(1:Np, 1:ND); 
sim.misc.X_uA = Ximp_temp((Np+1):end, 1:ND);
sim.misc.X_AI = Ximp_temp(1:Np, (ND+1):end); 
sim.misc.X_uI = Ximp_temp((Np+1):end, (ND+1):end);

sim.misc.RS_AA_T = sim.misc.X_AA' * ( S*sim.misc.X_AA - 1/dims.leff*C*sim.misc.X_uA );    
sim.misc.RS_AI_T = sim.misc.X_AA'*( S*sim.misc.X_AI - 1/dims.leff*C*sim.misc.X_uI );
sim.misc.RM_AA_T = sim.misc.X_AA' * M*sim.misc.X_AA;    
sim.misc.RM_AI_T = sim.misc.X_AA' * M*sim.misc.X_AI;
sim.misc.R_UA_T = sim.misc.X_uA; 
sim.misc.R_UI_T = sim.misc.X_uI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decay term computation

Xh = sim.results.Xh;
Xh = reshape(Xh, [], 2);
w = 2*pi*pars.f;
L_s = sim.matrices.Ls;

indA = 1:sim.Np;
indI = (sim.Np + sim.results.Nu_r) + (1:sim.results.Ni_s);
Qs_sector = sim.dims.Qs / sim.msh.symmetrySectors;

D0 = P_m2D*Xh(indA,:); D0 = reshape(D0(:,1) + 1i*D0(:,2), ND, []);
I0 = Xh(indI,:); I0 = I0(:,1) + 1i*I0(:,2);

Qfree_slave0 = [S(nfree, nfree)+1i*w*M(nfree, nfree) -1/dims.leff*C(nfree,:);
    1i*w*DR*transpose(C(nfree,:)) -speye(Nu, Nu)];
Qbnd_slave0 = [S(nfree, nd) + 1i*w*M(nfree, nd);
    1i*w*DR*transpose(C(nd,:))];

%initial conditions
DX0 = zeros(Np+Nu, Qs_sector);

DX0(ind_freeVars, :) = Qfree_slave0 \ ( -Qbnd_slave0*P_D2s*D0 + [zeros(Np_free, Qs_sector); 
    -DR*reshape(L_s*I0, Nu, [])] );
DX0(nd, :) = P_D2s * D0;

% time-stepping
N_decay = 50;
h_decay = zeros(Qs_sector*(ND+Nu), N_decay);
X_decay = zeros(Qs_sector*(Np+Nu), N_decay);

DXprev = real( DX0 );

%computing initial values
DDX = imag(DX0);
h_imp_temp = [-P_D2s'*L*(S*DXprev(1:Np,:) + w*M*DDX(1:Np,:) ...
    - 1/dims.leff*C*DXprev((Np+1):end,:));
    DXprev((Np+1):end,:)];
h_decay(:, 1) = reshape(h_imp_temp, [], 1);

h_decay2 = zeros(Qs_sector*(Nh+Nu), N_decay);


DXnew = zeros( size( DXprev ) );
for kt = 2:N_decay
    if kt < 2
        %zero boundary
        DFL = Mprev_slave*DXprev(ind_freeVars,:);
        DXnew(ind_freeVars_slave, :) = Qfree_slave \ DFL;
    else
        %boundary values dropping back to zero
        bndTerm = -(Sbnd + Mbnd) * DXnew(nd, :) + ...
            Mbnd * DXprev(nd, :);

        DXnew(ind_freeVars, :) = Qfree_slave \ ...
            (Mprev_slave * DXprev(ind_freeVars,:) + bndTerm  );
    end
    %extracting Lagrance multipliers and voltages
    h_imp_temp = [-P_D2s'*L*(S*DXnew(1:Np,:) + M*(DXnew(1:Np,:)-DXprev(1:Np,:))/dt_base ...
        - 1/dims.leff*C*DXnew((Np+1):end,:));
        DXnew((Np+1):end,:)];
    h_decay(:, kt) = reshape(h_imp_temp, [], 1);
    X_decay(:, kt) = reshape(DXnew, [], 1);
    
        %extracting Lagrance multipliers and voltages
    h_imp_temp = [-L*(S*DXnew(1:Np,:) + M*(DXnew(1:Np,:)-DXprev(1:Np,:))/dt_base ...
        - 1/dims.leff*C*DXnew((Np+1):end,:));
        DXnew((Np+1):end,:)];
    h_decay2(:, kt) = reshape(h_imp_temp, [], 1);
    
    DXprev = DXnew;    

    if kt == 2
        % maxes for checking convergence
        temp_maxes_comp = max(abs(DXnew), [], 1);            
    end
    if kt >= 2
        disp( max( max(abs(DXnew), [], 1)./temp_maxes_comp ) )

        if (max( max(abs(DXnew), [], 1)./temp_maxes_comp)) < 1e-14
            break;
        end
    end
end
sim.misc.h_decay = h_decay(:, 1:kt);
sim.misc.X_decay = X_decay(:, 1:kt);
sim.misc.h_decay2 = h_decay2(:, 1:kt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial conditions computation for CN-stepping


end