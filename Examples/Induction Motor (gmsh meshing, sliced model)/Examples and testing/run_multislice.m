%Running a multi-slice simulation.
%
% Not a standard feature of SMEKlib yet (Sep 2018), so very much manual
% labor.
%
% Some info:
% The AuI formulation is used, see e.g. the SMEKlib publication referred in
% the Readme. Each slice has its own vector potential (a) and voltage (u)
% variable. The same currents (i) travel trough all slices. Rotor currents
% are not eliminated, but kept as independent variables.
%
% Consequently, the Au part of each slice k is modelled by the system of eqs
%   (S+M*d/dt)a_k + C_J*u_k = 0
%   C_E*d/dt a_k - u_k = D_R*L*i.
% The slices are then coupled together by the voltage equation
%   (Sum_k L^T*u_k) + R_ew*i = U_supply.
%
% For laziness reasons, each slice is modelled as equally long to the full
% machine, i.e. dimsc.leff. The final results are then divided by the
% number of slices when appropriate.
%
% Initial conditions are obtained from the single-slice simulation run in
% setup.m
%
% (c) 2018 Antti Lehikoinen / Smeklab Ltd

%%{
Nslices = 5; %number of slices
angle_skew = 0.8 * 2*pi/dim.Qs;
skew_angles = linspace(-angle_skew/2, angle_skew/2, Nslices);

%setting parameters
w = 2*pi*pars.f;
wm = w/sim.dims.p * (1-pars.slip);
tsamples = pars.ts;
Nsamples = numel(tsamples);
dt = tsamples(2) - tsamples(1);

%non-FE stiffness and mass matrices (circuit matrices)
[Sc, Mc] = get_circuitMatrices_2(sim); Ni_r = size(sim.matrices.Lr, 2);
Mc = Mc/dt;

%numbers of variables, indices
Np = sim.Np;
Ntot = size(Mc, 1);
Nui = Ntot - size(sim.matrices.P,1);
Nu = sim.results.Nu_s + sim.results.Nu_r;
Ni = sim.results.Ni_s + Ni_r;
indI = (1:Ni) + sim.Np+Nu; %indices to currents
indAU = 1:(indI(1) - 1); %indices to vector pot. and voltages

%voltage supply; ugly hack for considering the larger number of slices
Uampl = pars.U / sim.msh.symmetrySectors * sim.dims.a * sqrt(2);
Ufun = @(t)(  Nslices*Uampl*[cos(w*t); cos(w*t - 2*pi/3); cos(w*t - 4*pi/3); zeros(Ni_r,1)] );

%constant non-ag stiffness matrix and mass matrix
Scell = cell(Nslices+1); [Scell{:}] = deal(sparse(Np+Nu, Np+Nu));
Mcell = cell(Nslices+1); [Mcell{:}] = deal(sparse(Np+Nu, Np+Nu));
Sagcell = cell(Nslices+1, 1); Jcell = cell(Nslices+1,1);
for kslice = 1:Nslices
    Scell{kslice,kslice} = Sc(indAU, indAU) ;
    Scell{kslice,end} = Sc(indAU, indI) ;
    Scell{end, kslice} = Sc(indI, indAU) ;
    
    Mcell{kslice,kslice} = Mc(indAU, indAU) ;
    Mcell{kslice,end} = Mc(indAU, indI) ;
    Mcell{end, kslice} = Mc(indI, indAU) ;
end
Scell{end,end} = Sc(indI,indI) * Nslices; 
Mcell{end,end} = Mc(indI,indI) * Nslices;
Sagcell{end} = sparse(Ni,Ni);
Jcell{end} = sparse(Ni,Ni);
Stot = cell2mat(Scell); 
Mtot = cell2mat(Mcell);

Jc = JacobianConstructor(sim.msh, Nodal2D(Operators.curl), Nodal2D(Operators.curl), false);

%parameters for Crank-Nicolson stepping
alpha2 = pars.alpha2;
alpha1 = 2 - alpha2;

%boundary condition matrix
P_AU = blkdiag(sim.matrices.P, eye(Nu));
Ptot = blkdiag(...
    kron(eye(Nslices), P_AU), eye(Ni));

%allocating vectors for unknowns
Xslice = zeros((Np+Nu)*Nslices + Ni, Nsamples);

%initial conditions
X0 = sim.results.X0;
Xslice(:, 1) = [repmat(X0(1:(Np+Nu)), Nslices, 1);
    X0((Np+Nu+1):(Np+Nu+Ni))];

FL = zeros((Np+Nu)*Nslices + Ni, 1);
Ntot_all = (Np+Nu)*Nslices + Ni;

%for plotting
indsI = ((Np+Nu)*Nslices+1):((Np+Nu)*Nslices + sim.results.Ni_s);
Iref = sim.Is;
tic
for kt = 1:Nsamples
    disp(['Time step ' num2str(kt) '...']);
    
    disp('Assembling slice matrices');
    for kslice = 1:Nslices
        %S_ag = get_MovingBandMatrix(wm*tsamples(kt) + skew_angles(kslice), msh, Np+Nu);
        S_ag = sim.msh.get_AGmatrix(wm*tsamples(kt)+ skew_angles(kslice), Np+Nu);
        Sagcell{kslice} = S_ag;
    end
    
    %first step? Initializing previous residual term for Crank-Nicolson
    %stepping
    if kt == 1
        res_t = zeros(Ntot_all, 1);
        for kslice = 1:Nslices
            %initializing previous residual term
            inds_slice = (1:(Np+Nu)) + (kslice-1)*(Np+Nu);
            [J, res] = Jc.eval(Xslice(inds_slice,kt), sim.nu_fun);
            res_t(inds_slice) = res;
        end
        Sag_here = blkdiag(kron(eye(Nslices), sim.msh.get_AGmatrix(0, Np+Nu)), zeros(Ni,Ni));
        res_prev = -res_t - (Stot + Sag_here)*Xslice(:,1) + [zeros((Np+Nu)*Nslices,1); Ufun(tsamples(kt))];
        continue;
    end
    
    FL = (2/alpha2)*Mtot*Xslice(:,kt-1) + [zeros((Np+Nu)*Nslices,1); Ufun(tsamples(kt))];
    Qconst = Stot + (2/alpha2)*Mtot + blkdiag(Sagcell{:});
    
    %Newton iteration
    Xslice(:,kt) = Xslice(:,kt-1); %initial condition
    for kiter = 1:15
        res_t = zeros(Ntot_all, 1);
        
        %assembling slice matrices
        for kslice = 1:Nslices
            inds_slice = (1:(Np+Nu)) + (kslice-1)*(Np+Nu);
            [J, res] = Jc.eval(Xslice(inds_slice,kt), sim.nu_fun); 
            Jcell{kslice} = J;
            res_t(inds_slice) = res;
        end
        
        %computing final residual
        res_tot = Ptot'*(res_t + Qconst*Xslice(:,kt) - FL - (alpha1/alpha2)*res_prev);
        
        %checking convergence
        resNorm = norm(res_tot) / norm(FL);
        disp([9 'Newton step ' num2str(kiter) ', relative residual ' num2str(resNorm) '.']);
        if resNorm < 1e-6
            break;
        end
        %solving
        Jtot = Ptot'*(blkdiag(Jcell{:}) + Qconst)*Ptot;        
        Xslice(:,kt) = Xslice(:,kt) - Ptot*( Jtot\res_tot );
    end
    
    %updating previous residual term
    res_prev = -res_t - (Stot + blkdiag(Sagcell{:}))*Xslice(:,kt) + ...
        [zeros((Np+Nu)*Nslices,1); Ufun(tsamples(kt))];
    
    figure(10); clf; hold on;
    plot( Xslice(indsI,1:kt)' );
    plot( Iref(:,1:kt)', 'Linestyle', '--');
    drawnow;
end
toc
Isliced = Xslice(indsI,:);
%}

%plotting currents (and the currents from harmonic approximation
figure(5); clf; hold on; box on;
hs = plot( tsamples(1:kt)*1e3, Iref(:,1:kt)', 'b' );
hm = plot( tsamples(1:kt)*1e3, Isliced(:,1:kt)', 'r' );
xlabel('Time (ms)')
ylabel('Phase current (A)')
%legend([hs(1) hm(1)], 'Single-slice', '11-slice')


figure(11); clf; hold on; box on;
kslice = 2; inds_slice = (1:(Np+Nu)) + (kslice-1)*(Np+Nu);
A = Xslice(inds_slice, 2);
drawFluxDensity(mshc, A, 'LineStyle', 'none'); colormap('jet'); colorbar; caxis([0 2]);
drawFluxLines(mshc, A, 16, 'k');
box on; axis tight; daspect([1 1 1]); drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computing torque (workaround)
Tsliced = zeros(1, Nsamples);
Xt_orig = sim.results.Xt; %saving original standard-sim solution
for kslice = 1:Nslices
    inds_slice = (1:(Np+Nu)) + (kslice-1)*(Np+Nu);
    sim.results.Xt = Xslice(inds_slice, :); %overwriting with a slice solution
    Tsliced = Tsliced + sim_compute_torque(sim, pars, 'stepping')/Nslices;
end
sim.results.Xt = Xt_orig; %reverting back to original

figure(6); clf; hold on; box on;
plot(tsamples(1:kt)*1e3, T(:,1:kt)', 'b' );
plot(tsamples(1:kt)*1e3, Tsliced(:,1:kt)', 'r' );
xlabel('Time (ms)');
ylabel('Torque (Nm)');
legend('1 slice', '5 slices');