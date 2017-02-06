%Time-stepping simulation, continued from the time-harmonic solution.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

dt = 1/f / 200; %time-step length
wm = w/2 * (1-slip);

Nsamples = ceil( 1/f / dt );
tsamples = (0:(Nsamples-1))*dt;

%updating necessary matrices (slip no longer as multiplier)
M_struct = assemble_matrix('', 'nodal', '', 'nodal', dims.sigma_rotor, horzcat(rotorConductors{1,:}), msh, []);
M_struct = assemble_matrix('', 'nodal', '', 'nodal', 6e6, shaft, msh, M_struct); %shaft
M = sparseFinalize(M_struct, Np, Np);

MJ = [-1/dims.leff*JF_r -JF_s*Ls]; %DC-current density matrix
MU = [1/dt*DRr*transpose(JF_r) real(Mur) sparse(Nu_r, Ni_s)]; %voltage equations
ME = [1/dt*dims.leff*Ls'*JF_s' sparse(Ni_s, Nu_r)]; %circuit equations

Mtemp = [MU;ME Z];

%pure mass-matrix
Mprev = 1/dt*[M sparse(Np, Nu_r+Ni_s);
    DRr*transpose(JF_r) sparse(Nu_r, Nu_r+Ni_s);
    dims.leff*Ls'*JF_s' sparse(Ni_s, Nu_r+Ni_s)];

%voltage function
phi0 = 0;
Ufun = @(t)( kron(eye(3), ones(1,1)) * U*[cos(w*t-phi0); cos(w*t - 2*pi/3-phi0); cos(w*t - 4*pi/3-phi0)] );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time-stepping

S_ag = get_MovingBandMatrix(0, msh);

Xsamples = zeros(Nvars, numel(tsamples)); 
Xsamples(:, 1) = Xtot(1:Nvars); 

tic;
Xfree_T = Xprev(1:(size(Xprev,1)/2));
for kt = 2:Nsamples
    disp(['Time step ' num2str(kt) '...']);
    
    S_ag = get_MovingBandMatrix(wm*tsamples(kt), msh);
    
    Qconst = [S_ag+M/dt MJ;
        Mtemp];
    
    FL = Mprev*Xsamples(:,kt-1) + [sparse(Np+Nu, 1); Ufun(tsamples(kt))];
    
    for kiter = 1:15
        [J, res] = assemble_Jacobian(nu_fun, PT*Xfree_T, [], msh);
        
        res_tot = PT'*(res + Qconst*(PT*Xfree_T) - FL);
        
        resNorm = norm(res_tot) / norm(FL);
        disp([9 'Newton step ' num2str(kiter) ', relative residual ' num2str(resNorm) '.']);
        %disp(['    Newton step ' num2str(kiter) ', relative residual ' num2str(resNorm) '.']); fflush(stdout); uncomment in Octave
        if resNorm < 1e-6
            break;
        end
        
        Jtot = PT' * (J + Qconst) * PT;
        
        Xfree_T = Xfree_T - Jtot\res_tot;
    end
    
    Xsamples(:,kt) = PT * Xfree_T;
end
toc


%plotting currents (and the currents from harmonic approximation
figure(5); clf; hold on; box on;
Iharm = bsxfun(@times, Xtot(indI), cos(w*tsamples)) + bsxfun(@times, Xtot(Nvars + indI), sin(w*tsamples));
h_td = plot( tsamples*1e3, 2*Xsamples(indI, :)' );
h_h = plot( tsamples*1e3, 2*Iharm', 'Linestyle', '--');
legend([h_td(1) h_h(1)], 'Time-stepping', 'Harmonic approximation');
title('Phase currents')
xlabel('Time (ms)');
ylabel('Current (A)');


%plotting line currents
figure(6); clf; hold on; box on;
Mline = [1 0 -1;-1 1 0;0 -1 1];
plot( tsamples*1e3, (Mline*Xsamples(indI, :))' )
title('Line currents');
xlabel('Time (ms)');
ylabel('Current (A)');

%ANOTHER AWESOME MIND-BOGGLING COLORDREAM
figure(7); clf; box on;
kt = 50;

drawFluxDensity(msh, Xsamples(indA, kt), wm*tsamples(kt), 'LineStyle', 'none'); 
colormap('jet'); colorbar; caxis([0 2]);
drawFluxLines(msh, Xsamples(indA, kt), 16, wm*tsamples(kt), 'k');
axis(dims.D_so/2*[-1 1 0 1]); axis tight; daspect([1 1 1]); 

title(['Flux lines and flux density (T) at t = ' num2str(1e3*tsamples(kt)) ' ms']);