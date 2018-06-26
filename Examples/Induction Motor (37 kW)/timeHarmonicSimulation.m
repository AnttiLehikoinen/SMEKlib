%Time-harmonic simulation using peak values
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

f = 50;
U = 400/sqrt(3)/2 *sqrt(2); %taking the second-turn side outside symmetry sector into account

w = 2*pi*f;

%assembling matrices
S_struct = assemble_matrix('grad', 'nodal', 'grad', 'nodal', nu_fun(0), [], msh, []);
S = sparseFinalize(S_struct, Np, Np);

%mass matrix for rotor
M_struct = assemble_matrix('', 'nodal', '', 'nodal', dims.slip*dims.sigma_rotor, horzcat(rotorConductors{1,:}), msh, []);
%M_struct = assemble_matrix('', 'nodal', '', 'nodal', dims.slip*6e6, shaft, msh, M_struct); %shaft
M = sparseFinalize(M_struct, Np, Np);

%rotor winding matrix
JF_struct = []; Nc_s = numel(statorConductors); Nc_r = numel(rotorConductors);
for k = 1:Nc_r
    JF_struct = assemble_vector('', 'nodal', dims.sigma_rotor, k, rotorConductors{k}, msh, JF_struct);
end
JF_r = sparseFinalize(JF_struct, Np, Nc_r);

%conductor areas for rotor
cAr = sum(JF_r*eye(Nc_r)/(dims.sigma_rotor),1);
DRr = sparsediag(dims.leff ./(dims.sigma_rotor*cAr)); %resistance matrix

%stator winding matrix
JF_struct = [];
for k = 1:Nc_s
    JF_struct = assemble_vector('', 'nodal', 1, k, statorConductors{k}, msh, JF_struct);
end
JF_s = sparseFinalize(JF_struct, Np, Nc_s);
cAs = sum(JF_s*eye(Nc_s),1);
DRs = sparsediag(dims.leff ./(dims.sigma_stator*cAs)) / dims.fillingFactor; %resistance matrix
JF_s = bsxfun(@times, JF_s, 1./cAs);

%number of voltages and currents
%(rotor currents eliminated)
Nu_s = size(Ls,1); Nu_r = size(Lr,1); Nu = Nu_r;
Ni_s = size(Ls,2); Ni_r = 0*size(Lr,2); Ni = Ni_s + Ni_r;
indA = (1:Np); indU = (Np+1):(Np+Nu); indI = (Np+Nu+1):(Np+Nu+Ni);
Nvars = Np + Nu + Ni;

%final winding matrices
MJ = [-1/dims.leff*JF_r -JF_s*Ls]; %DC-current density matrix

%eliminating rotor currents; assembling coefficient matrix of rotor bar
%voltages
Mur =  -speye(Nu_r)-DRr*Lr*(Zr\Lr');

MU = [DRr*1i*w*dims.slip*transpose(JF_r) Mur sparse(Nu_r, Ni_s)]; %voltage equations
ME = [dims.leff*1i*w*Ls'*JF_s' sparse(Ni_s, Nu_r)]; %circuit equations
Z = sparse(Ls'*DRs*Ls); %stator impedances

%total interpolation matrix for boundary conditions
PT = assemble_TotalMasterSlaveMatrix(Np + Nu + Ni, P_data, []);
PTT = blkdiag(PT, PT); %for harmonic analysis

S_ag = get_MovingBandMatrix(0, msh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonlinear time-harmonic simulation
%NOTE: due to nonlinearity, real and imaginary components have to be
%separated --> we have twice as large purely-real problem

%voltage vector
t0 = 0;
FI = kron(eye(3), ones(1,1)) * U*[exp(1i*w*t0); exp(1i*w*t0-1i*2*pi/3); exp(1i*w*t0-1i*4*pi/3)];

%assembling load vector
Ftemp = [zeros(Np+Nu,1); FI];
Ftot = [real(Ftemp); -imag(Ftemp)];

%assembling some constant matrices
Mtemp = [MU; ME Z]; %temporary matrix for shorter syntax
Sconst = [[S_ag real(MJ);real(Mtemp)] [w*M imag(MJ);imag(Mtemp)];
        -[w*M imag(MJ);imag(Mtemp)] [S_ag real(MJ);real(Mtemp)]];


% Newton iteration
Xprev = zeros( size(PTT, 2), 1); 
Xtot = PTT*Xprev;
for kiter = 1:15
    % assembling Jacobian and residual blocks
    [J11, J12, J21, J22, res11, res22] = assemble_ComplexJacobian(nu_fun, Xtot, [], msh);    
    
    %finalizing
    Jtot = PTT'*( [J11 J12; J21 J22] + Sconst )*PTT;    
    res_tot = PTT'*( Sconst*Xtot - Ftot + [res11;res22] );
    
    %checking convergence
    disp( norm(res_tot) / norm(Ftot) )
    if (norm(res_tot) / norm(Ftot)) < 1e-6
        break;
    end
    
    Xprev = Xprev - Jtot \ res_tot;
    Xtot = PTT*Xprev;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post-processing

%obtaining stator current phasor
Is = Xtot(indI) + 1i*Xtot(Nvars+indI);

%surface plot of vector potential
figure(3); clf;
msh_trimesh(msh, Xtot(indA), []);

%AWESOME color-y plottie-thingy
figure(4); clf; box on;
drawFluxDensity(msh, Xtot(indA), 'LineStyle', 'none'); colormap('jet'); colorbar; caxis([0 2]);
drawFluxLines(msh, Xtot(indA), 16, 'k');
axis(dims.D_so/2*[-1 1 0 1]); box on; axis tight; daspect([1 1 1]);
title('Flux lines and flux density (T)');
%set(gca,'xtick',[], 'ytick', []);