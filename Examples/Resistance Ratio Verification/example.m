%AC resistance verification for slot-bound conductors
%
% (c) 2017 Antti Lehikoinen / Aalto University
 
addpath(genpath('..\..\SMEKlib'));
 
%condType = 'rect';
condType = 'circ';
 
w_cond = 4e-3; %conductor width
w_gap = 0.0e-3; %vertical gap between conductors
w_slot = 5e-3; %slot width
N = 6; %number of conductors
 
sigma = 5.96e7; %copper conductivity
mur_iron = 5000; %iron relative permeability
f = 50; %supply frequency
 
h_yoke = 2e-2; %yoke height
h_free = 1e-2; %free space at slot opening
h_air = 1.5e-2;
 
w_domain = 3e-2; %domain width
 
leff = 1; %length in z-direction
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating geometry and meshing
 
%generating geometry description for conductors
c_centers = linspace(0, (N-1)*(w_cond+w_gap), N);
 
figure(1); clf; hold on; box on;
 
%conductors
if strcmp(condType, 'rect')
    G_cond_r = zeros(2+4*2, N);
elseif strcmp(condType, 'circ')
    G_cond_r = zeros(4, N);
end
for k = 1:N
    if strcmp(condType, 'rect')
        G_cond_r(:, k) = [3 4 w_cond/2*[-1 1 1 -1] c_centers(k)+w_cond/2*[-1 -1 1 1]]';
    elseif strcmp(condType, 'circ')
        G_cond_r(:, k) = [1 0 c_centers(k) w_cond/2]';
    else
        error('Invalid conductor shape.')
    end
end
mPlotG(G_cond_r, 'r'); axis equal;
 
%slot box
yb = (-w_cond/2-w_gap-h_yoke);
yb_s = -w_cond/2-w_gap;
yt = c_centers(end) + w_cond/2+w_gap + h_free;
G_slot = [3 0 w_domain/2*[-1 -1] w_slot/2*[-1 -1 1 1] w_domain/2*[1 1] ...
    yb yt*[1 1] yb_s*[1 1] yt*[1 1] yb]';
G_slot(2) = (numel(G_slot)-2)/2;
 
mPlotG(G_slot, 'c');
 
%full domain box
G_box = [3 0 w_domain/2*[-1 -1 1 1] ...
    yb (yt+h_air)*[1 1] yb]';
G_box(2) = (numel(G_box)-2)/2;
mPlotG(G_box, 'k:')
 
Gtot = zeropadcat(G_box, G_slot, G_cond_r);
 
%namespace
slot_ns = ['b' 'i' char('s'*ones(1,N))];
 
%set formula
slot_sf = 'b';
 
%creating empty model and adding geometry
slot_model = createpde(1);
%adding geometry to the model
[slot_dl, slot_bt] = decsg(Gtot, slot_sf, slot_ns);
geometryFromEdges(slot_model, slot_dl);
 
%meshing
generateMesh(slot_model, 'GeometricOrder', 'linear');
 
%extracting nicer data and refining
[p, e, t] = meshToPet(slot_model.Mesh);
[p, e, t] = refinemesh(slot_dl, p, e, t);
%[p, e, t] = refinemesh(slot_dl, p, e, t, 3:(3+N-1));
 
%getting conductor elements
el_iron = get_elementsInDomain(slot_bt(:,2), t(4,:)); el_iron = el_iron{1};
el_conductors = get_elementsInDomain(slot_bt(:,3:end), t(4,:));
 
msh = inittri(p, t(1:3,:));
 
figure(2); clf; hold on; box on; axis equal;
msh_triplot(msh, [], 'g');
msh_triplot(msh, el_iron, 'k');
msh_triplot(msh, horzcat(el_conductors{:}), 'r');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembling matrices and solving
 
Np = size(msh.p, 2); Ne = size(msh.t, 2);
Dnodes = unique(msh.edges(:, ~msh.e2t(2,:)) ); msh_plot(msh, Dnodes, 'ko');
freeVars = setdiff(1:(Np+N+1), Dnodes);
 
mu0 = pi*4e-7;
nu = ones(1, Ne)/mu0; nu(el_iron) = 1/(mur_iron*mu0);
 
%stiffness matrix
S_struct = assemble_matrix('grad', 'nodal', 'grad', 'nodal', nu, [], msh, []);
S = sparseFinalize(S_struct, Np, Np); clear S_struct;
 
%mass matrix
M_struct = assemble_matrix('', 'nodal', '', 'nodal', sigma, horzcat(el_conductors{:}), msh, []);
M = sparseFinalize(M_struct, Np, Np); clear M_struct;
 
%coupling matrix between voltages and the vector potential
CF_struct = [];
for kc = 1:N
    CF_struct = assemble_vector('', 'nodal', sigma, kc, el_conductors{kc}, msh, CF_struct);
end
CF = sparseFinalize(CF_struct, Np, N);
cA = sum(CF*speye(N, N), 1);
DR = sparsediag( leff ./ cA );
 
L = ones(N, 1); %loop matrix
 
%assembling final matrix and solving
w = 2*pi*f;
Qtot = [S+1i*w*M -1/leff*CF sparse(Np, 1);
    DR*1i*w*transpose(CF) -speye(N,N) DR*L;
    sparse(1, Np) transpose(L) sparse(1, 1)];
 
X = zeros(Np+N+1, 1);
X(freeVars) = Qtot(freeVars, freeVars) \ [zeros(numel(freeVars)-1,1);1];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post-processing
 
% plotting flux lines and current density
figure(3); clf; hold on; box on;
%current density:
J = zeros(Np, 1);
for k = 1:N
    n_c = unique( msh.t(:, el_conductors{k}) );
    J(n_c) = -1i*w*X(n_c)*sigma + sigma/leff*repmat(X(Np+k), numel(n_c,1), 1);
   
    subplot(1, 2, 1); hold on;
    drawCurrentDensity(msh, 1e6*J, horzcat(el_conductors{k}));
   
    subplot(1, 2, 2); hold on; box on;
    drawCurrentDensity(msh, 1e6*J, horzcat(el_conductors{k}));
end
subplot(1, 2, 1);
colorbar;
%flux lines:
drawFluxLines(msh, real(X(1:Np)), 20, 'k');
%plotting geometry objects for clarity
mPlotG(G_cond_r, 'k');
mPlotG(G_slot, 'k');
mPlotG(G_box, 'k');
axis([w_domain/2*[-1 1] yb (yt+h_air)]);
daspect([1 1 1]);
title('Flux lines and current density');
 
 
subplot(1, 2, 2); hold on; box on;
%drawCurrentDensity(msh, J, horzcat(el_conductors{:}));
mPlotG(G_cond_r, 'k');
axis([w_slot/2*[-1 1] yb_s yt]);
daspect([1 1 1]);
title('Zoom-in of current density')
 
 
%computing conductor impedances
I = X(end); U = X((Np+1):(Np+N));
 
Z = U/I
R = full(diag(DR));

real(Z)./R %impedance ratio