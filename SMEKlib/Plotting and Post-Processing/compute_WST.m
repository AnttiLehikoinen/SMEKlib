function [Force, Torque] = compute_WST(A, msh, arg1, varargin)
%compute_WST computes torque and force with Weighted Stress Tensor
% 
% Call syntax:
% [Force, Torque] = compute_WST(A, msh, arg1, varargin)
% where 
%   arg1 = list of elements to compute the forces from OR
%   arg1 = a triangulation for the force computation
% 
% For article reference, see:
% "The eggshell method for the computation of electromagnetic forces on 
% rigid bodies in 2D and 3D"
% 
% (c) 2017 Antti Lehikoinen / Aalto University

mu0 = pi*4e-7;

if numel(varargin)
    %computing torque for a time-dependent problem
    [Force, Torque] = internal_timedependentTorque(A, msh, varargin{1});
    return
end

%stationary problem (or time-harmonic); with the air-gap triangulation
%fully encompassing the rotor

if min(size(arg1)) == 1
    el_ag = arg1;
    msh_comp = msh;
else
    msh_comp = struct('t', arg1, 'p', msh.p);
    el_ag = 1:size(arg1, 2);
end


%get boundary (and free=middle) nodes of the air-gap triangulation
[~, nb_stator, n_free] = get_boundaryNodes(msh, msh_comp, el_ag);

%determining weighting function
Np = size(msh_comp.p, 2);
Ne = size(msh_comp.t, 2);

g = zeros(Np, 1);
g(nb_stator) = 1;
if any(n_free)
    S_s = assemble_matrix('grad', 'nodal', 'grad', 'nodal', 1, [], msh_comp, []);
    S = sparseFinalize(S_s, size(msh_comp.p,2), size(msh_comp.p,2));
    g(n_free) = S(n_free, n_free) \ (-S(n_free, nb_stator)*g(nb_stator));
end

%computing torque and force
[F, F0] = mappingTerms(msh_comp, el_ag);
detF = mappingDeterminant(F);

[X_quad, W_quad] = get_2DtriangleIntegrationPoints(1);
N_quad = numel(W_quad);

Torque = zeros(1, size(A,2));
Force = zeros(2, size(A,2));

for kcol = 1:size(A,2)
    for k_quad = 1:N_quad
        B = zeros(2, Ne);
        g_grad = zeros(2, Ne);
        x_global = mappingTimesVector(X_quad(:,k_quad), false, false, F) + F0;

        %computing flux density
        for kn = 1:3
            B = B + bsxfun(@times, evaluate_ShapeFunction('curl', 'nodal', kn, F, detF), ...
                transpose(A(msh_comp.t(kn, el_ag), kcol)));
            g_grad = g_grad + bsxfun(@times, evaluate_ShapeFunction('grad', 'nodal', kn, F, detF), ...
                transpose(g(msh_comp.t(kn, el_ag), :)));
        end
        Babs2 = sum(conj(B).*B, 1);
        
        %assembling stress tensor
        Tau = [conj(B(1,:)).*B(1,:)-0.5*Babs2;
            conj(B(2,:)).*B(1,:);
            conj(B(1,:)).*B(2,:);
            conj(B(2,:)).*B(2,:)-0.5*Babs2];

        Force(:,kcol) = Force(:,kcol) + W_quad(k_quad)*sum( ...
            bsxfun(@times, mappingTimesVector(g_grad, false, false, Tau), abs(detF)) , 2);
        
        Torque(:,kcol) = Torque(:,kcol) + W_quad(k_quad)*sum( ...
            crossProduct(x_global, mappingTimesVector(g_grad, false, false, Tau)).*abs(detF), 2);
    end
end

Force = real(Force) / mu0;
Torque = real(Torque) / mu0;

end

function [nb_rotor, nb_stator, n_free] = get_boundaryNodes(msh, msh_comp, el_ag)
% try to determine boundary nodes of the air-gap triangulation

%determining boundary nodes of the air-gap triangulation
if ~isfield(msh, 'rotel')
    error('A list of rotor elements must be listed in msh.rotel');
else
    rotorElements = msh.rotel;
end
nb_rotor = toRow(intersect(msh.t(:, rotorElements), msh_comp.t(:, el_ag)));
if isfield(msh, 'statel')
    nb_stator = toRow(intersect(msh.t(:, msh.statel), msh_comp.t(:, el_ag)));
else
    warning('No list of stator elements supplied in msh.statel. For PERIODIC problems, results may be incorrect.');
    [edges, e2t, ~] = getEdges(msh_comp.t(:, el_ag));
    nb_all = toRow( unique(edges(:, ~e2t(2,:))) ); %all boundary nodes
    nb_stator = toRow(setdiff(nb_all, nb_rotor));
end

n_free = toRow( setdiff( msh_comp.t(:, el_ag), [nb_rotor nb_stator] ) );

end

function [Force, Torque] = internal_timedependentTorque(A, msh, rotorAngles)

mu0 = pi*4e-7;

%determining boundary nodes based on the initial position
[t_ag, ~, pnew] = updateRotorPosition(rotorAngles(1), msh);
msh_comp.t = t_ag;
msh_comp.p = pnew;
el_ag = 1:size(t_ag,2);

Np = size(msh.p, 2);
Ne = size(t_ag, 2);

%getting boundary nodes
[nb_rotor, nb_stator, n_free] = get_boundaryNodes(msh, msh_comp, el_ag);

%determining weighting function
g = zeros(Np, 1);
g(nb_rotor) = -1;

Nsamples = numel(rotorAngles);


Torque = zeros(1, Nsamples);
Force = zeros(2, Nsamples);
[X_quad, W_quad] = get_2DtriangleIntegrationPoints(1);
N_quad = numel(W_quad);

for ksample = 1:Nsamples
    if isobject(msh.bandData)
        %[t_local, pnew, t_global] = msh.bandData.t_ag(rotorAngles(ksample));
        t_ag = msh.bandData.t_const;
        pnew = msh.bandData.p_virt;
        %msh_comp.t_global = t;
    else
        [t_ag, ~, pnew] = updateRotorPosition(rotorAngles(ksample), msh);
    end
    msh_comp.t = t_ag;
    msh_comp.p = pnew;
    
    %{
    [t_ag, ~, pnew] = updateRotorPosition(rotorAngles(ksample), msh);
    msh_comp.t = t_ag;
    msh_comp.p = pnew;
    %}
    
    [F, F0] = mappingTerms(msh_comp);
    detF = mappingDeterminant(F);
    
    %updating weighting function if required
    if any(n_free)
        S_s = assemble_matrix('grad', 'nodal', 'grad', 'nodal', 1, [], msh_comp, []);
        %handling periodicity if needed
        if isfield(msh, 'symmetrySectors')
            inds = 1:(S_s.ri-1);
            S_s.I(inds) = msh.bandData.el_table(2, S_s.I(inds));
            S_s.J(inds) = msh.bandData.el_table(2, S_s.J(inds));
        end        
        S = sparseFinalize(S_s, size(msh_comp.p,2), size(msh_comp.p,2));
        g(n_free) = S(n_free, n_free) \ (-S(n_free, nb_stator)*g(nb_stator));
    end
    
    for k_quad = 1:N_quad
        B = zeros(2, Ne);
        g_grad = zeros(2, Ne);
        x_global = mappingTimesVector(X_quad(:,k_quad), false, false, F) + F0;
        
        %computing flux density
        for kn = 1:3
            if isfield(msh, 'symmetrySectors')
                A_true = transpose(A( msh.bandData.el_table(2, msh_comp.t_ag(kn, :)), ksample)) ...
                    .* msh.bandData.el_table(3, :);
                g_true = transpose(g( msh.bandData.el_table(2, msh_comp.t_ag(kn, :)) ));
            else
                A_true = transpose( A(msh_comp.t(kn, :), ksample ) );
                g_true = transpose( g(msh_comp.t(kn, :) ) );
            end
            B = B + bsxfun(@times, evaluate_ShapeFunction('curl', 'nodal', kn, F, detF), ...
                A_true);
            g_grad = g_grad + bsxfun(@times,  evaluate_ShapeFunction('grad', 'nodal', kn, F, detF), g_true);
        end
        Babs2 = sum(conj(B).*B, 1);
        
        %assembling stress tensor
        Tau = [conj(B(1,:)).*B(1,:)-0.5*Babs2;
            conj(B(2,:)).*B(1,:);
            conj(B(1,:)).*B(2,:);
            conj(B(2,:)).*B(2,:)-0.5*Babs2];

        Force(:,ksample) = Force(:,ksample) + W_quad(k_quad)*sum( ...
            bsxfun(@times, mappingTimesVector(g_grad, false, false, Tau), abs(detF)) , 2);
        
        Torque(:,ksample) = Torque(:,ksample) + W_quad(k_quad)*sum( ...
            crossProduct(x_global, mappingTimesVector(g_grad, false, false, Tau)).*abs(detF), 2);
        
    end
end
Torque = Torque / mu0;
Force = Force / mu0;

end
