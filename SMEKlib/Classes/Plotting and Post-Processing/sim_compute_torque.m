function [Torque, Force] = sim_compute_torque(sim, pars, simtype)
%sim_compute_torque Torque computation.
%
% Computing forces/torque with the weighted stress tensor method, in the
% non-moving part of the AirgapTriangulation mesh.
%
% Limitations: Only works for two-layer triangulations now with non-curved
% elements.
%
% (c) 2018 Antti Lehikoinen / Aalto University

%initializing
if strcmp(simtype, 'static')
    A = sim.results.Xs(1:sim.Np, :);
elseif strcmp(simtype, 'stepping');
    A = sim.results.Xt(1:sim.Np,:);
elseif strcmp(simtype, 'harmonic');
    A = sim.results.Xh(1:sim.Np,:);
end
Nsamples = size(A,2);

rotorDisplacements = pars.rotorDisplacement;

Torque = zeros(1, Nsamples);
Force = zeros(2, Nsamples);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msh_ag = SimpleMesh(sim.msh.bandData.p_virt, sim.msh.bandData.t_const); %mesh for integration
n_bnd_moving = intersect( sim.msh.bandData.t_const, sim.msh.bandData.t_moving ); %boundary nodes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computing torque
for ks = 1:Nsamples
    %initializing/recomputing stuff
    if ks == 1 || ~isempty(rotorDisplacements)
        if ~isempty(rotorDisplacements)
            sim.msh.bandData.setEccentricity(rotorDisplacements(:,ks));
            msh_ag = SimpleMesh(sim.msh.bandData.p_virt, sim.msh.bandData.t_const);
        end
        
        % computing shape function values
        Nea = size(msh_ag.t, 2);
        [F,F0] = mappingTerms(msh_ag); detFa = mappingDeterminant(F);

        %pre-computing shape functions and indices
        Nf = size(msh_ag.t,1);
        [x_quad, w_quad] = get_2DtriangleIntegrationPoints(Nf+1);
        N_quad = numel(w_quad);
        Ngrada = cell(N_quad, Nf); N = Nodal2D(Operators.grad);
        x_global = zeros(2, Nea, N_quad);

        for k_quad = 1:N_quad
            x_global(:,:,k_quad) = bsxfun(@plus, mappingTimesVector(x_quad(:,k_quad), false, false, F), F0); %global coord for quad. point
            for k_shape = 1:Nf
                Ngrada{k_quad, k_shape} = N.eval(k_shape, x_quad(:,k_quad), msh_ag, F, detFa);
            end
        end

        %determining weighting function
        Np = size(msh_ag.p, 2);
        g = zeros(Np, 1);
        g(n_bnd_moving) = 1;
        g_grad = cell(N_quad, Nf); [g_grad{:}] = deal(zeros(2, size(msh_ag.t,2)));
        for k_quad = 1:N_quad
            for k_shape = 1:Nf
                g_grad{k_quad, k_shape} = g_grad{k_quad, k_shape} + ...
                    bsxfun(@times, Ngrada{k_quad, k_shape}, ...
                        transpose(g(msh_ag.t(k_shape, :), :)));
            end
        end
    end

    Ahere = sim.msh.bandData.tag_solution(A(:,ks));

    
    %integrating torque
    for k_quad = 1:numel(w_quad)
        %computing B
        B = zeros(2, Nea);
        for k = 1:Nf
            %inds = mshh.bandData.el_table(2, msh_ag.t(k,:) );
            %taking periodicity of solution into account
            %coeffs = mshh.bandData.el_table(3, msh_ag.t(k,:) );
            %B = B + bsxfun(@times, [0 1;-1 0]*Ngrada{k_quad, k}, ...
            %    coeffs.*transpose(A(inds, ks)));
            B = B + bsxfun(@times, [0 1;-1 0]*Ngrada{k_quad, k}, ...
                transpose(Ahere(msh_ag.t(k,:), :)) );
        end
        Babs2 = sum(conj(B).*B, 1);
        
        Tau = [conj(B(1,:)).*B(1,:)-0.5*Babs2;
            conj(B(2,:)).*B(1,:);
            conj(B(1,:)).*B(2,:);
            conj(B(2,:)).*B(2,:)-0.5*Babs2];
        
        for kf = 1:Nf
            Force(:,ks) = Force(:,ks) + w_quad(k_quad)*sum( ...
                bsxfun(@times, mappingTimesVector(g_grad{k_quad, kf}, false, false, Tau), abs(detFa)) , 2);
            
            Torque(:,ks) = Torque(:,ks) + w_quad(k_quad)*sum( ...
                crossProduct(x_global(:,:,k_quad), ...
                mappingTimesVector(g_grad{k_quad, kf}, false, false, Tau)).*abs(detFa), 2);
        end
    end
    
end

mu0 = pi*4e-7;
Torque = Torque / mu0  * sim.dims.leff * sim.msh.symmetrySectors;
Force = Force / mu0 * sim.dims.leff * sim.msh.symmetrySectors;

%reseting eccentricity
if ~isempty(rotorDisplacements)
    sim.msh.bandData.setEccentricity( [0;0] );
end

end