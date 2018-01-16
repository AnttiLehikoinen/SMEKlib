function J = assemble_Jacobian_old(nu_fun, A, elements, msh)
% J = assemble_Jacobian(nu_fun, A, elements, msh)
%
% assembles the Jacobian matrix for the 2D vector potential formulation,
% with the entries
% J_ij = ...
%

%no elements listed --> going over all
if ~any(elements)
    elements = 1:size(msh.t,2);
end

F = mappingTerms(msh, elements);
detF = mappingDeterminant(F);
DETF = abs(detF);

A = transpose(A);

%numbers of something
Np = size(msh.p, 2);
Ne = numel(elements);
NO_f = 3;

%reference first-order shape functions
phi_ref = {  [-1 -1 1]; [1 0 0]; [0 1 0] }; %nodal

%function handle to shape function
F1 = @(ind_n, X)( bsxfun(@times, 1 ./ detF, mappingTimesVector(phi_ref{ind_n}(1,1:2)', 1, 1, F) ) );

%indices
INDM = msh.t(:, elements);

%setting up quadrature points
[X_quad, W_quad] = get_2DtriangleIntegrationPoints(4);
N_quad = numel(W_quad); X_quad = [X_quad; ones(1,N_quad)];

E = zeros(1, Ne*NO_f^2); I = E; J = E;
Fcell = cell(NO_f,1);

%assembling matrix
for k_quad = 1:N_quad
    % pre-computing shape function values at the quadrature points
    for kb = 1:NO_f
        Fcell{kb} = F1(kb, X_quad(:,k_quad));
    end
    
    % computing B2 at the quadrature point
    B = zeros(2, Ne);
    for kb = 1:NO_f
        B = B + bsxfun(@times, Fcell{kb}, A(INDM(kb,:)) );
    end
    B2 = sum(B.^2,1);
    
    %computing reluctivity at the quadrature point
    [nu, dnu] = nu_fun(B2);
    nu = nu.*DETF; dnu = 2*dnu.*DETF;
    
    % computing the stiffness-matrix term of the Jacobian
    ri = 1;
    for k_test = 1:NO_f
        for k_shape = k_test:NO_f
            inds = ((ri-1)*Ne + 1):(ri*Ne);
            
            Etemp = W_quad(k_quad) * nu.*dotProduct( Fcell{k_test}, Fcell{k_shape} );
            
            %saving contribution
            E(inds) = E(inds) + Etemp;            
            %symmetric component
            if k_test ~= k_shape
                ri = ri + 1; inds = ((ri-1)*Ne + 1):(ri*Ne);
                E(inds) = E(inds) + Etemp;
            end
            ri = ri + 1;
        end
    end
    
    % computing the pure-Jacobian part
    ri = 1;
    for k_test = 1:NO_f
        for k_shape = k_test:NO_f
            
            inds = ((ri-1)*Ne + 1):(ri*Ne);
            
            Etemp = zeros(1, Ne);
            for l = 1:NO_f
                for m = 1:NO_f
                    Etemp = Etemp + dotProduct( Fcell{k_test}, Fcell{l} ) .* dotProduct( Fcell{k_shape}, Fcell{m} ) .* A(INDM(l,:)).*A(INDM(m,:));
                end
            end
            Etemp = W_quad(k_quad) * dnu.*Etemp;
            
            %saving contribution
            E(inds) = E(inds) + Etemp;
            %symmetric component
            if k_test ~= k_shape
                ri = ri + 1; inds = ((ri-1)*Ne + 1):(ri*Ne);
                E(inds) = E(inds) + Etemp;
            end
            ri = ri + 1;
        end
    end
end

%assembling indices
ri = 1;
for k_test = 1:NO_f
    for k_shape = k_test:NO_f
        inds = ((ri-1)*Ne + 1):(ri*Ne);
        
        I(inds) = INDM(k_test,:);
        J(inds) = INDM(k_shape,:);
        %symmetric component
        if k_test ~= k_shape
            ri = ri + 1; inds = ((ri-1)*Ne + 1):(ri*Ne);
            I(inds) = INDM(k_shape,:);
            J(inds) = INDM(k_test,:);
        end
        ri = ri + 1;
    end
end

J = sparse(I,J,E, Np, Np);

end