function [E, varargout] = assemble_Jacobian_data(nu_fun, X, elements, msh, evaluateResidual)
%assemble_Jacobian_data assembles data struct for the Jacobian matrix.
%
% Call syntax
% [E, varargout] = assemble_Jacobian_data(nu_fun, X, elements, msh, evaluateResidual)
% 
% with the acceptable combinations:
% 
% E = assemble_Jacobian_data(nu_fun, X, elements, msh, false)
% [E, I, J] = assemble_Jacobian_data(nu_fun, X, elements, msh, false)
% [E, E_res] = assemble_Jacobian_data(nu_fun, X, elements, msh, true)
% [E, I, J, E_res, I_res] = assemble_Jacobian_data(nu_fun, X, elements, msh, true)
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

%no elements listed --> going over all
if ~any(elements)
    elements = 1:size(msh.t,2);
end

Nout = nargout;

if evaluateResidual && (Nout==5)
    assembleIndices = true;
elseif ~evaluateResidual && (Nout == 3)
    assembleIndices = true;
else
    assembleIndices = false;
end

%default assumption = Jacobian symmetric
isSymmetric = true;

F = mappingTerms(msh, elements);
detF = mappingDeterminant(F);

DETF = abs(detF);

%numbers of something
Np = size(msh.p, 2);
Ne = numel(elements);
NO_f = 3;

%vector potential
A = transpose(X(1:Np));

%reference first-order shape functions
phi_ref = {  [-1 -1 1]; [1 0 0]; [0 1 0] }; %nodal

%function handle to shape function;
curlPhi = @(ind_n, X)( [0 1;-1 0]*mappingTimesVector(phi_ref{ind_n}(1,1:2)', 1, 1, F, [], detF) );

%indices
INDM = msh.t(:, elements);

%setting up quadrature points
[X_quad, W_quad] = get_2DtriangleIntegrationPoints(1);
N_quad = numel(W_quad); X_quad = [X_quad; ones(1,N_quad)];

E = zeros(1, Ne*NO_f^2);
Fcell = cell(NO_f,1);

if evaluateResidual
    E_res = zeros(1, Ne*NO_f);
end

%assembling matrix
for k_quad = 1:N_quad
    % pre-computing shape function values at the quadrature points
    for kb = 1:NO_f
        Fcell{kb} = curlPhi(kb, X_quad(:,k_quad));
    end
    
    % computing B2 at the quadrature point
    B = zeros(2, Ne);
    for kb = 1:NO_f
        B = B + bsxfun(@times, Fcell{kb}, A(INDM(kb,:)) );
    end
    
    %computing reluctivity at the quadrature point
    [nu, dnu] = nu_fun(B);
    %[nu, dnu] = nu_fun(B*0); dnu = 0*dnu; %enforced linearity
    
    if size(nu, 1) == 1
        %computing H and dH/dB
        H = bsxfun(@times, B, nu);

        dH = [nu+2*dnu.*B(1,:).*B(1,:);
            2*dnu.*B(1,:).*B(2,:);
            2*dnu.*B(2,:).*B(1,:);
            nu+2*dnu.*B(2,:).*B(2,:)];
    else
        %function returns H and dH/dB
        H = nu;
        dH = dnu;
        
        isSymmetric = false; %have to assume this for generality
    end
    
    % computing the stiffness-matrix term of the Jacobian
    ri = 1;
    for k_test = 1:NO_f
        for k_shape = k_test:NO_f
            inds = ((ri-1)*Ne + 1):(ri*Ne);
            
            Etemp = W_quad(k_quad) * ...
                dotProduct( Fcell{k_test}, mappingTimesVector(Fcell{k_shape}, false, false, dH) );
            
            %saving contribution
            E(inds) = E(inds) + Etemp;            
            
            %symmetric component?
            if k_test ~= k_shape
                if isSymmetric
                    ri = ri + 1; inds = ((ri-1)*Ne + 1):(ri*Ne);
                    E(inds) = E(inds) + Etemp;
                else
                    ri = ri + 1; inds = ((ri-1)*Ne + 1):(ri*Ne);
                    E(inds) = E(inds) + W_quad(k_quad) * ...
                        dotProduct( Fcell{k_shape}, mappingTimesVector(Fcell{k_test}, false, false, dH) );
                end
            end
            
            ri = ri + 1;
        end
        
        %computing residual contribution if necessary
        if evaluateResidual
            E_res( ((k_test-1)*Ne+1):(k_test*Ne) ) = E_res( ((k_test-1)*Ne+1):(k_test*Ne) ) + ...
                W_quad(k_quad) * dotProduct( Fcell{k_test}, H) .* DETF;
        end
    end
end

E = E .* repmat(DETF, 1, NO_f^2);

%assigning residual
if evaluateResidual && assembleIndices
    varargout{3} = E_res;
elseif evaluateResidual && ~assembleIndices
    varargout{1} = E_res;
end

if ~assembleIndices
    return;
end

%assembling indices
I = zeros(1, Ne*NO_f^2); J = I;
I_res = zeros(1, Ne*NO_f);
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
    
    if evaluateResidual
        I_res( ((k_test-1)*Ne+1):(k_test*Ne) ) = INDM(k_test, :);
    end
end

varargout{1} = I; varargout{2} = J;

if evaluateResidual
    varargout{4} = I_res;
end

end