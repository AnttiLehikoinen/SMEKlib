function S = assemble_matrix(op1, f1, op2, f2, v, elements, msh, S)
% S = assemble_matrix returns the matrix struct S
%   
%   S = assemble_matrix(op1, f1, op2, f2, v, elements, msh, S) assembles 
%   the sparse matrix struct S with the entries
%   S_ij = Int( v (op1 f1).(op2 f2) ),
%   integrated over the "elements".
%   
%   Possible input:
%   op1,op2 = '', 'grad', 'curl' or 'div'
%   f1,f2 = 'nodal', 'vector1', 'Nedelec', 'RT'
%   v = either a separate value for each element, or a constant to be used
%   with all elements
%   elements = either a list, or [] to evalute over the entire mesh
%   msh = the mesh struct
%   S = existing matrix struct, or [] to create new ones
%
% Copyright (c) 2013-2016 Antti Lehikoinen / Aalto University

%no elements listed --> going over all
if ~any(elements)
    elements = 1:size(msh.t,2);
end

%size of v(x) too small --> assuming constant value
if numel(v) == 1
    v = v(1) * ones(1, numel(elements));
elseif numel(v) ~= numel(elements);
    v = v(elements);
end

F = mappingTerms(msh, elements);
detF = mappingDeterminant(F);
DETF = abs(detF);

%numbers of something
Np = size(msh.p, 2);
Ne = numel(elements);
NO_f1 = 3;
NO_f2 = 3;

%restricting to the specified elements
%F = F(:, elements); detF = detF(1,elements); DETF = DETF(1,elements);

%reference first-order shape functions
n_ref = { [0 -1 1;1 0 0] ; [0 -1 0;1 0 0] ; [0 -1 0;1 0 -1] }; %Nedelec
RT_ref = { [1 0 0;0 1 0]; [1 0 -1;0 1 0]; [1 0 0;0 1 -1] }; RT_ref = {RT_ref{[3 1 2]}}'; %Raviart-Thomas
phi_ref = {  [-1 -1 1]; [1 0 0]; [0 1 0] }; %nodal
V_phi_ref = { [phi_ref{1}; 0 0 0]; [phi_ref{2};0 0 0]; [phi_ref{3};0 0 0]; %2D nodal
    [0 0 0;phi_ref{1}]; [0 0 0;phi_ref{2}]; [0 0 0;phi_ref{3}] };


%defining function handle for (op1 f1)
if strcmp(f1, 'Nedelec')
    if strcmp(op1, '')
        F1 = @(ind_n, X)( bsxfun(@times, msh.edgeDirections(ind_n, elements) ./ detF, ...
            mappingTimesVector(n_ref{ind_n}*X, 1, 1, F) ) );
    elseif strcmp(op1, 'curl')
        F1 = @(ind_n, X)( 2 ./ detF .* msh.edgeDirections(ind_n, elements) );
    elseif strcmp(op1, 'div')
        if strcomp(f2, 'nodal') && strcmp(op2, '')
            error('This seems to have been a dead end')
            S = - assemble_matrix('', 'Nedelec', 'grad', 'nodal', v, elements, msh, S) + ...
                assemble_boundary_matrix('', 'nodal', '', 'Nedelec', v, elements, msh, S);
            return;
        else
            error('Invalid operator 1.');
        end
    end
    INDM_1 = msh.t2e(:, elements);
elseif strcmp(f1, 'RT')
    if strcmp(op1, '')
        F1 = @(ind_n, X)( bsxfun(@times, msh.edgeDirections(ind_n, elements) ./ detF, ...
            mappingTimesVector(RT_ref{ind_n}*X, 0, 0, F) ) );
    elseif strcmp(op1, 'curl')
        error('Curl not implemented for Raviart-Thomas elements.')
    elseif strcmp(op1, 'div')
        F1 = @(ind_n, X)( 2 ./ detF .* msh.edgeDirections(ind_n, elements) );
    end
    INDM_1 = msh.t2e(:, elements);    
elseif strcmp(f1, 'nodal')
    if strcmp(op1, '')
        F1 = @(ind_n, X)( repmat(phi_ref{ind_n}*X, 1, Ne) );
    elseif strcmp(op1, 'grad')
        F1 = @(ind_n, X)( bsxfun(@times, 1 ./ detF, ...
            mappingTimesVector(phi_ref{ind_n}(1,1:2)', 1, 1, F) ) );
        %F1 = @(ind_n, X)( mappingTimesVector(phi_ref{ind_n}(1,1:2)', 1, 1, F, [], detF) );
    elseif strcmp(op1, 'curl')
        F1 = @(ind_n, X)( [0 1;-1 0] * bsxfun(@times, 1 ./ detF, ...
            mappingTimesVector(phi_ref{ind_n}(1,1:2)', 1, 1, F) ) );
    else
        error('Invalid operator 1.')
    end
    INDM_1 = msh.t(:, elements);
elseif strcmp(f1, 'vector1')
    NO_f1 = 6;
    if strcmp(op1, '')
        F1 = @(ind_n, X)( repmat(V_phi_ref{ind_n}*X, 1, Ne) );
    elseif strcmp(op1, 'grad')
        error('Invalid operator 1.');
    elseif strcmp(op1, 'div')
        F1 = @(ind_n, X)( [1 0]*bsxfun(@times, 1 ./ detF, mappingTimesVector(V_phi_ref{ind_n}(1,1:2)', 1, 1, F) ) + ...
            [0 1]*bsxfun(@times, 1 ./ detF, mappingTimesVector(V_phi_ref{ind_n}(2,1:2)', 1, 1, F) ) );
    elseif strcmp(op1, 'curl')
        F1 = @(ind_n, X)( -[0 1]*bsxfun(@times, 1 ./ detF, mappingTimesVector(V_phi_ref{ind_n}(1,1:2)', 1, 1, F) ) + ...
            [1 0]*bsxfun(@times, 1 ./ detF, mappingTimesVector(V_phi_ref{ind_n}(2,1:2)', 1, 1, F) ) );
    else
        error('Invalid operator 1.')
    end
    INDM_1 = [msh.t(:, elements); bsxfun(@plus, msh.t(:, elements), Np)];
else
    error(['Undefined function type ' f1]);
end

%defining function handle for (op2 f2)
if strcmp(f2, 'Nedelec')
    if strcmp(op2, '')
        F2 = @(ind_n, X)( bsxfun(@times, msh.edgeDirections(ind_n, elements) ./ detF, ...
            mappingTimesVector(n_ref{ind_n}*X, 1, 1, F) ) );
    elseif strcmp(op2, 'curl')
        F2 = @(ind_n, X)( 2 ./ detF .* msh.edgeDirections(ind_n, elements) );
    else
        error('Invalid operator 2.');
    end
    INDM_2 = msh.t2e(:, elements);
elseif strcmp(f2, 'RT')
    if strcmp(op2, '')
        F2 = @(ind_n, X)( bsxfun(@times, msh.edgeDirections(ind_n, elements) ./ detF, ...
            mappingTimesVector(RT_ref{ind_n}*X, 0, 0, F) ) );
    elseif strcmp(op2, 'curl')
        error('Curl not implemented for Raviart-Thomas elements.')
    elseif strcmp(op2, 'div')
        F2 = @(ind_n, X)( 2 ./ detF .* msh.edgeDirections(ind_n, elements) );
    end
    INDM_2 = msh.t2e(:, elements);    
elseif strcmp(f2, 'nodal')
    if strcmp(op2, '')
        F2 = @(ind_n, X)( repmat(phi_ref{ind_n}*X, 1, Ne) );
    elseif strcmp(op2, 'grad')
        F2 = @(ind_n, X)( bsxfun(@times, 1 ./ detF, ...
            mappingTimesVector(phi_ref{ind_n}(1,1:2)', 1, 1, F) ) );
        %F2 = @(ind_n, X)( mappingTimesVector(phi_ref{ind_n}(1,1:2)', 1, 1, F, [], detF) );
    elseif strcmp(op2, 'curl')
        F2 = @(ind_n, X)( [0 1;-1 0] * bsxfun(@times, 1 ./ detF, ...
            mappingTimesVector(phi_ref{ind_n}(1,1:2)', 1, 1, F) ) );
    else
        error('Invalid operator 2.');
    end
    INDM_2 = msh.t(:, elements);
elseif strcmp(f2, 'vector1')
    NO_f2 = 6;
    if strcmp(op2, '')
        F2 = @(ind_n, X)( repmat(V_phi_ref{ind_n}*X, 1, Ne) );
    elseif strcmp(op2, 'grad')
        error('Invalid operator 2.');
    elseif strcmp(op2, 'div')
        F2 = @(ind_n, X)( [1 0]*bsxfun(@times, 1 ./ detF, mappingTimesVector(V_phi_ref{ind_n}(1,1:2)', 1, 1, F) ) + ...
            [0 1]*bsxfun(@times, 1 ./ detF, mappingTimesVector(V_phi_ref{ind_n}(2,1:2)', 1, 1, F) ) );
    elseif strcmp(op2, 'curl')
        F2 = @(ind_n, X)( -[0 1]*bsxfun(@times, 1 ./ detF, mappingTimesVector(V_phi_ref{ind_n}(1,1:2)', 1, 1, F) ) + ...
            [1 0]*bsxfun(@times, 1 ./ detF, mappingTimesVector(V_phi_ref{ind_n}(2,1:2)', 1, 1, F) ) );
    else
        error('Invalid operator 2.')
    end
    INDM_2 = [msh.t(:, elements); bsxfun(@plus, msh.t(:, elements), Np)];
else
    error(['Undefined function type ' f1]);
end

%initializing quadrature calculation
E = cell(NO_f1, NO_f2); [E{:}] = deal( zeros(1, numel(elements)) );
[X_quad, W_quad] = get_2DtriangleIntegrationPoints(3);
N_quad = numel(W_quad); X_quad = [X_quad; ones(1,N_quad)];

%calculating entries
for k_quad = 1:N_quad
    for ind_1 = 1:NO_f1
        for ind_2 = 1:NO_f2
            E{ind_1, ind_2} = E{ind_1, ind_2} + W_quad(k_quad) * ...
                dotProduct( F1(ind_1, X_quad(:,k_quad)), F2(ind_2, X_quad(:,k_quad)) );
        end
    end
end

%adding entries to matrix struct
for ind_1 = 1:NO_f1
    for ind_2 = 1:NO_f2
        S = sparseAdd( INDM_1(ind_1,:), INDM_2(ind_2,:), ...
            v .* E{ind_1,ind_2} .* DETF, S);
    end
end

end