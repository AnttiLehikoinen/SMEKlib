function S = assemble_PS_vector(op1, f1, v, column, pointSources, msh, S)
% S = assemble_vector(op1, f1, v, column, pointSources, msh, S) assembles
% the point-source load vector

elements = pointSources(3,:);

if numel(v) == 1
    v = v(1) * ones(1, numel(elements));
end

[F, F0] = mappingTerms(msh);
detF = mappingDeterminant(F);

%numbers of something
Np = size(msh.p, 2);
Ne = numel(elements);
NO_f1 = 3;

%restricting to the specified elements
F = F(:, elements); F0 = F0(:,elements); detF = detF(1,elements);

%coordinates of point-sources
X_PS = pointSources(1:2, :); %global coordinates
x_ps = [bsxfun(@times, mappingTimesVector(X_PS - F0, 1, 0, F), 1./detF);
    ones(1,Ne)];

%reference first-order shape functions
n_ref = { [0 -1 1;1 0 0] ; [0 -1 0;1 0 0] ; [0 -1 0;1 0 -1] }; %Nedelec
RT_ref = { [1 0 0;0 1 0]; [1 0 -1;0 1 0]; [1 0 0;0 1 -1] }; RT_ref = {RT_ref{[3 1 2]}}'; %Raviart-Thomas
phi_ref = {  [-1 -1 1]; [1 0 0]; [0 1 0] }; %nodal
V_phi_ref = { [phi_ref{1}; 0 0 0]; [phi_ref{2};0 0 0]; [phi_ref{3};0 0 0]; %2D nodal
    [0 0 0;phi_ref{1}]; [0 0 0;phi_ref{2}]; [0 0 0;phi_ref{3}] };


%defining function handle for (op1 f1)
if strcmp(f1, 'Nedelec')
    error('Not implemented')
elseif strcmp(f1, 'RT')
    error('Not implemented')
elseif strcmp(f1, 'nodal')
    if strcmp(op1, '')
        F1 = @(ind_n, X)( phi_ref{ind_n}*X );
    elseif strcmp(op1, 'grad')
        F1 = @(ind_n, X)( bsxfun(@times, 1 ./ detF, ...
            mappingTimesVector(phi_ref{ind_n}(1,1:2)', 1, 1, F) ) );
    elseif strcmp(op1, 'curl')
        F1 = @(ind_n, X)( [0 1;-1 0] * bsxfun(@times, 1 ./ detF, ...
            mappingTimesVector(phi_ref{ind_n}(1,1:2)', 1, 1, F) ) );
    else
        error('Invalid operator 1.')
    end
    INDM_1 = msh.t(:, elements);
elseif strcmp(f1, 'vector1')
    error('Not implemented')
else
    error(['Undefined function type ' f1]);
end

INDM_2 = column * ones(1, numel(elements));
NO_f2 = 1;

%initializing quadrature calculation
E = cell(NO_f1, NO_f2); [E{:}] = deal( zeros(1, numel(elements)) );

%calculating entries
for ind_1 = 1:NO_f1
    for ind_2 = 1:NO_f2
        E{ind_1, ind_2} = E{ind_1, ind_2} + F1(ind_1, x_ps);
    end
end

%adding entries to matrix struct
for ind_1 = 1:NO_f1
    for ind_2 = 1:NO_f2
        S = sparseAdd( INDM_1(ind_1,:), INDM_2(ind_2,:), ...
            v .* E{ind_1,ind_2}, S);
    end
end


end