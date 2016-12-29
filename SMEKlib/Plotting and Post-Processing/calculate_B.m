function [Babs, B] = calculate_B(X, msh)
%calculate_B calculates flux density.
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

if size(msh.t, 1) > 3
    % higher-order mesh
    t = aux_mesh(msh);
    msh = struct('t', t, 'p', msh.p);
end

A = X(1:size(msh.p,2));

Ne = size(msh.t, 2);
phiGrad = [-1 -1; 1 0; 0 1]';   %gradients of the reference shape functions (1st order)

F = mappingTerms(msh);
detF = mappingDeterminant(F);

dPhi = @(ind_n)( mappingTimesVector(phiGrad(:,ind_n), 1, 1, F, [], detF) );

B = zeros(2,Ne);
for kn = 1:3
    B = B + bsxfun(@times, dPhi(kn), transpose(A(msh.t(kn,:))));        
end

B = [0 1;-1 0] * B;
Babs = sqrt( dotProduct(B,B) );


end