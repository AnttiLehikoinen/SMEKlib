function [Babs, B] = calculate_B(X, msh)
%calculate_B calculates flux density.
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

%{
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
%}

if size(msh.t,1) == 3
    % first-order elements
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
elseif size(msh.t,1) == 6
    xref = [0 0;1 0;0 1;0.5 0;0.5 0.5;0 0.5]';
    B = cell(6, 1); Babs = zeros(6, size(msh.t,2));
    [B{:}] = deal(zeros(2, size(msh.t,2)));
    N = Nodal2D(Operators.curl);
    F = msh.getMappingMatrix;
    detF = matrixDeterminant(F);
    for kp = 1:6
        for kf = 1:6
            B{kp} = B{kp} + bsxfun(@times, N.eval(kf, xref(:,kp), msh, F, detF), ...
                transpose(X(msh.t(kf,:))) );
        end
        Babs(kp,:) = sum(B{kp}.^2, 1).^0.5;
    end
end
        


end