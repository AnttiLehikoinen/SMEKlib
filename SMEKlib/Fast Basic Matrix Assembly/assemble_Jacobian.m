function [J, varargout] = assemble_Jacobian(nu_fun, X, elements, msh)
%assemble_Jacobian assembles the Jacobian matrix.
%
% J = assemble_Jacobian(nu_fun, X, elements, msh)
% assembles the Jacobian matrix for the 2D vector potential formulation,
% with the entries
% J_ij = ...
%
% [J, r] = assemble_Jacobian(nu_fun, X, elements, msh) also computes the
% stiffness matrix contribution to the residual
% r = S*X(inds_A)
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

Nvars = size(X, 1);
%Nvars = size(msh.p, 2);

if nargout > 0
    % also computing residual
    [E, I, J, E_res, I_res] = assemble_Jacobian_data(nu_fun, X, elements, msh, true);
    
    J = sparse(I, J, E, Nvars, Nvars);
    varargout{1} = sparse(I_res, ones(size(I_res)), E_res, size(X, 1), 1);
else
    %only the Jacobian itself computed
    [E, I, J] = assemble_Jacobian_data(nu_fun, X, elements, msh, false);
    J = sparse(I, J, E, Nvars, Nvars);
end

end
        
        