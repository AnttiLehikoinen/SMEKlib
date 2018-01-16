function [nu, varargout] = calculate_reluctivity(varargin)
%calculate_reluctivity updates reluctivity in each element.
%
% nu = calculate_reluctivity(B2, nu_struct)
% calculates the reluctivity nu in each element of the mesh, based on which
% the supplied nu_struct was initialized. B2 is either a 1xNe vector of squared
% flux densities in each element, or a 2xNe vector of flux density
% components.
%
% Use [nu, dnu] = calculate_reluctivity(B2, nu_struct) to also calculate
% dnu/dB^2 in the similar fashion.
% 
% Call with calculate_reluctivity(nu_struct) to assume 0 flux.

if nargin == 1
    nu_struct = varargin{1};
    B2 = zeros(1, nu_struct.Ne);
else
    B2 = varargin{1};
    nu_struct = varargin{2};
end

if (size(B2, 2) == 1) && (nu_struct.Ne > 1)
    B2 = zeros(1, nu_struct.Ne);
elseif size(B2,1) > 1
    B2 = sum(B2.*conj(B2),1);
end

if isfield(nu_struct, 'Bnu')
    %[nu, varargout] = internal_directInterPolation(B2, nu_struct);
    if nargout > 1
        [nu, varargout{1}] = internal_directInterPolation(B2, nu_struct);
        return;
    else
        nu = internal_directInterPolation(B2, nu_struct);
        return;
    end
end

if nargout > 1
    alsoDer = true;
else
    alsoDer = false;
end

nu = zeros(size(B2,1), nu_struct.Ne);
if alsoDer
    dnu = zeros(size(B2,1), nu_struct.Ne);
end

for kmat = 1:nu_struct.N_mat
    nu(:, nu_struct.nu_cell{kmat, 3} ) = ppval( nu_struct.nu_cell{kmat, 1}, B2(:, nu_struct.nu_cell{kmat, 3}) );
    if alsoDer
        dnu(:, nu_struct.nu_cell{kmat, 3} ) = ppval( nu_struct.nu_cell{kmat, 2}, B2(:, nu_struct.nu_cell{kmat, 3}) );
    end
end

if alsoDer
    varargout{1} = dnu;
end

end

function [nu, varargout] = internal_directInterPolation(B2, nu_struct)

mu0 = pi*4e-7;

if nargout > 1
    alsoDer = true;
else
    alsoDer = false;
end

if isa(B2, 'gpuArray')
    nu = zeros(size(B2,1), nu_struct.Ne, 'gpuArray');
else
    nu = zeros(size(B2,1), nu_struct.Ne);
end
    
if alsoDer
    if isa(B2, 'gpuArray')
        dnu = zeros(size(B2,1), nu_struct.Ne, 'gpuArray');
    else
        dnu = zeros(size(B2,1), nu_struct.Ne);
    end
end

for kmat = 1:nu_struct.N_mat
    nu(:, nu_struct.mats{kmat} ) = interp1( nu_struct.Bnu{kmat}(:,1), nu_struct.Bnu{kmat}(:,2), B2(:, nu_struct.mats{kmat}), 'linear', 1/mu0);
    if alsoDer
        dnu(:, nu_struct.mats{kmat} ) = interp1( nu_struct.Bdnu{kmat}(:,1), nu_struct.Bdnu{kmat}(:,2), B2(:, nu_struct.mats{kmat}), 'linear', 0);
    end
end

if alsoDer
    varargout{1} = dnu;
end

end
