function nu_struct = initialize_reluctivityStruct(msh, useDefaultMaterials, varargin)
%initialize_reluctivityStruct initializes a struct for reluctivity
% computation.
%
% initializes a struct of reluctivity curve data (B^2, nu) used for later
% calculations. Both the reluctivity curves and their derivatives (w.r.t.
% B^2) are initialized. Spline interpolation will be used. The mesh struct
% msh must have the field matel containing the material index for each
% element. The useDefaultMaterials flag indicates whether the default
% material curves are used.
%
% Use the syntax
% initialize_reluctivityStruct(msh, useDefaultMaterials, BH_n1, BH_n2, ...)
% to supply extra material data, with n1 < n2 ...
% If useDefaultMaterials = false, only the supplied data is used.
% If useDefaultMaterials = true, only materials with their index exceeding
% the number of default materials (=get_defaultMaterials(0)) are
% read from the supplied data, with the one with the smallest index
% corresponding to BH_n1 etc.
%
% NOTE: In general assumes very "nice" material data; BH curves extending
% to high-enough values of B, monotonicity (or close) of the (B^2,
% nu)-curves etc, etc.
% Use initialize_reluctivityStruct_interp1 for somewhat increased
% robustness.
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

%materials present in the model
mats = unique(msh.matel);
Nmats = numel(mats);

if useDefaultMaterials
    N_defaultMats = get_defaultMaterials(0); %number of default materials available
    N_defMatsUsed = numel(intersect(mats, 1:N_defaultMats)); %number of default materials used
else
    N_defMatsUsed = 0;
end

nu_cell = cell(Nmats, 3);
for kmat = 1:Nmats
    matIndex = mats(kmat);
    if kmat <= N_defMatsUsed
        BH = get_defaultMaterials(matIndex);
    else 
        BH = varargin{kmat - N_defMatsUsed};
    end

    %computing reluctivity curve
    nu = BH(:,2) ./ BH(:,1); nu(1) = nu(2);

    %initializing interpolation splines
    pp = spline(BH(:,1).^2, nu); %spline for nu(B^2)
    ppder = derivate_pp(pp); %spline for d/dB^2 nu(B^2)
    matInElements = find( msh.matel == matIndex ); %which elements have which material

    nu_cell{kmat, 1} = pp; nu_cell{kmat, 2} = ppder; nu_cell{kmat, 3} = matInElements;
end

nu_struct = struct('Ne', size(msh.t,2), 'N_mat', Nmats);
nu_struct.nu_cell = nu_cell;

end