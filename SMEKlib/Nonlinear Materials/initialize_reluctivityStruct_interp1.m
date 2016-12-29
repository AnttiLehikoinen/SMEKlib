function nu_struct = initialize_reluctivityStruct_interp1(msh, useDefaultMaterials, varargin)
%initialize_reluctivityStruct_interp1 initializes a struct for reluctivity
% computation.
% 
% Call syntax identical to initialize_reluctivityStruct, but linear
% interpolation is used instead of splines.
% Robustness is also somewhat improved.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

N_INTERPOLATIONPOINTS_1 = 50; %number of interpolation points used for the available BH curve
N_INTERPOLATIONPOINTS_2 = 100; %number of interpolation points AFTER the BH data, up to B2_MAX T^2
B2_MAX = 1000; %max value of B^2 considered

%materials present in the model
mats = unique(msh.matel);
Nmats = numel(mats);

if useDefaultMaterials
    N_defaultMats = get_defaultMaterials(0); %number of default materials available
    N_defMatsUsed = numel(intersect(mats, 1:N_defaultMats)); %number of default materials used
else
    N_defMatsUsed = 0;
end

Bnu = cell(Nmats, 1);
Bdnu = cell(Nmats, 1);
mats_cell = cell(Nmats, 1);
for kmat = 1:Nmats
    matIndex = mats(kmat);
    if kmat <= N_defMatsUsed
        BH = get_defaultMaterials(matIndex);
    else 
        BH = varargin{kmat - N_defMatsUsed};
    end
    
    %getting nu and extended (extrapolated) range of BH data
    [nu, BH2] = internal_getNu(BH);
    Bmax = max(BH(:,1));
    
    
    %initializing temporary interpolation splines (for easy derivation)
    pp = spline(BH2(:,1).^2, nu); %spline for nu(B^2)
    %pp = pchip(BH2(:,1).^2, nu); %these might also be worth a try
    %pp = csaps(BH2(:,1).^2, nu);

    ppder = derivate_pp(pp); %spline for d/dB^2 nu(B^2)
    
    %interpolation points (T^2)
    B2samples = [linspace(0, Bmax.^2, N_INTERPOLATIONPOINTS_1)  internal_linspace(Bmax^2, 10, ceil(N_INTERPOLATIONPOINTS_2)/2)...
        aux_intlogspace(10, B2_MAX, ceil(N_INTERPOLATIONPOINTS_2)/2)];
    
    %interpolated values for nu and dnu/dB^2
    Bnu{kmat} = [B2samples;  ppval(pp, B2samples)]';
    Bdnu{kmat} = [B2samples;  max(ppval(ppder, B2samples), 0)]';
    
    mats_cell{kmat} = find( msh.matel == matIndex );
end

nu_struct = struct('Ne', size(msh.t,2), 'N_mat', Nmats);
nu_struct.Bnu = Bnu;
nu_struct.Bdnu = Bdnu;
nu_struct.mats = mats_cell;

end

function y = internal_linspace(x0, x1, N)
%internal_linspace returns N points in ]x0, x1[

if (x1-x0)<=0
    y = [];
else
    y = linspace(x0 + (x1-x0)/N, x1, N);
end

end


function [nu, BH2] = internal_getNu(BH)
%returns reluctivity

mu0 = pi*4e-7;

BH2 = internal_extendBH(BH);

nu_init = BH(:,2) ./ BH(:,1); nu_init(1) = nu_init(2);

nu = min(interp1(BH(:,1), nu_init, BH2(:,1), 'linear', 'extrap'), 1/mu0);

end

function BH = internal_extendBH(BH)
% linearly extends the BH curve to higher flux densities

dB = BH(end, 1) - BH(end-1, 1);

mu0 = pi*4e-7;

%number of extra points
Next = 100;

B_ext = [BH(:,1)' BH(end,1)+(1:Next)*dB];
H_ext = min(interp1(BH(:,1), BH(:,2), B_ext, 'linear', 'extrap'), B_ext/mu0);

BH2 = [B_ext' H_ext'];

BH = BH2;
end