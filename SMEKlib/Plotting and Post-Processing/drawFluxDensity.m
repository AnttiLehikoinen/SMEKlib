function h = drawFluxDensity(msh, A, varargin)
%drawFluxDensity plots the flux density
%
% drawFluxDensity(msh, A, args) plots the flux density using "fill" and the
% arguments args
%
% If the mesh msh has a field "rotel" listing the rotor elements, the call 
% syntax can be drawFluxDensity(msh, A, , rotorAngle, args)
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

[p_plot, plotArgs] = aux_Plotting_parseInput(msh, varargin{:});

if size(msh.t, 1) > 3
    % higher-order mesh
    t = aux_mesh(msh);
else
    t = msh.t;
end

X = zeros(3, size(t,2));
Y = X;

for kn = 1:3
    X(kn,:) = p_plot(1, t(kn,:));
    Y(kn,:) = p_plot(2, t(kn,:));
end

%calculating flux density
Babs = transpose(calculate_B(A, msh));

%for backwards compatibility (up to Matlab 2014 and earlier at least)
Babs = repmat(Babs', 3, 1);

h = fill(X,Y, Babs, plotArgs{:});


end