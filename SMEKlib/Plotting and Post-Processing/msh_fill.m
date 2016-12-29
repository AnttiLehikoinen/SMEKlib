function h = msh_fill(msh, elements, varargin)
%msh_patch patch plot for the mesh
% h = msh_fill(msh, elements, plotArgs)
% plots filled elements defined on the mesh msh.
%
% If the mesh msh has a field "rotel" listing the rotor elements, the call 
% syntax can be msh_fill(msh, elements, rotorAngle, plotArgs)
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

[p_plot, plotArgs] = aux_Plotting_parseInput(msh, varargin{:});

if any(elements)
else
    elements = 1:size(msh.t,2);
end
X = [p_plot(1, msh.t(1, elements)); p_plot(1, msh.t(2, elements)); p_plot(1, msh.t(3, elements))];
Y = [p_plot(2, msh.t(1, elements)); p_plot(2, msh.t(2, elements)); p_plot(2, msh.t(3, elements))];

h = patch(X, Y, plotArgs{:});

end