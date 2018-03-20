function h = msh_trimesh(msh, A, elements, varargin)
%msh_trimesh 3D surface plot.
% 
% Call syntax:
% msh_trimesh(msh, A, []): plot all elements in mesh
% msh_trimesh(msh, A, elements): plot only particular elements
% msh_trimesh(msh, A, elements, plot_args) incluce e.g. plot-color
%
% If the mesh msh has a field "rotel" listing the rotor elements, the call 
% syntax can be msh_trimesh(msh, A, elements, rotorAngle, plot_args):
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

[p_plot, plotArgs] = aux_Plotting_parseInput(msh, varargin{:});

if size(msh.t, 1) > 3
    % higher-order mesh
    t = aux_mesh(msh);
else
    t = msh.t;
end
    
    
Np = size(msh.p,2);
if any(elements)
    h = trimesh( t(:,elements)', p_plot(1,:), p_plot(2,:), A(1:Np), plotArgs{:} );
else
    h = trimesh( t', p_plot(1,:), p_plot(2,:), A(1:Np), plotArgs{:} );
end

end