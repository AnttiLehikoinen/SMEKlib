function [p_plot, plotArgs] = aux_Plotting_parseInput(msh, varargin)
%aux_Plotting_parseInput Auxiliary function for plotting.
%
% Parses input data depending on if the mesh msh has a list of rotor
% elements.
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

if isfield(msh, 'rotel') && numel(varargin) && isnumeric(varargin{1}) && isscalar(varargin{1})
    %taking rotation into account if necessary
    rotorAngle = varargin{1};
    
    if numel(varargin) > 1
        plotArgs = {varargin{2:end}};
    else
        plotArgs = cell(0,0);
    end
    
    p_plot = msh.p;
    rotorNodes = unique( [msh.t(1,msh.rotel) msh.t(2,msh.rotel) msh.t(3,msh.rotel) ] );
    p_plot(:,rotorNodes) = [cos(rotorAngle) -sin(rotorAngle);sin(rotorAngle) cos(rotorAngle)] * p_plot(:,rotorNodes);
else
    p_plot = msh.p;
    plotArgs = varargin;
end

end