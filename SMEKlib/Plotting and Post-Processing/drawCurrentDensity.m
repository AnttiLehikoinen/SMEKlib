function h = drawCurrentDensity(msh, J, elements, varargin)
%drawCurrentDensity draws current density.
% 
% h = drawCurrentDensity(msh, J, elements, plotArgs)
% draws the current density J on the specified elements ([] = entire mesh).
%
% If the mesh msh has a field "rotel" listing the rotor elements, the call 
% syntax can be drawCurrentDensity(msh, J, elements, rotorAngle, plotArgs)
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University



[p_plot, plotArgs] = aux_Plotting_parseInput(msh, varargin{:});

if isempty(elements)
    elements = 1:size(msh.t,2);
end

X = zeros(3, size(msh.t(:,elements),2));
Y = X;

for kn = 1:3
    X(kn,:) = p_plot(1, msh.t(kn,elements));
    Y(kn,:) = p_plot(2, msh.t(kn,elements));
end

Jplot = transpose(J);

size(X)
size(Y)
size(Jplot(msh.t(:,elements)))

colormap('jet');
if numel(varargin) == 0
   %h = fill(X,Y, Jplot(msh.t(:,elements)), 'LineStyle', 'none');
   h = patch(X,Y, Jplot(msh.t(:,elements)), 'LineStyle', 'none');
else
    h = patch(X,Y, Jplot(msh.t(:,elements)), plotArgs{:});
end


end