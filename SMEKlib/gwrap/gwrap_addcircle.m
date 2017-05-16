function [] = gwrap_addcircle(xc, yc, radius, Npoints, varargin)
%gwrap_addCircle adds a circle to the geometry.
% 
% gwrap_addCircle(xc, yc, radius, Npoints)
% gwrap_addCircle(xc, yc, radius, Npoints, circleName)
% 
% (c) 2017 Antti Lehikoinen / Aalto University

angles = linspace(0, 2*pi*(1-1/Npoints), Npoints);

gwrap_addpolygon(xc+radius*cos(angles), yc+radius*sin(angles), varargin{:});

end