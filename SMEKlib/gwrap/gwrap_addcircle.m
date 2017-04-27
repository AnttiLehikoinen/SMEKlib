function [] = gwrap_addcircle(xc, yc, radius, Npoints, surfaceName)
%gwrap_addCircle adds a named circular surface to the geometry.
% 
% [] = gwrap_addCircle(xc, yc, radius, Npoints, surfaceName)
% adds a circle centered at (xc, yc) with a given radius, drawn with
% Npoints line segments, belonging to the named surface surfaceName.
% 
% See "help gwrap_addpolygon" for permissible surfaceNames.
% 
% (c) 2017 Antti Lehikoinen / Aalto University

angles = linspace(0, 2*pi*(1-1/Npoints), Npoints);
gwrap_addpolygon(xc+radius*cos(angles), yc+radius*sin(angles), surfaceName);

end