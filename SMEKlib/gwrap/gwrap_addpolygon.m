function [] = gwrap_addpolygon(x, y, surfaceName)
%gwrap_addpolygon adds a polygon to the geometry.
% 
% [] = gwrap_addpolygon(x, y, surfaceName)
% adds a named polygon-surface to the geometry under construction.
% 
% In general, the polygons must NOT overlap (or even touch). However, there are two
% exceptions, indicated by surfaceName:
%   surfaceName = 'Hole' --> a hole in the geometry, remains unmeshed.
%   surfaceName = 'OuterBoundary' --> outer boundary for the domain to be
%       meshed. Must be added LAST.
%
% (c) 2017 Antti Lehikoinen / Aalto University

global gwrap_fid gwrap_Npoint gwrap_Nline gwrap_Nlineloop gwrap_Nsurface ...
    gwrap_PhysicalSurfaces

sName = surfaceName;
fprintf(gwrap_fid, '// Describing polygon in the surface "%s"...\n\n', sName);


Np = numel(x);
nbias_point = gwrap_Npoint - 1;

%printing points
for k = 1:Np
    fprintf(gwrap_fid, 'Point (%d) = {%f, %f, 0};\n', k+nbias_point, x(k), y(k));
end
gwrap_Npoint = gwrap_Npoint + Np;
fprintf(gwrap_fid, '\n');

%printing lines
nbias_line = gwrap_Nline - 1;
for k = 1:Np
    is = k + nbias_point; %index of start node
    ie = mod(k, Np) + 1 + nbias_point; %index of end node
    fprintf(gwrap_fid, 'Line (%d) = {%d, %d};\n', k+nbias_line, is, ie);
end
gwrap_Nline = gwrap_Nline + Np;
fprintf(gwrap_fid, '\n');

%printing line loops
fprintf(gwrap_fid, 'Line Loop (%d) = {', gwrap_Nlineloop);
for k = 1:Np
    if k < Np
        fprintf(gwrap_fid, ' % d,', k + nbias_line);
    else
        fprintf(gwrap_fid, ' % d', k + nbias_line);
    end
end
fprintf(gwrap_fid, '};\n\n');

%printing surface
if strcmp(sName, 'OuterBoundary')
    %OuterBoundary presumably contains all other named surfaces. Thus, it
    %has to be defined as
    %Plane Surface (number) = {line_loop_of_OuterBoundary, other_loops};
    fprintf(gwrap_fid, 'Plane Surface (%d) = { %d', gwrap_Nsurface, gwrap_Nlineloop);
    for k = 1:(gwrap_Nlineloop-1)
        fprintf(gwrap_fid, ', %d', k);
    end
    fprintf(gwrap_fid, '};');    
elseif strcmp(sName, 'Hole')
    %pass
else
    fprintf(gwrap_fid, 'Plane Surface (%d) = { %d };\n', gwrap_Nsurface, gwrap_Nlineloop);
end

%assigning surface to a named group if necessary
if ~strcmp(sName, 'Hole')
    if gwrap_PhysicalSurfaces.isKey(sName)
        gwrap_PhysicalSurfaces(sName) = [gwrap_PhysicalSurfaces(sName) gwrap_Nsurface];
    else
        gwrap_PhysicalSurfaces(sName) = gwrap_Nsurface;
    end
end

gwrap_Nlineloop = gwrap_Nlineloop + 1;
gwrap_Nsurface = gwrap_Nsurface + 1;
end