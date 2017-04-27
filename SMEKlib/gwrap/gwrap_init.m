function [] = gwrap_init(gmspath, varargin)
%gmshWrap_init initializes the gmsh wrapper tool - a minimal wrapper for
%gmsh.
% 
% gwrap_init(gmspath) performs basic initialization. The input gmspath
% points to the location of your gmsh installation. 
% Output files are
%  saved under 
%   gmshpath/gmsh_temp_geometry.geo and
%   gmshpath/gmsh_temp_geometry.msh.
%
% Alternatively, you can specify the output file name by calling
% gwrap_init(gmspath, filename)
% 
% (c) 2017 Antti Lehikoinen / Aalto University

global gwrap_fid gwrap_Npoint gwrap_Nline gwrap_Nlineloop gwrap_Nsurface ...
    gwrap_PhysicalSurfaces gwrap_filename

global gwrap_gmshpath
gwrap_gmshpath = gmspath;
if gwrap_gmshpath(end)~='\'
    gwrap_gmshpath = [gwrap_gmshpath '\'];
end

if numel(varargin)
    gwrap_filename = [gwrap_gmshpath varargin{1}];
else
    gwrap_filename = [gwrap_gmshpath 'gmsh_temp_geometry.geo'];
end

gwrap_fid = fopen(gwrap_filename, 'w');
gwrap_Npoint = 1;
gwrap_Nline = 1;
gwrap_Nlineloop = 1;
gwrap_Nsurface = 1;

gwrap_PhysicalSurfaces = containers.Map();

end