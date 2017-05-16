function [] = gwrap_mesh(varargin)
%gwrap_mesh
%
% gwrap_mesh() : calling after gwrap_init()
% gwrap_mesh(gmsh_path, geo_path)
% (c) 2017 Antti Lehikoinen / Aalto University

global gwrap_gmshpath gwrap_filename
if numel(varargin)
   gmsh_path = varargin{1};
   if gmsh_path(end)~='\'
        gmsh_path = [gmsh_path '\'];
    end
   
   geo_path = varargin{2};
   system(['"' gmsh_path '"' 'gmsh ' '"' geo_path '"' ' -2']);
else
    system(['"' gwrap_gmshpath '"' 'gmsh ' '"' gwrap_filename '"' ' -2']);
end

end