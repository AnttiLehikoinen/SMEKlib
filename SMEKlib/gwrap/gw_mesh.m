function [] = gw_mesh(gm, varargin)

if numel(varargin)
    filename = varargin{1};
else
    filename = [gm.gpath 'gm_geo.geo'];
end

system(['"' gm.gpath '"' 'gmsh ' '"' filename '"' ' -2']);

end