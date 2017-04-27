function [t, p, t_inEntity, Name2id] = gwrap_finalize()
%gwrap_finalize finalizes and meshes the geometry.
%
% See "help gwrap_loadmesh" for description of the output files.
% 
% (c) 2017 Antti Lehikoinen / Aalto University

global gwrap_fid gwrap_PhysicalSurfaces gwrap_gmshpath gwrap_filename

%adding descriptions of physical surfaces
allSurfaces = gwrap_PhysicalSurfaces.keys;
Nsurfaces = numel(allSurfaces);
if Nsurfaces
    fprintf(gwrap_fid, '\n\n// Describing physical surfaces...\n');
    for k = 1:Nsurfaces
        surfName = allSurfaces{k};
        planeSurfaces = gwrap_PhysicalSurfaces(surfName);
        %fprintf(gwrap_fid, '\nPhysical Surface("%s") = {%d};\n\n', sName, gwrap_Nsurface);
        fprintf(gwrap_fid, 'Physical Surface("%s") = {%d', surfName, planeSurfaces(1));
        for kp = 2:numel(planeSurfaces)
            fprintf(gwrap_fid, ', %d', planeSurfaces(kp));
        end
        fprintf(gwrap_fid, '};\n');
    end
end
fclose(gwrap_fid);

%meshing with gmsh
disp('Meshing with gmsh...');
system(['"' gwrap_gmshpath '"' 'gmsh ' '"' gwrap_filename '"' ' -2']);

%loading mesh
fprintf('Done. Loading mesh...');
[t, p, t_inEntity, Name2id] = gwrap_loadmesh();
fprintf('Done\n');

end