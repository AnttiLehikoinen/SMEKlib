function [] = gw_writeFile(gm, varargin)

%opening file
if numel(varargin)
    filename = varargin{1};
else
    filename = [gm.gpath 'gm_geo.geo'];
end
fid = fopen(filename, 'w');

%printing points
for k = 1:gm.N_points
    if gm.p(3, k) > 0
        fprintf(fid, 'Point (%d) = {%f, %f, 0, %f};\n', k, gm.p(1, k), gm.p(2,k), gm.p(3,k));
    else
        fprintf(fid, 'Point (%d) = {%f, %f, 0};\n', k, gm.p(1, k), gm.p(2,k));
    end
end

fprintf(fid, '\n');

%printing lines
for k = 1:gm.N_lines
    fprintf(fid, 'Line (%d) = {%d, %d};\n', k, gm.l(1,k), gm.l(2,k));
end

fprintf(fid, '\n');

%printing line loops
for kll = 1:gm.N_lineloops
    fprintf(fid, 'Line Loop (%d) = {', kll);
    Np = numel(gm.ll{kll});
    for k = 1:Np
        if k < Np
            fprintf(fid, ' % d,', gm.ll{kll}(k));
        else
            fprintf(fid, ' % d', gm.ll{kll}(k));
        end
    end
    fprintf(fid, '};\n');
end

fprintf(fid, '\n');

%printing plane surfaces
names = gm.surfaces.keys;
for kn = 1:numel(names)
    if strcmp(names{kn}, 'OuterBoundary')
        surfaceLoops = gm.surfaces.get(names{kn}); surfaceLoops = surfaceLoops(1);
        fprintf(fid, 'Plane Surface (%d) = { %d', surfaceLoops, surfaceLoops);
        holeLoops = setdiff(1:gm.N_lineloops, surfaceLoops);
        for kh = 1:numel(holeLoops)
            fprintf(fid, ', %d', abs(holeLoops(kh)));
        end
        fprintf(fid, '};\n');
    else
        surfaceLoops = gm.surfaces.get(names{kn});
        for kl = 1:numel(surfaceLoops)
            surfloop = surfaceLoops(kl);
            if surfloop < 0
                % surface is a hole --> not printing anything
                continue;
            end
            %surface is not a hole, continuing normally
            
            fprintf(fid, 'Plane Surface (%d) = { %d', surfloop, surfloop);
            
            holeLoops = gm.holesInSurfaces.get( ['ll_' num2str(surfloop)] );
            for kh = 1:numel(holeLoops)
                fprintf(fid, ', %d', abs(holeLoops(kh)));
            end
            fprintf(fid, '};\n');
        end
    end
end

fprintf(fid, '\n');

%printing physical surfaces
for kn = 1:numel(names)
    surfaceLoops = gm.surfaces.get(names{kn});
    
    if any(surfaceLoops<0)
        %surface is a hole
        continue;
    end
    
    fprintf(fid, 'Physical Surface("%s") = {%d', names{kn}, surfaceLoops(1));
    for kp = 2:numel(surfaceLoops)
        fprintf(fid, ', %d', surfaceLoops(kp));
    end
    fprintf(fid, '};\n');
end

%printing physical lines
fprintf(fid, '\n');
names = gm.physicalLines.keys();
for kn = 1:numel(names)
    lines = gm.physicalLines.get(names{kn});
    
    fprintf(fid, 'Physical Line("%s") = {%d', names{kn}, lines(1));
    for kp = 2:numel(lines)
        fprintf(fid, ', %d', lines(kp));
    end
    fprintf(fid, '};\n');
end

fclose(fid);

end
    