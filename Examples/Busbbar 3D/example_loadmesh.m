%Load busbar mesh.
%
% Ugly and non-general implementation, to be fixed.
% 
% (c) 2018 Antti Lehikoinen / Aalto University

addpath(genpath('..\..\SMEKlib'));

filename = 'busbar.msh';

fid = fopen(filename);

%skip unnecessary stuff until the line '$PhysicalNames'
while true
    dline = fgetl(fid);
    if any(dline) && dline(2) == 'P'
        break;
    end
end

%scanning physical names
N_physicalNames = fscanf(fid, '%d', [1 1]);

Surfaces = SLContainer();

fgetl(fid);
for k = 1:N_physicalNames
    %temp = fscanf(fid, '%d %d %s', [1 3])
    dline = fgetl(fid);
    temp = sscanf(dline, '%d %d %999999c');
    Name = char(temp(4:(end-1)))';
    id = temp(2);
    Surfaces.add(Name, id);
end

%skipping lines again
while true
    dline = fgetl(fid);
    if any(dline) && dline(1) ~= '$'
        break;
    end
end
Np = sscanf(dline, '%d');
temp = fscanf(fid, '%d %f %f %f', [4 Np])';
p = temp(:, 2:4)';

%skipping lines again
while true
    dline = fgetl(fid);
    if any(dline) && dline(1) ~= '$'
        break;
    end
end

%scanning elements
Ne = sscanf(dline, '%d');

%scanning named surface triangles, if any
fpos = ftell(fid);
triangles = zeros(3, Ne);
l_inEntity = zeros(1, Ne);
for ke = 1:Ne
    ldata = sscanf(fgetl(fid), '%d');
    if ldata(2) == 2
        triangles(:,ke) = ldata((end-2):end);
        l_inEntity(ke) = ldata(4);
        fpos = ftell(fid);
    else
        fseek(fid, fpos, 'bof');
        break;
    end
end
triangles = triangles(:, 1:(ke-1));
l_inEntity = l_inEntity(1:(ke-1));
l_entities = unique(l_inEntity);
Ne = Ne - (ke-1);

%named-something
names = Surfaces.keys();
for kn = 1:numel(names)
    if ~ismember(Surfaces.get(names{kn}), l_entities)
        %physical group is not a line
        continue;
    end
    
    els = find( l_inEntity == Surfaces.get(names{kn}) );
    if any(els)
        Surfaces.set(names{kn}, triangles(:, els));
    end
    names{kn} = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scanning tetras

%scanning elements
try
    %assuming a fixed number (=2) of tags
    temp = fscanf(fid, '%d %d %d %d %d %d %d %d %d', [9 Ne]);
    t_inEntity = temp(4,:);
    t = temp(6:9, :);
catch
    %something failed; looping over elements 1-by-1
    t = zeros(3, Ne);
    t_inEntity = zeros(1, Ne);    
    for ke = 1:Ne
        dline = fgetl(fid);
        temp = strsplit(dline);
        t_inEntity(ke) = str2double(temp{4});
        t(:,ke) = [str2double(temp{end-2}); str2double(temp{end-1}); str2double(temp{end})];
    end
end
fclose(fid);

for kn = 1:numel(names)
    if isempty(names{kn})
        continue;
    end

    els = find( t_inEntity == Surfaces.get(names{kn}) );
    if any(els)
        Surfaces.set(names{kn}, els);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_out = unique(Surfaces.get('Outer'));

figure(1); clf; hold on;
%tetramesh(t', p'); %plotting mesh, SLOW
plot3(p(1,n_out), p(2,n_out), p(3,n_out), 'ko'); %plotting boundary nodes for verification