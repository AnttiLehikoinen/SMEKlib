function msh = msh_removeElements(msh, els)
%msh_removeElements removes elements from mesh.
%
% TODO: update named edges
%
% (c) 2018 Antti Lehikoinen / Aalto University

%original sizes
Ne_full = size(msh.t, 2);
Np_full = size(msh.p, 2);

%removing elements from master mesh
[el_master, IA] = setdiff(1:Ne_full, els);
Ne_master = numel(el_master);
el_inds_new = zeros(1, Ne_full); el_inds_new(IA) = 1:Ne_master;

%removing in-slot nodes from master mesh
n_slave = unique(msh.t(:,els));
n_slave = toRow( setdiff( n_slave, msh.t(:, el_master) ) );

%[n_master, IA] = intersect(1:Np_full,  msh.t(:, el_master) );
[n_master, IA] = setdiff(1:Np_full, n_slave);
Np_master = numel(n_master);
n_inds_new = zeros(1, Np_full); n_inds_new(IA) = 1:Np_master;

%initializing master mesh
msh.p = msh.p(:, n_master);
msh.t = reshape( n_inds_new( msh.t(:, el_master) ), size(msh.t,1), Ne_master);
msh.matel = msh.matel(el_master);

%setting edges
[medges, me2t, mt2e] = getEdges(msh.t);
msh.edges = medges; msh.e2t = me2t;
msh.t2e = mt2e;

% going through named elements and nodes
msh = msh_updateNamedElements(msh, el_inds_new);

%nodes
keys = msh.namedNodes.keys;
for k = 1:numel(keys)
    key = keys{k};
    msh.namedNodes.set(key, n_inds_new( msh.namedNodes.get(key) ) );
end

%updating airgap triangulation
msh.bandData.agNodes_global = n_inds_new( msh.bandData.agNodes_global );
msh.bandData.el_table(2,:) = n_inds_new( msh.bandData.el_table(2,:) );

try
    %necessary for the AirgapTriangulation class
    msh.bandData.setConstantAGmatrix(Np_master);
catch
    %does not work for the struct-based bandData approach.
end
    
%edges: TODO

end
