function msh = msh_updateNamedNodes(msh)
%msh_updateNamedNodes attempting to update named nodes after moving to
%2nd-order mesh.
%
% (c) 2018 Antti Lehikoinen / Aalto University

keys = {'Dirichlet', 'Periodic_master', 'Periodic_slave', 'n_ag_s', 'n_ag_r'};

edges = int_getEdges(msh, 'Dirichlet');
n_dir = toRow( edges );
n_dir_extra = setdiff(msh.namedNodes.get('Dirichlet'), n_dir);
n_dir_all = union(n_dir_extra, n_dir);
msh.namedNodes.set('Dirichlet', [n_dir_extra sortRadialEdges(msh.p, edges)]);

edges = int_getEdges(msh, 'Periodic_master');
msh.namedNodes.set('Periodic_master', setdiff(sortRadialEdges(msh.p, edges), n_dir_all));

edges = int_getEdges(msh, 'Periodic_slave');
msh.namedNodes.set('Periodic_slave', setdiff(sortRadialEdges(msh.p, edges), n_dir_all));

edges = int_getEdges(msh, 'n_ag_s');
msh.namedNodes.set('n_ag_s', sortSegmentEdges(msh.p, edges));

edges = int_getEdges(msh, 'n_ag_r');
msh.namedNodes.set('n_ag_r', sortSegmentEdges(msh.p, edges));
    
end

function e = int_getEdges(msh, key)

n = msh.namedNodes.get(key);
edge_defs = [n(1:(end-1)); n(2:end)];
edge_inds = edgeDefs2edges(msh, edge_defs);
edge_inds = edge_inds(edge_inds>0);
e = msh.edges(:, edge_inds);

end
