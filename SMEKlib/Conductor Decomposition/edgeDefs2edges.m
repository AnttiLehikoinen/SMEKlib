function edgeList = edgeDefs2edges(msh, edgeDefs)
% given a definition of edges in the format [start_nodes; end_nodes]
% this function returns indices to corresponding edges as defined in the
% mesh struct msh

[~, edgeList] = ismember(sort(edgeDefs,1)', msh.edges', 'rows');

end