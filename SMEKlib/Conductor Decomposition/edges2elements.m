function elems = edges2elements(msh, edgeDefs)
% given a definition of edges in the format [start_nodes; end_nodes]
% this function returns indices to the first element each of these edges
% belongs to, as defined in the mesh struct msh


edgeList = edgeDefs2edges(msh, edgeDefs);
elems = msh.e2t(1, edgeList);


end