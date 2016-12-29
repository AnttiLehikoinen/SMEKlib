function elems = edges2elements(msh, edgeDefs)
%edges2elements element that an edge belongs to.
% 
% Function
% elems = edges2elements(msh, edgeDefs)
% when given a definition of edges in the format [start_nodes; end_nodes]
% returns indices to the first element each of these edges
% belongs to, as defined in the mesh struct msh
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University


edgeList = edgeDefs2edges(msh, edgeDefs);
elems = msh.e2t(1, edgeList);


end