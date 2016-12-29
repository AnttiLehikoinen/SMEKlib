function edgeList = edgeDefs2edges(msh, edgeDefs)
%edgeDefs2edges returns edge indices corresponding to edge definitions
%
% Function
% edgeDefs2edges(msh, edgeDefs)
% when given a definition of edges in the format [start_nodes; end_nodes]
% returns indices to corresponding edges as defined in the
% mesh struct msh
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

[~, edgeList] = ismember(sort(edgeDefs,1)', msh.edges', 'rows');

end