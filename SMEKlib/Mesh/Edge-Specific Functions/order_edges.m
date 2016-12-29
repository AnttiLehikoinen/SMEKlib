function sorted_edges = order_edges(bnd_edges)
%order_edges orders edges.
% 
% When given a 2xN array of edges (defined by end-nodes), function
% sorted_edges = order_edges(bnd_edges)
% tries to order the edges in such a way that a polygon is formed.
% Generally, it only works then the given edges can be formed into a closed
% loop.
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University


Ne = size(bnd_edges, 2);
bnd_nodes = unique(bnd_edges(:));

%checking for odd elements (not used, prolly not working)
odd_node = -1;
for k_n = 1:numel(bnd_nodes)*0
    inds = find( bnd_edges(:) == bnd_nodes(k_n) );
    if numel(inds) == 1
        odd_node = bnd_edges(inds);
        odd_column = floor( (inds-1)/2 ) + 1;
        break;
    end
end
if odd_node > 0
    prev_node = odd_node; prev_column = odd_column;
else
    prev_node = bnd_edges(1,1); prev_column = 1;
end


sorted_edges = zeros(2, Ne);
for ke = 1:Ne   
    new_node = bnd_edges(:, prev_column);    
    sorted_edges(:,ke) = [prev_node new_node( new_node ~= prev_node)];
    
    cands = find( bnd_edges(:) == sorted_edges(2,ke) );    
    cands = cands( (floor((cands-1)/2)+1) ~= prev_column );
    
    prev_node = bnd_edges(cands);
    prev_column = floor((cands-1)/2)+1;
end

end