function [n_out, e_out] = findOuterNodes(t)
%findOuterNodes returns the outer nodes of the triangulation.
% 
% [n_out, e_out] = findOuterNodes(t) returns the outer nodes n_out of the
% triangulation t, along with the outer edges e_out
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

[edges, e2t, ~] = getEdges(t);

e_out = edges(:, ~e2t(2,:));

e_out = order_edges(e_out);
n_out = e_out(1,:);

end