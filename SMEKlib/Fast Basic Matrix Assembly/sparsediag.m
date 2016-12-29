function D = sparsediag(V)
%sparsediag sparse diagonal matrix.
% 
% D = sparsediag(V)
% is equivalent to sparse(diag(V))
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

Np = numel(V);
D = sparse(1:Np, 1:Np, V, Np, Np);

end