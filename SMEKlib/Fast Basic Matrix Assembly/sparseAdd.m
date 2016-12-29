function S = sparseAdd(I, J, E, S)
%sparseAdd add entries to the sparse matrix struct.
% 
% sparseAdd(I, J, E, S) adds the triplets (I,J,E) to the sparse struct S.
% Use S = [] to initialize a new struct.
%
% Copyright (c) 2013-2016 Antti Lehikoinen / Aalto University

if ~isstruct(S)
    SIZEFACTOR = 1000;
    S = struct('I', zeros(SIZEFACTOR, 1), 'J', zeros(SIZEFACTOR, 1), 'E', zeros(SIZEFACTOR, 1), ...
        'ri', 1, 'maxSize', SIZEFACTOR);
end

ri = S.ri;
addSize = numel(E);
if (ri + addSize -1) > S.maxSize
    %incrSize = 2*addSize;
    incrSize = 2*numel(S.I);
    S.I = [S.I; zeros(incrSize,1)];
    S.J = [S.J; zeros(incrSize,1)];
    S.E = [S.E; zeros(incrSize,1)];
    S.maxSize = size(S.E,1);
end

S.I(ri:(ri+addSize-1)) = I;
S.J(ri:(ri+addSize-1)) = J;
S.E(ri:(ri+addSize-1)) = E;
S.ri = ri + addSize;

end