function S = sparseAdd2(varargin)
%sparseAdd2 sparseAdd with more diverse input; yet slower.
%
% ...typically not worth using.
% 
% sparseAdd2(I, J, E, S)
% 
% sparseAdd2(S_to_add, S_existing)

if (numel(varargin) == 4) && ~isa(varargin{1}, 'struct')
    % call syntax: sparseAdd(I, J, E, S)
    I = varargin{1};
    J = varargin{2};
    E = varargin{3};
    S = varargin{4};
elseif numel(varargin) == 2
    % call syntax: sparseAdd(S_to_add, S_existing)
    S = varargin{2};
    ri = varargin{1}.ri - 1;
    I = varargin{1}.I(1:ri);
    J = varargin{1}.J(1:ri);
    E = varargin{1}.E(1:ri);
end

if ~isstruct(S)
    SIZEFACTOR = 1000;
    S = struct('I', zeros(SIZEFACTOR, 1), 'J', zeros(SIZEFACTOR, 1), 'E', zeros(SIZEFACTOR, 1), ...
        'ri', 1, 'maxSize', SIZEFACTOR);
end

ri = S.ri;
addSize = numel(E);
if (ri + addSize -1) > S.maxSize
    incrSize = 2*addSize;
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