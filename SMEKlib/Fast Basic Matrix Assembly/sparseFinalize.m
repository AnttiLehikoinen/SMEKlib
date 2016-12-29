function Sf = sparseFinalize(S, varargin)
%sparseFinalize assemble sparse matrix from sparse struct.
% 
% Call syntax
% S = sparseFinalize(Sf, args)
% where Sf = sparse struct and args is a list of arguments (see Matlab's
% "sparse" for examples).
%
% Copyright (c) 2013-2016 Antti Lehikoinen / Aalto University

ri = S.ri-1;
 
Sf = sparse(S.I(1:ri), S.J(1:ri), S.E(1:ri), varargin{:});

end