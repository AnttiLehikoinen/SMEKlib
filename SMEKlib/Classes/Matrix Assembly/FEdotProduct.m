function E = FEdotProduct(u, v, varargin)
%dotProduct vectorized dot product for FE matrix assembly.
%
% Call syntax:
%   dotProduct(N1, N2)
%   dotProduct(N1, v, N2)

if ~numel(varargin)
    E = sum( u.*v, 1);
elseif size(v,1)<=3
    %E = sum(u.*(varargin{1}*v), 1); %assuming broadcasting here
    E = sum(u.*bsxfun(@times, v, varargin{1}), 1);
else
    %tensor material supplied
    E = sum(u.*matrixTimesVector(v, varargin{1}, false, false), 1);   
    %error('Not implemented or something.');
end

end