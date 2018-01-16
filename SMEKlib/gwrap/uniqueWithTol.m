function varargout = uniqueWithTol(varargin)
%uniqueWithTol wrapper for uniquetol
%
% (c) 2017 Antti Lehikoinen / Aalto University

[varargout{1:nargout}] = uniquetol(varargin{:});

return
try
    [varargout{1:nargout}] = uniquetol(varargin{:});
catch
    error('Matlab built-in function uniquewithtol not available. Replacement not implemented.');
end

end